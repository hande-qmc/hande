module annihilation

use const
use fciqmc_data

implicit none

contains

    subroutine direct_annihilation(sys, rng, tinitiator, determ_flags)

        ! Annihilation algorithm.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        ! In:
        !    sys: system being studied.
        !    tinitiator: true if the initiator approximation is being used.
        ! In/Out:
        !    rng: random number generator.
        !    determ_flags: A list of flags specifying whether determinants in
        !        walker_dets are deterministic or not.

        use parallel, only: nthreads, nprocs
        use spawn_data, only: annihilate_wrapper_spawn_t
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        logical, intent(in) :: tinitiator
        integer, intent(inout), optional :: determ_flags(:)

        integer, parameter :: thread_id = 0

        call annihilate_wrapper_spawn_t(qmc_spawn, tinitiator)

        if (qmc_spawn%head(thread_id,0) > 0) then
            ! Have spawned walkers on this processor.

            if (tinitiator) then 
                call annihilate_main_list_initiator()
            else
                call annihilate_main_list()
            end if

            ! Remove determinants with zero walkers on them from the main
            ! walker list.
            call remove_unoccupied_dets(rng, determ_flags)

            ! Remove low-population spawned walkers by stochastically
            ! rounding their population up to one or down to zero.
            if (real_amplitudes) call round_low_population_spawns(rng)

            ! Insert new walkers into main walker list.
            call insert_new_walkers_wrapper(sys, determ_flags)

        else

            ! No spawned walkers so we only have to check to see if death has
            ! killed the entire population on a determinant.
            call remove_unoccupied_dets(rng, determ_flags)

        end if

    end subroutine direct_annihilation

    subroutine annihilate_main_list()

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        use basis, only: total_basis_length
        use search, only: binary_search

        integer :: i, pos, k, istart, iend, nannihilate
        integer(int_p) :: old_pop(sampling_size)
        integer(i0) :: f(total_basis_length)
        logical :: hit
        integer, parameter :: thread_id = 0

        nannihilate = 0
        istart = 1
        iend = tot_walkers
        do i = 1, qmc_spawn%head(thread_id,0)
            f = int(qmc_spawn%sdata(:total_basis_length,i), i0)
            call binary_search(walker_dets, f, istart, iend, hit, pos)
            if (hit) then
                ! Annihilate!
                old_pop = walker_population(:,pos)
                walker_population(:,pos) = walker_population(:,pos) + &
                    int(qmc_spawn%sdata(qmc_spawn%bit_str_len+1:qmc_spawn%bit_str_len+qmc_spawn%ntypes,i), int_p)
                nannihilate = nannihilate + 1
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                nparticles = nparticles + real(abs(walker_population(:,pos)) - abs(old_pop),dp)/real_factor
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                k = i - nannihilate
                qmc_spawn%sdata(:,k) = qmc_spawn%sdata(:,i)
            end if
        end do

        qmc_spawn%head(thread_id,0) = qmc_spawn%head(thread_id,0) - nannihilate

    end subroutine annihilate_main_list

    subroutine annihilate_main_list_initiator()

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        ! This version is for the initiator algorithm, where we also need to
        ! discard spawned walkers which are on previously unoccupied determinants
        ! and which are from non-initiator or non-sign-coherent events.

        use basis, only: total_basis_length
        use search, only: binary_search

        integer :: i, ipart, pos, k, istart, iend, nannihilate
        integer(int_p) :: old_pop(sampling_size)
        integer(i0) :: f(total_basis_length)
        logical :: hit, discard
        integer, parameter :: thread_id = 0

        nannihilate = 0
        istart = 1
        iend = tot_walkers
        do i = 1, qmc_spawn%head(thread_id,0)
            f = int(qmc_spawn%sdata(:total_basis_length,i), i0)
            call binary_search(walker_dets, f, istart, iend, hit, pos)
            if (hit) then
                old_pop = walker_population(:,pos)
                ! Need to take into account that the determinant might not have
                ! a non-zero population for all particle types.
                do ipart = 1, sampling_size
                    if (walker_population(ipart,pos) /= 0_int_p) then
                        ! Annihilate!
                        walker_population(ipart,pos) = walker_population(ipart,pos) + &
                                                        int(qmc_spawn%sdata(ipart+qmc_spawn%bit_str_len,i), int_p)
                    else if (.not.btest(qmc_spawn%sdata(qmc_spawn%flag_indx,i),ipart+qmc_spawn%bit_str_len)) then
                        ! Keep only if from a multiple spawning event or an
                        ! initiator.
                        ! If this is the case, then sdata(flag_indx,i) 
                        ! does not have a bit set in corresponding to 2**pop_indx, 
                        ! where pop_indx is the index of this walker type in the
                        ! qmc_spawn%sdata array (i.e. ipart+bit_str_len).
                        walker_population(ipart,pos) = int(qmc_spawn%sdata(ipart+qmc_spawn%bit_str_len,i), int_p)
                    end if
                end do
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                nparticles = nparticles + real(abs(walker_population(:,pos)) - abs(old_pop), dp)/real_factor
                ! One more entry to be removed from the qmc_spawn%sdata array.
                nannihilate = nannihilate + 1
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                ! Keep only progeny spawned by initiator determinants
                ! or multiple sign-coherent events.  If neither of these
                ! conditions are met then the j-th bit of sdata(flag_indx,i) is set,
                ! where j is the particle index in qmc_spawn%sdata.
                discard = .true.
                do ipart = qmc_spawn%bit_str_len+1, qmc_spawn%bit_str_len+qmc_spawn%ntypes
                    if (btest(qmc_spawn%sdata(qmc_spawn%flag_indx,i),ipart)) then
                        ! discard attempting spawnings from non-initiator walkers
                        ! onto unoccupied determinants.
                        ! note that the number of particles (nparticles) was not
                        ! updated at the time of spawning, so doesn't change.
                        qmc_spawn%sdata(ipart,i-nannihilate) = 0_int_s
                    else
                        ! keep!
                        qmc_spawn%sdata(ipart,i-nannihilate) = qmc_spawn%sdata(ipart,i)
                        discard = .false.
                    end if
                end do
                if (discard) then
                    ! Don't need to keep any particles from the current slot so can
                    ! just overwrite them...
                    nannihilate = nannihilate + 1
                else
                    ! Need to copy the bit string across...
                    qmc_spawn%sdata(:total_basis_length,i-nannihilate) = qmc_spawn%sdata(:total_basis_length,i)
                end if
            end if
        end do

        qmc_spawn%head(thread_id,0) = qmc_spawn%head(thread_id,0) - nannihilate

    end subroutine annihilate_main_list_initiator

    subroutine remove_unoccupied_dets(rng, determ_flags)

        ! Remove any determinants with 0 population.
        ! This can be done in a more efficient manner by doing it only when
        ! necessary...
        ! Also, before doing this, stochastically round up or down any
        ! populations which are less than one.
        ! The above steps are not performed for deterministic states which are
        ! kept in walker_dets whatever their sign.

        ! In/Out:
        !    rng: random number generator.
        !    determ_flags: A list of flags specifying whether determinants in
        !        walker_dets are deterministic or not.

        use basis, only: total_basis_length
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use qmc_common, only: stochastic_round

        type(dSFMT_t), intent(inout) :: rng
        integer, intent(inout), optional :: determ_flags(:)

        integer :: nzero, i, k, itype
        integer(int_p) :: old_pop(sampling_size)
        real(dp) :: r
        logical :: determ_det

        nzero = 0
        do i = 1, tot_walkers

            determ_det = .false.
            if (present(determ_flags)) determ_det = determ_flags(i) == 0

            ! Stochastically round the walker populations up to real_factor
            ! (which is equal to 1 in the decoded representation) or down to
            ! zero. This is not done for deterministic states.
            ! [review] - JSS: should we only call stochastic_round (as an opimisation) if real_population?
            ! [review] - JSS: (similarly throughout annihilation).  Only just noticed this...
            ! [reply] - NSB: OK, although it is a deoptimisation if reals are being used! It is
            ! [reply] - NSB: also using more global data. I haven't added an if-statement in the
            ! [reply] - NSB: other place that stochastic_round is called, round_low_population_spawns,
            ! [reply] - NSB: as this is only called for reals anyway.
            if (real_amplitudes .and. (.not. determ_det)) then
                old_pop = walker_population(:,i)
                call stochastic_round(rng, walker_population(:,i), real_factor, qmc_spawn%ntypes)
                nparticles = nparticles + real(abs(walker_population(:,i)) - abs(old_pop),dp)/real_factor
            end if

            if (all(walker_population(:,i) == 0_int_p) .and. (.not. determ_det)) then
                nzero = nzero + 1
            else if (nzero > 0) then
                k = i - nzero
                walker_dets(:,k) = walker_dets(:,i)
                walker_population(:,k) = walker_population(:,i)
                walker_data(:,k) = walker_data(:,i)
                determ_flags(k) = determ_flags(i)
            end if
        end do
        tot_walkers = tot_walkers - nzero

    end subroutine remove_unoccupied_dets

    subroutine round_low_population_spawns(rng)

        ! Loop over all spawned walkers. For each walker with a population of
        ! less than one, round it up to one with a probability equal to its
        ! population, otherwise down to zero. This will ensure that the
        ! expectation value of the amplitudes on all determinants are the same
        ! as before, but will prevent many low-weight walkers remaining in the
        ! simulation.

        ! Note that all deterministic states are always kept in the main list.
        ! The walkers in the spawned list at this point only contain states
        ! not in the main list, so cannot contain any deterministic states.
        ! Therefore, as an optimisation, we don't need to check if determinants
        ! are deterministic or not before any rounding.

        ! In/Out:
        !    rng: random number generator.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use qmc_common, only: stochastic_round

        type(dSFMT_t), intent(inout) :: rng

        integer :: i, k, itype, nremoved
        integer(int_s) :: real_factor_s
        real(dp) :: r
        integer, parameter :: thread_id = 0

        real_factor_s = int(real_factor, int_s)

        nremoved = 0
        ! [note] - It might be more efficient to combine this with insert_new_walkers_wrapper.
        ! [note] - The number of particles to insert should be small by this point though...
        do i = 1, qmc_spawn%head(thread_id,0)

            ! spawned_population holds the spawned population in its encoded
            ! form (see comments for walker_population).
            associate(spawned_population => qmc_spawn%sdata(qmc_spawn%bit_str_len+1:qmc_spawn%bit_str_len+qmc_spawn%ntypes, i))

                ! Stochastically round the walker populations up or down to
                ! real_factor (which is equal to 1 in the decoded representation).
                call stochastic_round(rng, spawned_population, real_factor_s, qmc_spawn%ntypes)

                ! If all the amplitudes for this determinant were zeroed then we
                ! don't want to add it to the main list.
                if (all(spawned_population == 0_int_s)) then
                    nremoved = nremoved + 1
                else
                    ! Shuffle this determinant down to fill in any newly opened
                    ! slots.
                    k = i - nremoved
                    qmc_spawn%sdata(:,k) = qmc_spawn%sdata(:,i)
                end if

            end associate

        end do

        qmc_spawn%head(thread_id,0) = qmc_spawn%head(thread_id,0) - nremoved

    end subroutine round_low_population_spawns

    subroutine insert_new_walkers_wrapper(sys, determ_flags)

        ! Insert new walkers into the main walker list from the spawned list.
        ! This is done after all particles have been annihilated, so the spawned
        ! list contains only new walkers.

        ! In:
        !    sys: system being studied.
        !    determ_flags: A list of flags specifying whether determinants in
        !        walker_dets are deterministic or not.

        use basis, only: total_basis_length
        use search, only: binary_search
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(inout), optional :: determ_flags(:)

        integer :: i, istart, iend, j, k, pos
        integer(int_p) :: spawned_population(sampling_size)
        real(dp) :: real_population(sampling_size)
        logical :: hit
        integer, parameter :: thread_id = 0

        ! Merge new walkers into the main list.

        ! Both the main list and the spawned list are sorted: the spawned list
        ! is sorted explicitly and the main list is sorted by construction via
        ! merging.

        ! 1. Find the position where the spawned walker should go.
        ! 2. Move all walkers above it to create the vacant slot for the new
        ! walker.  As we know how many elements we are inserting, we only need
        ! move a given walker at most once.
        ! 3. Insert the new walker at the bottom of the shifted block so it
        ! doesn't have to be moved again to accommodate other new walkers.

        ! We can make the search faster as we iterate through the spawned
        ! walkers in descending order, so once we know where one walker goes, we
        ! know that the next new walker has to go below it, allowing us to
        ! search through an ever-decreasing number of elements.

        istart = 1
        iend = tot_walkers
        do i = qmc_spawn%head(thread_id,0), 1, -1

            ! spawned det is not in the main walker list.
            call binary_search(walker_dets, int(qmc_spawn%sdata(:total_basis_length,i), i0), istart, iend, hit, pos)
            ! f should be in slot pos.  Move all determinants above it.
            do j = iend, pos, -1
                ! i is the number of determinants that will be inserted below j.
                k = j + i
                walker_dets(:,k) = walker_dets(:,j)
                walker_population(:,k) = walker_population(:,j)
                walker_data(:,k) = walker_data(:,j)
                if (present(determ_flags)) determ_flags(k) = determ_flags(j)
            end do

            ! Insert new walker into pos and shift it to accommodate the number
            ! of elements that are still to be inserted below it.
            k = pos + i - 1

            ! The encoded spawned walker sign.
            associate(spawned_population => qmc_spawn%sdata(qmc_spawn%bit_str_len+1:qmc_spawn%bit_str_len+qmc_spawn%ntypes, i))
                call insert_new_walker(sys, k, int(qmc_spawn%sdata(:total_basis_length,i), i0), int(spawned_population, int_p))
                ! Extract the real sign from the encoded sign.
                real_population = real(spawned_population,dp)/real_factor
                nparticles = nparticles + abs(real_population)
            end associate

            ! A deterministic state can never leave the main list so cannot be
            ! in the spawned list at this point. So set the flag to specify
            ! that this state is not deterministic.
            if (present(determ_flags)) determ_flags(k) = 1

            ! Next walker will be inserted below this one.
            iend = pos - 1
        end do

        ! Update tot_walkers
        tot_walkers = tot_walkers + qmc_spawn%head(thread_id,0)

    end subroutine insert_new_walkers_wrapper

    subroutine insert_new_walker(sys, pos, det, population)

        ! Insert a new determinant, det, at position pos in walker_dets. Also
        ! insert a new population at position pos in walker_population and
        ! calculate and insert all new values of walker_data for the given
        ! walker. Note that all data at position pos in these arrays will be
        ! overwritten.

        ! In:
        !    sys: system being studied.
        !    pos: The position in the walker arrays in which to insert the new
        !        data.
        !    det: The determinant to insert into walker_dets.
        !    population: The population to insert into walker_population.

        use basis, only: basis_length, total_basis_length
        use calc, only: doing_calc, hfs_fciqmc_calc, dmqmc_calc
        use calc, only: trial_function, neel_singlet
        use heisenberg_estimators, only: neel_singlet_data
        use hfs_data, only: O00
        use proc_pointers, only: sc0_ptr, op0_ptr
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: pos
        integer(i0), intent(in) :: det(total_basis_length)
        integer(int_p), intent(in) :: population(sampling_size)

        ! Insert the new determinant.
        walker_dets(:,pos) = det
        ! Insert the new population.
        walker_population(:,pos) = population
        ! Calculate and insert all new components of walker_data.
        if (.not. doing_calc(dmqmc_calc)) walker_data(1,pos) = sc0_ptr(sys, det) - H00
        if (trial_function == neel_singlet) walker_data(sampling_size+1:sampling_size+2,pos) = neel_singlet_data(det)
        if (doing_calc(hfs_fciqmc_calc)) then
            ! Set walker_data(2:,k) = <D_i|O|D_i> - <D_0|O|D_0>.
            walker_data(2,pos) = op0_ptr(sys, det) - O00
        else if (doing_calc(dmqmc_calc)) then
            ! Set the energy to be the average of the two induvidual energies.
            walker_data(1,pos) = (walker_data(1,pos) + sc0_ptr(sys, walker_dets((basis_length+1):(2*basis_length),pos)) - H00)/2
            if (replica_tricks) then
                walker_data(2:sampling_size,pos) = walker_data(1,pos)
            end if
        end if

    end subroutine insert_new_walker

end module annihilation
