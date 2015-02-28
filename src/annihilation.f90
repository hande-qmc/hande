module annihilation

use const
use fciqmc_data

implicit none

contains

    subroutine direct_annihilation(sys, rng, tinitiator, nspawn_events, determ)

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
        !    determ (optional): Derived type containing information on the
        !       semi-stochastic part of the simulation.
        ! Out:
        !    nspawn_events (optional): number of successful spawning events on
        !       the processor.

        use parallel, only: nthreads, nprocs, iproc
        use spawn_data, only: annihilate_wrapper_spawn_t, calc_events_spawn_t
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        logical, intent(in) :: tinitiator
        integer, optional, intent(out) :: nspawn_events
        type(semi_stoch_t), intent(inout), optional :: determ

        integer, parameter :: thread_id = 0

        if (present(nspawn_events)) nspawn_events = calc_events_spawn_t(qmc_spawn)

        ! If performing a semi-stochastic calculation then the annihilation
        ! process is slightly different, so call the correct routines depending
        ! on the situation.
        if (present(determ)) then
            if (determ%separate_annihilation) then
                call annihilate_wrapper_spawn_t(qmc_spawn, tinitiator)
            else
                call annihilate_wrapper_spawn_t(qmc_spawn, tinitiator, determ%sizes(iproc))
            end if

            call annihilate_main_list_wrapper(sys, rng, tinitiator, qmc_spawn, determ_flags=determ%flags)
        else
            call annihilate_wrapper_spawn_t(qmc_spawn, tinitiator)
            call annihilate_main_list_wrapper(sys, rng, tinitiator, qmc_spawn)
        end if

    end subroutine direct_annihilation

    subroutine direct_annihilation_received_list(sys, rng, tinitiator)

        ! Annihilation algorithm for non-blocking communications.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! If doing load balancing as well the annihilation procedure changes a
        ! bit. If load balancing has been decided upon then proc_map will have been
        ! updated by now. This means any walkers spawned during the current iteration
        ! will have been added to the modified section of the spawned walker array.
        ! To take care of annihilation we first merge the received list into the
        ! main list and then redistribute sections of the main list according
        ! the the updated proc map, as in normal load balancing.
        ! These walkers will be added to the spawned list and then communicated,
        ! using non-blocking comms instead of the normal MPI_AlltoAll. So they
        ! will need to be evolved upon receipt, as is the case for nomal
        ! non-blocking communications.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        ! In:
        !    sys: system being studied.
        !    tinitiator: true if the initiator approximation is being used.
        ! In/Out:
        !    rng: random number generator.

        use parallel, only: nthreads, nprocs, iproc
        use spawn_data, only: annihilate_wrapper_non_blocking_spawn, calculate_displacements, &
                              non_blocking_send
        use sort, only: qsort
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        logical, intent(in) :: tinitiator

        integer, parameter :: thread_id = 0

        ! Perform annihilation inside received list. This involves annihilating
        ! walkers which were spawned onto this processor from other processors
        ! (not including the current processor) from  the previous iteration.
        ! They have since been evolved so they can be annihilated with the main list.
        ! First annihilate within the received_list.
        call annihilate_wrapper_non_blocking_spawn(received_list, tinitiator)
        ! Annihilate with main list.
        call annihilate_main_list_wrapper(sys, rng, tinitiator, received_list)

    end subroutine direct_annihilation_received_list

    subroutine direct_annihilation_spawned_list(sys, rng, tinitiator, send_counts, req_data_s, non_block_spawn,&
                                                nspawn_events)

        ! Annihilation algorithm for non-blocking communications.
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
        !    send_counts: array of messages sizes. Will be allocated in
        !       calculate_displacements and sent in non_blocking_send.
        !    req_data_s: array of requests for non-blocking send of walkers.
        ! Out:
        !    non_block_spawn: number of spawned particles on current processor
        !       during current MC cycle.
        !    nspawn_events (optional): number of successful spawning events on
        !       the processor.

        use parallel, only: nthreads, nprocs, iproc
        use spawn_data, only: annihilate_wrapper_non_blocking_spawn, calculate_displacements, &
                              non_blocking_send
        use sort, only: qsort
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        logical, intent(in) :: tinitiator
        integer, intent(inout) :: send_counts(0:)
        integer, intent(inout) :: req_data_s(0:)
        integer, intent(out) :: non_block_spawn(:)
        integer, optional, intent(out) :: nspawn_events

        integer, parameter :: thread_id = 0

        ! Need to calculate how many walkers we are going to send to all other
        ! processors. Need to do it now as spawn%head changes meaning upon annihilation.
        call calculate_displacements(qmc_spawn, send_counts, non_block_spawn)
        if (present(nspawn_events)) nspawn_events = non_block_spawn(1)

        ! Perform annihilation within the spawned walker list.
        ! This involves locating, compressing and sorting the section of the spawned
        ! list which needs to be annihilated with the main list on this processor.
        call annihilate_wrapper_non_blocking_spawn(qmc_spawn, tinitiator, iproc)
        ! Annihilate portion of spawned list with main list.
        call annihilate_main_list_wrapper(sys, rng, tinitiator, qmc_spawn, qmc_spawn%head_start(thread_id, iproc)+nthreads)
        ! Communicate walkers spawned onto other processors during this
        ! evolution step to their new processors.
        call non_blocking_send(qmc_spawn, send_counts, req_data_s)

    end subroutine direct_annihilation_spawned_list

    subroutine annihilate_main_list_wrapper(sys, rng, tinitiator, spawn, lower_bound, determ_flags)

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process during non-blocking
        ! communications.

        ! In:
        !    sys: system being studied.
        !    tinitiator: true if the initiator approximation is being used.
        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object containing spawned particles. For non-blocking
        !       communications a subsection of the spawned walker list will be annihilated
        !       with the main list, otherwise the entire list will be annihilated and merged.
        !    determ_flags (optional): A list of flags specifying whether determinants in
        !        walker_dets are deterministic or not.
        ! In (optional):
        !     lower_bound: starting point we annihiliate from in spawn_t object.

        use system, only: sys_t
        use spawn_data, only: spawn_t
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        logical, intent(in) :: tinitiator
        integer, optional, intent(in) :: lower_bound
        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout), optional :: determ_flags(:)

        integer, parameter :: thread_id = 0
        integer :: spawn_start

        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if

        if (spawn%head(thread_id,0) >= spawn_start) then
            ! Have spawned walkers on this processor.

            if (tinitiator) then
                call annihilate_main_list_initiator(spawn, sys%basis%tensor_label_len, lower_bound)
            else
                call annihilate_main_list(spawn, sys%basis%tensor_label_len, lower_bound)
            end if

            ! Remove determinants with zero walkers on them from the main
            ! walker list.
            call remove_unoccupied_dets(rng, determ_flags)

            ! Remove low-population spawned walkers by stochastically
            ! rounding their population up to one or down to zero.
            if (real_amplitudes) call round_low_population_spawns(rng, lower_bound)

            ! Insert new walkers into main walker list.
            call insert_new_walkers(sys, spawn, determ_flags, lower_bound)

        else

            ! No spawned walkers so we only have to check to see if death has
            ! killed the entire population on a determinant.
            call remove_unoccupied_dets(rng, determ_flags)

        end if

    end subroutine annihilate_main_list_wrapper

    subroutine annihilate_main_list(spawn, tensor_label_len, lower_bound)

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        ! In:
        !    tensor_label_len: number of elements in the bit array describing the position
        !       of the particle in the space (i.e.  determinant label in vector/pair of
        !       determinants label in array).
        ! In/Out:
        !    spawn: spawn_t obeject containing spawned particles to be annihilated with main
        !       list.
        ! In (optional):
        !    lower_bound: starting point we annihiliate from in spawn_t object.
        !       Default: 1.

        use search, only: binary_search
        use spawn_data, only: spawn_t

        integer, intent(in) :: tensor_label_len
        type(spawn_t), intent(inout) :: spawn
        integer, intent(in), optional :: lower_bound

        integer :: i, pos, k, istart, iend, nannihilate, spawn_start
        integer(int_p) :: old_pop(sampling_size)
        integer(i0) :: f(tensor_label_len)

        logical :: hit
        integer, parameter :: thread_id = 0

        nannihilate = 0
        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if
        istart = 1
        iend = tot_walkers

        do i = spawn_start, spawn%head(thread_id,0)
            f = int(spawn%sdata(:tensor_label_len,i), i0)
            call binary_search(walker_dets, f, istart, iend, hit, pos)
            if (hit) then
                ! Annihilate!
                old_pop = walker_population(:,pos)
                walker_population(:,pos) = walker_population(:,pos) + &
                    int(spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes,i), int_p)
                nannihilate = nannihilate + 1
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                nparticles = nparticles + real(abs(walker_population(:,pos)) - abs(old_pop),p)/real_factor
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                k = i - nannihilate
                spawn%sdata(:,k) = spawn%sdata(:,i)
            end if
        end do

        spawn%head(thread_id,0) = spawn%head(thread_id,0) - nannihilate

    end subroutine annihilate_main_list

    subroutine annihilate_main_list_initiator(spawn, tensor_label_len, lower_bound)

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list starting from lower bound in spawn.

        ! This version is for the initiator algorithm, where we also need to
        ! discard spawned walkers which are on previously unoccupied determinants
        ! and which are from non-initiator or non-sign-coherent events.

        ! In:
        !    tensor_label_len: number of elements in the bit array describing the position
        !       of the particle in the space (i.e.  determinant label in vector/pair of
        !       determinants label in array).
        ! In/Out:
        !    spawn: spawn_t object we wish to annihilate with main list.
        ! In (Optional):
        !    lower_bound: starting point we annihiliate from in spawn_t object.

        use search, only: binary_search

        integer, intent(in) :: tensor_label_len
        type(spawn_t), intent(inout) :: spawn
        integer, intent(in), optional :: lower_bound

        integer :: i, ipart, pos, k, istart, iend, nannihilate, spawn_start
        integer(int_p) :: old_pop(sampling_size)
        integer(i0) :: f(tensor_label_len)
        logical :: hit, discard
        integer, parameter :: thread_id = 0

        nannihilate = 0
        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if
        istart = 1
        iend = tot_walkers
        do i = spawn_start, spawn%head(thread_id,0)
            f = int(spawn%sdata(:tensor_label_len,i), i0)
            call binary_search(walker_dets, f, istart, iend, hit, pos)
            if (hit) then
                old_pop = walker_population(:,pos)
                ! Need to take into account that the determinant might not have
                ! a non-zero population for all particle types.
                do ipart = 1, sampling_size
                    if (walker_population(ipart,pos) /= 0_int_p) then
                        ! Annihilate!
                        walker_population(ipart,pos) = walker_population(ipart,pos) + &
                                                        int(spawn%sdata(ipart+spawn%bit_str_len,i), int_p)
                    else if (.not.btest(spawn%sdata(spawn%flag_indx,i),ipart-1)) then
                        ! Keep only if from a multiple spawning event or an
                        ! initiator.
                        ! If this is the case, then sdata(flag_indx,i)
                        ! does not have a bit set in corresponding to 2**(ipart-1)
                        ! where ipart+bit_str_len is the index of this walker type
                        ! in the spawn%sdata array.
                        walker_population(ipart,pos) = int(spawn%sdata(ipart+spawn%bit_str_len,i), int_p)
                    end if
                end do
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                nparticles = nparticles + real(abs(walker_population(:,pos)) - abs(old_pop),p)/real_factor
                ! One more entry to be removed from the spawn%sdata array.
                nannihilate = nannihilate + 1
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                ! Keep only progeny spawned by initiator determinants
                ! or multiple sign-coherent events.  If neither of these
                ! conditions are met then the (j-1)-th bit of sdata(flag_indx,i) is set,
                ! where j is the particle index in spawn%sdata(bit_str_len+1:,:).
                discard = .true.
                do ipart = 1, spawn%ntypes
                    if (btest(spawn%sdata(spawn%flag_indx,i),ipart-1)) then
                        ! discard attempting spawnings from non-initiator walkers
                        ! onto unoccupied determinants.
                        ! note that the number of particles (nparticles) was not
                        ! updated at the time of spawning, so doesn't change.
                        spawn%sdata(spawn%bit_str_len+ipart,i-nannihilate) = 0_int_s
                    else
                        ! keep!
                        spawn%sdata(spawn%bit_str_len+ipart,i-nannihilate) = spawn%sdata(spawn%bit_str_len+ipart,i)
                        discard = .false.
                    end if
                end do
                if (discard) then
                    ! Don't need to keep any particles from the current slot so can
                    ! just overwrite them...
                    nannihilate = nannihilate + 1
                else
                    ! Need to copy the bit string across...
                    spawn%sdata(:tensor_label_len,i-nannihilate) = spawn%sdata(:tensor_label_len,i)
                end if
            end if
        end do

        spawn%head(thread_id,0) = spawn%head(thread_id,0) - nannihilate

    end subroutine annihilate_main_list_initiator

    subroutine deterministic_annihilation(sys, rng, determ)

        ! Add in the deterministic spawnings to the main list.

        ! In:
        !    sys: system being studied.
        !    determ: Derived type containing information on the semi-stochastic
        !       part of the simulation.
        ! In/Out:
        !    rng: random number generator.

        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        type(semi_stoch_t), intent(in), optional :: determ

        integer :: i, ind
        integer(int_p) :: nspawn, old_pop(sampling_size)
        real(p) :: scaled_amp, spawn_sign

        do i = 1, size(determ%vector)
            ind = determ%indices(i)

            scaled_amp = determ%vector(i)*real_factor
            spawn_sign = sign(1.0_p, scaled_amp)
            ! Stochastically round the scaled amplitude to the nearest integer
            ! in order to encode it.
            scaled_amp = abs(scaled_amp)
            nspawn = int(scaled_amp, int_p)
            scaled_amp = scaled_amp - nspawn
            if (scaled_amp > get_rand_close_open(rng)) nspawn = nspawn + 1_int_p

            ! Add in the now-encoded deterministic spawning amplitude.
            old_pop = walker_population(:,ind)
            walker_population(1,ind) = walker_population(1,ind) + spawn_sign*nspawn
            nparticles = nparticles + real(abs(walker_population(:,ind)) - abs(old_pop),p)/real_factor
        end do

    end subroutine deterministic_annihilation

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

        use dSFMT_interface, only: dSFMT_t
        use stoch_utils, only: stochastic_round

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
            if (real_amplitudes .and. (.not. determ_det)) then
                old_pop = walker_population(:,i)
                call stochastic_round(rng, walker_population(:,i), real_factor, qmc_spawn%ntypes)
                nparticles = nparticles + real(abs(walker_population(:,i)) - abs(old_pop),p)/real_factor
            end if

            if (all(walker_population(:,i) == 0_int_p) .and. (.not. determ_det)) then
                nzero = nzero + 1
            else if (nzero > 0) then
                k = i - nzero
                walker_dets(:,k) = walker_dets(:,i)
                walker_population(:,k) = walker_population(:,i)
                walker_data(:,k) = walker_data(:,i)
                if (present(determ_flags)) determ_flags(k) = determ_flags(i)
            end if
        end do
        tot_walkers = tot_walkers - nzero

    end subroutine remove_unoccupied_dets

    subroutine round_low_population_spawns(rng, lower_bound)

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
        ! In (optional):
        !    lower_bound: starting point we annihiliate from in spawn_t object.
        !       Default: 1.

        use dSFMT_interface, only: dSFMT_t
        use stoch_utils, only: stochastic_round

        type(dSFMT_t), intent(inout) :: rng
        integer, optional, intent(in) :: lower_bound

        integer :: i, k, itype, nremoved, spawn_start
        integer(int_s) :: real_factor_s
        real(dp) :: r
        integer, parameter :: thread_id = 0

        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if

        real_factor_s = int(real_factor, int_s)

        nremoved = 0
        ! [note] - It might be more efficient to combine this with insert_new_walkers.
        ! [note] - The number of particles to insert should be small by this point though...
        do i = spawn_start, qmc_spawn%head(thread_id,0)

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

    subroutine insert_new_walkers(sys, spawn, determ_flags, lower_bound)

        ! Insert new walkers into the main walker list from the spawned list.
        ! This is done after all particles have been annihilated, so the spawned
        ! list contains only new walkers.

        ! In:
        !    sys: system being studied.
        ! In/Out:
        !    spawn: spawn_t object containing list of spawned particles.
        ! In (optional):
        !    determ_flags: A list of flags specifying whether determinants in
        !        walker_dets are deterministic or not.
        !    lower_bound: starting point we annihiliate from in spawn_t object.
 
        use search, only: binary_search
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout), optional :: determ_flags(:)
        integer, intent(in), optional :: lower_bound

        integer :: i, istart, iend, j, k, pos, spawn_start, disp
        integer(int_p) :: spawned_population(sampling_size)
        real(p) :: real_population(sampling_size)

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

        if (present(lower_bound)) then
            spawn_start = lower_bound
            disp = lower_bound - 1
        else
            spawn_start = 1
            disp = 0
        end if
        istart = 1
        iend = tot_walkers
        do i = spawn%head(thread_id,0), spawn_start, -1

            ! spawned det is not in the main walker list.
            call binary_search(walker_dets, int(spawn%sdata(:sys%basis%tensor_label_len,i), i0), istart, iend, hit, pos)
            ! f should be in slot pos.  Move all determinants above it.
            do j = iend, pos, -1
                ! i is the number of determinants that will be inserted below j.
                k = j + i - disp
                walker_dets(:,k) = walker_dets(:,j)
                walker_population(:,k) = walker_population(:,j)
                walker_data(:,k) = walker_data(:,j)
                if (present(determ_flags)) determ_flags(k) = determ_flags(j)
            end do

            ! Insert new walker into pos and shift it to accommodate the number
            ! of elements that are still to be inserted below it.
            k = pos + i - 1 - disp

            ! The encoded spawned walker sign.
            associate(spawned_population => spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes, i), &
                    tbl=>sys%basis%tensor_label_len)
                call insert_new_walker(sys, k, int(spawn%sdata(:tbl,i), i0), int(spawned_population, int_p))
                ! Extract the real sign from the encoded sign.
                real_population = real(spawned_population,p)/real_factor
                nparticles = nparticles + abs(real_population)
            end associate

            ! A deterministic state can never leave the main list so cannot be
            ! in the spawned list at this point. So set the flag to specify
            ! that this state is not deterministic.
            if (present(determ_flags)) determ_flags(k) = 1

            ! Next walker will be inserted below this one.
            iend = pos - 1
        end do

        ! Update tot_walkers.
        tot_walkers = tot_walkers + spawn%head(thread_id,0) - disp

    end subroutine insert_new_walkers

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

        use calc, only: doing_calc, hfs_fciqmc_calc, dmqmc_calc
        use calc, only: trial_function, neel_singlet, propagate_to_beta
        use heisenberg_estimators, only: neel_singlet_data
        use hfs_data, only: O00
        use proc_pointers, only: sc0_ptr, op0_ptr, trial_dm_ptr
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: pos
        integer(i0), intent(in) :: det(sys%basis%tensor_label_len)
        integer(int_p), intent(in) :: population(sampling_size)

        ! Insert the new determinant.
        walker_dets(:,pos) = det
        ! Insert the new population.
        walker_population(:,pos) = population
        ! Calculate and insert all new components of walker_data.
        walker_data(1,pos) = sc0_ptr(sys, det) - H00
        if (trial_function == neel_singlet) walker_data(sampling_size+1:sampling_size+2,pos) = neel_singlet_data(sys, det)
        if (doing_calc(hfs_fciqmc_calc)) then
            ! Set walker_data(2:,k) = <D_i|O|D_i> - <D_0|O|D_0>.
            walker_data(2,pos) = op0_ptr(sys, det) - O00
        else if (doing_calc(dmqmc_calc)) then
            if (propagate_to_beta) then
                ! Store H^T_ii-H_jj so we can propagate with ~ 1 + \Delta\beta(H^T_ii - H_jj),
                ! where H^T_ii is the "trial" Hamiltonian.
                associate(bl=>sys%basis%string_len)
                    walker_data(1,pos) = -trial_dm_ptr(sys, walker_dets((bl+1):(2*bl),pos)) + (walker_data(1,pos) + H00)
                end associate
            else
                ! Set the energy to be the average of the two induvidual energies.
                associate(bl=>sys%basis%string_len)
                    walker_data(1,pos) = (walker_data(1,pos) + sc0_ptr(sys, walker_dets((bl+1):(2*bl),pos)) - H00)/2
                end associate
                if (replica_tricks) then
                    walker_data(2:sampling_size,pos) = walker_data(1,pos)
                end if
            end if
        end if

    end subroutine insert_new_walker

end module annihilation
