module annihilation

use const
use fciqmc_data

implicit none

real(dp) :: annihilation_comms_time = 0.0_dp

contains

    subroutine direct_annihilation(sys, tinitiator)

        ! Annihilation algorithm.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        ! In:
        !    sys: system being studied.
        !    tinitiator: true if the initiator approximation is being used.

        use parallel, only: nthreads, nprocs
        use spawn_data, only: annihilate_wrapper_spawn_t
        use sort, only: qsort
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        logical, intent(in) :: tinitiator

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
            call remove_unoccupied_dets()

            ! Insert new walkers into main walker list.
            call insert_new_walkers(sys)

        else

            ! No spawned walkers so we only have to check to see if death has
            ! killed the entire population on a determinant.
            call remove_unoccupied_dets()

        end if

    end subroutine direct_annihilation

    subroutine annihilate_main_list()

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        use basis, only: total_basis_length
        use search, only: binary_search

        integer :: i, pos, k, nannihilate, istart, iend, old_pop(sampling_size)
        integer(i0) :: f(total_basis_length)
        logical :: hit
        integer, parameter :: thread_id = 0

        nannihilate = 0
        istart = 1
        iend = tot_walkers
        do i = 1, qmc_spawn%head(thread_id,0)
            f = qmc_spawn%sdata(:total_basis_length,i)
            call binary_search(walker_dets, f, istart, iend, hit, pos)
            if (hit) then
                ! Annihilate!
                old_pop = walker_population(:,pos)
                walker_population(:,pos) = walker_population(:,pos) + &
                    qmc_spawn%sdata(qmc_spawn%bit_str_len+1:qmc_spawn%bit_str_len+qmc_spawn%ntypes,i)
                nannihilate = nannihilate + 1
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                nparticles = nparticles + abs(walker_population(:,pos)) - abs(old_pop(:))
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

        integer :: i, ipart, pos, k, nannihilate, istart, iend, old_pop(sampling_size)
        integer(i0) :: f(total_basis_length)
        logical :: hit, discard
        integer, parameter :: thread_id = 0

        nannihilate = 0
        istart = 1
        iend = tot_walkers
        do i = 1, qmc_spawn%head(thread_id,0)
            f = qmc_spawn%sdata(:total_basis_length,i)
            call binary_search(walker_dets, f, istart, iend, hit, pos)
            if (hit) then
                old_pop = walker_population(:,pos)
                ! Need to take into account that the determinant might not have
                ! a non-zero population for all particle types.
                do ipart = 1, sampling_size
                    if (walker_population(ipart,pos) /= 0) then
                        ! Annihilate!
                        walker_population(ipart,pos) = walker_population(ipart,pos) + &
                                                        qmc_spawn%sdata(ipart+qmc_spawn%bit_str_len,i)
                    else if (.not.btest(qmc_spawn%sdata(qmc_spawn%flag_indx,i),ipart+qmc_spawn%bit_str_len)) then
                        ! Keep only if from a multiple spawning event or an
                        ! initiator.
                        ! If this is the case, then sdata(flag_indx,i) 
                        ! does not have a bit set in corresponding to 2**pop_indx, 
                        ! where pop_indx is the index of this walker type in the
                        ! qmc_spawn%sdata array (i.e. ipart+bit_str_len).
                        walker_population(ipart,pos) = qmc_spawn%sdata(ipart+qmc_spawn%bit_str_len,i)
                    end if
                end do
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                nparticles = nparticles + abs(walker_population(:,pos)) - abs(old_pop)
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
                        qmc_spawn%sdata(ipart,i-nannihilate) = 0
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

    subroutine remove_unoccupied_dets()

        ! Remove any determinants with 0 population.
        ! This can be done in a more efficient manner by doing it only when necessary...

        use basis, only: total_basis_length

        integer :: nzero, i, k

        nzero = 0
        do i = 1, tot_walkers
            if (all(walker_population(:,i) == 0)) then
                nzero = nzero + 1
            else if (nzero > 0) then
                k = i - nzero
                walker_dets(:,k) = walker_dets(:,i)
                walker_population(:,k) = walker_population(:,i)
                walker_data(:,k) = walker_data(:,i)
            end if
        end do
        tot_walkers = tot_walkers - nzero

    end subroutine remove_unoccupied_dets

    subroutine insert_new_walkers(sys)

        ! Insert new walkers into the main walker list from the spawned list.
        ! This is done after all particles have been annihilated, so the spawned
        ! list contains only new walkers.

        ! In:
        !    sys: system being studied.

        use basis, only: basis_length, total_basis_length
        use calc, only: doing_calc, hfs_fciqmc_calc, dmqmc_calc
        use determinants, only: decode_det
        use search, only: binary_search
        use system, only: sys_t
        use calc, only: trial_function, neel_singlet
        use hfs_data, only: O00
        use proc_pointers, only: sc0_ptr, op0_ptr
        use heisenberg_estimators, only: neel_singlet_data

        type(sys_t), intent(in) :: sys

        integer :: i, istart, iend, j, k, pos
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
            ! spawned det is not in the main walker list
            call binary_search(walker_dets, qmc_spawn%sdata(:total_basis_length,i), istart, iend, hit, pos)
            ! f should be in slot pos.  Move all determinants above it.
            do j = iend, pos, -1
                ! i is the number of determinants that will be inserted below j.
                k = j + i
                walker_dets(:,k) = walker_dets(:,j)
                walker_population(:,k) = walker_population(:,j)
                walker_data(:,k) = walker_data(:,j)
            end do
            ! Insert new walker into pos and shift it to accommodate the number
            ! of elements that are still to be inserted below it.
            k = pos + i - 1
            walker_dets(:,k) = qmc_spawn%sdata(:total_basis_length,i)
            walker_population(:,k) = qmc_spawn%sdata(qmc_spawn%bit_str_len+1:qmc_spawn%bit_str_len+qmc_spawn%ntypes,i)
            nparticles = nparticles + abs(qmc_spawn%sdata(qmc_spawn%bit_str_len+1:qmc_spawn%bit_str_len+qmc_spawn%ntypes,i))
            walker_data(1,k) = sc0_ptr(sys, walker_dets(:,k)) - H00
            if (trial_function == neel_singlet) walker_data(sampling_size+1:sampling_size+2,k) = neel_singlet_data(walker_dets(:,k))
            if (doing_calc(hfs_fciqmc_calc)) then
                ! Set walker_data(2:,k) = <D_i|O|D_i> - <D_0|O|D_0>.
                walker_data(2,k) = op0_ptr(sys, walker_dets(:,k)) - O00
            else if (doing_calc(dmqmc_calc)) then
                ! Set the energy to be the average of the two induvidual energies.
                walker_data(1,k) = (walker_data(1,k) + sc0_ptr(sys, walker_dets((basis_length+1):(2*basis_length),k)) - H00)/2
                if (replica_tricks) then
                    walker_data(2:sampling_size,k) = walker_data(1,k)
                end if
            end if
            ! Next walker will be inserted below this one.
            iend = pos - 1
        end do

        ! Update tot_walkers
        tot_walkers = tot_walkers + qmc_spawn%head(thread_id,0)

    end subroutine insert_new_walkers

end module annihilation
