module annihilation

use const
use fciqmc_data

implicit none

contains

    subroutine direct_annihilation(sc0)

        ! Annihilation algorithm.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        interface
            function sc0(f) result(hmatel)
                use basis, only: basis_length
                use const, only: i0, p
                implicit none
                real(p) :: hmatel
                integer(i0), intent(in) :: f(basis_length)
            end function sc0
        end interface

#ifdef PARALLEL
        ! 0. Send spawned walkers to the processor which "owns" them and receive
        ! the walkers "owned" by this processor.
        call distribute_walkers()
#endif

        if (spawning_head(0) > 0) then
            ! Have spawned walkers on this processor.

            ! 1. Sort spawned walkers list.
            call sort_spawned_lists()

            ! 2. Annihilate within spawned walkers list.
            ! Compress the remaining spawned walkers list.
            call annihilate_spawned_list()

            ! 3. Annihilate main list.
            call annihilate_main_list()

            ! 4. Remove determinants with zero walkers on them from the main
            ! walker list.
            call remove_unoccupied_dets()

            ! 5. Insert new walkers into main walker list.
            call insert_new_walkers(sc0)

        else

            ! No spawned walkers so we only have to check to see if death has
            ! killed the entire population on a determinant.
            call remove_unoccupied_dets()

        end if

    end subroutine direct_annihilation

    subroutine direct_annihilation_initiator(sc0)

        ! Annihilation algorithm.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        ! This version is for use with the initiator-FCIQMC algorithm.

        interface
            function sc0(f) result(hmatel)
                use basis, only: basis_length
                use const, only: i0, p
                implicit none
                real(p) :: hmatel
                integer(i0), intent(in) :: f(basis_length)
            end function sc0
        end interface

#ifdef PARALLEL
        ! 0. Send spawned walkers to the processor which "owns" them and receive
        ! the walkers "owned" by this processor.
        call distribute_walkers()
#endif

        if (spawning_head(0) > 0) then
            ! Have spawned walkers on this processor.

            ! 1. Sort spawned walkers list.
            call sort_spawned_lists()

            ! 2. Annihilate within spawned walkers list.
            ! Compress the remaining spawned walkers list and update the parent
            ! flag.
            call annihilate_spawned_list_initiator()

            ! 3. Annihilate main list.
            ! This also removes spawned walkers that don't come from initiators
            ! or sign-coherent events and are on unoccupied determinants.
            call annihilate_main_list_initiator()

            ! 4. Remove determinants with zero walkers on them from the main
            ! walker list.
            call remove_unoccupied_dets()

            ! 5. Insert new walkers into main walker list.
            call insert_new_walkers(sc0)

        else

            ! No spawned walkers so we only have to check to see if death has
            ! killed the entire population on a determinant.
            call remove_unoccupied_dets()

        end if

    end subroutine direct_annihilation_initiator

    subroutine distribute_walkers()

        ! Send spawned walkers to the pre-designated processor which hosts the
        ! determinant upon which the walker has been spawned.

#ifdef PARALLEL

        use parallel

        use basis, only: basis_length

        integer :: send_counts(0:nprocs-1), send_displacements(0:nprocs-1)
        integer :: receive_counts(0:nprocs-1), receive_displacements(0:nprocs-1)
        integer :: send_counts_info(0:nprocs-1), send_displacements_info(0:nprocs-1)
        integer :: receive_counts_info(0:nprocs-1), receive_displacements_info(0:nprocs-1)
        integer :: s(2,0:nprocs-1)
        integer :: r(2,0:nprocs-1)
        integer :: i, step, ierr
        integer(i0), pointer :: tmp_dets(:,:)
        integer, pointer :: tmp_info(:,:)

        ! Send spawned walkers to the processor which "owns" them and receive
        ! the walkers "owned" by this processor.

        ! The walkers are already stored in the spawned walker arrays in blocks,
        ! where each block corresponds to determinants owned by a given
        ! processor.

        ! Tests on cx2 indicate that there is not much difference between
        ! sending messages of the same size using MPI_AlltoAll and
        ! MPI_AlltoAllv (though MPI_AlltoAllv is very slightly slower, by a few
        ! percent).  Therefore it is likely that using MPI_AlltoAllv will be
        ! more efficient as it allows us to only send spawned walkers rather
        ! than the entire spawned lists.  It does require an additional
        ! communication to set up however, so for calculations with large
        ! numbers of walkers maybe MPI_AlltoAll would be more efficient?

        ! Find out how many walkers we are going to send and receive.
        step = spawning_block_start(1)
        forall (i=0:nprocs-1)
            s(1,i) = spawning_head(i) - spawning_block_start(i)
            s(2,i) = D0_population
            send_displacements(i) = i*step
        end forall

        call MPI_AlltoAll(s, 2, MPI_INTEGER, r, 2, MPI_INTEGER, MPI_COMM_WORLD, ierr)

        send_counts = s(1,:)
        receive_counts = r(1,:)
        D0_population = r(2, D0_proc)

        ! Want spawning data to be continuous after move, hence need to find the
        ! receive displacements.
        receive_displacements(0) = 0
        do i=1, nprocs-1
            receive_displacements(i) = receive_displacements(i-1) + receive_counts(i-1)
        end do

        ! Send spawning populations.
        send_counts_info = send_counts*spawned_info_size
        receive_counts_info = receive_counts*spawned_info_size
        send_displacements_info = send_displacements*spawned_info_size
        receive_displacements_info = receive_displacements*spawned_info_size
        call MPI_AlltoAllv(spawned_walker_info, send_counts_info, send_displacements_info, MPI_INTEGER, &
                           spawned_walker_info_recvd, receive_counts_info, receive_displacements_info, MPI_INTEGER, &
                           MPI_COMM_WORLD, ierr)
        ! Send spawning determinants.
        ! Each element contains basis_length integers (of type
        ! i0/mpi_det_integer) so we need to change the counts and
        ! displacements accordingly:
        send_counts = send_counts*basis_length
        receive_counts = receive_counts*basis_length
        send_displacements = send_displacements*basis_length
        receive_displacements = receive_displacements*basis_length
        call MPI_AlltoAllv(spawned_walker_dets, send_counts, send_displacements, mpi_det_integer, &
                           spawned_walker_dets_recvd, receive_counts, receive_displacements, mpi_det_integer, &
                           MPI_COMM_WORLD, ierr)

        ! Swap pointers so that spawned_walker_dets and
        ! spawned_walker_info point to the received data.
        tmp_dets => spawned_walker_dets
        spawned_walker_dets => spawned_walker_dets_recvd
        spawned_walker_dets_recvd => tmp_dets
        tmp_info => spawned_walker_info
        spawned_walker_info => spawned_walker_info_recvd
        spawned_walker_info_recvd => tmp_info

        ! Set spawning_head(0) to be the number of walkers now on this
        ! processor.
        spawning_head(0) = receive_displacements(nprocs-1) + receive_counts(nprocs-1)

#endif

    end subroutine distribute_walkers

    subroutine annihilate_spawned_list()

        ! Annihilate the spawned walker list and compress the remaining
        ! elements.

        ! The spawned walker list is already sorted, so annihilation amounts to
        ! looping through the list and adding consective walker populations
        ! together if they're the same walker.

        integer :: islot, k

        ! islot is the current element in the spawned walkers lists.
        islot = 1
        ! k is the current element which is being compressed into islot (if
        ! k and islot refer to the same determinants).
        k = 1
        self_annihilate: do
            ! Set the current free slot to be the next unique spawned walker.
            spawned_walker_dets(:,islot) = spawned_walker_dets(:,k) 
            spawned_walker_info(1,islot) = spawned_walker_info(1,k) 
            compress: do
                k = k + 1
                if (k > spawning_head(0)) exit self_annihilate
                if (all(spawned_walker_dets(:,k) == spawned_walker_dets(:,islot))) then
                    ! Add the populations of the subsequent identical walkers.
                    spawned_walker_info(1,islot) = spawned_walker_info(1,islot) + spawned_walker_info(1,k)
                else
                    ! Found the next unique spawned walker.
                    exit compress
                end if
            end do compress
            ! All done?
            if (islot == spawning_head(0)) exit self_annihilate
            ! go to the next slot if the current determinant wasn't completed
            ! annihilated.
            if (spawned_walker_info(1,islot) /= 0) islot = islot + 1
        end do self_annihilate

        ! We didn't check if the population on the last determinant is
        ! completely annihilated or not.
        if (spawned_walker_info(1, islot) == 0) islot = islot - 1

        ! update spawning_head(0)
        spawning_head(0) = islot

    end subroutine annihilate_spawned_list

    subroutine annihilate_spawned_list_initiator()

        ! Annihilate the spawned walker list and compress the remaining
        ! elements.

        ! The spawned walker list is already sorted, so annihilation amounts to
        ! looping through the list and adding consective walker populations
        ! together if they're the same walker.

        ! This version is for the initiator algorithm, whereby we also need to
        ! take care of the parent flag (ie handle the origin of the spawned
        ! walkers).

        integer :: islot, k, pop_sign

        ! islot is the current element in the spawned walkers lists.
        islot = 1
        ! k is the current element which is being compressed into islot (if
        ! k and islot refer to the same determinants).
        k = 1
        self_annihilate: do
            ! Set the current free slot to be the next unique spawned walker.
            spawned_walker_dets(:,islot) = spawned_walker_dets(:,k) 
            spawned_walker_info(:,islot) = spawned_walker_info(:,k) 
            compress: do
                k = k + 1
                if (k > spawning_head(0)) exit self_annihilate
                if (all(spawned_walker_dets(:,k) == spawned_walker_dets(:,islot))) then
                    ! Update the parent flag.
                    ! Note we ignore the possibility of multiple spawning events
                    ! onto the same unoccupied determinant.  Such events become
                    ! vanishingly unlikely with the size of the determinant
                    ! space (according to Ali at least!).
                    pop_sign = spawned_walker_info(1,islot)*spawned_walker_info(1,k)
                    if (pop_sign > 0) then
                        ! Sign coherent event.
                        ! Set parent_flag to 2 (indicating multiple
                        ! sign-coherent spawning events).
                        spawned_walker_info(2,islot) = 2
                    else
                        ! Keep the parent_flag of the largest spawning event.
                        if (spawned_walker_info(1,k) > spawned_walker_info(1,islot)) then
                            spawned_walker_info(2,islot) = spawned_walker_info(2,k)
                        end if
                    end if
                    ! Add the populations of the subsequent identical walkers.
                    spawned_walker_info(1,islot) = spawned_walker_info(1,islot) + spawned_walker_info(1,k)
                else
                    ! Found the next unique spawned walker.
                    exit compress
                end if
            end do compress
            ! All done?
            if (islot == spawning_head(0)) exit self_annihilate
            ! go to the next slot if the current determinant wasn't completed
            ! annihilated.
            if (spawned_walker_info(1,islot) /= 0) islot = islot + 1
        end do self_annihilate

        ! We didn't check if the population on the last determinant is
        ! completely annihilated or not.
        if (spawned_walker_info(1, islot) == 0) islot = islot - 1

        ! update spawning_head(0)
        spawning_head(0) = islot

    end subroutine annihilate_spawned_list_initiator

    subroutine annihilate_main_list()

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        use basis, only: basis_length

        integer :: i, pos, k, nannihilate, istart, iend, old_pop
        integer(i0) :: f(basis_length)
        logical :: hit

        nannihilate = 0
        istart = 1
        iend = tot_walkers
        do i = 1, spawning_head(0)
            f = spawned_walker_dets(:,i)
            call search_walker_list(f, istart, iend, hit, pos)
            if (hit) then
                ! Annihilate!
                old_pop = walker_population(pos)
                walker_population(pos) = walker_population(pos) + spawned_walker_info(1,i)
                nannihilate = nannihilate + 1
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                nparticles = nparticles + abs(walker_population(pos)) - abs(old_pop)
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                k = i - nannihilate
                spawned_walker_dets(:,k) = spawned_walker_dets(:,i)
                spawned_walker_info(:,k) = spawned_walker_info(:,i)
            end if
        end do

        spawning_head(0) = spawning_head(0) - nannihilate

    end subroutine annihilate_main_list

    subroutine annihilate_main_list_initiator()

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        ! This version is for the initiator algorithm, where we also need to
        ! discard spawned walkers which are on previously unoccupied determinants
        ! and which are from non-initiator or non-sign-coherent events.

        use basis, only: basis_length

        integer :: i, pos, k, nannihilate, istart, iend, old_pop
        integer(i0) :: f(basis_length)
        logical :: hit

        nannihilate = 0
        istart = 1
        iend = tot_walkers
        do i = 1, spawning_head(0)
            f = spawned_walker_dets(:,i)
            call search_walker_list(f, istart, iend, hit, pos)
            if (hit) then
                ! Annihilate!
                old_pop = walker_population(pos)
                walker_population(pos) = walker_population(pos) + spawned_walker_info(1,i)
                nannihilate = nannihilate + 1
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                nparticles = nparticles + abs(walker_population(pos)) - abs(old_pop)
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                ! Keep only progeny spawned by initiator determinants
                ! (parent_flag=0) or multiple sign-coherent events
                ! (parent_flag=2).
                if (spawned_walker_info(2,i) == 1) then
                    ! discard attempting spawnings from non-initiator walkers
                    ! onto unoccupied determinants.
                    nannihilate = nannihilate + 1
                    nparticles = nparticles - abs(spawned_walker_info(1,i))
                else
                    ! keep!
                    k = i - nannihilate
                    spawned_walker_dets(:,k) = spawned_walker_dets(:,i)
                    spawned_walker_info(:,k) = spawned_walker_info(:,i)
                end if
            end if
        end do

        spawning_head(0) = spawning_head(0) - nannihilate

    end subroutine annihilate_main_list_initiator

    subroutine remove_unoccupied_dets()

        ! Remove any determinants with 0 population.
        ! This can be done in a more efficient manner by doing it only when necessary...

        integer :: nzero, i, k

        nzero = 0
        do i = 1, tot_walkers
            if (walker_population(i) == 0) then
                nzero = nzero + 1
            else if (nzero > 0) then
                k = i - nzero
                walker_dets(:,k) = walker_dets(:,i)
                walker_population(k) = walker_population(i)
                walker_energies(k) = walker_energies(i)
            end if
        end do
        tot_walkers = tot_walkers - nzero

    end subroutine remove_unoccupied_dets

    subroutine insert_new_walkers(sc0)

        ! Insert new walkers into the main walker list from the spawned list.
        ! This is done after all particles have been annihilated, so the spawned
        ! list contains only new walkers.

        use basis, only: basis_length
        use determinants, only: decode_det
        use system, only: nel
        use hamiltonian, only: slater_condon0_hub_real

        interface
            function sc0(f) result(hmatel)
                use basis, only: basis_length
                use const, only: i0, p
                implicit none
                real(p) :: hmatel
                integer(i0), intent(in) :: f(basis_length)
            end function sc0
        end interface

        integer :: i, istart, iend, j, k, pos
        logical :: hit

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
        do i = spawning_head(0), 1, -1
            ! spawned det is not in the main walker list
            call search_walker_list(spawned_walker_dets(:,i), istart, iend, hit, pos)
            ! f should be in slot pos.  Move all determinants above it.
            do j = iend, pos, -1
                ! i is the number of determinants that will be inserted below j.
                k = j + i
                walker_dets(:,k) = walker_dets(:,j)
                walker_population(k) = walker_population(j)
                walker_energies(k) = walker_energies(j)
            end do
            ! Insert new walker into pos and shift it to accommodate the number
            ! of elements that are still to be inserted below it.
            k = pos + i - 1
            walker_dets(:,k) = spawned_walker_dets(:,i)
            walker_population(k) = spawned_walker_info(1,i)
            nparticles = nparticles + abs(spawned_walker_info(1,i))
            walker_energies(k) = sc0(walker_dets(:,k)) - H00
            ! Next walker will be inserted below this one.
            iend = pos - 1
        end do
        
        ! Update tot_walkers
        tot_walkers = tot_walkers + spawning_head(0)

    end subroutine insert_new_walkers

end module annihilation
