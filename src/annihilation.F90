module annihilation

use const
use fciqmc_data

implicit none

contains

    subroutine direct_annihilation(sc0)

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

        ! 1. Sort spawned walkers list.
        call sort_spawned_lists()

        ! 2. Annihilate within spawned walkers list.
        ! Compress the remaining spawned walkers list.
        call annihilate_spawned_list()

        ! 4. Annilate main list.
        call annihilate_main_list()

        ! 5. Insert new walkers into main walker list.
        call insert_new_walkers(sc0)

    end subroutine direct_annihilation

    subroutine distribute_walkers()

        use parallel

        use basis, only: basis_length

        integer :: send_counts(0:nprocs-1), send_displacements(0:nprocs-1)
        integer :: receive_counts(0:nprocs-1), receive_displacements(0:nprocs-1)
        integer :: s(2,0:nprocs-1)
        integer :: r(2,0:nprocs-1)
        integer :: i, step, ierr
        integer(i0), pointer :: tmp_dets(:,:)
        integer, pointer :: tmp_population(:)

        if (nprocs == 1) then
            ! No need to communicate!
        else
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
#ifdef PARALLEL
            call MPI_AlltoAll(s, 1, MPI_INTEGER, r, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
#endif
            send_counts = s(1,:)
            receive_counts = r(1,:)
            D0_population = r(2, D0_proc)

            ! Want spawning data to be continuous after move, hence need to find the
            ! receive displacements.
            receive_displacements(0) = 0
            do i=1, nprocs-1
                receive_displacements(i) = receive_displacements(i-1) + receive_counts(i-1)
            end do

#ifdef PARALLEL
            ! Send spawning populations.
            call MPI_AlltoAllv(spawned_walker_population, send_counts, send_displacements, MPI_INTEGER, &
                               spawned_walker_population_recvd, receive_counts, receive_displacements, MPI_INTEGER, &
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
                               spawned_walker_dets, receive_counts, receive_displacements, mpi_det_integer, &
                               MPI_COMM_WORLD, ierr)
#endif

            ! Swap pointers so that spawned_walker_dets and
            ! spawned_walker_population point to the received data.
            tmp_dets => spawned_walker_dets
            spawned_walker_dets => spawned_walker_dets_recvd
            spawned_walker_dets_recvd => tmp_dets
            tmp_population => spawned_walker_population
            spawned_walker_population => spawned_walker_population_recvd
            spawned_walker_population_recvd => tmp_population

            ! Set spawning_head(0) to be the number of walkers now on this
            ! processor.
            spawning_head(0) = receive_displacements(nprocs-1) + receive_counts(nprocs-1)
        end if

    end subroutine distribute_walkers

    subroutine annihilate_spawned_list()

        ! Annihilate the spawned walker list and compress the remaining
        ! elements.

        ! The spawned walker list is already sorted, so annihilation amounts to
        ! looping through the list and adding consective walker populations
        ! together if they're the same walker.

        integer :: islot, k, nremoved

        ! islot is the current element in the spawned walkers lists.
        islot = 1
        ! k is the current element which is being compressed into islot (if
        ! k and islot refer to the same determinants).
        k = 1
        nremoved = 0
        self_annihilate: do
            ! Set the current free slot to be the next unique spawned walker.
            spawned_walker_dets(:,islot) = spawned_walker_dets(:,k) 
            spawned_walker_population(islot) = spawned_walker_population(k) 
            compress: do
                k = k + 1
                if (k > spawning_head(0)) exit self_annihilate
                if (all(spawned_walker_dets(:,k) == spawned_walker_dets(:,islot))) then
                    ! Add the populations of the subsequent identical walkers.
                    spawned_walker_population(islot) = spawned_walker_population(islot) + spawned_walker_population(k)
                    nremoved = nremoved + 1
                else
                    ! Found the next unique spawned walker.
                    exit compress
                end if
            end do compress
            ! go to the next slot.
            islot = islot + 1
            if (islot > spawning_head(0)) exit self_annihilate
        end do self_annihilate

        ! update spawning_head(0)
        spawning_head(0) = spawning_head(0) - nremoved

    end subroutine annihilate_spawned_list

    subroutine annihilate_main_list()

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        use basis, only: basis_length

        integer :: i, pos, k, nannihilate, nzero, istart, iend
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
                walker_population(pos) = walker_population(pos) + spawned_walker_population(i)
                nannihilate = nannihilate + 1
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                k = i - nannihilate
                spawned_walker_dets(:,k) = spawned_walker_dets(:,i)
                spawned_walker_population(k) = spawned_walker_population(i)
            end if
        end do

        spawning_head(0) = spawning_head(0) - nannihilate

        ! Remove any determinants with 0 population.
        ! This can be done in a more efficient manner by doing it only when necessary...
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

    end subroutine annihilate_main_list

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
            walker_population(k) = spawned_walker_population(i)
            walker_energies(k) = sc0(walker_dets(:,k)) - H00
            ! Next walker will be inserted below this one.
            iend = pos - 1
        end do
        
        ! Update tot_walkers
        tot_walkers = tot_walkers + spawning_head(0)

    end subroutine insert_new_walkers

end module annihilation
