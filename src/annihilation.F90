module annihilation

use const
use fciqmc_data

implicit none

contains

    subroutine direct_annihilation()

        ! Annihilation algorithm.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

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
            call insert_new_walkers()

        else

            ! No spawned walkers so we only have to check to see if death has
            ! killed the entire population on a determinant.
            call remove_unoccupied_dets()

        end if

    end subroutine direct_annihilation

    subroutine direct_annihilation_initiator()

        ! Annihilation algorithm.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        ! This version is for use with the initiator-FCIQMC algorithm.

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
            call insert_new_walkers()

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
        integer :: i, ierr
        integer(i0), pointer :: tmp_walkers(:,:)

        ! Send spawned walkers to the processor which "owns" them and receive
        ! the walkers "owned" by this processor.

        ! The walkers are already stored in the spawned walkers array in blocks,
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
        forall (i=0:nprocs-1)
            send_counts(i) = spawning_head(i) - spawning_block_start(i)
        end forall

        call MPI_AlltoAll(send_counts, 1, MPI_INTEGER, receive_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

        ! Want spawning data to be continuous after move, hence need to find the
        ! receive displacements.
        receive_displacements(0) = 0
        do i=1, nprocs-1
            receive_displacements(i) = receive_displacements(i-1) + receive_counts(i-1)
        end do

        ! Set spawning_head(0) to be the number of walkers now on this
        ! processor.
        ! This is given by the number of items sent by the last processor plus
        ! the displacement used by the last processor (which is the number of
        ! items sent by all other processors).
        spawning_head(0) = receive_displacements(nprocs-1) + receive_counts(nprocs-1)

        ! Send spawned_walkers.
        ! Each element contains spawned_size integers (of type
        ! i0/mpi_det_integer) so we need to change the counts and
        ! displacements accordingly:
        send_counts = send_counts*spawned_size
        receive_counts = receive_counts*spawned_size
        send_displacements = spawning_block_start(:nprocs-1)*spawned_size
        receive_displacements = receive_displacements*spawned_size
        call MPI_AlltoAllv(spawned_walkers, send_counts, send_displacements, mpi_det_integer, &
                           spawned_walkers_recvd, receive_counts, receive_displacements, mpi_det_integer, &
                           MPI_COMM_WORLD, ierr)

        ! Swap pointers so that spawned_walkers points to the received data.
        tmp_walkers => spawned_walkers
        spawned_walkers => spawned_walkers_recvd
        spawned_walkers_recvd => tmp_walkers

#endif

    end subroutine distribute_walkers

    subroutine annihilate_spawned_list()

        ! Annihilate the spawned walker list and compress the remaining
        ! elements.

        ! The spawned walker list is already sorted, so annihilation amounts to
        ! looping through the list and adding consective walker populations
        ! together if they're the same walker.

        use basis, only: basis_length

        integer :: islot, k

        ! islot is the current element in the spawned walkers lists.
        islot = 1
        ! k is the current element which is being compressed into islot (if
        ! k and islot refer to the same determinants).
        k = 1
        self_annihilate: do
            ! Set the current free slot to be the next unique spawned walker.
            spawned_walkers(:,islot) = spawned_walkers(:,k) 
            compress: do
                k = k + 1
                if (k > spawning_head(0)) exit self_annihilate
                if (all(spawned_walkers(:basis_length,k) == spawned_walkers(:basis_length,islot))) then
                    ! Add the populations of the subsequent identical walkers.
                    spawned_walkers(spawned_pop:spawned_hf_pop,islot) =    &
                         spawned_walkers(spawned_pop:spawned_hf_pop,islot) &
                       + spawned_walkers(spawned_pop:spawned_hf_pop,k)
                else
                    ! Found the next unique spawned walker.
                    exit compress
                end if
            end do compress
            ! All done?
            if (islot == spawning_head(0)) exit self_annihilate
            ! go to the next slot if the current determinant wasn't completed
            ! annihilated.
            if (any(spawned_walkers(spawned_pop:spawned_hf_pop,islot) /= 0)) islot = islot + 1
        end do self_annihilate

        ! We didn't check if the population on the last determinant is
        ! completely annihilated or not.
        if (all(spawned_walkers(spawned_pop:spawned_hf_pop, islot) == 0)) islot = islot - 1

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

        use basis, only: basis_length

        integer :: islot, k, pop_sign

        ! islot is the current element in the spawned walkers lists.
        islot = 1
        ! k is the current element which is being compressed into islot (if
        ! k and islot refer to the same determinants).
        k = 1
        self_annihilate: do
            ! Set the current free slot to be the next unique spawned walker.
            spawned_walkers(:,islot) = spawned_walkers(:,k) 
            compress: do
                k = k + 1
                if (k > spawning_head(0)) exit self_annihilate
                if (all(spawned_walkers(:basis_length,k) == spawned_walkers(:basis_length,islot))) then
                    ! Update the parent flag.
                    ! Note we ignore the possibility of multiple spawning events
                    ! onto the same unoccupied determinant.  Such events become
                    ! vanishingly unlikely with the size of the determinant
                    ! space (according to Ali at least!).
                    pop_sign = spawned_walkers(spawned_pop,islot)*spawned_walkers(spawned_pop,k)
                    if (pop_sign > 0) then
                        ! Sign coherent event.
                        ! Set parent_flag to 2 (indicating multiple
                        ! sign-coherent spawning events).
                        spawned_walkers(spawned_parent,islot) = 2
                    else
                        ! Keep the parent_flag of the largest spawning event.
                        if (spawned_walkers(spawned_pop,k) > spawned_walkers(spawned_pop,islot)) then
                            spawned_walkers(spawned_parent,islot) = spawned_walkers(spawned_parent,k)
                        end if
                    end if
                    ! Add the populations of the subsequent identical walkers.
                    spawned_walkers(spawned_pop,islot) = &
                                   spawned_walkers(spawned_pop,islot) + spawned_walkers(spawned_pop,k)
                else
                    ! Found the next unique spawned walker.
                    exit compress
                end if
            end do compress
            ! All done?
            if (islot == spawning_head(0)) exit self_annihilate
            ! go to the next slot if the current determinant wasn't completed
            ! annihilated.
            if (spawned_walkers(spawned_pop,islot) /= 0) islot = islot + 1
        end do self_annihilate

        ! We didn't check if the population on the last determinant is
        ! completely annihilated or not.
        if (spawned_walkers(spawned_pop, islot) == 0) islot = islot - 1

        ! update spawning_head(0)
        spawning_head(0) = islot

    end subroutine annihilate_spawned_list_initiator

    subroutine annihilate_main_list()

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        use basis, only: basis_length

        integer :: i, pos, k, nannihilate, istart, iend, old_pop(sampling_size)
        integer(i0) :: f(basis_length)
        logical :: hit

        nannihilate = 0
        istart = 1
        iend = tot_walkers
        do i = 1, spawning_head(0)
            f = spawned_walkers(:basis_length,i)
            call search_walker_list(f, istart, iend, hit, pos)
            if (hit) then
                ! Annihilate!
                old_pop = walker_population(:,pos)
                walker_population(:,pos) = walker_population(:,pos) + spawned_walkers(spawned_pop:spawned_hf_pop,i)
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
                spawned_walkers(:,k) = spawned_walkers(:,i)
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
            f = spawned_walkers(:basis_length,i)
            call search_walker_list(f, istart, iend, hit, pos)
            if (hit) then
                ! Annihilate!
                old_pop = walker_population(1,pos)
                walker_population(1,pos) = walker_population(1,pos) + spawned_walkers(spawned_pop,i)
                nannihilate = nannihilate + 1
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                nparticles = nparticles + abs(walker_population(1,pos)) - abs(old_pop)
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                ! Keep only progeny spawned by initiator determinants
                ! (parent_flag=0) or multiple sign-coherent events
                ! (parent_flag=2).
                if (spawned_walkers(spawned_parent,i) == 1) then
                    ! discard attempting spawnings from non-initiator walkers
                    ! onto unoccupied determinants.
                    ! note that the number of particles (nparticles) was not
                    ! updated at the time of spawning, so doesn't change.
                    nannihilate = nannihilate + 1
                else
                    ! keep!
                    k = i - nannihilate
                    spawned_walkers(:,k) = spawned_walkers(:,i)
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
            if (all(walker_population(:,i) == 0)) then
                nzero = nzero + 1
            else if (nzero > 0) then
                k = i - nzero
                walker_dets(:,k) = walker_dets(:,i)
                walker_population(:,k) = walker_population(:,i)
                walker_energies(:,k) = walker_energies(:,i)
            end if
        end do
        tot_walkers = tot_walkers - nzero

    end subroutine remove_unoccupied_dets

    subroutine insert_new_walkers()

        ! Insert new walkers into the main walker list from the spawned list.
        ! This is done after all particles have been annihilated, so the spawned
        ! list contains only new walkers.

        use basis, only: basis_length
        use calc, only: doing_calc, hfs_fciqmc_calc
        use determinants, only: decode_det
        use system, only: nel
        use hamiltonian, only: slater_condon0_hub_real
        use hfs_data, only: lmask, O00
        use operators, only: calc_orb_occ
        use proc_pointers, only: sc0

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
            call search_walker_list(spawned_walkers(:basis_length,i), istart, iend, hit, pos)
            ! f should be in slot pos.  Move all determinants above it.
            do j = iend, pos, -1
                ! i is the number of determinants that will be inserted below j.
                k = j + i
                walker_dets(:,k) = walker_dets(:,j)
                walker_population(:,k) = walker_population(:,j)
                walker_energies(:,k) = walker_energies(:,j)
            end do
            ! Insert new walker into pos and shift it to accommodate the number
            ! of elements that are still to be inserted below it.
            k = pos + i - 1
            walker_dets(:,k) = spawned_walkers(:basis_length,i)
            walker_population(:,k) = spawned_walkers(spawned_pop:spawned_hf_pop,i)
            nparticles = nparticles + abs(spawned_walkers(spawned_pop:spawned_hf_pop,i))
            walker_energies(1,k) = sc0(walker_dets(:,k)) - H00
            if (doing_calc(hfs_fciqmc_calc)) then
                ! Set walker_energies(2:,k) = <D_i|O|D_i>.
                walker_energies(2,k) = calc_orb_occ(walker_dets(:,k), lmask) - O00
            end if
            ! Next walker will be inserted below this one.
            iend = pos - 1
        end do
        
        ! Update tot_walkers
        tot_walkers = tot_walkers + spawning_head(0)

    end subroutine insert_new_walkers

end module annihilation
