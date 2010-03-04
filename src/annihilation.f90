module annihilation

use const
use fciqmc_data

implicit none

contains

    subroutine direct_annihilation()

        ! 1. Sort spawned walkers list.
        call sort_spawned_lists()

        ! 2. Annihilate within spawned walkers list.
        ! Compress the remaining spawned walkers list.
        call annihilate_spawned_list()

        ! 4. Annilate main list.
        call annihilate_main_list()

        ! 5. Insert new walkers into main walker list.
        call insert_new_walkers()

    end subroutine direct_annihilation

    subroutine annihilate_spawned_list()

        ! Annihilate the spawned walker list and compress the remaining
        ! elements.

        ! The spawned walker list is already sorted, so annihilation amounts to
        ! looping through the list and adding consective walker populations
        ! together if they're the same walker.

        integer :: i, k, nremoved

        nremoved = 0
        i = 0
        k = 0
        annihilate: do
            i = i + 1
            if (i > spawning_head) exit annihilate
            k = k + 1
            ! Compress.
            spawned_walker_dets(:,k) = spawned_walker_dets(:,i)
            spawned_walker_population(k) = spawned_walker_population(i)
            i = i + 1
            do while (all(spawned_walker_dets(:,k) == spawned_walker_dets(:,i)))
                ! Annihilate.
                spawned_walker_population(k) = spawned_walker_population(k) + spawned_walker_population(i)
                nremoved = nremoved + 1
                i = i + 1
                if (i > spawning_head) exit annihilate
            end do
        end do annihilate

        ! The last element in the spawned list is now k.
        spawning_head = k

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
        do i = 1, spawning_head
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

        spawning_head = spawning_head - nannihilate

        if (nannihilate > 0) then
            ! Remove any determinants with 0 population.
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
        end if

    end subroutine annihilate_main_list

    subroutine insert_new_walkers()

        ! Insert new walkers into the main walker list from the spawned list.
        ! This is done after all particles have been annihilated, so the spawned
        ! list contains only new walkers.

        use basis, only: basis_length
        use determinants, only: decode_det
        use system, only: nel
        use hamiltonian, only: slater_condon0_hub_k

        integer :: i, istart, iend, j, k, pos, occ_list(nel)
        integer(i0) :: f(basis_length)
        logical :: hit

        ! Merge new walkers into the main list.
        istart = 1
        iend = tot_walkers
        do i = spawning_head, 1, -1
            f = spawned_walker_dets(:,i)
            ! f is not in the main walker list
            call search_walker_list(f, istart, iend, hit, pos)
            ! f should be in slot pos.  Move all determinants above it.
            do j = iend, pos+1, -1
                ! i is the number of determinants that will be inserted below j.
                k = pos + i
                walker_dets(:,k) = walker_dets(:,i)
                walker_population(k) = walker_population(i)
                walker_energies(k) = walker_energies(i)
            end do
            ! Insert new walker
            walker_dets(:,pos) = f
            walker_population(pos) = spawned_walker_population(i)
            occ_list = decode_det(f)
            walker_energies(pos) = slater_condon0_hub_k(occ_list)
            ! Next walker will be inserted below this one.
            iend = pos
        end do
        
        ! Update tot_walkers
        tot_walkers = tot_walkers + spawning_head

    end subroutine insert_new_walkers

end module annihilation
