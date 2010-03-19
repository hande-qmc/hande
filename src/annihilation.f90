module annihilation

use const
use fciqmc_data

implicit none

contains

    subroutine direct_annihilation(sc0)

        interface
            function sc0(f) result(hmatel)
                use basis, only: basis_length
                use const, only: i0, dp
                implicit none
                real(dp) :: hmatel
                integer(i0), intent(in) :: f(basis_length)
            end function sc0
        end interface

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

    subroutine annihilate_spawned_list()

        ! Annihilate the spawned walker list and compress the remaining
        ! elements.

        ! The spawned walker list is already sorted, so annihilation amounts to
        ! looping through the list and adding consective walker populations
        ! together if they're the same walker.

        integer :: i, k, nremoved

        nremoved = 0
        i = 1
        k = 0

        ! Go through the sorted spawned_walker lists and compress it by adding
        ! up the populations of duplicate walkers and removing the duplicate
        ! determinants.

        ! Treat the first entry separately to make the subsequent removal of
        ! walkers with 0 population easy.
        if (spawning_head == 1) then
            ! Nothing to do bar ensure spawning_head remains 1.
            k = 1
            i = i + 1
        else if (spawning_head > 1) then
            k = k + 1
            ! Compress.
            i = i + 1
            do while (all(spawned_walker_dets(:,k) == spawned_walker_dets(:,i)))
                ! Annihilate.
                spawned_walker_population(k) = spawned_walker_population(k) + spawned_walker_population(i)
                nremoved = nremoved + 1
                i = i + 1
                if (i > spawning_head) exit
            end do
        end if

        annihilate: do
            if (i > spawning_head) exit annihilate
            if (spawned_walker_population(k) /= 0) then
                ! Go to the next slot.
                k = k + 1
            end if
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
                use const, only: i0, dp
                implicit none
                real(dp) :: hmatel
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
        do i = spawning_head, 1, -1
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
        tot_walkers = tot_walkers + spawning_head

    end subroutine insert_new_walkers

end module annihilation
