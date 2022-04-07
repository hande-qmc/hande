module mrcc_utils

! Module containing utilities to extract a model space (set of secondary references) 
! for a multi-reference CCMC calculation from a short FCIQMC run.

! Parallelisation strategy: for each MPI rank, there exists a local running total array, which is update after each
! report loop; these arrays are merged and truncated at the end of the information-gathering period (set in input).

! This mirrors the local iteration-by-iteration update strategy, where the local running list is updated by the same 
! merge-and-truncate algorithm.

use const
use qmc_data, only: qmc_state_t, highest_dets_t

implicit none

contains

    subroutine init_highest_dets_t(qs, sys, fciqmc_in)

        ! Initialises the derived type storing highest (cumulative) occupied determinants

        use parallel
        use errors, only: stop_all
        use qmc_data, only: fciqmc_in_t
        use system, only: sys_t
        use checking, only: check_allocate

        type(fciqmc_in_t), intent(in) :: fciqmc_in
        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(inout) :: qs

        integer :: ierr
        

        ! Brian: [todo] - These should be sorted before review!
        if (nprocs > 1) then
            call stop_all('init_highest_dets_t', 'FCIQMC-tailored mr-CCMC does not currently support MPI')
        end if

        if (sys%basis%bit_string_len > 1) then
            call stop_all('init_highest_dets_t', 'FCIQMC-tailored mr-CCMC does not currently handle more than 64 basis functions.')
        end if

        qs%doing_fciqmc_tailored_mrccmc = .true.
        qs%highest_dets%nreports = fciqmc_in%mrcc_selection_nreports
        qs%highest_dets%target_ndets = fciqmc_in%mrcc_selection_target_ndets
        qs%highest_dets%max_ndets = fciqmc_in%mrcc_selection_max_ndets
        qs%highest_dets%det_file = fciqmc_in%mrcc_selection_dets_file

        allocate(qs%highest_dets%det_list_main((1+sys%basis%bit_string_len),qs%highest_dets%max_ndets), stat=ierr)
        call check_allocate('qs%highest_dets%det_list_main', (1+sys%basis%bit_string_len)*qs%highest_dets%max_ndets, ierr)

        ! We only collect the highest target_ndets determinants from each iteration, and then sort and merge it into the main
        ! determinant list, which has more than target_ndets slots, so that we have some buffer space in the scenario where
        ! some determinants are initially less occupied but becomes more occupied later on
        allocate(qs%highest_dets%det_list_new((1+sys%basis%bit_string_len),qs%highest_dets%target_ndets), stat=ierr)
        call check_allocate('qs%highest_dets%det_list_new', (1+sys%basis%bit_string_len)*qs%highest_dets%target_ndets, ierr)

        qs%highest_dets%ndets = 0

    end subroutine init_highest_dets_t

    subroutine update_highest_dets(qs)

        ! Called after each report loop to update the main list

        type(qmc_state_t), intent(inout) :: qs

        call get_highest_occupied_dets(qs%highest_dets, qs%psip_list)
        call merge_highest_dets(qs%highest_dets)

    end subroutine update_highest_dets

    subroutine get_highest_occupied_dets(hd, psip_list)

        ! Extract from a psip_list the (absolute) highest occupied determinants.
        ! This just wraps around semi_stoch::find_most_populated_dets, which does the same thing.

        use qmc_data, only: particle_t
        use semi_stoch, only: find_most_populated_dets
        use sort, only: qsort

        type(highest_dets_t), intent(inout) :: hd
        type(particle_t), intent(in) :: psip_list

        ! The second last argument is for the array holding determinants, and the last for the corresponding population
        call find_most_populated_dets(psip_list%states, psip_list%pops, psip_list%nstates, hd%target_ndets, &
                                      hd%det_list_new(2:,:), hd%det_list_new(1,:))
        ! find_most_populated_dets returns the absolute population, so we don't need to worry about sorting

        ! Only use the population to sort (nsort=1)
        call qsort(hd%det_list_new, nsort=1)

    end subroutine get_highest_occupied_dets

    subroutine merge_highest_dets(hd)

        ! Merge det_list_new with det_list_main, the former containing target_ndets determinants, and the latter containing 
        ! (max_ndets >=) ndets (>= target_ndets) determinants, which is larger.

        use search, only: linear_search, binary_search
        use sort, only: qsort

        type(highest_dets_t), intent(inout) :: hd

        integer :: i, pos
        logical :: hit

        if (hd%ndets == 0) then
            ! This is the first time we're merging, just copy everything over
            hd%ndets = hd%target_ndets
            ! qsort puts the largest element last, we follow that ordering
            hd%det_list_main(:,hd%ndets+1:) = hd%det_list_new
        else
            ! We need to merge the new list and the main list
            do i = hd%target_ndets, 1, -1
                ! Search through the (unsorted since det_list_main is sorted by population) determinant (sub)list
                call linear_search(hd%det_list_main(2:,:), hd%det_list_new(2:,i), hd%max_ndets-hd%ndets+1, hd%max_ndets, hit, pos)
                if (hit) then
                    ! Augment the population
                    hd%det_list_main(1, pos) = hd%det_list_main(1, pos) + hd%det_list_new(1, i)
                else 
                    ! Now that it's not a hit, we just prepend it, and do a sort afterwards

                    if (hd%ndets < hd%max_ndets) then
                        ! Come on in we have room
                        hd%det_list_main(:, hd%max_ndets-hd%ndets) = hd%det_list_new(:, i)
                        hd%ndets = hd%ndets + 1
                    else
                        ! The main list is full, we only keep it if it's not the smallest element
                        if (hd%det_list_new(1, i) > hd%det_list_main(1, 1)) then
                            ! We need to insert it, and discard the last element
                            hd%det_list_main(:, 1) = hd%det_list_new(:, i)
                        end if
                    end if
                    call qsort(hd%det_list_main, nsort=1)
                end if
            end do
        end if

    end subroutine merge_highest_dets

    subroutine write_out_highest_dets(sys, hd, filename)

        use parallel
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(highest_dets_t), intent(in) :: hd
        character(*), intent(in) :: filename

        integer :: fileunit, i

        if (parent) then
            open(file=filename, newunit=fileunit, action='write', status='replace')
            do i = hd%max_ndets, hd%max_ndets-hd%target_ndets+1, -1
                write(fileunit, '(I0, 1X, I0)') hd%det_list_main(2:,i), hd%det_list_main(1,i)
            end do
            close(fileunit)
        end if
        
    end subroutine

end module mrcc_utils