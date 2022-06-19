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

    subroutine init_highest_dets_t(qs, sys, fciqmc_in, ss_in)

        ! Initialises the derived type storing highest (cumulative) occupied determinants

        use parallel
        use errors, only: stop_all
        use qmc_data, only: fciqmc_in_t, semi_stoch_in_t
        use system, only: sys_t
        use checking, only: check_allocate

        type(fciqmc_in_t), intent(in), optional :: fciqmc_in
        type(semi_stoch_in_t), intent(in), optional :: ss_in
        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(inout) :: qs

        integer :: ierr, nelements
        logical :: fciqmc_set, ss_set

        fciqmc_set = .false.
        ss_set = .false.
        if (present(fciqmc_in)) fciqmc_set = .true.
        if (present(ss_in)) ss_set = .true.
        
        ! Brian: [todo] - These should be sorted before review!
        if (nprocs > 1) then
            call stop_all('init_highest_dets_t', 'FCIQMC-tailored mr-CCMC does not currently support MPI')
        end if

        if (sys%basis%bit_string_len > 1) then
            call stop_all('init_highest_dets_t', 'FCIQMC-tailored mr-CCMC does not currently handle more than 64 basis functions.')
        end if

        qs%averaging = .true.

        if (fciqmc_set) then
            qs%highest_dets%nreports = fciqmc_in%mrcc_selection_nreports

            if (fciqmc_in%mrcc_selection_target_ndets > 0) then
                qs%highest_dets%target_ndets = fciqmc_in%mrcc_selection_target_ndets
                qs%highest_dets%using_target = .true.
            else if (fciqmc_in%mrcc_selection_wfn_weightage > 0.0_p) then
                qs%highest_dets%target_ndets = fciqmc_in%mrcc_selection_det_list_new_size
                allocate(qs%highest_dets%nw_cumul_proc(nprocs), source=0.0_dp, stat=ierr)
                call check_allocate('nw_cumul_proc', nprocs, ierr)
                qs%highest_dets%wfn_weightage = fciqmc_in%mrcc_selection_wfn_weightage
                qs%highest_dets%using_weightage = .true.
            end if
            
            qs%highest_dets%max_ndets = fciqmc_in%mrcc_selection_max_ndets
            qs%highest_dets%det_file = fciqmc_in%mrcc_selection_dets_file
        else if (ss_set) then
            qs%highest_dets%nreports = ss_in%nreports
            if (ss_in%target_size > 0) then
                qs%highest_dets%target_ndets = ss_in%target_size
                qs%highest_dets%using_target = .true.
                qs%highest_dets%max_ndets = ss_in%max_ndets
            else if (ss_in%wfn_weightage > 0.0_p) then
                ! Just to be safe..
                qs%highest_dets%target_ndets = 2*qs%psip_list%nstates
                allocate(qs%highest_dets%nw_cumul_proc(nprocs), source=0.0_dp, stat=ierr)
                call check_allocate('nw_cumul_proc', nprocs, ierr)
                qs%highest_dets%wfn_weightage = ss_in%wfn_weightage
                qs%highest_dets%using_weightage = .true.

                ! [todo] - document why this is different than fciqmc_set, or better yet make them consistent
                qs%highest_dets%max_ndets = 2*qs%psip_list%nstates
            end if
        end if

        if (qs%highest_dets%using_target) then
            ! Both the absolute and the signed populations need to be stored
            nelements = 2 + sys%basis%bit_string_len
        else
            ! No truncation will be done, so no absolute population need to be stored
            nelements = 1 + sys%basis%bit_string_len
        end if


        allocate(qs%highest_dets%det_list_main(nelements,qs%highest_dets%max_ndets), stat=ierr)
        call check_allocate('qs%highest_dets%det_list_main', nelements*qs%highest_dets%max_ndets, ierr)

        ! We only collect the highest target_ndets determinants from each report, and then sort and merge it into the main
        ! determinant list, which has more than target_ndets slots, so that we have some buffer space in the scenario where
        ! some determinants are initially less occupied but becomes more occupied later on
        allocate(qs%highest_dets%det_list_new(nelements,qs%highest_dets%target_ndets), stat=ierr)
        call check_allocate('qs%highest_dets%det_list_new', nelements*qs%highest_dets%target_ndets, ierr)

        qs%highest_dets%ndets = 0

    end subroutine init_highest_dets_t

    subroutine update_highest_dets(qs, sys)

        ! Called after each report loop to update the main list

        use, intrinsic :: iso_fortran_env, only: error_unit
        use utils, only: int_fmt
        use parallel, only: iproc
        use system, only: sys_t

        type(qmc_state_t), intent(inout) :: qs
        type(sys_t), intent(in) :: sys
        real(dp) :: norm

        norm = sqrt(sum((real(qs%psip_list%pops(1,:qs%psip_list%nstates),p)/qs%psip_list%pop_real_factor)**2))

        ! iproc is 0-indexed
        call get_highest_occupied_dets(sys, qs%highest_dets, qs%psip_list, norm)
        call merge_highest_dets(sys, qs%highest_dets)

        if (qs%highest_dets%error) then
            write (error_unit,'(1X,"# Error: No space left in highest_det%det_list_main on processor",'//int_fmt(iproc,1)//',".")')&
                              iproc
            write (error_unit,'(1X,"# Error: HANDE will exit at the end of this report loop.")')
        end if

    end subroutine update_highest_dets

    subroutine get_highest_occupied_dets(sys, hd, psip_list, norm)

        ! Extract from a psip_list the (absolute) highest occupied determinants.
        ! This just wraps around semi_stoch::find_most_populated_dets, which does the same thing.

        use system, only: sys_t
        use qmc_data, only: particle_t
        use semi_stoch, only: find_most_populated_dets
        use sort, only: qsort

        type(highest_dets_t), intent(inout) :: hd
        type(particle_t), intent(in) :: psip_list
        real(dp), intent(in) :: norm
        type(sys_t), intent(in) :: sys

        integer :: target_ndets

        ! Careful not to collect meaningless data
        target_ndets = min(hd%target_ndets, psip_list%nstates)

        if (hd%using_target) then
            ! The second last argument is for the array holding determinants, and the last for the corresponding population
            call find_most_populated_dets(psip_list%states, psip_list%pops, psip_list%nstates, target_ndets, &
                                          hd%det_list_new(3:,:), hd%det_list_new(2,:), abs_pop = .false.)

            ! The first element is absolute pop, second is the signed population
            hd%det_list_new(1,:) = abs(hd%det_list_new(2,:))

            ! Only use the absolute population to sort (nsort=1)
            hd%ndets_new = target_ndets
            call qsort(hd%det_list_new, head=hd%ndets_new, nsort=1)
            
        else if (hd%using_weightage) then
            call truncate_wfn(psip_list, psip_list%nstates, norm, hd%wfn_weightage, &
                              sys%basis%bit_string_len, hd%ndets_new, hd%det_list_new)

            ! Sort this (hopefully short) list by determinants before merging with main list
            call qsort(hd%det_list_new, hd%ndets_new, sys%basis%bit_string_len)
        end if

    end subroutine get_highest_occupied_dets

    subroutine merge_highest_dets(sys, hd)

        ! Merge det_list_new with det_list_main, the former containing target_ndets determinants, and the latter containing 
        ! (max_ndets >=) ndets (>= target_ndets) determinants, which is larger.

        use system, only: sys_t
        use search, only: linear_search, binary_search
        use sort, only: qsort

        type(sys_t), intent(in) :: sys
        type(highest_dets_t), intent(inout) :: hd

        integer :: i, j, k, pos, istart, iend, nmerge
        logical :: hit
        integer(i0) :: f(sys%basis%bit_string_len)

        if (hd%ndets == 0) then
            ! This is the first time we're merging, just copy everything over
            ! qsort puts the largest element last
            hd%ndets = hd%ndets_new
            hd%det_list_main(:,1:hd%ndets) = hd%det_list_new(:,1:hd%ndets)
        else
            ! We need to merge the new list and the main list
            if (hd%using_target) then
                do i = hd%ndets_new, 1, -1
                    ! Search through the (unsorted since det_list_main is sorted by population) determinant (sub)list
                    call linear_search(hd%det_list_main(3:,:), hd%det_list_new(3:,i), 1, hd%ndets, hit, pos)
                    if (hit) then
                        ! Augment the (signed) population
                        hd%det_list_main(2, pos) = hd%det_list_main(2, pos) + hd%det_list_new(2, i)
                        hd%det_list_main(1, pos) = abs(hd%det_list_main(2, pos))
                    else 
                        ! Now that it's not a hit, we just append it, and do a sort afterwards

                        if (hd%ndets < hd%max_ndets) then
                            ! Come on in we have room
                            hd%det_list_main(:, hd%ndets+1) = hd%det_list_new(:, i)
                            hd%ndets = hd%ndets + 1
                        else
                            ! The main list is full, we only keep it if it's not the smallest element
                            if (hd%det_list_new(1, i) > hd%det_list_main(1, 1)) then
                                ! We need to insert it, and discard the last element
                                hd%det_list_main(:, 1) = hd%det_list_new(:, i)
                            end if
                        end if
                        call qsort(hd%det_list_main, head=hd%ndets, nsort=1)
                    end if
                end do
            else if (hd%using_weightage) then
                ! This is like the spawn to main merge
                nmerge = 0
                istart = 1
                iend = hd%ndets
                do i = 1, hd%ndets_new
                    f = int(hd%det_list_new(:sys%basis%bit_string_len,i), i0)
                    call binary_search(hd%det_list_main(:sys%basis%bit_string_len,:), f, istart, iend, hit, pos)
                    if (hit) then
                        ! Merge populations
                        hd%det_list_main(sys%basis%bit_string_len+1, pos) = hd%det_list_main(sys%basis%bit_string_len+1, pos) + &
                                                                      hd%det_list_new(sys%basis%bit_string_len+1, i)
                        nmerge = nmerge + 1
                        istart = pos + 1
                    else
                        ! Compress new list.
                        k = i - nmerge
                        hd%det_list_new(:,k) = hd%det_list_new(:,i)
                    end if
                end do
                hd%ndets_new = hd%ndets_new - nmerge
                ! Now merge
                if (hd%ndets_new + hd%ndets > hd%max_ndets) then
                    ! Handle the error in update_highest_dets
                    hd%error = .true.
                else
                    istart = 1
                    iend = hd%ndets

                    do i = hd%ndets_new, 1, -1
                        f = int(hd%det_list_new(:sys%basis%bit_string_len,i), i0)
                        call binary_search(hd%det_list_main(:sys%basis%bit_string_len,:), f, istart, iend, hit, pos)
                        ! f should be in slot pos. Move all determinants above it.
                        do j = iend, pos, -1
                            ! i is the number of determinants that will be inserted below j.
                            k = j + i
                            hd%det_list_main(:, k) = hd%det_list_main(:, j)
                        end do

                        ! Insert new determinant into pos
                        k = pos + i - 1
                        hd%det_list_main(:,k) = hd%det_list_new(:,i)

                        ! Next determinant will be inserted below this one.
                        iend = pos - 1
                    end do
                    hd%ndets = hd%ndets + hd%ndets_new
                end if
            end if
        end if

    end subroutine merge_highest_dets

    subroutine write_out_highest_dets(qs, sys, hd, filename)

        use parallel
        use errors, only: warning
        use system, only: sys_t
        use bit_utils, only: bit_str_cmp
        use search, only: linear_search, binary_search
        use sort, only: argsort
        use checking, only: check_allocate, check_deallocate

        type(qmc_state_t), intent(in) :: qs
        type(sys_t), intent(in) :: sys
        type(highest_dets_t), intent(inout) :: hd
        character(*), intent(in) :: filename

        integer :: det_unit, det_pop_unit, i, pos, ierr
        integer, allocatable :: args(:)
        logical :: done_ref, hit
        integer(i0) :: bs_flipped(sys%basis%bit_string_len)
        real(p) :: weightage
        real(dp) :: norm

        done_ref = .false.

        deallocate(hd%det_list_new, stat=ierr)
        call check_deallocate('hd%det_list_new', ierr)

        if (parent) then
            if (hd%using_target) then
                open(file=filename, newunit=det_unit, action='write', status='replace')
                open(file=filename//'.pop', newunit=det_pop_unit, action='write', status='replace')
                
                do i = hd%ndets, hd%ndets-hd%target_ndets+1, -1
                    if (.not. done_ref) then
                        ! We can't have the reference determinant in the list of secondary references, 
                        ! so actively check if the current determinant is the reference, and do not write out if it is
                        done_ref = (bit_str_cmp(hd%det_list_main(3:,i), qs%ref%f0) == 0)
                        if (i == hd%max_ndets .and. .not. done_ref) then
                            ! Warn the user if the highest occupied determinant is not the reference
                            call warning('write_out_highest_dets','The reference is not the highest occupied determinant!',1)
                        end if
                        if (.not. done_ref) then
                            ! Not reference and haven't seen reference yet, write out current list and keep an eye out still
                            write(det_unit, '(I0, 1X, I0)') hd%det_list_main(3:,i)
                            write(det_pop_unit, '(I0, 1X, I0)') hd%det_list_main(3:,i), hd%det_list_main(1,i)
                        end if
                    else
                        ! Don't need to check if this is the reference as we've seen it already.
                        write(det_unit, '(I0, 1X, I0)') hd%det_list_main(3:,i)
                        write(det_pop_unit, '(I0, 1X, I0)') hd%det_list_main(3:,i), hd%det_list_main(1,i)
                    end if

                    if (sys%Ms == 0) then
                        ! Ms=0 configurations that are true eigenfunctions of S^2 should be formed by a linear combination 
                        ! of ab and ba (i.e., singlet: (ab-ba)/sqrt(2), triplet Ms=0: (ab+ba)/sqrt(2)), so should always include
                        ! a spin-flipped version of the determinant, if not already present
                        ! Triplet Ms = +-1 states are already eigenfunctions of S^2 so do not need augmenting
                        ! Higher spin multiplicities are not supported yet.
                        call get_spin_flipped_bstring(hd%det_list_main(3:,i), bs_flipped, sys)
                        call linear_search(hd%det_list_main(3:,:), bs_flipped, hd%max_ndets-hd%target_ndets+1, &
                                           hd%max_ndets, hit, pos)
                        if (.not. hit) then
                            write(det_unit, '(I0, 1X, I0)') bs_flipped
                        end if
                    end if
                end do
                close(det_pop_unit)
                close(det_unit)
            else

                ! for now print out everything
                open(file=filename, newunit=det_unit, action='write', status='replace')
                open(file=filename//'.pop', newunit=det_pop_unit, action='write', status='replace')
                
                do i = 1, hd%ndets
                    if (.not. (bit_str_cmp(hd%det_list_main(:sys%basis%bit_string_len,i), qs%ref%f0) == 0)) then
                        write(det_unit, '(I0, 1X, I0)') hd%det_list_main(:sys%basis%bit_string_len,i)
                        write(det_pop_unit, '(I0, 1X, I0)') hd%det_list_main(:sys%basis%bit_string_len,i), &
                                                            hd%det_list_main(sys%basis%bit_string_len+1,i)
                    end if
                    
                    if (sys%Ms == 0) then
                        ! Ms=0 configurations that are true eigenfunctions of S^2 should be formed by a linear combination 
                        ! of ab and ba (i.e., singlet: (ab-ba)/sqrt(2), triplet Ms=0: (ab+ba)/sqrt(2)), so should always include
                        ! a spin-flipped version of the determinant, if not already present
                        ! Triplet Ms = +-1 states are already eigenfunctions of S^2 so do not need augmenting
                        ! Higher spin multiplicities are not supported yet.
                        call get_spin_flipped_bstring(hd%det_list_main(:sys%basis%bit_string_len,i), bs_flipped, sys)
                        call binary_search(hd%det_list_main(:sys%basis%bit_string_len,:), bs_flipped, 1, hd%ndets, hit, pos)
                        if (.not. hit) then
                            write(det_unit, '(I0, 1X, I0)') bs_flipped
                        end if
                    end if
                end do
                close(det_pop_unit)
                close(det_unit)
            end if
        end if

    end subroutine

    pure subroutine get_spin_flipped_bstring(bs_in, bs_out, sys)

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: bs_in(sys%basis%bit_string_len)
        integer(i0), intent(out) :: bs_out(sys%basis%bit_string_len)

        integer(i0), dimension(sys%basis%bit_string_len) :: bs_alpha, bs_beta

        bs_alpha = iand(bs_in, sys%basis%alpha_mask)
        bs_beta = iand(bs_in,sys%basis%beta_mask)

        bs_out = ior(ishft(bs_alpha, 1), ishft(bs_beta, -1))

    end subroutine get_spin_flipped_bstring

    subroutine truncate_wfn(psip_list_in, ndets_in, norm, wfn_weight, bit_string_len, ndets_out, wfn_out)

        ! On output dets_out contains the determinants whose combined weightage makes up the top wfn_weight of the local
        ! (i.e. on-processor) slice of the wavefunction, pops_out the corresponding signed populations, and these will be sorted 
        ! by the absolute population, finally ndets_out is the number of states included in the lists.

        ! NOTE: It is assumed that size(dets_in,2) >= size(dets_out,2) and
        ! similarly for pops_in and pops_out.

        ! In:
        !    dets_in: List of determinants to consider.
        !    pops_in: Populations corresponding to the determinants in dets_in.
        !    ndets_in: The number of determinants in dets_in to consider.
        !    ndets_out: The number of determinants to keep when finding the
        !        most populated determinants. Anything beyond ndets_out in
        !        dets_out and pops_out may be garbage.
        !    norm: the normalisation of the current (slice of if using MPI) wavefunction
        ! Out:
        !    dets_out: List of most populated determinants.
        !    pops_out: Absolute populations corresponding to the determinants in dets_out.

        use sort, only: argsort
        use checking, only: check_allocate, check_deallocate
        use qmc_data, only: particle_t
        use errors, only: stop_all

        type(particle_t), intent(in) :: psip_list_in
        integer, intent(in) :: ndets_in
        real(dp), intent(in) :: norm
        real(p), intent(in) :: wfn_weight
        integer, intent(in) :: bit_string_len
        
        integer, intent(out) :: ndets_out
        ! wfn_out(:,i) is the truncated version of wfn_in, holds determinant i (in the 1:bit_string_len elements) 
        ! and the signed population on it (in the bit_string_len+1 element)
        integer(int_s), intent(out) :: wfn_out(:,:)
        
        integer :: i, j, min_ind, ierr
        integer(i0) :: min_pop
        integer, allocatable :: args(:)
        real(p) :: weight

        allocate(args(ndets_in), source=0, stat=ierr)
        call check_allocate('args', ndets_in, ierr)

        ! You're out of luck if you're doing Hellmann-Feynman sampling
        call argsort(psip_list_in%pops(:,1:ndets_in), args, ndets_in, 1, sort_abs=.true.)

        weight = 0.0_p
        i = ndets_in
        j = 1

        do while (weight < wfn_weight)
            ! Note that if we use MPI, this will be taking wfn_weight of each of the local slice of the wavefunction,
            ! to be combined at the end.
            weight = weight + ((real(psip_list_in%pops(1, args(i)),p)/psip_list_in%pop_real_factor)/norm)**2
            ! Bitstring
            wfn_out(1:bit_string_len,j) = int(psip_list_in%states(1:bit_string_len,args(i)), int_s)
            ! Population
            wfn_out(bit_string_len+1,j) = int(psip_list_in%pops(1,args(i)), int_s)
            i = i - 1
            j = j + 1
            if (j > ubound(wfn_out,dim=2)) call stop_all('truncate_wfn', & 
            'det_list_new is not big enough to hold the specified percentage of the tailoring wavefunction!')
        end do

        ndets_out = j - 1

    end subroutine truncate_wfn

end module mrcc_utils