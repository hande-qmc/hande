module uccmc_utils

use const, only: i0, int_p, int_64, p, dp, debug

implicit none

contains

    subroutine allocate_time_average_lists(qs, states, pops, nstates)

        use qmc_data, only: qmc_state_t

        type(qmc_state_t), intent(in) :: qs
        integer(i0), intent(out), allocatable :: states(:,:)
        real(p), intent(out), allocatable :: pops(:)
        integer, intent(out) :: nstates

        allocate(states(size(qs%psip_list%states(:,1)),size(qs%psip_list%states(1,:))))
        allocate(pops(size(qs%psip_list%states(1,:))))
        states(:,:qs%psip_list%nstates)  = qs%psip_list%states(:,:qs%psip_list%nstates)
        pops(:) = (real(qs%psip_list%pops(1,:))/qs%psip_list%pop_real_factor)
        nstates = qs%psip_list%nstates
    end subroutine allocate_time_average_lists

    subroutine add_ci_contribution(cluster, cdet, states, pops, nstates)

        use determinant_data, only: det_info_t
        use ccmc_data, only: cluster_t

        type(det_info_t), intent(in) :: cdet
        type(cluster_t), intent(in) :: cluster
        integer(i0), intent(inout) :: states(:,:)
        integer, intent(inout) :: nstates
        real(p), intent(inout) :: pops(:)

        integer(i0) :: state(size(states(:,1)))
        logical :: hit
        integer :: pos
        real(p) :: population

        state = cdet%f 
        hit = .false.
        do pos = 1, nstates
            if (all(state == states(:,pos))) then
                hit = .true.
                exit
            end if
        end do
        population = cluster%amplitude*cluster%cluster_to_det_sign/cluster%pselect
        if (hit) then
           pops(pos) = pops(pos) + population 
        else
            states(:,nstates+1) = state
            pops(nstates+1) = population
            nstates = nstates + 1
        end if
    end subroutine add_ci_contribution

    subroutine add_t_contributions(qs, states, pops, sq, nstates)

        use qmc_data, only: qmc_state_t
        use search, only: binary_search

        type(qmc_state_t), intent(in) :: qs
        integer(i0), intent(inout) :: states(:,:)
        real(p), intent(inout) :: pops(:), sq(:,:)
        integer, intent(inout) :: nstates

        integer :: i, j, k, pos
        logical :: hit
        integer(i0) :: state(size(states(:,1)))

        do i = 1, qs%psip_list%nstates
            state = qs%psip_list%states(:,i) 
            call binary_search(states, state, 1, nstates, hit, pos)
            if (hit) then
                  pops(pos) = pops(pos) + (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)
                  sq(2,pos) = sq(2,pos) + (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)**2 
               else
                   do j = nstates, pos, -1
                       k = j + 1 
                       states(:,k) = states(:,j)
                       pops(k) = pops(j)
                       sq(1,k) = sq(1,j)
                       sq(2,k) = sq(2,j)
                   end do
                   states(:,pos) = qs%psip_list%states(:,i)
                   pops(pos) = (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)
                   sq(1,pos) = qs%psip_list%states(1,i)
                   sq(2,pos) = (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)**2
                   nstates = nstates + 1
               end if
        end do
    end subroutine add_t_contributions
    
    pure function earliest_unset(f, f0, nel, basis) result (early)
        
         ! Function to find earliest unset bit in a determinant bit string.
! [review] - AJWT:  Not really what this does - it finds the earliest different to f0.
! [review] - AJWT:  Perhaps split into two functions, earliest_unset and earliest_diff?

         ! In:
         !    f: bit_string encoding determinant
! [review] - AJWT:  Explain what f0 is doing here
         !    nel: number of electrons in the system
         !    basis: basis_t object with information on one-electron basis in use.
         !
         ! Returns:
         !    Zero-based index of the earliest unset bit.
     
         use basis_types, only: basis_t
         use bit_utils, only: count_set_bits

         type(basis_t), intent(in) :: basis
! [review] - AJWT: Do we need the whole tot_string_len when we actually only care about the bit_string_len?
         integer(i0), intent(in) :: f(basis%tot_string_len), f0(basis%tot_string_len)
         integer, intent(in) :: nel
         integer :: i, early
         integer(i0) :: diff
         integer(i0) :: f0_loc
         
         early = 0
         i = 0
          
! [review] - AJWT: This is very dangerous to only take the least significant word - it should be extended to the full bit_string_len
         diff = ieor(f(1), f0(1))
         if (diff /= 0) then
             do 
                 if (btest(diff, i)) then
                     early = i
                     exit
                 end if
                 i = i + 1
             end do
         else
! [review] - AJWT: I think this is finding the the highest set bit in f0, but I don't know why.
             f0_loc = f0(1)
             do while (f0_loc/=0)
                 f0_loc = ibclr(f0_loc, early)
                 early = early + 1
             end do
         end if
    end function

    pure function latest_unset(f, f0, nel, basis) result (late)
        
         ! Function to find earliest unset bit in a determinant bit string.
! [review] - AJWT:  Not really what this does - it finds the latest different to f0 in the occupied in f0.
! [review] - AJWT:  Perhaps split into two functions, latest_unset and latest_diff?

         ! In:
         !    f: bit_string encoding determinant
! [review] - AJWT:  Explain what f0 is doing here
         !    nel: number of electrons in the system
         !    basis: basis_t object with information on one-electron basis in use.
     
         use basis_types, only: basis_t
         use bit_utils, only: count_set_bits

         type(basis_t), intent(in) :: basis
         integer(i0), intent(in) :: f(basis%tot_string_len), f0(basis%tot_string_len)
         integer, intent(in) :: nel
         integer :: late
         integer(i0) :: diff
         integer(i0) :: f0_loc
         
         late = 0
          
! [review] - AJWT: This is very dangerous to only take the least significant word - it should be extended to the full bit_string_len
         diff = ieor(f(1),f0(1))
         diff = iand(diff, f0(1))
         if (diff /= 0) then
             do 
                 diff = ibclr(diff, late)
                 if (diff /= 0) then
                     late = late + 1 
                 else
                     exit
                 end if
             end do
         else
! [review] - AJWT: I think this is finding the the highest set bit in f0, but I don't know why.
             f0_loc = f0(1)
             do while (f0_loc/=0)
                 f0_loc = ibclr(f0_loc, late)
                 late = late + 1
             end do
         end if
    end function

    subroutine add_info_str_trot(basis, f0, nel, f)

! [review] - AJWT: This also stores the 'latest' unset bit.
        ! Sets bits within bit string to give excitation level at end of bit strings.
        ! This routine sets ex level from provided reference.

        use basis_types, only: basis_t
        use excitations, only: get_excitation_level
        type(basis_t), intent(in) :: basis
        integer(i0), intent(inout) :: f(:)
        integer(i0), intent(in) :: f0(:)
        integer, intent(in) :: nel

        integer(i0) :: counter(basis%tot_string_len)
       


! [review] - AJWT: Presumable one needs info_string_len >=2.
        if (basis%info_string_len/=0) then

            f(basis%bit_string_len+2) = latest_unset(f, f0, nel, basis)     
            f(basis%bit_string_len+1) = nel - get_excitation_level(f(:basis%bit_string_len), f0(:basis%bit_string_len)) 

        end if
 
    end subroutine add_info_str_trot

end module
