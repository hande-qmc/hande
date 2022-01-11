module uccmc_utils

! Module containing utility functions only used in UCCMC.

use const, only: i0, int_p, int_64, p, dp, debug

implicit none

contains

    subroutine allocate_time_average_lists(psip_list, states, pops, nstates)

        ! Allocate arrays to store average wavefunction in UCCMC. 

        ! In:
        !    psip_list: particle_t object storing current state of the QMC system. 
        ! Out:
        !    states: array for excip labels, initialised to current qs%psip_list%states
        !    pops: array for excip populations, initialised to current qs%psip_list%pops(1,:)
        !    nstates: number of excips

        use qmc_data, only: particle_t

        type(particle_t), intent(in) :: psip_list 
        integer(i0), intent(out), allocatable :: states(:,:)
        real(p), intent(out), allocatable :: pops(:)
        integer, intent(out) :: nstates

        allocate(states(size(psip_list%states(:,1)),size(psip_list%states(1,:))))
        allocate(pops(size(psip_list%states(1,:))))
        states(:,:psip_list%nstates)  = psip_list%states(:,:psip_list%nstates)
        pops(:) = (real(psip_list%pops(1,:))/psip_list%pop_real_factor)
        nstates = psip_list%nstates
    end subroutine allocate_time_average_lists

    subroutine add_ci_contribution(cluster, cdet, states, pops, nstates)

        ! Add contribution to CI wavefunction - used for variational energy estimator. 

        ! In:
        !    cdet: det_info_t object describing current location of excip.
        !    cluster: cluster_t object describing excip
        ! In/Out:
        !    states: array of determinant labels in the wavefunction.
        !    nstates: number of states in the wavefunction.
        !    pops: excip populations on determinant.

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

    subroutine add_t_contributions(psip_list, states, pops, sq, nstates)

        ! Add contributions to UCC wavefunction.

        ! In:
        !    psip_list: particle_t object storing current state of the QMC system. 

        ! In/Out:
        !    states: array of determinant labels in the wavefunction.
        !    nstates: number of states in the wavefunction.
        !    pops: excip populations.
        !    sq: squared excip populations.

        use qmc_data, only: particle_t
        use search, only: binary_search

        type(particle_t), intent(in) :: psip_list
        integer(i0), intent(inout) :: states(:,:)
        real(p), intent(inout) :: pops(:), sq(:,:)
        integer, intent(inout) :: nstates

        integer :: i, j, k, pos
        logical :: hit
        integer(i0) :: state(size(states(:,1)))

        do i = 1, psip_list%nstates
            state = psip_list%states(:,i) 
            call binary_search(states, state, 1, nstates, hit, pos)
            if (hit) then
                  pops(pos) = pops(pos) + (real(psip_list%pops(1,i))/psip_list%pop_real_factor)
                  sq(2,pos) = sq(2,pos) + (real(psip_list%pops(1,i))/psip_list%pop_real_factor)**2 
               else
                   do j = nstates, pos, -1
                       k = j + 1 
                       states(:,k) = states(:,j)
                       pops(k) = pops(j)
                       sq(1,k) = sq(1,j)
                       sq(2,k) = sq(2,j)
                   end do
                   states(:,pos) = psip_list%states(:,i)
                   pops(pos) = (real(psip_list%pops(1,i))/psip_list%pop_real_factor)
                   sq(1,pos) = psip_list%states(1,i)
                   sq(2,pos) = (real(psip_list%pops(1,i))/psip_list%pop_real_factor)**2
                   nstates = nstates + 1
               end if
        end do
    end subroutine add_t_contributions
    
    pure function earliest_unset(f, f0, nel, basis) result (early)
        
         ! Function to find earliest unset bit different to reference 
         ! Currently not in use. Reverses the order of some parameters in the trotter expansion of the 
         ! UCCMC function relative to using earliest_unset. [todo] implement switch?

         ! f0 in a determinant bit string.

         ! In:
         !    f: bit string encoding determinant.
         !    f0: bit string encoding reference determinant.
         !    nel: number of electrons in the system.
         !    basis: basis_t object with information on one-electron basis in use.
         !
         ! Returns:
         !    Zero-based index of the earliest unset bit.
     
         use basis_types, only: basis_t
         use bit_utils, only: count_set_bits
         use const, only: i0_end

         type(basis_t), intent(in) :: basis
         integer(i0), intent(in) :: f(basis%bit_string_len), f0(basis%bit_string_len)
         integer, intent(in) :: nel
         integer :: i, early, ibasis
         integer(i0) :: diff(basis%bit_string_len)
         integer(i0) :: f0_loc(basis%bit_string_len)
         
         early = 0
         i = 0
          
         diff = ieor(f, f0)
         ! If f is different from f0, find earliest different unset bit.
         if (any(diff /= 0)) then
             basis_loop: do ibasis = 1, basis%bit_string_len
                 do i = 0, i0_end
                    if (btest(diff(ibasis), i)) then
                        early = i + (i0_end + 1) * (ibasis - 1)
                        exit basis_loop
                    end if
                 end do
             end do basis_loop
         else
         ! Else find latest set bit in the reference + 1. This is not necessarily the first
         ! unset bit in f0. This will always be higher than the value for any other f.
             f0_loc = f0
             do while (any(f0_loc/=0))
                 ibasis = early/(i0_end + 1) + 1
                 i = early - (i0_end + 1)*(ibasis - 1)
                 f0_loc = ibclr(f0_loc(ibasis), i)
                 early = early + 1
             end do
         end if
    end function

    pure function latest_unset(f, f0, nel, basis) result (late)
        
         ! Function to find latest unset bit in a determinant bit string which is set in the reference f0.

         ! In:
         !    f: bit string encoding determinant
         !    f0: bit string encoding reference determinant.
         !    nel: number of electrons in the system
         !    basis: basis_t object with information on one-electron basis in use.
     
         use basis_types, only: basis_t
         use bit_utils, only: count_set_bits
         use const, only: i0_end

         type(basis_t), intent(in) :: basis
         integer(i0), intent(in) :: f(basis%bit_string_len), f0(basis%bit_string_len)
         integer, intent(in) :: nel
         integer :: late, ibasis, i
         integer(i0) :: diff(basis%bit_string_len)
         integer(i0) :: f0_loc(basis%bit_string_len)
         
         late = 0
          
         diff = ieor(f, f0)
         diff = iand(diff, f0)
         ! If f is different from f0, find latest different unset bit.
         if (any(diff /= 0)) then
             basis_loop: do ibasis = 1, basis%bit_string_len
                 do i = 0, i0_end
                    late = i + (i0_end + 1) * (ibasis - 1)
                    diff = ibclr(diff(ibasis), i)
                    if (all(diff == 0))  exit basis_loop
                 end do
             end do basis_loop
         else
         ! Else find latest set bit in the reference + 1. This is not necessarily the first
         ! unset bit in f0. This will always be higher than the value for any other f.
             f0_loc = f0
             do while (any(f0_loc/=0))
                 ibasis = late/(i0_end + 1) + 1
                 i = late - (i0_end + 1)*(ibasis - 1)
                 f0_loc = ibclr(f0_loc(ibasis), i)
                 late = late + 1
             end do
         end if
    end function

    subroutine add_info_str_trot(basis, f0, nel, f)

        ! Sets bits within info string to give excitation level and the highest energy orbital excited from for a determinant.
    
        ! In:
        !
        !    basis: basis_t object with information on one-electron basis in use.
        !    f0: bit string encoding reference determinant.
        !    nel: number of electrons in the system
        ! In/Out:
        !    f: bit string encoding determinant

        use basis_types, only: basis_t
        use errors, only: stop_all
        use excitations, only: get_excitation_level

        type(basis_t), intent(in) :: basis
        integer(i0), intent(inout) :: f(:)
        integer(i0), intent(in) :: f0(:)
        integer, intent(in) :: nel

        integer(i0) :: counter(basis%tot_string_len)
       
        if (basis%info_string_len>=0) then

            f(basis%bit_string_len+2) = latest_unset(f, f0, nel, basis)     
            f(basis%bit_string_len+1) = nel - get_excitation_level(f(:basis%bit_string_len), f0(:basis%bit_string_len)) 
        else
            call stop_all('add_info_str_trot', 'Determinant bit string does not have room for information bits.')
        end if
 
    end subroutine add_info_str_trot

end module
