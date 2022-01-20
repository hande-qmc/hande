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
                  sq(size(state) + 1,pos) = sq(size(state) + 1,pos) + (real(psip_list%pops(1,i))/psip_list%pop_real_factor)**2 
               else
                   do j = nstates, pos, -1
                       k = j + 1 
                       states(:,k) = states(:,j)
                       pops(k) = pops(j)
                       sq(:,k) = sq(:,j)
                   end do
                   states(:,pos) = psip_list%states(:,i)
                   pops(pos) = (real(psip_list%pops(1,i))/psip_list%pop_real_factor)
                   sq(:size(state),pos) = psip_list%states(:,i)
                   sq(size(state) + 1,pos) = (real(psip_list%pops(1,i))/psip_list%pop_real_factor)**2
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

    pure subroutine collapse_deexcitor_onto_cluster(basis, excitor_excitation, f0, cluster_excitor, &
                                                  cluster_annihilation, cluster_creation, cluster_population, &
                                                  excitor_annihilation, cluster_excitation)

        use basis_types, only: basis_t, reset_extra_info_bit_string

        use bit_utils, only: count_set_bits
        use const, only: i0_end

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f0(basis%tot_string_len)
        integer(i0), intent(inout) :: cluster_excitor(basis%tot_string_len)
        complex(p), intent(inout) :: cluster_population
        integer(i0), intent(in) :: excitor_excitation(basis%tot_string_len)
        integer(i0), intent(inout) :: cluster_annihilation(basis%tot_string_len)
        integer(i0), intent(inout) :: cluster_creation(basis%tot_string_len)
        integer(i0), intent(inout) :: cluster_excitation(basis%tot_string_len)
        integer(i0), intent(in) :: excitor_annihilation(basis%tot_string_len)
    

        integer :: ibasis, ibit
        integer(i0) :: permute_operators(basis%tot_string_len)
        integer(i0) :: cluster_loc(basis%tot_string_len)
        integer(i0) :: excit_swap(basis%tot_string_len)

        do ibasis = basis%bit_string_len, 1, -1
            do ibit = i0_end, 0, -1
                if (btest(excitor_excitation(ibasis),ibit)) then
                    if (.not. btest(f0(ibasis),ibit)) then
                        ! Exciting from this orbital.
                        cluster_excitor(ibasis) = ibclr(cluster_excitor(ibasis),ibit)
                        ! We must permute it with all creation operators of lower index in the cluster.
                        permute_operators = iand(not(basis%excit_mask(:,basis%basis_lookup(ibit,ibasis))), cluster_creation)
                        ! Also permutes with all annihilation operators of lower index in the EXCITOR.
                        ! e.g. in a_b a_a a_a^+ a_b^+ a_c^+, a_b must permute with a_a in the excitor and a_a^+ in the cluster
                        ! to reach the correct position to cancel out with a_b^+.
                        excit_swap = iand(not(basis%excit_mask(:,basis%basis_lookup(ibit,ibasis))), excitor_annihilation)
                        permute_operators = ieor(permute_operators, excit_swap)
                        ! Exclude swapping with itself.
                        permute_operators(ibasis) = ibclr(permute_operators(ibasis),ibit)
                    else
                        ! Exciting into this orbital.
                        cluster_excitor(ibasis) = ibset(cluster_excitor(ibasis),ibit)
                        ! Need to swap it with every remaining creation operator and 
                        ! annihilation operators with a higher index already in the cluster.
                        permute_operators = iand(basis%excit_mask(:,basis%basis_lookup(ibit,ibasis)),cluster_annihilation)
                        permute_operators = ior(permute_operators, cluster_creation)
                    end if
                    if (mod(sum(count_set_bits(permute_operators)),2) == 1) &
                        cluster_population = -cluster_population

                    cluster_loc = cluster_excitor
                    call reset_extra_info_bit_string(basis, cluster_loc)
                    ! Unlike for excitation, must update the cluster at every step to get the correct number of permutations.
                    cluster_excitation = ieor(f0, cluster_loc)
                    cluster_annihilation = iand(cluster_excitation, f0)
                    cluster_creation = iand(cluster_excitation, cluster_loc)
                end if
            end do
        end do
    end subroutine collapse_deexcitor_onto_cluster

    subroutine var_energy_uccmc(sys, states, pops, nstates, var_energy, D0_pop)

       ! Computes the variational energy of a wavefunction expressed in CI coefficients.
       !
       ! IN:
       !    sys: sys_t object encoding the system
       !    states: list of determinant labels encoded as integers
       !    pops: CI population on each determinant
       !    nstates: number of determinants in the wavefunction
       !    D0_pop: population on D0 in cluster expansion (different from pops(D0_pos) in
       !    unitary CC
       !
       ! OUT:
       !    var_energy: total variational energy estimator (NOTE: not just correlation energy)


       use excitations, only: excit_t, get_excitation
       use hamiltonian, only: get_hmatel
       use energy_evaluation, only: hmatel_t
       use system, only: sys_t
       use qmc_data, only: particle_t
       use read_in_symmetry, only: cross_product_basis_read_in
       use determinants, only: decode_det

       type(sys_t), intent(in) :: sys
       integer(i0), intent(in) :: states(:,:)
       integer, intent(in) :: nstates
       real(p), intent(in) :: pops(:), D0_pop
       real(p), intent(out) :: var_energy
       real(p) :: normalisation

       type(excit_t) :: excitation
       type(hmatel_t) :: hmatel
       integer :: occ_list(sys%nel)

       integer :: i, j
       integer :: ij_sym, ab_sym

       normalisation = 0.0_p
       var_energy = 0.0_p
       do i = 1, nstates
           normalisation = normalisation + (pops(i)/D0_pop)**2
           do j = 1, nstates
               hmatel = get_hmatel(sys, states(:,i), states(:,j))
               if (i>=j) print*, states(1,i), states(1, j), hmatel%r
               var_energy = var_energy + hmatel%r*(pops(i)/D0_pop)*(pops(j)/D0_pop)
           end do
       end do
       var_energy = var_energy/normalisation
   end subroutine var_energy_uccmc
end module
