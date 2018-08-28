module excit_gen_utils
use const
implicit none

!Routines for the power_pitzer/heat bath excit gens

contains

    subroutine select_ij_heat_bath(rng, nel, ij_weights_precalc, cdet, i, j, i_ind, j_ind, &
            ij_weights_occ, ij_weights_occ_tot, ji_weights_occ, ji_weights_occ_tot, allowed_excitation)
        ! Routine to select i and j according to the heat bath algorithm. Used by heat_bath_uniform, heat_bath_single,
        ! power_pitzer_orderM_ij

        ! In:
        !   nel: number of electrons (= sys%nel)
        !   ij_weights_precalc: precalculated weights for j given i
        !   cdet: current determinant to attempt spawning from.
        ! In/Out:
        !   rng: random number generator
        ! Out:
        !   allowed_excitation: true if excitation with ij is possible.
        !   ij_weights_occ: weights of j for all occupied spinorbitals, given i.
        !   ji_weights_occ: weights of i for all occupied spinorbitals, given j (reverse selection for pgen calculation).
        !   ij_weights_occ_tot: sum of weights of j for all occupied spinorbitals, given i.
        !   ji_weights_occ_tot: sum of weights of i for all occupied spinorbitals, given j.

        use determinant_data, only: det_info_t
        use dSFMT_interface, only: dSFMT_t
        use alias, only: select_weighted_value

        integer, intent(in) :: nel
        real(p), intent(in) :: ij_weights_precalc(:,:)
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: ij_weights_occ(:), ji_weights_occ(:)
        real(p), intent(out) :: ij_weights_occ_tot, ji_weights_occ_tot
        logical, intent(out) :: allowed_excitation
        
        integer :: pos_occ, i_ind, j_ind, i, j

        i_ind = select_weighted_value(rng, nel, cdet%i_d_occ%weights, cdet%i_d_occ%weights_tot)
        i = cdet%occ_list(i_ind)

        ij_weights_occ_tot = 0.0_p
        do pos_occ = 1, nel
            ij_weights_occ(pos_occ) = ij_weights_precalc(cdet%occ_list(pos_occ),i)
            ij_weights_occ_tot = ij_weights_occ_tot + ij_weights_occ(pos_occ)
        end do

        if (ij_weights_occ_tot > 0.0_p) then
            ! There is a j for this i.
            j_ind = select_weighted_value(rng, nel, ij_weights_occ, ij_weights_occ_tot)
            j = cdet%occ_list(j_ind)

            ! Pre-compute the other direction (first selecting j then i) as well as that is required for pgen.
            ji_weights_occ_tot = 0.0_p
            do pos_occ = 1, nel
                ji_weights_occ(pos_occ) = ij_weights_precalc(cdet%occ_list(pos_occ),j)
                ji_weights_occ_tot = ji_weights_occ_tot + ji_weights_occ(pos_occ)
            end do
            allowed_excitation = .true.
        else
            allowed_excitation = .false.
        end if

    end subroutine select_ij_heat_bath

    subroutine init_double_weights_ab(sys, i, j, weight)
        ! WARNING: this routine assumes that i /= j!
        ! Routine that helps set up weights in a heat bath manner for the part where it loops over a and b.

        ! In:
        !   sys: system information
        !   i,j: orbitals i and j. i/=j.
        ! In/Out:
        !   weight: Hijab summed over a and b in this function. Can be an ongoing sum over i and j as well.

        use system, only: sys_t
        use read_in_symmetry, only: cross_product_basis_read_in
        use proc_pointers, only: slater_condon2_excit_ptr, abs_hmatel_ptr
        use hamiltonian_data, only: hmatel_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j
        real(p), intent(inout) :: weight

        integer :: ij_sym, isymb, i_tmp, j_tmp, a, b, a_tmp, b_tmp
        type(hmatel_t) :: hmatel

        if (j < i) then
            i_tmp = j
            j_tmp = i
        else
            i_tmp = i
            j_tmp = j
        end if 
                    
        ! The symmetry of b (=b_cdet), isymb, is given by
        ! (sym_i_cdet* x sym_j_cdet* x sym_a_cdet)* = sym_b_cdet
        ! (at least for Abelian point groups)
        ! ij_sym: symmetry conjugate of the irreducible representation spanned by the codensity
        !        \phi_i_cdet*\phi_j_cdet. (We assume that ij is going to be in the bra of the excitation.)
        ! [todo] - Check whether order of i and j matters here.

        ij_sym = sys%read_in%sym_conj_ptr(sys%read_in, cross_product_basis_read_in(sys, i_tmp, j_tmp))
        
        !$omp parallel do default(none) &
        !$omp shared(sys,i,j,i_tmp,j_tmp,ij_sym,slater_condon2_excit_ptr,abs_hmatel_ptr) &
        !$omp private(a,b,a_tmp,b_tmp,isymb,hmatel) reduction(+:weight)
        do a = 1, sys%basis%nbasis
            if ((a /= i_tmp) .and. (a /= j_tmp)) then
                isymb = sys%read_in%sym_conj_ptr(sys%read_in, &
                    sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, sys%basis%basis_fns(a)%sym))
                do b = 1, sys%basis%nbasis 
                    if ((((sys%basis%basis_fns(i_tmp)%Ms == sys%basis%basis_fns(a)%Ms) .and. &
                        (sys%basis%basis_fns(j_tmp)%Ms == sys%basis%basis_fns(b)%Ms)) .or. &
                        ((sys%basis%basis_fns(i_tmp)%Ms == sys%basis%basis_fns(b)%Ms) .and. &
                        (sys%basis%basis_fns(j_tmp)%Ms == sys%basis%basis_fns(a)%Ms))) .and. &
                        (sys%basis%basis_fns(b)%sym == isymb) .and. (a /= b) .and. (b /= i_tmp) .and. &
                        (b /= j_tmp)) then
                        if (b < a) then
                            a_tmp = b
                            b_tmp = a
                        else
                            a_tmp = a
                            b_tmp = b
                        end if
                        ! slater_condon2 does not check whether ij -> ab is allowed by symmetry/spin
                        ! but we have checked for that here so it is ok.
                        hmatel = slater_condon2_excit_ptr(sys, i_tmp, j_tmp, a_tmp, b_tmp, .false.)
                        weight = weight + abs_hmatel_ptr(hmatel)
                    end if
                end do
            end if
        end do
        !$omp end parallel do

    end subroutine init_double_weights_ab

    !--- Helper Functions for calculating weights in decoder ---

    pure subroutine find_i_d_weights(nel, i_weights_precalc, d)

        ! Find weights to select i in a double excitation using pre-calculed heat bath weights.

        ! In:
        !   nel: number of electrons
        !   i_weights_precalc: precalculated weights that i could have. Need to reduce them to the ones
        !        that are actually occupied.
        ! In/Out:
        !   d: det_code_t object. Information about the determinant to decode.
        
        use determinant_data, only: det_info_t
        
        integer, intent(in) :: nel
        real(p), intent(in) :: i_weights_precalc(:)
        type(det_info_t), intent(inout) :: d

        integer :: pos_occ

        d%i_d_occ%weights_tot = 0.0_p
        do pos_occ = 1, nel
            d%i_d_occ%weights(pos_occ) = i_weights_precalc(d%occ_list(pos_occ))
            d%i_d_occ%weights_tot = d%i_d_occ%weights_tot + d%i_d_occ%weights(pos_occ)
        end do

    end subroutine find_i_d_weights
    
    pure subroutine find_ia_single_weights(sys, d)

        ! Create a random single excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! This calculates the weights for i and a as exact as possibly.
        ! For i, the weight is \sum_a H_ia, for a given i (a|i) it is H_ia.

        ! In:
        !    sys: system object being studied.
        ! In/Out:
        !    d: information about current determinant to decode

        use determinant_data, only: det_info_t
        use proc_pointers, only: abs_hmatel_ptr
        use system, only: sys_t
        use hamiltonian_data, only: hmatel_t
        use hamiltonian_molecular, only: slater_condon1_mol
        use hamiltonian_periodic_complex, only: slater_condon1_periodic_complex

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(inout) :: d
        
        type(hmatel_t) :: hmatel
        integer :: i_ind, a_ind
        
        d%i_s_occ%weights_tot = 0.0_p
        do i_ind = 1, sys%nel
            d%i_s_occ%weights(i_ind) = 0.0_p
            do a_ind = 1, sys%nvirt
                ! [todo] - this reduces the computational time (by calculation ia_weights here)
                ! but more memory costs as after we have selected i, we don't need all elements in ia_weights. Need to balance
                ! these costs.
                ! [todo] - could consider rewritting with proc_pointer here but that would imply having to rewrite slater
                ! condon 1 mol / periodic complex. slater_condon1_excit_mol is not safe enough (no checks).
                if (sys%read_in%comp) then
                    hmatel%r = 0.0_p
                    hmatel%c = slater_condon1_periodic_complex(sys, d%occ_list, d%occ_list(i_ind), d%unocc_list(a_ind), &
                        .false.)
                else
                    hmatel%r = slater_condon1_mol(sys, d%occ_list, d%occ_list(i_ind), d%unocc_list(a_ind), .false.)
                    hmatel%c = cmplx(0.0_p, 0.0_p, p)
                end if
                d%ia_s_weights_occ(a_ind, i_ind) = abs_hmatel_ptr(hmatel)
                d%i_s_occ%weights(i_ind) = d%i_s_occ%weights(i_ind) + d%ia_s_weights_occ(a_ind, i_ind)
            end do
            d%i_s_occ%weights_tot = d%i_s_occ%weights_tot + d%i_s_occ%weights(i_ind)
        end do

    end subroutine find_ia_single_weights
    
    pure subroutine find_diff_ref_cdet(sys, d, pp)

        ! Pre-calculate the reordered list of occupied list that is equal to the
        ! occupied list of orbitals of the reference where orbitals that differ between
        ! the current determinant and the reference are substituted by occupied orbitals
        ! from the current determinant. These substituted orbitals are sorted by spin.
        !
        ! In:
        !    sys: system being studied (contains required basis information).
        !    pp: information for pp excitation generators
        ! In/Out:
        !    d: det_info_t variable.  The following components are set here:
        !        ref_cdet_occ_list: reordered list of occ. spin orbitals.

        use system, only: sys_t
        use excit_gens, only: excit_gen_power_pitzer_t
        use excitations, only: get_excitation_locations
        use determinant_data, only: det_info_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(inout) :: d
        type(excit_gen_power_pitzer_t), intent(in) :: pp

        integer :: ref_store(sys%nel), det_store(sys%nel)
        integer :: ii, jj, nex, t

        call get_excitation_locations(pp%occ_list, d%occ_list, ref_store, det_store, sys%nel, nex)
        ! These orbitals might not be aligned in the most efficient way:
        !  They may not match in spin, so first deal with this

        d%ref_cdet_occ_list = pp%occ_list(:sys%nel)
        ! ref store (e.g.) contains the indices within excit_gen_pp%occ_list of the orbitals
        ! which have been excited from.
        do ii=1, nex
            associate(bfns=>sys%basis%basis_fns)
                if (bfns(pp%occ_list(ref_store(ii)))%Ms /= bfns(d%occ_list(det_store(ii)))%Ms) then
                    jj = ii + 1
                    do while (bfns(pp%occ_list(ref_store(ii)))%Ms /= bfns(d%occ_list(det_store(jj)))%Ms)
                        jj = jj + 1
                    end do
                    ! det's jj now points to an orb of the same spin as ref's ii, so swap cdet_store's ii and jj.
                    t = det_store(ii)
                    det_store(ii) = det_store(jj)
                    det_store(jj) = t
                end if
            end associate
            d%ref_cdet_occ_list(ref_store(ii)) = d%occ_list(det_store(ii))
        end do

    end subroutine find_diff_ref_cdet

end module excit_gen_utils
