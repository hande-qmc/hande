module excit_gen_heat_bath_mol

! A module containing excitations generators for molecules which weight excitations according to the heat bath algorithm
! by Holmes et al. (Holmes, A. A.; Changlani, H. J.; Umrigar, C. J. J. Chem. Theory Comput. 2016, 12, 1561–1571).

use const, only: i0, p, depsilon

implicit none

contains

    subroutine init_excit_mol_heat_bath(sys, hb, original)

        ! Generate excitation tables for all spinorbitals for the gen_excit_mol_heat_bath
        ! excitation generator.

        ! In:
        !    sys: system object being studied.
        !    original: boolean, true if original heat bath, false if single excitations are sampled uniformly/exactly.
        ! In/Out:
        !    hb: an empty excit_gen_heat_bath_t object which gets filled with
        !           the alias tables required to generate excitations.

        use system, only: sys_t
        use sort, only: qsort
        use proc_pointers, only: create_weighted_excitation_list_ptr, slater_condon2_excit_ptr
        use proc_pointers, only: abs_hmatel_ptr
        use excit_gens, only: excit_gen_heat_bath_t
        use alias, only: generate_alias_tables
        use read_in_symmetry, only: cross_product_basis_read_in
        use hamiltonian_data, only: hmatel_t
        use errors, only: stop_all, warning
        type(sys_t), intent(in) :: sys
        type(excit_gen_heat_bath_t), intent(inout) :: hb
        logical, intent(in) :: original
        type(hmatel_t) :: hmatel

        integer :: i, j, a, b, ij_sym, isymb, ims, isyma
        integer :: i_tmp, j_tmp, a_tmp, b_tmp
        real(p) :: i_weight, ij_weight, ija_weight, ijab_weight
        
        integer :: j_nonzero(sys%basis%nbasis,sys%basis%nbasis) 
        ! Temp storage
        allocate(hb%i_weights(sys%basis%nbasis)) ! This will hold S_p in the Holmes JCTC paper.
        allocate(hb%ij_weights(sys%basis%nbasis,sys%basis%nbasis)) ! This stores D_pq.
        allocate(hb%ija_weights(sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis))
        allocate(hb%ija_aliasU(sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis))
        allocate(hb%ija_aliasK(sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis))
        allocate(hb%ija_weights_tot(sys%basis%nbasis,sys%basis%nbasis))
        allocate(hb%ijab_weights(sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis))
        allocate(hb%ijab_aliasU(sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis))
        allocate(hb%ijab_aliasK(sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis))
        allocate(hb%ijab_weights_tot(sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis))
        allocate(hb%ia_weights(sys%basis%nbasis,sys%basis%nbasis)) ! implicit sum over j

        ! [todo] - consider setting hmatel%r and hmatel%c to zero to avoid undefined behaviour.
        ! [todo] - although it is probably fine, abs_hmatel_ptr ignores the irrelevant bit.

        if (original) then
            ! check whether we are allowed to use the original heat bath algorithm. Exit if not.
            call check_heat_bath_bias(sys)
        end if

        j_nonzero = 0

        hb%ia_weights = 0.0_p

        do i = 1, sys%basis%nbasis
            i_weight = 0.0_p
            do j = 1, sys%basis%nbasis
                ij_weight = 0.0_p
                hb%ija_weights_tot(j,i) = 0.0_p
                if (i /= j) then
                    if (j < i) then
                        i_tmp = j
                        j_tmp = i
                    else
                        i_tmp = i
                        j_tmp = j
                    end if
                    ! The symmetry of b, isymb, is given by
                    ! (sym_i* x sym_j* x sym_a)* = sym_b
                    ! (at least for Abelian point groups)
                    ! ij_sym: symmetry conjugate of the irreducible representation spanned by the codensity
                    !        \phi_i*\phi_j. (We assume that ij is going to be in the bra of the excitation.)
                    ! [todo] - Check whether order of i and j matters here.
                    ij_sym = sys%read_in%sym_conj_ptr(sys%read_in, cross_product_basis_read_in(sys, i_tmp, j_tmp))
                    do a = 1, sys%basis%nbasis
                        ija_weight = 0.0_p
                        hb%ijab_weights_tot(a,j,i) = 0.0_p
                        if ((a /= i) .and. (a /= j)) then
                            isymb = sys%read_in%sym_conj_ptr(sys%read_in, &
                                        sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, sys%basis%basis_fns(a)%sym))
                            do b = 1, sys%basis%nbasis
                                ijab_weight = 0.0_p
                                ! Check spin conservation and symmetry conservation.
                                if ((((sys%basis%basis_fns(i_tmp)%Ms == sys%basis%basis_fns(a)%Ms) .and. &
                                    (sys%basis%basis_fns(j_tmp)%Ms == sys%basis%basis_fns(b)%Ms)) .or. &
                                    ((sys%basis%basis_fns(i_tmp)%Ms == sys%basis%basis_fns(b)%Ms) .and. &
                                    (sys%basis%basis_fns(j_tmp)%Ms == sys%basis%basis_fns(a)%Ms))) .and. &
                                    (sys%basis%basis_fns(b)%sym == isymb) .and. ((b /= a) .and. (b /= i) .and. &
                                    (b /= j))) then
                                    if (b < a) then
                                        a_tmp = b
                                        b_tmp = a
                                    else
                                        a_tmp = a
                                        b_tmp = b
                                    end if
                                    hmatel = slater_condon2_excit_ptr(sys, i_tmp, j_tmp, a_tmp, b_tmp, .false.)
                                    i_weight = i_weight + abs_hmatel_ptr(hmatel)
                                    ij_weight = ij_weight + abs_hmatel_ptr(hmatel)
                                    ija_weight = ija_weight + abs_hmatel_ptr(hmatel)
                                    ijab_weight = ijab_weight + abs_hmatel_ptr(hmatel)
                                end if
                                hb%ijab_weights(b,a,j,i) = ijab_weight
                                hb%ijab_weights_tot(a,j,i) = hb%ijab_weights_tot(a,j,i) + ijab_weight
                            end do
                        else
                            do b = 1, sys%basis%nbasis
                                hb%ijab_weights(b,a,j,i) = 0.0_p
                            end do
                        end if
                        hb%ija_weights(a,j,i) = ija_weight
                        hb%ija_weights_tot(j,i) = hb%ija_weights_tot(j,i) + ija_weight
                        hb%ia_weights(a,i) = hb%ia_weights(a,i) + ija_weight
                        if (ija_weight > depsilon) then
                            j_nonzero(a,i) = j_nonzero(a,i) + 1
                        end if
                    end do
                else
                    do a = 1, sys%basis%nbasis
                        hb%ija_weights(a,j,i) = 0.0_p
                        hb%ijab_weights_tot(a,j,i) = 0.0_p
                        do b = 1, sys%basis%nbasis
                            hb%ijab_weights(b,a,j,i) = 0.0_p
                        end do
                    end do
                end if
                hb%ij_weights(j,i) = ij_weight
            end do
            hb%i_weights(i) = i_weight

            ! Old test for bias. Replaced by subroutine check_heat_bath_bias
            ! Test that all single excitation i -> a that are allowed have some non zero weight ija for some j.
            !do a = 1, sys%basis%nbasis
            !    if ((j_nonzero(a,i) < (sys%basis%nbasis - sys%nel)) .and. (original == .true.) .and. (i /= a)) then
            !        ! no non zero weight ija for some j and (i /= a). Now check that i -> a is valid.
            !        ims = sys%basis%basis_fns(i)%ms
            !        isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(i)%sym, &
            !                                        sys%read_in%pg_sym%gamma_sym)
            !        if ((sys%basis%basis_fns(a)%sym == isyma) .and. (sys%basis%basis_fns(a)%ms == ims)) then
            !            call stop_all('init_excit_mol_heat_bath','Not all possible single excitations can be accounted for. &
            !                Use another excitation generator.')
            !        end if
            !    end if
            !end do
        end do

        do i = 1, sys%basis%nbasis
            do j = 1, sys%basis%nbasis
                if (abs(hb%ija_weights_tot(j,i)) > 0.0_p) then
                    call generate_alias_tables(sys%basis%nbasis, hb%ija_weights(:,j,i), hb%ija_weights_tot(j,i), &
                                                hb%ija_aliasU(:,j,i), hb%ija_aliasK(:,j,i))
                    do a = 1, sys%basis%nbasis
                        if (abs(hb%ijab_weights_tot(a,j,i)) > 0.0_p) then
                            call generate_alias_tables(sys%basis%nbasis, hb%ijab_weights(:,a,j,i), hb%ijab_weights_tot(a,j,i), &
                                                        hb%ijab_aliasU(:,a,j,i), hb%ijab_aliasK(:,a,j,i))
                        end if
                    end do
                end if
            end do
        end do

    end subroutine init_excit_mol_heat_bath

    subroutine check_heat_bath_bias(sys)
        
        ! This routine checks for bias in the heat_bath excitation generator where a single excitation
        ! i -> a is only chosen after i, j and a have been selected. This means that if there is no
        ! j that allows us to select i, j and a but i -> a is valid, i -> a cannot be selected even
        ! though it should have that chance. If a bias is found, stop the calculation run.

        ! The bias follows the check by Holmes et al. (Section D, p. 1565)
        ! (Holmes, A. A.; Changlani, H. J.; Umrigar, C. J. J. Chem. Theory Comput. 2016, 12, 1561–1571)
        ! The number of irreducible representations in the full symmetry group are found by the
        ! number of i for which Hij == 0 for any j /= i for the system we study except that we only
        ! have one electron instead of N electrons in M spinorbitals.
        ! If the number of irreducible representations is bigger than either the number of up spins or
        ! down spins (whatever is bigger), then the heat bath algorithm might be biased. If not, then
        ! it is not biased (a sufficient but not necessary condition).

        ! In:
        !   sys: information about the system to be studied.

        use errors, only: stop_all
        use system, only: sys_t
        use molecular_integrals, only: get_one_body_int_mol_real, get_one_body_int_mol_complex
        
        type(sys_t), intent(in) :: sys

        integer :: i, j, counter
        complex(p) :: hij_comp
        real(p) :: hij_real

        counter = 0 ! number of irreps.
        hij_comp = cmplx(0.0_p, 0.0_p, p)
        hij_real = 0.0_p

        do i = 1, sys%basis%nbasis
            do j = 1, sys%basis%nbasis
                if (i /= j) then
                    if (sys%read_in%comp) then
                        ! add the abs value so that minus/plus cancellations do not happen.
                        hij_comp = hij_comp + abs(get_one_body_int_mol_complex(sys%read_in%one_e_h_integrals, &
                                                sys%read_in%one_e_h_integrals_imag, i, j, sys))
                    else
                        hij_real = hij_real + abs(get_one_body_int_mol_real(sys%read_in%one_e_h_integrals, i, j, sys))
                    end if
                end if
            end do
            ! [todo] compare as > 0 or > depsilon?
            if (((sys%read_in%comp) .and. (abs(real(hij_comp)) > 0.0_p) .and. (abs(aimag(hij_comp)) > 0.0_p)) .or. &
                ((sys%read_in%comp == .false.) .and. (abs(hij_real) > 0.0_p))) then
                counter = counter + 1
            end if
        end do

        if ((counter < max(sys%nalpha, sys%nbeta)) == .false.) then
            ! [todo] - should only the parent process deal with this?
            call stop_all('check_heat_bath_bias','Maybe not all possible single excitations can be accounted for. &
                            Use another excitation generator.')
        end if

    end subroutine check_heat_bath_bias

    subroutine gen_excit_mol_heat_bath(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! This is the main heat bath algorithm.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Excitation generation data.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant in molecular systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2, get_excitation_locations
        use proc_pointers, only: slater_condon1_excit_ptr, slater_condon2_excit_ptr, abs_hmatel_ptr
        use proc_pointers, only: create_weighted_excitation_list_ptr
        use checking, only: check_allocate, check_deallocate
        use system, only: sys_t
        use excit_gen_mol, only: gen_single_excit_mol_no_renorm
        use excit_gens, only: excit_gen_heat_bath_t, excit_gen_data_t
        use alias, only: select_weighted_value_precalc, select_weighted_value
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use hamiltonian_data, only: hmatel_t
        use read_in_symmetry, only: cross_product_basis_read_in
        use search, only: binary_search
        use sort, only: qsort

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation
        real(p), allocatable :: i_weights_occ(:)
        real(p), allocatable :: ij_weights_occ(:)
        real(p), allocatable :: ji_weights_occ(:)
        integer :: occ_list(sys%nel)

        integer :: pos_bas, pos_occ, pos_q, i_ind, j_ind, i, j, a, b, ierr, ims, isyma
        
        real(p) :: i_weights_occ_tot, ij_weights_occ_tot, ji_weights_occ_tot, pgen_ij, hmatel_mod_ia, x, psingle
        real(p) :: pgen_ija, pgen_ijb, pgen_jia, pgen_jib, psingle_q
        logical :: double
        type(hmatel_t) :: hmatel_single

        occ_list = cdet%occ_list
        ! call qsort(occ_list,sys%nel)
        associate( hb => excit_gen_data%excit_gen_hb )
            ! 2b. Select orbitals to excite from
            allocate(i_weights_occ(1:sys%nel), stat=ierr)
            call check_allocate('i_weights_occ', sys%nel, ierr)

            ! [todo] can move iido to double.
            i_weights_occ_tot = 0.0_p
            do pos_occ = 1, sys%nel
                i_weights_occ(pos_occ) = hb%i_weights(cdet%occ_list(pos_occ))
                i_weights_occ_tot = i_weights_occ_tot + i_weights_occ(pos_occ)
            end do
                
            i_ind = select_weighted_value(rng, sys%nel, i_weights_occ, i_weights_occ_tot)
            i = occ_list(i_ind)

            allocate(ij_weights_occ(1:sys%nel), stat=ierr)
            call check_allocate('ij_weights_occ', sys%nel, ierr)

            ij_weights_occ_tot = 0.0_p
            do pos_occ = 1, sys%nel
                ij_weights_occ(pos_occ) = hb%ij_weights(cdet%occ_list(pos_occ),i)
                ij_weights_occ_tot = ij_weights_occ_tot + ij_weights_occ(pos_occ)
            end do

            j_ind = select_weighted_value(rng, sys%nel, ij_weights_occ, ij_weights_occ_tot)
            j = occ_list(j_ind)

            allocate(ji_weights_occ(1:sys%nel), stat=ierr)
            call check_allocate('ji_weights_occ', sys%nel, ierr)
                
            ji_weights_occ_tot = 0.0_p
            do pos_occ = 1, sys%nel
                ji_weights_occ(pos_occ) = hb%ij_weights(cdet%occ_list(pos_occ),j)
                ji_weights_occ_tot = ji_weights_occ_tot + ji_weights_occ(pos_occ)
            end do

            allowed_excitation = .true.
            if (abs(hb%ija_weights_tot(j, i)) > 0.0_p) then
                a = select_weighted_value_precalc(rng, sys%basis%nbasis, hb%ija_aliasU(:, j, i), hb%ija_aliasK(:, j, i))
                if (.not.btest(cdet%f(sys%basis%bit_lookup(2,a)), sys%basis%bit_lookup(1,a))) then
                    ims = sys%basis%basis_fns(i)%ms
                    isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(i)%sym, &
                                                    sys%read_in%pg_sym%gamma_sym)    
                    if ((sys%basis%basis_fns(a)%sym == isyma) .and. (sys%basis%basis_fns(a)%ms == ims)) then
                        connection%from_orb(1) = i
                        connection%to_orb(1) = a
                        connection%nexcit = 1
                        call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)
                        
                        hmatel_single = slater_condon1_excit_ptr(sys, cdet%occ_list, i, a, connection%perm)
                        hmatel_mod_ia = abs_hmatel_ptr(hmatel_single)
                        x = get_rand_close_open(rng)                        
                        if (hmatel_mod_ia < hb%ijab_weights_tot(a,j,i)) then
                            psingle = hmatel_mod_ia/(hb%ijab_weights_tot(a,j,i) + hmatel_mod_ia)
                        else
                            psingle = 0.5_p
                        end if
                        if (x < psingle) then
                            double = .false.
                        else
                            double = .true.
                        end if
                    else
                        ! i -> a not possible without j and b.
                        double = .true.
                        psingle = 0.0_p
                    end if
                else
                    allowed_excitation = .false.
                end if
            else
                allowed_excitation = .false.
            end if

            if (allowed_excitation) then
                if (double) then
                    b = select_weighted_value_precalc(rng, sys%basis%nbasis, hb%ijab_aliasU(:, a, j, i), &
                                                        hb%ijab_aliasK(:, a, j, i))
                    if (.not.btest(cdet%f(sys%basis%bit_lookup(2,b)), sys%basis%bit_lookup(1,b))) then
                        allowed_excitation = .true.
                        
                        pgen_ija = ((i_weights_occ(i_ind)/i_weights_occ_tot) * (ij_weights_occ(j_ind)/ij_weights_occ_tot)) * &
                            (hb%ija_weights(a,j,i)/hb%ija_weights_tot(j,i)) * (1.0_p - psingle) * &
                            (hb%ijab_weights(b,a,j,i)/hb%ijab_weights_tot(a,j,i))

                        ims = sys%basis%basis_fns(i)%ms
                        isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(i)%sym, &
                                                                            sys%read_in%pg_sym%gamma_sym)

                        connection%from_orb(1) = i
                        connection%to_orb(1) = b
                        connection%nexcit = 1

                        if ((sys%basis%basis_fns(b)%sym == isyma) .and. (sys%basis%basis_fns(b)%ms == ims)) then
                            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)
                            hmatel_single = slater_condon1_excit_ptr(sys, cdet%occ_list, i, b, connection%perm)
                            hmatel_mod_ia = abs_hmatel_ptr(hmatel_single)
                            if (hmatel_mod_ia < hb%ijab_weights_tot(b,j,i)) then
                                psingle = hmatel_mod_ia/(hb%ijab_weights_tot(b,j,i) + hmatel_mod_ia)
                            else
                                psingle = 0.5_p
                            end if
                        else
                            psingle = 0.0_p
                        end if

                        pgen_ijb = ((i_weights_occ(i_ind)/i_weights_occ_tot) * (ij_weights_occ(j_ind)/ij_weights_occ_tot)) * &
                            (hb%ija_weights(b,j,i)/hb%ija_weights_tot(j,i)) * (1.0_p - psingle) * &
                            (hb%ijab_weights(a,b,j,i)/hb%ijab_weights_tot(b,j,i))

                        ims = sys%basis%basis_fns(j)%ms
                        isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(j)%sym, &
                                                    sys%read_in%pg_sym%gamma_sym)

                        connection%from_orb(1) = j
                        connection%to_orb(1) = a
                        connection%nexcit = 1

                        if ((sys%basis%basis_fns(a)%sym == isyma) .and. (sys%basis%basis_fns(a)%ms == ims)) then
                            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)
                            hmatel_single = slater_condon1_excit_ptr(sys, cdet%occ_list, j, a, connection%perm)
                            hmatel_mod_ia = abs_hmatel_ptr(hmatel_single)
                            if (hmatel_mod_ia < hb%ijab_weights_tot(a,i,j)) then
                                psingle = hmatel_mod_ia/(hb%ijab_weights_tot(a,i,j) + hmatel_mod_ia)
                            else
                                psingle = 0.5_p
                            end if
                        else
                            psingle = 0.0_p
                        end if

                        pgen_jia = ((i_weights_occ(j_ind)/i_weights_occ_tot) * (ji_weights_occ(i_ind)/ji_weights_occ_tot)) * &
                            (hb%ija_weights(a,i,j)/hb%ija_weights_tot(i,j)) * (1.0_p - psingle) * &
                            (hb%ijab_weights(b,a,i,j)/hb%ijab_weights_tot(a,i,j))

                        connection%from_orb(1) = j
                        connection%to_orb(1) = b
                        connection%nexcit = 1
                        
                        if ((sys%basis%basis_fns(b)%sym == isyma) .and. (sys%basis%basis_fns(b)%ms == ims)) then
                            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)
                            hmatel_single = slater_condon1_excit_ptr(sys, cdet%occ_list, j, b, connection%perm)
                            hmatel_mod_ia = abs_hmatel_ptr(hmatel_single)
                            if (hmatel_mod_ia < hb%ijab_weights_tot(b,i,j)) then
                                psingle = hmatel_mod_ia/(hb%ijab_weights_tot(b,i,j) + hmatel_mod_ia)
                            else
                                psingle = 0.5_p
                            end if
                        else
                            psingle = 0.0_p
                        end if

                        pgen_jib = ((i_weights_occ(j_ind)/i_weights_occ_tot) * (ji_weights_occ(i_ind)/ji_weights_occ_tot)) * &
                            (hb%ija_weights(b,i,j)/hb%ija_weights_tot(i,j)) * (1.0_p - psingle) * &
                            (hb%ijab_weights(a,b,i,j)/hb%ijab_weights_tot(b,i,j))

                        pgen = pgen_ija + pgen_ijb + pgen_jia + pgen_jib

                        if (i < j) then
                            connection%from_orb(1) = i
                            connection%from_orb(2) = j
                        else
                            connection%from_orb(1) = j
                            connection%from_orb(2) = i
                        end if
                        if (a < b) then
                            connection%to_orb(1) = a
                            connection%to_orb(2) = b
                        else
                            connection%to_orb(2) = a
                            connection%to_orb(1) = b
                        end if
                        connection%nexcit = 2

                        call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)
                        hmatel = slater_condon2_excit_ptr(sys, connection%from_orb(1), connection%from_orb(2), &
                                                    connection%to_orb(1), connection%to_orb(2), connection%perm)
                    else
                        allowed_excitation = .false.
                        ! We have not found a valid excitation.
                        ! Return a null excitation.
                        hmatel%c = cmplx(0.0_p, 0.0_p, p)
                        hmatel%r = 0.0_p
                        pgen = 1.0_p
                    end if
                else
                    connection%from_orb(1) = i
                    connection%to_orb(1) = a
                    connection%nexcit = 1
                    hmatel = slater_condon1_excit_ptr(sys, cdet%occ_list, i, a, connection%perm)
                    pgen = 0.0_p
                    do pos_q = 1, sys%nel
                        if ((i /= cdet%occ_list(pos_q)) .and. (a /= cdet%occ_list(pos_q))) then
                            if (hmatel_mod_ia < hb%ijab_weights_tot(a,cdet%occ_list(pos_q),i)) then
                                psingle_q = hmatel_mod_ia/(hb%ijab_weights_tot(a,cdet%occ_list(pos_q),i) + hmatel_mod_ia)
                            else
                                psingle_q = 0.5_p
                            end if
                            pgen = pgen + (psingle_q * (ij_weights_occ(pos_q)/ij_weights_occ_tot) * &
                                (hb%ija_weights(a, cdet%occ_list(pos_q), i)/hb%ija_weights_tot(cdet%occ_list(pos_q), i)))
                        end if
                    end do
                    pgen = pgen * (i_weights_occ(i_ind)/i_weights_occ_tot)
                end if
            else
                ! We have not found a valid excitation.
                ! Return a null excitation.
                hmatel%c = cmplx(0.0_p, 0.0_p, p)
                hmatel%r = 0.0_p
                pgen = 1.0_p
            end if

            ! deallocate weight arrays if allocated
            if (allocated(i_weights_occ)) deallocate(i_weights_occ)
            if (allocated(ij_weights_occ)) deallocate(ij_weights_occ)
            if (allocated(ji_weights_occ)) deallocate(ji_weights_occ)
        end associate

    end subroutine gen_excit_mol_heat_bath
    
    subroutine gen_excit_mol_heat_bath_uniform(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! This is the heat bath algorithm that selects single excitations uniformly.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Excitation generation data.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant in molecular systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2, get_excitation_locations
        use proc_pointers, only: slater_condon2_excit_ptr
        use proc_pointers, only: create_weighted_excitation_list_ptr
        use checking, only: check_allocate, check_deallocate
        use system, only: sys_t
        use excit_gen_mol, only: gen_single_excit_mol
        use excit_gens, only: excit_gen_heat_bath_t, excit_gen_data_t
        use alias, only: select_weighted_value_precalc, select_weighted_value
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use hamiltonian_data, only: hmatel_t
        use qmc_data, only: excit_gen_heat_bath_uniform, excit_gen_heat_bath_single
        use read_in_symmetry, only: cross_product_basis_read_in
        use search, only: binary_search
        use sort, only: qsort

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation
        real(p), allocatable :: i_weights_occ(:)
        real(p), allocatable :: ij_weights_occ(:)
        real(p), allocatable :: ji_weights_occ(:)
        integer :: occ_list(sys%nel)

        integer :: pos_bas, pos_occ, i_ind, j_ind, i, j, a, b, ierr
        
        real(p) :: i_weights_occ_tot, ij_weights_occ_tot, ji_weights_occ_tot

        ! 1. Select single or double.
        if (get_rand_close_open(rng) < excit_gen_data%pattempt_single) then  
            ! We have a single
            select case(excit_gen_data%excit_gen)
            case (excit_gen_heat_bath_uniform)
                call gen_single_excit_mol(rng, sys, excit_gen_data%pattempt_single, cdet, pgen, connection, hmatel, &
                                                allowed_excitation)
            case (excit_gen_heat_bath_single)
                call gen_single_excit_heat_bath_exact(rng, sys, excit_gen_data%pattempt_single, cdet, pgen, connection, hmatel,&
                                                    allowed_excitation)
            case default
                ! [todo] debug. remove.
                print *, "something went wrong"
            end select
        else
            occ_list = cdet%occ_list
            ! call qsort(occ_list,sys%nel)
            ! We have a double
            associate( hb => excit_gen_data%excit_gen_hb )
                ! 2b. Select orbitals to excite from
                allocate(i_weights_occ(1:sys%nel), stat=ierr)
                call check_allocate('i_weights_occ', sys%nel, ierr)
                
                i_weights_occ_tot = 0.0_p
                do pos_occ = 1, sys%nel
                    i_weights_occ(pos_occ) = hb%i_weights(cdet%occ_list(pos_occ))
                    i_weights_occ_tot = i_weights_occ_tot + i_weights_occ(pos_occ)
                end do
                
                i_ind = select_weighted_value(rng, sys%nel, i_weights_occ, i_weights_occ_tot)
                i = occ_list(i_ind)

                allocate(ij_weights_occ(1:sys%nel), stat=ierr)
                call check_allocate('ij_weights_occ', sys%nel, ierr)

                ij_weights_occ_tot = 0.0_p
                do pos_occ = 1, sys%nel
                    ij_weights_occ(pos_occ) = hb%ij_weights(cdet%occ_list(pos_occ),i)
                    ij_weights_occ_tot = ij_weights_occ_tot + ij_weights_occ(pos_occ)
                end do

                j_ind = select_weighted_value(rng, sys%nel, ij_weights_occ, ij_weights_occ_tot)
                j = occ_list(j_ind)

                allocate(ji_weights_occ(1:sys%nel), stat=ierr)
                call check_allocate('ji_weights_occ', sys%nel, ierr)
                
                ji_weights_occ_tot = 0.0_p
                do pos_occ = 1, sys%nel
                    ji_weights_occ(pos_occ) = hb%ij_weights(cdet%occ_list(pos_occ),j)
                    ji_weights_occ_tot = ji_weights_occ_tot + ji_weights_occ(pos_occ)
                end do

                pgen = ((i_weights_occ(i_ind)/i_weights_occ_tot) * (ij_weights_occ(j_ind)/ij_weights_occ_tot)) + &
                        ((i_weights_occ(j_ind)/i_weights_occ_tot) * (ji_weights_occ(i_ind)/ji_weights_occ_tot))

                if (abs(hb%ija_weights_tot(j, i)) > 0.0_p) then
                    a = select_weighted_value_precalc(rng, sys%basis%nbasis, hb%ija_aliasU(:, j, i), hb%ija_aliasK(:, j, i))
                    if ((abs(hb%ijab_weights_tot(a, j, i)) > 0.0_p) .and. (.not.btest(cdet%f(sys%basis%bit_lookup(2,a)), &
                            sys%basis%bit_lookup(1,a)))) then
                        b = select_weighted_value_precalc(rng, sys%basis%nbasis, hb%ijab_aliasU(:, a, j, i), &
                                                            hb%ijab_aliasK(:, a, j, i))
                        if (.not.btest(cdet%f(sys%basis%bit_lookup(2,b)), sys%basis%bit_lookup(1,b))) then
                            allowed_excitation = .true.               
                        else
                            allowed_excitation = .false.
                        end if
                    else
                        allowed_excitation = .false.
                    end if
                else
                    allowed_excitation = .false.
                end if

                if (allowed_excitation) then
                    if (i < j) then
                        connection%from_orb(1) = i
                        connection%from_orb(2) = j
                    else
                        connection%from_orb(1) = j
                        connection%from_orb(2) = i
                    end if
                    if (a < b) then
                        connection%to_orb(1) = a
                        connection%to_orb(2) = b
                    else
                        connection%to_orb(2) = a
                        connection%to_orb(1) = b
                    end if
                    connection%nexcit = 2
                    ! 4b. Parity of permutation required to line up determinants.
                    ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                    call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

                    ! 5b. Find the connecting matrix element.
                    hmatel = slater_condon2_excit_ptr(sys, connection%from_orb(1), connection%from_orb(2), &
                                                      connection%to_orb(1), connection%to_orb(2), connection%perm)
                    pgen = (1.0_p - excit_gen_data%pattempt_single) * pgen * &
                        (((hb%ija_weights(a, j, i)/hb%ija_weights_tot(j, i)) * &
                        (hb%ijab_weights(b, a, j, i)/hb%ijab_weights_tot(a, j, i))) + &
                        ((hb%ija_weights(b, j, i)/hb%ija_weights_tot(j, i)) * &
                        (hb%ijab_weights(a, b, j, i)/hb%ijab_weights_tot(b, j, i))))
                else
                    ! We have not found a valid excitation ij -> ab.  Such
                    ! events are not worth the cost of renormalising the generation
                    ! probabilities.
                    ! Return a null excitation.
                    hmatel%c = cmplx(0.0_p, 0.0_p, p)
                    hmatel%r = 0.0_p
                    pgen = 1.0_p
                end if
                ! deallocate weight arrays if allocated
                if (allocated(i_weights_occ)) deallocate(i_weights_occ)
                if (allocated(ij_weights_occ)) deallocate(ij_weights_occ)
                if (allocated(ji_weights_occ)) deallocate(ji_weights_occ)
            end associate
        end if

    end subroutine gen_excit_mol_heat_bath_uniform

    subroutine gen_single_excit_heat_bath_exact(rng, sys, pattempt_single, cdet, pgen, connection, hmatel, &
            allowed_excitation)

        ! Create a random single excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! This calculates the weights for i and a as exact as possibly.
        ! For i, the weight is \sum_a H_ia, for a given i it is H_ia.

        ! In:
        !    sys: system object being studied.
        !    pattempt_single: probability of having chosen to attempt a single and not a double
        !                    excitation.
        !    cdet: determinant to attempt spawning from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !       determinant and a connected determinant in molecular systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation1
        use proc_pointers, only: abs_hmatel_ptr, slater_condon1_excit_ptr
        use system, only: sys_t
        use hamiltonian_data, only: hmatel_t
        use hamiltonian_molecular, only: slater_condon1_mol
        use hamiltonian_periodic_complex, only: slater_condon1_periodic_complex
        use dSFMT_interface, only: dSFMT_t
        use alias, only: select_weighted_value

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: pattempt_single
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation

        integer :: virt_list(sys%basis%nbasis - sys%nel)
        real(p) :: i_weights(sys%nel), ia_weights(sys%basis%nbasis - sys%nel, sys%nel), ia_weights_tot(sys%nel)
        real(p) :: i_weights_tot
        integer :: pos, virt_pos, occ_pos, i_ind, a_ind, i, a, nvirt

        nvirt = sys%basis%nbasis - sys%nel
        allowed_excitation = .true.
        connection%nexcit = 1
        
        ! Make list of virtual orbitals, virt_list.
        virt_pos = 1
        occ_pos = 1
        do pos = 1, sys%basis%nbasis
            if (occ_pos > sys%nel) then
                virt_list(virt_pos) = pos
                virt_pos = virt_pos + 1
            else
                if (cdet%occ_list(occ_pos) == pos) then
                    occ_pos = occ_pos + 1
                else
                    virt_list(virt_pos) = pos
                    virt_pos = virt_pos + 1
                end if
            end if
        end do
        
        ! Calculate weights.
        i_weights_tot = 0.0_p
        ia_weights_tot = 0.0_p
        do i_ind = 1, sys%nel
            i_weights(i_ind) = 0.0_p
            do a_ind = 1, nvirt
                ! [todo] - this reduces the computational time (by calculation ia_weights here)
                ! but more memory costs as after we have selected i, we don't need all elements in ia_weights. Need to balance these
                ! costs.
                ! [todo] - could consider rewritting with proc_pointer here but that would imply having to rewrite slater
                ! condon 1 mol / periodic complex. slater_condon1_excit_mol is not safe enough (no checks).
                if (sys%read_in%comp) then
                    hmatel%r = 0.0_p
                    hmatel%c = slater_condon1_periodic_complex(sys, cdet%occ_list, cdet%occ_list(i_ind), virt_list(a_ind), &
                        .false.)
                else
                    hmatel%r = slater_condon1_mol(sys, cdet%occ_list, cdet%occ_list(i_ind), virt_list(a_ind), .false.)
                    hmatel%c = cmplx(0.0_p, 0.0_p, p)
                end if
                ia_weights(a_ind, i_ind) = abs_hmatel_ptr(hmatel)
                ia_weights_tot(i_ind) = ia_weights_tot(i_ind) + ia_weights(a_ind, i_ind)
                i_weights(i_ind) = i_weights(i_ind) + ia_weights(a_ind, i_ind)
            end do
            i_weights_tot = i_weights_tot + i_weights(i_ind)
        end do

        if (i_weights_tot < depsilon) then
            ! no allowed single excitations.
            allowed_excitation = .false.
            hmatel%c = cmplx(0.0_p, 0.0_p, p)
            hmatel%r = 0.0_p
            pgen = 1.0_p ! Avoid any dangerous division by pgen by returning a sane (but cheap) value.
        else
            ! Select i.
            i_ind = select_weighted_value(rng, sys%nel, i_weights, i_weights_tot)
            i = cdet%occ_list(i_ind)

            ! Select a.
            a_ind = select_weighted_value(rng, nvirt, ia_weights(:, i_ind), ia_weights_tot(i_ind))
            a = virt_list(a_ind)

            connection%from_orb(1) = i
            connection%to_orb(1) = a

            ! Calculate pgen and hmatel.
            pgen = pattempt_single * (i_weights(i_ind)/i_weights_tot) * (ia_weights(a_ind, i_ind)/ia_weights_tot(i_ind))

            ! Parity of permutation required to line up determinants.
            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)

            ! Find the connecting matrix element.
            hmatel = slater_condon1_excit_ptr(sys, cdet%occ_list, connection%from_orb(1), connection%to_orb(1), connection%perm)
        end if

    end subroutine gen_single_excit_heat_bath_exact
    
end module excit_gen_heat_bath_mol
