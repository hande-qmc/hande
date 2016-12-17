module excit_gen_cauchy_schwarz_mol

! A module containing excitations generators for molecules which weight excitations according to the exchange matrix elements.

use const, only: i0, p

implicit none

! [review] - JSS: some code-style uniformity please.  e.g. vertical space before/after procedures, around the procedure comments,
! [review] - JSS: indented comments, end do/end if instead of end do and end if, spaces around binary operators (=, <, ...), ...

contains

    subroutine init_excit_mol_cauchy_schwarz_occ_ref(sys, ref, cs)
        ! [review] - JSS: How good is this really?  I suspect it's okish for single-reference truncated calculations and quickly
        ! [review] - JSS: becomes bad for multi-reference/as the truncation level increases.
        ! [review] - JSS: Is this an experiment?  Ready for use in production calculations?

        ! Generate excitation tables from the reference for the 
        ! gen_excit_mol_cauchy_schwarz_occ_ref excitation generator.
        ! This creates a random excitation from a det and calculates both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        ! Weight the double excitations according the the Cauchy-Schwarz bound
        ! <ij|ab> <= Sqrt(<ia|ai><jb|bj>)
        ! This is an O(M/64) version which pretends the determinant excited from is the reference,
        ! then modifies the selected orbitals to be those of the determinant given.
        ! Each occupied and virtual not present in the det we're actually given will be
        ! mapped to the one of the equivalent free numerical index in the reference.

        ! In:
        !    sys: system object being studied.
        !    ref: the reference from which we are exciting.
        ! In/Out:
        !    cs: an empty excit_gen_cauchy_schwarz_t object which gets filled with
        !           the alias tables required to generate excitations.

        use system, only: sys_t
        use qmc_data, only: reference_t
        use sort, only: qsort
        use proc_pointers, only: create_weighted_excitation_list_ptr
        use excit_gens, only: excit_gen_cauchy_schwarz_t
        use alias, only: generate_alias_tables
        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: ref
        type(excit_gen_cauchy_schwarz_t), intent(inout) :: cs

        integer :: i, j, ind_a, ind_b, maxv, nv, jsym
        
        ! Temp storage
        maxv = max(sys%nvirt_alpha,sys%nvirt_beta)
        allocate(cs%ia_aliasP(maxv,sys%nel))
        allocate(cs%ia_aliasY(maxv,sys%nel))
        allocate(cs%ia_weights(maxv,sys%nel))
        allocate(cs%ia_weights_tot(sys%nel))
        allocate(cs%jb_aliasP(maxval(sys%read_in%pg_sym%nbasis_sym_spin), sys%sym0_tot:sys%sym_max_tot, sys%nel))
        allocate(cs%jb_aliasY(maxval(sys%read_in%pg_sym%nbasis_sym_spin), sys%sym0_tot:sys%sym_max_tot, sys%nel))
        allocate(cs%jb_weights(maxval(sys%read_in%pg_sym%nbasis_sym_spin), sys%sym0_tot:sys%sym_max_tot, sys%nel))
        allocate(cs%jb_weights_tot(sys%sym0_tot:sys%sym_max_tot, sys%nel))
        allocate(cs%occ_list(sys%nel+1))  ! The +1 is a pad to allow loops to look better
        allocate(cs%virt_list_alpha(sys%nvirt_alpha))
        allocate(cs%virt_list_beta(sys%nvirt_beta))
        
        cs%occ_list(:sys%nel) = ref%occ_list0(:sys%nel)
        cs%occ_list(sys%nel+1) = sys%basis%nbasis*2  ! A pad 
        ! Now sort this
        ! [review] - JSS: is this not always true?
        call qsort(cs%occ_list,sys%nel)

        ! Make the unocc list.
        j = 1     ! The next occ to look at
        ind_a = 0 ! The present position in the virt_list we're making
        ind_b = 0 ! The present position in the virt_list we're making
        
        do i = 1, sys%basis%nbasis
            if (i==cs%occ_list(j)) then ! Our basis fn is in the ref
                j = j + 1
            else ! Need to store it as a virt
                if (sys%basis%basis_fns(i)%Ms == -1) then ! beta
                    ind_b = ind_b + 1
                    cs%virt_list_beta(ind_b) = i
                else
                    ind_a = ind_a + 1
                    cs%virt_list_alpha(ind_a) = i
                end if
            end if
        end do

        do i = 1, sys%nel
            j = cs%occ_list(i)  ! The elec we're looking at
            if (sys%basis%basis_fns(j)%Ms == -1) then ! beta
                nv = sys%nvirt_beta
                call create_weighted_excitation_list_ptr(sys, j, 0, cs%virt_list_beta, nv, cs%ia_weights(:,i), cs%ia_weights_tot(i))
                call generate_alias_tables(nv, cs%ia_weights(:,i), cs%ia_weights_tot(i), cs%ia_aliasP(:,i), cs%ia_aliasY(:,i))
                do jsym = sys%sym0_tot, sys%sym_max_tot
                    call create_weighted_excitation_list_ptr(sys, j, 0, sys%read_in%pg_sym%sym_spin_basis_fns(:,1,jsym), &
                        sys%read_in%pg_sym%nbasis_sym_spin(1,jsym), cs%jb_weights(:,jsym,i), cs%jb_weights_tot(jsym,i))
                    call generate_alias_tables(sys%read_in%pg_sym%nbasis_sym_spin(1,jsym), cs%jb_weights(:,jsym,i), &
                        cs%jb_weights_tot(jsym,i), cs%jb_aliasP(:,jsym,i), cs%jb_aliasY(:,jsym,i))
                end do
            else ! alpha
                nv = sys%nvirt_alpha
                call create_weighted_excitation_list_ptr(sys, j, 0, cs%virt_list_alpha, nv, cs%ia_weights(:,i), cs%ia_weights_tot(i))
                call generate_alias_tables(nv, cs%ia_weights(:,i), cs%ia_weights_tot(i), cs%ia_aliasP(:,i), cs%ia_aliasY(:,i))
                do jsym = sys%sym0_tot, sys%sym_max_tot
                    call create_weighted_excitation_list_ptr(sys, j, 0, sys%read_in%pg_sym%sym_spin_basis_fns(:,2,jsym), &
                        sys%read_in%pg_sym%nbasis_sym_spin(2,jsym), cs%jb_weights(:,jsym,i), cs%jb_weights_tot(jsym,i))
                    call generate_alias_tables(sys%read_in%pg_sym%nbasis_sym_spin(2,jsym), cs%jb_weights(:,jsym,i), &
                        cs%jb_weights_tot(jsym,i), cs%jb_aliasP(:,jsym,i), cs%jb_aliasY(:,jsym,i))
                end do
            end if
        end do
    end subroutine init_excit_mol_cauchy_schwarz_occ_ref

    subroutine gen_excit_mol_cauchy_schwarz_occ_ref(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        ! Weight the double excitations according the the Cauchy-Schwarz bound
        ! <ij|ab> <= Sqrt(<ia|ai><jb|bj>)
        ! This is an O(M/64) version which pretends the determinant excited from is the reference,
        ! then modifies the selected orbitals to be those of the determinant given.
        ! Each occupied and virtual not present in the det we're actually given will be
        ! mapped to the one of the equivalent free numerical index in the reference.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Excitation generation data.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
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
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2,get_excitation_locations
        use proc_pointers, only: slater_condon1_excit_ptr, slater_condon2_excit_ptr
        use system, only: sys_t
        use excit_gen_mol, only: calc_pgen_single_mol_no_renorm,find_ia_mol
        use excit_gens, only: excit_gen_cauchy_schwarz_t, excit_gen_data_t
        use alias, only: select_weighted_value_prec
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use hamiltonian_data, only: hmatel_t


        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel 
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation


        integer ::  ij_spin
      
        integer :: a, b, i, j, a_ind, b_ind, i_ind, j_ind

        integer :: nex
        integer :: cdet_store(sys%nel)
        integer :: ref_store(sys%nel)
        integer :: ii, jj, t

        ! 1. Select single or double.
        if (get_rand_close_open(rng) < excit_gen_data%pattempt_single) then

            ! 2a. Select orbital to excite from and orbital to excite into.
            call find_ia_mol(rng, sys, sys%read_in%pg_sym%gamma_sym, cdet%f, cdet%occ_list, connection%from_orb(1), &
                             connection%to_orb(1), allowed_excitation)
            connection%nexcit = 1

            if (allowed_excitation) then
                ! 3a. Probability of generating this excitation.
                pgen = excit_gen_data%pattempt_single*calc_pgen_single_mol_no_renorm(sys, connection%to_orb(1))

                ! 4a. Parity of permutation required to line up determinants.
                call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)

                ! 5a. Find the connecting matrix element.
                hmatel = slater_condon1_excit_ptr(sys, cdet%occ_list, connection%from_orb(1), & 
                                            connection%to_orb(1), connection%perm)
            else
                ! Forbidden---connection%to_orb(1) is already occupied.
                hmatel%c = cmplx(0.0_p, 0.0_p, p)
                hmatel%r = 0.0_p
                pgen = 1.0_p ! Avoid any dangerous division by pgen by returning a sane (but cheap) value.
            end if

        else
            ! We have a double
            associate( cs => excit_gen_data%excit_gen_cs )

                ! 2b. Select orbitals to excite from
                
                call choose_ij_ind(rng, sys, cs%occ_list, i_ind, j_ind, ij_spin)

                ! At this point we pretend we're the reference, and fix up mapping ref's orbitals to cdet's orbitals later.
                i = cs%occ_list(i_ind)
                j = cs%occ_list(j_ind)

                ! We now need to select the orbitals to excite into which we do with weighting:
                ! p(ab|ij) = p(a|i) p(b|j) + p(a|j) p(b|i)
                ! We actually choose a|i then b|j, but since we could have also generated the excitation b from i and a from j, we
                ! need to include that prob too.

                ! Given i, use the alias table to select a
                if (sys%basis%basis_fns(i)%Ms < 0) then
                    a_ind = select_weighted_value_prec(rng, sys%nvirt_beta, cs%ia_aliasP(:,i_ind), cs%ia_aliasY(:,i_ind))
                    ! [review] - JSS: is 4 copies of the identical comment really necessary?  I don't think any of them are...(also
                    ! [review] - JSS: not clear if the comment is correct here!)
                    ! Use the alias method to select i with the appropriate probability
                    a = cs%virt_list_beta(a_ind) 
                else
                    a_ind = select_weighted_value_prec(rng, sys%nvirt_alpha, cs%ia_aliasP(:,i_ind), cs%ia_aliasY(:,i_ind))
                    ! Use the alias method to select i with the appropriate probability
                    a = cs%virt_list_alpha(a_ind) 
                end if 

                ! Given j use the alias table to select b
                if (sys%basis%basis_fns(j)%Ms < 0) then
                    b_ind = select_weighted_value_prec(rng, sys%nvirt_beta, cs%ia_aliasP(:,j_ind), cs%ia_aliasY(:,j_ind))
                    ! Use the alias method to select i with the appropriate probability
                    b = cs%virt_list_beta(b_ind) 
                else
                    b_ind = select_weighted_value_prec(rng, sys%nvirt_alpha, cs%ia_aliasP(:,j_ind), cs%ia_aliasY(:,j_ind))
                    ! Use the alias method to select i with the appropriate probability
                    b = cs%virt_list_alpha(b_ind) 
                end if 

                ! 3b. Probability of generating this excitation.

                ! Calculate p(ab|ij) = p(a|i) p(j|b) + p(b|i)p(a|j)
                if (ij_spin==0) then 
                    ! Not possible to have chosen the reversed excitation
                    pgen = cs%ia_weights(a_ind,i_ind) / cs%ia_weights_tot(i_ind) &
                            * cs%ia_weights(b_ind,j_ind) / cs%ia_weights_tot(j_ind)
                else
                    ! i and j have same spin, so could have been selected in the other order.
                    pgen = ( cs%ia_weights(a_ind, i_ind) * cs%ia_weights(b_ind, j_ind) &
                             + cs%ia_weights(b_ind, i_ind) * cs%ia_weights(a_ind, j_ind) ) &
                           / (cs%ia_weights_tot(i_ind)*cs%ia_weights_tot(j_ind))
                end if
                pgen = excit_gen_data%pattempt_double * pgen * 2.0_p/(sys%nel*(sys%nel-1)) ! pgen(ab)
                connection%nexcit = 2
                allowed_excitation = (a/=b)

                if (allowed_excitation) then
                    ! Now do the translation.  We need a list of the substitutio ns in cdet vs the ref.  i.e. for each orbital in
                    ! the ref-lined-up cdet which is not in ref, we need the location.
                    ! This is currently done with an O(N) step, but might be sped up at least.

                    call get_excitation_locations(cs%occ_list, cdet%occ_list, ref_store, cdet_store, sys%nel, nex)
                    ! These orbitals might not be aligned in the most efficient way:
                    !  They may not match in spin, so first deal with this

                    ! ref store (e.g.) contains the indices within cs%occ_list of the orbitals
                    ! which have been excited from.
                    ! [review] - JSS: this could/should be done once per determinant/excitor rather than once per excit gen.
                    do ii=1, nex
                        associate(bfns=>sys%basis%basis_fns, ref_orb=>cs%occ_list(ref_store(ii)), det=>cdet%occ_list)
                            if (bfns(ref_orb)%Ms /= bfns(det(cdet_store(ii)))%Ms) then
                                jj = ii + 1
                                do while (bfns(ref_orb)%Ms /= bfns(det(cdet_store(jj)))%Ms)
                                    jj = jj + 1
                                end do
                                ! det's jj now points to an orb of the same spin as ref's ii, so swap cdet_store's ii and jj.
                                t = cdet_store(ii)
                                cdet_store(ii) = cdet_store(jj)
                                cdet_store(jj) = t 
                            end if
                        end associate
                    end do
                    ! At this point we may want to align the orbitals even further to avoid strange cases where selection
                    ! probabilities are zero, but let's leave that for another day.

                    ! Now see if i and j are in ref_store or a and b are in cdet_store and map appropriately

                    connection%from_orb(1)=i
                    connection%from_orb(2)=j
                    connection%to_orb(1)=a
                    connection%to_orb(2)=b
                    ! ref store (e.g.) contains the indices within cs%occ_list of the orbitals which are excited out of ref into
                    ! cdet det store (e.g.) contains the indices within cdet%occ_list of the orbitals which are in cdet (excited out
                    ! of ref).  i_ind and j_ind are the indices of the orbitals in cs%occ_list which we're exciting from.
                    do ii=1, nex
                        if (ref_store(ii)==i_ind) then  ! from_orb(1) isn't actually in cdet, so we replace it with the orb that is
                            connection%from_orb(1) = cdet%occ_list(cdet_store(ii))
                        else if (ref_store(ii)==j_ind) then ! from_orb(2) isn't actually in cdet, so we replace it with the orb that is
                            connection%from_orb(2) = cdet%occ_list(cdet_store(ii))
                        end if
                        if (cdet%occ_list(cdet_store(ii))==a) then
                            connection%to_orb(1) = cs%occ_list(ref_store(ii))
                        else if (cdet%occ_list(cdet_store(ii))==b) then
                            connection%to_orb(2) = cs%occ_list(ref_store(ii))
                        end if
                    end do
                    ! Now order the excitation
                    if( connection%to_orb(1) > connection%to_orb(2) ) then
                        t = connection%to_orb(1)
                        connection%to_orb(1) = connection%to_orb(2)
                        connection%to_orb(2) = t
                    end if
                    if( connection%from_orb(1) > connection%from_orb(2) ) then
                        t = connection%from_orb(1)
                        connection%from_orb(1) = connection%from_orb(2)
                        connection%from_orb(2) = t
                    end if


                    ! 4b. Parity of permutation required to line up determinants.
                    ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                    call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

                    ! 5b. Find the connecting matrix element.
                    hmatel = slater_condon2_excit_ptr(sys, connection%from_orb(1), connection%from_orb(2), &
                                                      connection%to_orb(1), connection%to_orb(2), connection%perm)
                else
                    ! Carelessly selected ij with no possible excitations.  Such
                    ! events are not worth the cost of renormalising the generation
                    ! probabilities.
                    ! Return a null excitation.
                    hmatel%c = cmplx(0.0_p, 0.0_p, p)
                    hmatel%r = 0.0_p
                    pgen = 1.0_p
                end if
            end associate
        end if

    end subroutine gen_excit_mol_cauchy_schwarz_occ_ref

    ! [review] - JSS: close to choose_ij_mol but without symmetry?  If so, unnecessary code duplication.
    subroutine choose_ij_ind(rng, sys, occ_list, i_ind, j_ind, ij_spin)

        ! Randomly select two occupied orbitals in a determinant from which
        ! electrons are excited as part of a double excitation.
        !
        ! In:
        !    sys: system object being studied.
        !    occ_list: integer list of occupied spin-orbitals in the determinant.
        !        (min length: sys%nel.)
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    i, j: orbitals in determinant from which two electrons are excited.
        !        Note that i,j are ordered such that i<j.
        ! [review] JSS: ij_sym has vanished.
        !    ij_sym: symmetry conjugate of the irreducible representation spanned by the codensity
        !        \phi_i*\phi_j. (We assume that ij is going to be in the bra of the excitation.)
        !    ij_spin: spin label of the combined ij codensity.
        !        ij_spin = -2   i,j both down
        !                =  0   i up and j down or vice versa
        !                =  2   i,j both up

        use system, only: sys_t

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(:)
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: i_ind, j_ind, ij_spin

        integer :: ind, i, j

        ! See comments in choose_ij_k for how the occupied orbitals are indexed
        ! to allow one random number to decide the ij pair.

        ind = int(get_rand_close_open(rng)*sys%nel*(sys%nel-1)/2) + 1

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.
        ! Note that the indexing scheme for the strictly lower triangular array
        ! assumes j>i.  As occ_list is ordered, this means that we will return
        ! i,j (referring to spin-orbitals) where j>i.  This ordering is
        ! convenient subsequently, e.g. is assumed in the
        ! find_excitation_permutation2 routine.
        j_ind = int(1.50_p + sqrt(2*ind-1.750_p))
        i_ind = ind - ((j_ind-1)*(j_ind-2))/2

        i = occ_list(i_ind)
        j = occ_list(j_ind)

        ! ij_spin = -2 (down, down), 0 (up, down or down, up), +2 (up, up)
        ij_spin = sys%basis%basis_fns(i)%Ms + sys%basis%basis_fns(j)%Ms

    end subroutine choose_ij_ind

    subroutine gen_excit_mol_cauchy_schwarz_occ(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        ! Weight the double excitations according the the Cauchy-Schwarz bound
        ! <ij|ab> <= Sqrt(<ia|ai><jb|bj>)
        ! This requires a lookup of O(M) two-electron integrals in its setup.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Excitation generation data.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
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
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2
        use proc_pointers, only: slater_condon1_excit_ptr, slater_condon2_excit_ptr, create_weighted_excitation_list_ptr
        use system, only: sys_t
        use excit_gen_mol, only: calc_pgen_single_mol_no_renorm, find_ia_mol, choose_ij_mol
        use hamiltonian_data, only: hmatel_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use search, only: binary_search
        use checking, only: check_allocate, check_deallocate
        use excit_gens, only: excit_gen_cauchy_schwarz_t, excit_gen_data_t
        use alias, only: select_weighted_value

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation

        integer ::  ij_spin, ij_sym, ierr
        logical :: found
        real(p), allocatable :: ia_weights(:), ja_weights(:), jb_weights(:)

        !real(p) :: ia_weights(max(sys%nalpha + sys%nel,sys%nbeta+ sys%nel)), ja_weights(max(sys%nalpha+sys%nel,sys%nbeta+sys%nel)),
        !jb_weights(max(sys%nalpha+sys%nel,sys%nbeta+sys%nel))
      
        real(p) :: ia_weights_tot, ja_weights_tot, jb_weights_tot
        integer :: a, b, i, j, a_ind, b_ind, a_ind_rev, b_ind_rev, isymb, imsb, isyma

        ! 1. Select single or double.
        if (get_rand_close_open(rng) < excit_gen_data%pattempt_single) then

            ! [review] - JSS: this block is identical for renorm, no_renorm and the three Cauchy-Schwarz generators.  Should abstract it.
            ! 2a. Select orbital to excite from and orbital to excite into.
            call find_ia_mol(rng, sys, sys%read_in%pg_sym%gamma_sym, cdet%f, cdet%occ_list, connection%from_orb(1), &
                             connection%to_orb(1), allowed_excitation)
            connection%nexcit = 1

            if (allowed_excitation) then
                ! 3a. Probability of generating this excitation.
                pgen = excit_gen_data%pattempt_single*calc_pgen_single_mol_no_renorm(sys, connection%to_orb(1))

                ! 4a. Parity of permutation required to line up determinants.
                call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)

                ! 5a. Find the connecting matrix element.
                hmatel = slater_condon1_excit_ptr(sys, cdet%occ_list, connection%from_orb(1), connection%to_orb(1), connection%perm)
            else
                ! Forbidden---connection%to_orb(1) is already occupied.
                hmatel%c = cmplx(0.0_p, 0.0_p, p)
                hmatel%r = 0.0_p
                pgen = 1.0_p ! Avoid any dangerous division by pgen by returning a sane (but cheap) value.
            end if

        else
            ! We have a double

            ! 2b. Select orbitals to excite from
            
            call choose_ij_mol(rng, sys, cdet%occ_list, i, j, ij_sym, ij_spin)

            ! Now we've chosen i and j. 
 
            ! We now need to select the orbitals to excite into which we do with weighting:
            ! p(ab|ij) = p(a|i) p(b|j) + p(a|j) p(b|i)
            
            ! We actually choose a|i then b|j, but since we could have also generated the excitation b from i and a from j, we need to include that prob too.

            ! Given i, construct the weights of all possible a
            if (sys%basis%basis_fns(i)%Ms < 0) then
                ! [review] - JSS: is it really worth constructing this (O(N) time) rather than just going for a binary (or even
                ! [review] - JSS: linear) search on the cumulative table directly?
                allocate(ia_weights(1:sys%nvirt_beta), stat=ierr)
                call check_allocate('ia_weights', sys%nvirt_beta, ierr)
                call create_weighted_excitation_list_ptr(sys, i, 0, cdet%unocc_list_beta, sys%nvirt_beta, ia_weights, ia_weights_tot)
                ! Use the alias method to select i with the appropriate probability
                a_ind = select_weighted_value(rng, sys%nvirt_beta, ia_weights, ia_weights_tot)
                a = cdet%unocc_list_beta(a_ind) 
            else
                allocate(ia_weights(1:sys%nvirt_alpha), stat=ierr)
                call check_allocate('ia_weights', sys%nvirt_alpha, ierr)

                call create_weighted_excitation_list_ptr(sys, i, 0, cdet%unocc_list_alpha, sys%nvirt_alpha, ia_weights, ia_weights_tot)
                ! Use the alias method to select i with the appropriate probability
                a_ind = select_weighted_value(rng, sys%nvirt_alpha, ia_weights, ia_weights_tot)
                a = cdet%unocc_list_alpha(a_ind) 
            end if

            ! Given i,j,a construct the weights of all possible b
            ! This requires that total symmetry and spin are conserved.
            ! The symmetry of b (isymb) is given by 
            ! (sym_i* x sym_j* x sym_a)* = sym_b
            ! (at least for Abelian point groups)
            isymb = sys%read_in%sym_conj_ptr(sys%read_in, &
                    sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, sys%basis%basis_fns(a)%sym))
            ! Ms_i + Ms_j = Ms_a + Ms_b => Ms_b = Ms_i + Ms_j - Ms_a
            ! Ms_k is +1 if up and -1 if down but imsb is +2 if up and +1 if down, 
            ! therefore a conversion is necessary.
            imsb = (ij_spin-sys%basis%basis_fns(a)%Ms+3)/2
            
            if (sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb) > 0) then 
                allocate(jb_weights(1:sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb)), stat=ierr)
                call check_allocate('jb_weights', sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb), ierr)
                call create_weighted_excitation_list_ptr(sys, j, a, sys%read_in%pg_sym%sym_spin_basis_fns(:,imsb,isymb), &
                                        sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb), jb_weights, jb_weights_tot)
            else
                jb_weights_tot = 0.0_p
            end if

            ! Test whether at least one possible b given i,j,a exists.
            ! Note that we did not need a btest for orbital a because we only considered
            ! virtual orbitals there.

            if (jb_weights_tot > 0.0_p) then
                ! Use the alias method to select b with the appropriate probability
                b_ind = select_weighted_value(rng, sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb), jb_weights, jb_weights_tot)
                b = sys%read_in%pg_sym%sym_spin_basis_fns(b_ind,imsb,isymb)

                if (.not.btest(cdet%f(sys%basis%bit_lookup(2,b)), sys%basis%bit_lookup(1,b))) then
         
                    ! 3b. Probability of generating this excitation.

                    ! Calculate p(ab|ij) = p(a|i) p(j|b) + p(b|i)p(a|j)
                    if (ij_spin==0) then 
                        ! not possible to have chosen the reversed excitation
                        pgen=ia_weights(a_ind)/ia_weights_tot*jb_weights(b_ind)/jb_weights_tot
                    else 
                        ! i and j have same spin, so could have been selected in the other order.
                        if (imsb == 1) then 
                            ! find index b as if we had it selected first and as a from list of unoccupied virtual orbitals.
                            ! This uses the fact that the position of b in cdet%unocc_list_{beta,alpha} is the same
                            ! position index as in b's index in ia_weights is. 
                            call binary_search(cdet%unocc_list_beta, b, 1, sys%nvirt_beta, found, b_ind_rev)
                        else 
                            call binary_search(cdet%unocc_list_alpha, b, 1, sys%nvirt_alpha, found, b_ind_rev)
                        end if

                        isyma = sys%read_in%sym_conj_ptr(sys%read_in, &
                                sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, isymb))
                        ! imsa = imsb
                        allocate(ja_weights(1:sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma)), stat=ierr)
                        call check_allocate('ja_weights', sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma), ierr)
                        call create_weighted_excitation_list_ptr(sys, j, b, sys%read_in%pg_sym%sym_spin_basis_fns(:,imsb,isyma), &
                                            sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma), ja_weights, ja_weights_tot)
                        call binary_search(sys%read_in%pg_sym%sym_spin_basis_fns(:,imsb,isyma), a, 1, &
                                    sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma), found, a_ind_rev)
                        pgen=       (ia_weights(a_ind)*jb_weights(b_ind)) / (ia_weights_tot*jb_weights_tot) + &
                            (ia_weights(b_ind_rev)*ja_weights(a_ind_rev) ) / (ia_weights_tot*ja_weights_tot)
                    end if
                    pgen = excit_gen_data%pattempt_double * pgen *2.0_p/(sys%nel*(sys%nel-1)) ! pgen(ab)
                    connection%nexcit = 2

                    ! if create_weighted_excitation_list checks this then this check is redundant
                    ! allowed_excitation = a /= b

                    if (i<j) then
                        connection%from_orb(1) = i
                        connection%from_orb(2) = j
                    else
                        connection%from_orb(2) = i
                        connection%from_orb(1) = j
                    end if
                    if (a<b) then
                        connection%to_orb(1) = a
                        connection%to_orb(2) = b 
                    else
                        connection%to_orb(2) = a
                        connection%to_orb(1) = b 
                    end if
                    allowed_excitation = .true.
                else
                    allowed_excitation = .false.
                end if
            else
                allowed_excitation = .false. 
            end if
           
            if (allowed_excitation) then

                ! 4b. Parity of permutation required to line up determinants.
                ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

                ! 5b. Find the connecting matrix element.
                hmatel = slater_condon2_excit_ptr(sys, connection%from_orb(1), connection%from_orb(2), &
                                            connection%to_orb(1), connection%to_orb(2), connection%perm)
            else
                ! Carelessly selected ij with no possible excitations.  Such
                ! events are not worth the cost of renormalising the generation
                ! probabilities.
                ! Return a null excitation.
                hmatel%c = cmplx(0.0_p, 0.0_p, p)
                hmatel%r = 0.0_p
                pgen = 1.0_p
            end if

        end if
    end subroutine gen_excit_mol_cauchy_schwarz_occ

end module excit_gen_cauchy_schwarz_mol
