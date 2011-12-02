module spawning_mol_system

! Module for spawning routine(s) related to the molecular system (ie one read in
! from an FCIDUMP file).

! TODO:
!  * pgen evaluation
!  * test!
!  * comment

use const, only: i0, p

implicit none

contains

!--- Top level spawning routines ---

    subroutine spawn_mol(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant for the 
        ! molecular systems.
        !
        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use excitations, only: excit
        use fciqmc_data, only: tau

        use dSFMT_interface, only:  genrand_real2

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel, pspawn

        ! 1. Generate random excitation.
        call gen_excit_mol(cdet, pgen, connection, hmatel)

        ! 2. Calculate P_spawn.
        pspawn = tau*abs(hmatel)/pgen

        ! 3. Attempt spawning.
        nspawn = int(pspawn)
        pspawn = pspawn - nspawn
        if (pspawn > genrand_real2()) nspawn = nspawn + 1
        if (hmatel > 0) then
            nspawn = -sign(nspawn, parent_sign)
        else
            nspawn = sign(nspawn, parent_sign)
        end if

    end subroutine spawn_mol

    subroutine spawn_mol_no_renorm(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant for the 
        ! molecular systems.
        !
        ! This doesn't use excitation generators which exclude the case where,
        ! having selected (e.g.) 2 occupied orbitals (i and j) and the first virtual
        ! orbital (a), the final orbital is already occupied and so that
        ! excitation is impossible.  Whilst it is somewhat wasteful (generating
        ! excitations which can't be performed), there is a balance between the
        ! cost of generating forbidden excitations and the O(N) cost of
        ! renormalising the generation probabilities.
        !
        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use excitations, only: excit
        use fciqmc_data, only: tau

        use dSFMT_interface, only:  genrand_real2

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel, pspawn

        ! 1. Generate random excitation.
        call gen_excit_mol_no_renorm(cdet, pgen, connection, hmatel)

        ! 2. Calculate P_spawn.
        pspawn = tau*abs(hmatel)/pgen

        ! 3. Attempt spawning.
        nspawn = int(pspawn)
        pspawn = pspawn - nspawn
        if (pspawn > genrand_real2()) nspawn = nspawn + 1
        if (hmatel > 0) then
            nspawn = -sign(nspawn, parent_sign)
        else
            nspawn = sign(nspawn, parent_sign)
        end if

    end subroutine spawn_mol_no_renorm

!--- Excitation generation ---

    subroutine gen_excit_mol(cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate both the probability 
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a 
        !    determinant and a connected determinant in molecular systems.

        use determinants, only: det_info
        use excitations, only: excit
        use fciqmc_data, only: pattempt_single
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2
        use hamiltonian_molecular, only: slater_condon1_mol, slater_condon2_mol

        use dSFMT_interface, only: genrand_real2

        type(det_info), intent(in) :: cdet
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        integer :: ij_sym, ij_spin

        ! 1. Select single or double.
        
        if (genrand_real2() < pattempt_single) then
        
            ! 2a. Select orbital to excite from and orbital to excit into.
            call choose_ia_mol(cdet%f, cdet%occ_list, connection%from_orb(1), connection%to_orb(1))
            connection%nexcit = 1
            
            ! 3a. Probability of generating this excitation.

            ! 4a. Parity of permutation required to line up determinants.
            call find_excitation_permutation1(cdet%f, connection)

            ! 5a. Find the connecting matrix element.
            hmatel = slater_condon1_mol(cdet%occ_list, connection%from_orb(1), connection%to_orb(1), connection%perm)

        else

            ! 2b. Select orbitals to excite from and orbitals to excite into.
            call choose_ij_mol(cdet%occ_list, connection%from_orb(1), connection%from_orb(2), ij_sym, ij_spin)
            call choose_ab_mol(cdet%f, ij_sym, ij_spin, connection%to_orb(1), connection%to_orb(2))
            connection%nexcit = 2
            
            ! 3b. Probability of generating this excitation.

            ! 4b. Parity of permutation required to line up determinants.
            ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
            call find_excitation_permutation2(cdet%f, connection)

            ! 5b. Find the connecting matrix element.
            hmatel = slater_condon2_mol(connection%from_orb(1), connection%from_orb(2), &
                                        connection%to_orb(1), connection%to_orb(2), connection%perm)

        end if

    end subroutine gen_excit_mol

    subroutine gen_excit_mol_no_renorm(cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate both the probability 
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a 
        !    determinant and a connected determinant in molecular systems.

        use basis, only: basis_fns
        use determinants, only: det_info
        use excitations, only: excit
        use fciqmc_data, only: pattempt_single
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2
        use hamiltonian_molecular, only: slater_condon1_mol, slater_condon2_mol

        use dSFMT_interface, only: genrand_real2

        type(det_info), intent(in) :: cdet
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        logical :: allowed_excitation
        integer :: ij_sym, ij_spin

        ! 1. Select single or double.
        
        if (genrand_real2() < pattempt_single) then
        
            ! 2a. Select orbital to excite from and orbital to excit into.
            call find_ia_mol(cdet%f, cdet%occ_list, connection%from_orb(1), connection%to_orb(1), allowed_excitation)
            connection%nexcit = 1

            if (allowed_excitation) then
                ! 3a. Probability of generating this excitation.
                pgen = calc_pgen_single_mol_no_renorm(connection%to_orb(1))
            
                ! 4a. Parity of permutation required to line up determinants.
                call find_excitation_permutation1(cdet%f, connection)

                ! 5a. Find the connecting matrix element.
                hmatel = slater_condon1_mol(cdet%occ_list, connection%from_orb(1), connection%to_orb(1), connection%perm)
            else
                ! Forbidden---connection%to_orb(1) is already occupied.
                hmatel = 0.0_p
                pgen = 1.0_p ! Avoid any dangerous division by pgen by returning a sane (but cheap) value.
            end if

        else

            ! 2b. Select orbitals to excite from and orbitals to excite into.
            call choose_ij_mol(cdet%occ_list, connection%from_orb(1), connection%from_orb(2), ij_sym, ij_spin)
            call find_ab_mol(cdet%f, ij_sym, ij_spin, connection%to_orb(1), connection%to_orb(2), allowed_excitation)
            connection%nexcit = 2

            if (allowed_excitation) then
                ! 3b. Probability of generating this excitation.
                pgen = calc_pgen_double_mol_no_renorm(connection%to_orb(1), connection%to_orb(2), ij_spin)

                ! 4b. Parity of permutation required to line up determinants.
                ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                call find_excitation_permutation2(cdet%f, connection)

                ! 5b. Find the connecting matrix element.
                hmatel = slater_condon2_mol(connection%from_orb(1), connection%from_orb(2), &
                                            connection%to_orb(1), connection%to_orb(2), connection%perm)
            else
                ! Forbidden---connection%to_orb(2) is already occupied.
                hmatel = 0.0_p
                pgen = 1.0_p ! Avoid any dangerous division by pgen by returning a sane (but cheap) value.
            end if

        end if

    end subroutine gen_excit_mol_no_renorm

!--- Select random orbitals involved in a valid single excitation ---

    subroutine choose_ia_mol(f, occ_list, i, a)

        ! Randomly choose a single excitation, i->a, of a determinant for
        ! molecular systems.

        ! In:
        !    f: bit string representation of the Slater determinant from which
        !    an electron is excited.
        !    occ_list: integer list of occupied spin-orbitals in the determinant.
        ! Out:

        use basis, only: basis_length
        use system, only: nel

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: occ_list(nel)
        integer, intent(out) :: i, a

        logical :: allowed_excitation
        
        allowed_excitation = .false.
        do while (.not.allowed_excitation)
            call find_ia_mol(f, occ_list, i, a, allowed_excitation)
        end do

    end subroutine choose_ia_mol

!--- Select random orbitals involved in a valid double excitation ---

    subroutine choose_ij_mol(occ_list, i, j, ij_sym, ij_spin)

        use basis, only: basis_fns
        use system, only: nel
        use point_group_symmetry, only: cross_product_pg_basis

        use dSFMT_interface, only: genrand_real2

        integer, intent(in) :: occ_list(nel)
        integer, intent(out) :: i, j, ij_sym, ij_spin

        integer :: ind

        ! See comments in choose_ij_k for how the occupied orbitals are indexed
        ! to allow one random number to decide the ij pair.

        ind = int(genrand_real2()*nel*(nel-1)/2) + 1

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.
        ! Note that the indexing scheme for the strictly lower triangular array
        ! assumes j>i.  As occ_list is ordered, this means that we will return
        ! i,j (referring to spin-orbitals) where j>i.  This ordering is
        ! convenient subsequently, e.g. is assumed in the
        ! find_excitation_permutation2 routine.
        j = int(1.50_p + sqrt(2*ind-1.750_p))
        i = ind - ((j-1)*(j-2))/2

        i = occ_list(i)
        j = occ_list(j)

        ij_sym = cross_product_pg_basis(i,j)
        ! ij_spin = -2 (down, down), 0 (up, down or down, up), +2 (up, up)
        ij_spin = basis_fns(i)%Ms + basis_fns(j)%Ms

    end subroutine choose_ij_mol

    subroutine choose_ab_mol(f, sym, spin, a, b)

        use basis, only: basis_length
        use system, only: nel

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: sym, spin
        integer, intent(out) :: a, b

        logical :: allowed_excitation
        
        allowed_excitation = .false.
        do while (.not.allowed_excitation)
            call find_ab_mol(f, sym, spin, a, b, allowed_excitation)
        end do

    end subroutine choose_ab_mol

!--- Select random orbitals in single excitations ---

    subroutine find_ia_mol(f, occ_list, i, a, allowed_excitation)

        use basis, only: basis_length, basis_fns, bit_lookup
        use point_group_symmetry, only: nbasis_sym_spin, sym_spin_basis_fns
        use system, only: nel

        use dSFMT_interface, only: genrand_real2

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: occ_list(nel)
        integer, intent(out) :: i, a
        logical, intent(out) :: allowed_excitation

        integer :: ims, isym, ind

        ! Select an occupied orbital at random.
        i = occ_list(int(genrand_real2()*nel)+1)

        ! Conserve symmetry (spatial and spin) in selecting a.
        ims = (basis_fns(i)%Ms+3)/2 
        isym = basis_fns(i)%sym
        ind = int(nbasis_sym_spin(ims,isym)*genrand_real2())+1
        if (nbasis_sym_spin(ims,isym) == 0) then
            ! No orbitals with the correct symmetry.
            allowed_excitation = .false.
        else
            a = sym_spin_basis_fns(ind,ims,isym)
            ! Is a already occupied in the determinant f?  If so, the excitation is
            ! not permitted.
            allowed_excitation = .not.btest(f(bit_lookup(2,a)), bit_lookup(1,a))
        end if

    end subroutine find_ia_mol

!--- Select random orbitals in double excitations ---

    subroutine find_ab_mol(f, sym, spin, a, b, allowed_excitation)

        use basis, only: basis_length, basis_fns, bit_lookup, nbasis
        use point_group_symmetry, only: nbasis_sym_spin, sym_spin_basis_fns, cross_product_pg_sym
        use system, only: nel

        use dSFMT_interface, only: genrand_real2

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: sym, spin
        integer, intent(out) :: a, b
        logical, intent(out) :: allowed_excitation

        integer :: ims, isym, ind, shift, na, fac

        ! Select a virtual orbital at random.
        ! Ensure spin can be conserved...
        select case(spin)
        case(-2)
            ! a must be down.
            fac = 2
            shift = 1
            na = nbasis/2
        case(0)
            ! a must be down or we can arbitrarily require it to be down.
            fac = 1
            shift = 0
            na = nbasis
        case(2)
            ! a must be up.
            fac = 2
            shift = 0
            na = nbasis/2
        end select
        do
            ! We assume that the user is not so crazy that he/she is
            ! running a calculation where there exists no virtual 
            ! orbitals of a given spin.
            ! random integer between 1 and # possible a orbitals.
            a = int(genrand_real2()*na) + 1
            ! convert to down orbital (ie odd integer between 1 and
            ! nbasis-1) or up orbital (ie even integer between 2 and nbasis)
            ! or to any orbital.
            a = fac*a-shift
            ! If a is unoccupied, then have found first orbital to excite into.
            if (.not.btest(f(bit_lookup(2,a)), bit_lookup(1,a))) exit
        end do

        ! Conserve symmetry (spatial and spin) in selecting the fourth orbital.
        ! Ms_i + Ms_j = Ms_a + Ms_b (Ms_i = -1,+1)
        ! => Ms_b = Ms_i + Ms_j - Ms_a
        ims = (spin-basis_fns(a)%Ms+3)/2
        ! sym_i x sym_j x sym_a = sym_b
        ! (at least for Abelian point groups)
        isym = cross_product_pg_sym(sym, basis_fns(a)%sym)

        if (nbasis_sym_spin(ims,isym) == 0) then
            ! No orbitals with the correct symmetry.
            allowed_excitation = .false.
        else if (spin /= 0 .and. isym == basis_fns(a)%sym .and. nbasis_sym_spin(ims,isym) == 1) then
            allowed_excitation = .false.
        else
            do
                ind = int(nbasis_sym_spin(ims,isym)*genrand_real2())+1
                b = sym_spin_basis_fns(ind,ims,isym)
                if (b /= a) exit
            end do

            ! Is b already occupied in the determinant f?  If so, the excitation is
            ! not permitted.
            allowed_excitation = .not.btest(f(bit_lookup(2,b)), bit_lookup(1,b))

            ! It is useful to return a,b ordered (e.g. for the find_excitation_permutation2 routine).
            if (a > b) then
                ind = a
                a = b
                b = ind
            end if
        end if

    end subroutine find_ab_mol

!--- Excitation generation probabilities ---

    function calc_pgen_single_mol_no_renorm(a) result(pgen)

        use basis, only: basis_fns
        use fciqmc_data, only: pattempt_single
        use system, only: nel
        use point_group_symmetry, only: nbasis_sym_spin

        real(p) :: pgen
        integer, intent(in) :: a

        integer :: ims, isym

        ! We explicitly reject excitations i->a where a is already
        ! occupied, so the generation probability, pgen, is simple:
        !   pgen = p_single p(i) p(a|i)
        !        = p_single 1/nel 1/nbasis_sym_spin(ms_i, sym_i)

        ims = (basis_fns(a)%Ms+3)/2
        isym = basis_fns(a)%sym
        pgen = pattempt_single/(nel*nbasis_sym_spin(ims,isym))

    end function calc_pgen_single_mol_no_renorm

    function calc_pgen_double_mol_no_renorm(a, b, spin) result(pgen)

        ! WARNING: We assume that the excitation is actually valid. 
        ! This routine does *not* calculate the correct probability that
        ! a forbidden excitation (e.g. due to a or b actually being occupied) is
        ! generated.  The correct way to handle those excitations is to set
        ! hmatel to be zero.  The generation probabilites of allowed excitations
        ! correctly take into account events where a/b is occupied.  This
        ! problem arises because if a or b are actually occpied, then p(a|ijb)
        ! or p(b|ijb) = 0.  We do not handle such cases.

        use basis, only: basis_fns
        use fciqmc_data, only: pattempt_double
        use system, only: nel, nvirt, nvirt_alpha, nvirt_beta
        use point_group_symmetry, only: nbasis_sym_spin

        real(p) :: pgen
        integer, intent(in) :: a, b, spin

        integer :: imsa, isyma, imsb, isymb, n_aij
        real(p) :: p_bija, p_aijb

        ! We explicitly reject excitations i,j->a,b where b is already
        ! occupied, so the generation probability, pgen, is simple:
        ! pgen = p_double p(i,j) [ p(a|i,j) p(b|i,j,a) + p(b|i,j) p(a|i,j,b) ]
        !
        ! p(i,j) = 1/binom(nel,2) = 2/(nel*(nel-1))
        ! p(a|i,j) = | 1/(nbasis-nel) if i,j are (up,down) or (down,up)
        !            | 1/(nbasis/2-nalpha) if i,j are (up,up)
        !            | 1/(nbasis/2-nbeta) if i,j are (down,down)
        ! p(b|i,j,a) = 1/nbasis_sym_spin(ms_b, sym_b), where ms_b and sym_b are
        ! such that spin and spatial symmetry are conserved
        ! p(b|i,j) = p(a|i,j) by symmetry.
        ! p(a|i,j,b) = 1/nbasis_sym_spin(ms_a, sym_a), where ms_a and sym_a are
        ! such that spin and spatial symmetry are conserved.
        ! n.b. p(a|i,j,b) =/= p(b|i,j,a) in general.

        select case(spin)
        case(-2)
            n_aij = nvirt_beta
        case(0)
            n_aij = nvirt
        case(2)
            n_aij = nvirt_alpha
        end select

        imsa = (basis_fns(a)%ms+3)/2
        isyma = basis_fns(a)%sym
        imsb = (basis_fns(b)%ms+3)/2
        isymb = basis_fns(b)%sym

        if (isyma == isymb .and. imsa == imsb) then
            ! b cannot be the same as a.
            p_aijb = 1.0_p/(nbasis_sym_spin(imsa, isyma)-1)
            p_bija = 1.0_p/(nbasis_sym_spin(imsb, isymb)-1)
        else
            p_aijb = 1.0_p/nbasis_sym_spin(imsa, isyma)
            p_bija = 1.0_p/nbasis_sym_spin(imsb, isymb)
        end if

        pgen = 2*pattempt_double/(nel*(nel-1)*n_aij)*(p_bija+p_aijb)

    end function calc_pgen_double_mol_no_renorm
            
end module spawning_mol_system
