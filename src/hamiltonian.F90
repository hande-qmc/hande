module hamiltonian

use const
use basis
use determinants
use parallel

implicit none

contains

    elemental function get_hmatel_dets(d1, d2) result(hmatel)

        ! In:
        !    d1, d2: integer labels of two determinants, as stored in the
        !            dets array.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants, 
        !    < D1 | H | D2 >.

        ! This is just a wrapper around get_hmatel (which is itself a wrapper
        ! around system-specific functions) but is handy for computing matrix
        ! elements (slowly!) when we have the entire Hilbert space of
        ! determinants stored in dets_list.

        real(p) :: hmatel
        integer, intent(in) :: d1, d2

        hmatel = get_hmatel(dets_list(:,d1), dets_list(:,d2))

    end function get_hmatel_dets

    pure function get_hmatel(f1, f2) result(hmatel)

        ! In:
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants, 
        !    < D1 | H | D2 >.

        ! This is just a wrapper function around the system specific get_hmatel
        ! functions.

        ! Having separate functions for the different systems might seem
        ! somewhat redundant (a lot of the code in the functions is similar)
        ! but enables us to use only one test for the system type.  A small
        ! efficiency for not much effort. :-)

        use hamiltonian_molecular, only: get_hmatel_mol

        real(p) :: hmatel
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)

        select case(system_type)
        case(hub_k)
            hmatel = get_hmatel_k(f1, f2)
        case(hub_real)
            hmatel = get_hmatel_real(f1, f2)
        case(heisenberg)
            hmatel = get_hmatel_heisenberg(f1, f2)
        case(read_in)
            hmatel = get_hmatel_mol(f1, f2)
        end select

    end function get_hmatel

    pure function get_hmatel_k(f1, f2) result(hmatel)

        ! In:
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants, 
        !    < D1 | H | D2 >, where the determinants are formed from
        !    momentum space basis functions.

        ! Used in the momentum space formulation of the Hubbard model only.

        use excitations, only: excit, get_excitation

        real(p) :: hmatel
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        logical :: non_zero
        type(excit) :: excitation

        hmatel = 0.0_p
        non_zero = .false.

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! Assume D1 and D2 are of the same symmetry.  Namely:

            ! We assume Ms is conserved (ie has already been checked for).

            ! In the momentum space description the overall crystal 
            ! momentum must be conserved up to a reciprocal lattice
            ! vector (i.e. satisfy translational symmetry).
            ! We assume this is also already checked.

        excitation = get_excitation(f1,f2)
        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then
            non_zero = .true.
        end if

        if (non_zero) then
            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.
            case(0)

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
                hmatel = slater_condon0_hub_k(f1)

!            case(1)

                ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

                ! Single excitations are not connected in the momentum space
                ! basis.

                ! One electron operator
                ! The kinetic operator is diagonal in the momentum space basis.

                ! Two electron operator
                ! < ij | aj > = 0 only if crystal momentum is conserved up to
                ! a reciprocal lattice vector.
                ! As k_i /= k_j, this cannot be met.

            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >

                ! Two electron operator
                hmatel = slater_condon2_hub_k(excitation%from_orb(1), excitation%from_orb(2), &
                                            & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)

            end select
        end if

    end function get_hmatel_k

    pure function get_hmatel_real(f1, f2) result(hmatel)

        ! In:
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants, 
        !    < D1 | H | D2 >, where the determinants are formed from
        !    real space basis functions.

        ! Used in the real space formulation of the Hubbard model only.

        use excitations, only: excit, get_excitation

        real(p) :: hmatel
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        logical :: non_zero
        type(excit) :: excitation

        hmatel = 0.0_p
        non_zero = .false.

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! We assume Ms is conserved (ie has already been checked for).
        excitation = get_excitation(f1, f2)
        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then
            ! Space group symmetry not currently implemented.
            non_zero = .true.
        end if

        ! Matrix elements in the real space formulation are quite simple.

        ! 1. < i | T | i > = 0
        !    Thus the one-electron terms only occur between single excitation
        !    matrix elements.
        ! 2. < m,s1 n,s2 | U | p,s1 q,s2 > = U \delta_{m,n} \delta_{m,p} \delta_{m,q} , s1/=s2
        !    Thus the Coulomb integrals that occur in < D | H | D_i^a > and
        !    < D | H | D_{ij}^{ab} > are zero.

        if (non_zero) then
            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.
            case(0)

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
                hmatel = slater_condon0_hub_real(f1)
    
            case(1)

                hmatel = slater_condon1_hub_real(excitation%from_orb(1), excitation%to_orb(1), excitation%perm)

!            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >
                !                         = 0 within the real space formulation
                !                             of the Hubbard model.

            end select
        end if

    end function get_hmatel_real

    pure function get_hmatel_heisenberg(f1, f2) result(hmatel)

        ! In:
        !    f1, f2: bit string representation of the basis functions 
        !        D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two basis functions, 
        !    < D1 | H | D2 >, where the basis functions are formed from
        !    the list of spins on each site.

        ! Used in the Heisenberg model only.

        use excitations, only: excit, get_excitation

        real(p) :: hmatel
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        type(excit) :: excitation

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! We assume Ms is conserved (ie has already been checked for).
        excitation = get_excitation(f1, f2)
        ! Connected determinants can differ by (at most) 1 spin orbitals.
        ! Space group symmetry not currently implemented.

        select case(excitation%nexcit)
        case(0)

            ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
            if (abs(staggered_magnetic_field) > depsilon) then
                hmatel = diagonal_element_heisenberg_staggered(f1)
            else
                hmatel = diagonal_element_heisenberg(f1)
            end if

        case(1)

            hmatel = offdiagonal_element_heisenberg(excitation%from_orb(1), excitation%to_orb(1))

        case default

            hmatel = 0.0_p

        end select

    end function get_hmatel_heisenberg

    pure function slater_condon0_hub_k(f) result(hmatel)
        
        ! In:
        !    f: bit string representation of the Slater determinant.
        !    occ_list: integer list of occupied spin-orbitals in a determinant, D_i.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Hubbard model in momentum space.

        use system, only: nalpha, nbeta, hub_k_coulomb

        real(p) :: hmatel
        integer(i0), intent(in) :: f(basis_length)
        integer :: occ_list(nel)
        integer :: i

        call decode_det(f, occ_list)

        ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >

        ! Two electron operator
        ! 1/2 \sum_i \sum_j < ij || ij >
        ! Some points to note:
        !   a) Crystal momentum is conserved for all < ij | ij > and < ij | ji > integrals 
        !      by defintion.
        !   b) If i,j are of the same spin, then < ij | ij > = < ij | ji > and
        !      so < ij || ij > = 0.
        !   c) < ij | ij > = U/nsites for all i,j.
        !   d) The double sum has 2*nalpha*nbeta terms corresponding to i,j of
        !      different spins.
        !   e) Thus  1/2 \sum_i \sum_j < ij || ij > = nalpha*nbeta*U/nsites.
        hmatel = nalpha*nbeta*hub_k_coulomb

        ! One electron operator
        ! Get directly rather than incur the cost of the if test in get_one_e_int_k.
        do i = 1, nel
            hmatel = hmatel + basis_fns(occ_list(i))%sp_eigv
        end do

    end function slater_condon0_hub_k

    pure function slater_condon2_hub_k(i, j, a, b, perm) result(hmatel)

        ! In:
        !    i,j:  index of the spin-orbital from which an electron is excited in
        !          the reference determinant.
        !    a,b:  index of the spin-orbital into which an electron is excited in
        !          the excited determinant.
        !    perm: true if D and D_i^a are connected by an odd number of
        !          permutations.
        ! Returns:
        !    < D | H | D_ij^ab >, the Hamiltonian matrix element between a 
        !    determinant and a double excitation of it in the momemtum space
        !    formulation of the Hubbard model.

        use hubbard_k, only: get_two_e_int_k

        real(p) :: hmatel
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        hmatel = get_two_e_int_k(i, j, a, b)

        if (perm) hmatel = -hmatel

    end function slater_condon2_hub_k

    pure subroutine slater_condon2_hub_k_excit(f, connection, hmatel)

        ! Generate the matrix element between a determinant and a double
        ! excitation in the momentum space formulation of the Hubbard model.
        ! WARNING: this routine assumes that the excitation is allowed (i.e.
        ! conserves crystal momentum).  It is, however, faster as symmetry
        ! checking is skipped.

        ! In:
        !    f: bit string representation of the Slater determinant, D.
        ! In/Out:
        !    connection: excit type describing the excitation between |D> and
        !    |D_ij^ab>.  On entry, only the from_orb and to_orb fields must be
        !    set.  On exit the from_orb and to_orb fields will be ordered
        !    and the perm field will be set.
        ! Out:
        !    hmatel: < D | H | D_ij^ab >, the Hamiltonian matrix element between a 
        !    determinant and a double excitation of it in the momemtum space
        !    formulation of the Hubbard model.

        use excitations, only: excit, find_excitation_permutation2
        use basis, only: basis_length

        integer(i0), intent(in) :: f(basis_length)
        type(excit), intent(inout) :: connection
        real(p), intent(out) :: hmatel

        integer :: tmp

        ! The permuting algorithm works by lining up the min(i,j) with
        ! min(a,b) and max(i,j) with max(a,b) and hence we can find out
        ! whether the Coulomb or exchange integral is non-zero.  
        ! Thus (i,j) and (a,b) must be ordered.
        if (connection%from_orb(1) > connection%from_orb(2)) then
            ! Swap.
            tmp = connection%from_orb(1)
            connection%from_orb(1) = connection%from_orb(2)
            connection%from_orb(2) = tmp
        end if
        if (connection%to_orb(1) > connection%to_orb(2)) then
            ! Swap
            tmp = connection%to_orb(1)
            connection%to_orb(1) = connection%to_orb(2)
            connection%to_orb(2) = tmp
        end if

        ! a) Sign from value of U as well---U can be negative!
        hmatel = hub_k_coulomb

        ! b) Negative sign from permuting the determinants so that they line
        ! up?
        call find_excitation_permutation2(f, connection)
        if (connection%perm) then
            ! Matrix element gets a -sign from rearranging determinants so
            ! that they maximally line up.
            hmatel = -hmatel
        end if

        ! c) Because the only non-zero excitations are when i is alpha and
        ! j is beta or vice-versa, only the Coulomb integral or the exchange
        ! integral is non-zero.  If it's the exchange
        ! integral, then we obtain an additional minus sign.
        if (mod(connection%from_orb(1)-connection%to_orb(1),2) /= 0) then
            ! (i',a') are (alpha,beta) or (beta,alpha).
            ! Thus it is the exchange integral which contributes to the
            ! connecting matrix element.
            hmatel = -hmatel
        end if

    end subroutine slater_condon2_hub_k_excit

    pure function slater_condon0_hub_real(f) result(hmatel)

        ! In:
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Hubbard model in real space.

        use hubbard_real, only: t_self_images, get_one_e_int_real, get_coulomb_matel_real

        real(p) :: hmatel
        integer(i0), intent(in) :: f(basis_length)
        integer :: root_det(nel)
        integer :: i

        ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
        hmatel = 0.0_p

        ! < i | T | i > = 0 within the real space formulation of the
        ! Hubbard model, unless site i is its own periodic image, in
        ! which case it has a kinetic interaction with its self-image.
        ! This only arises if there is at least one crystal cell vector
        ! which is a unit cell vector.
        if (t_self_images) then
            call decode_det(f, root_det)
            do i = 1, nel
                hmatel = hmatel + get_one_e_int_real(root_det(i), root_det(i))
            end do
        end if

        ! Two electron operator
        hmatel = hmatel + get_coulomb_matel_real(f)

    end function slater_condon0_hub_real

    pure function slater_condon1_hub_real(i, a, perm) result(hmatel)

        ! In:
        !    i: index of the spin-orbital from which an electron is excited in
        !        the reference determinant.
        !    a: index of the spin-orbital into which an electron is excited in
        !        the excited determinant.
        !    perm: true if D and D_i^a are connected by an odd number of
        !        permutations.
        ! Returns:
        !    < D | H | D_i^a >, the Hamiltonian matrix element between a 
        !        determinant and a single excitation of it.

        use hubbard_real, only: get_one_e_int_real

        real(p) :: hmatel
        integer, intent(in) :: i, a
        logical, intent(in) :: perm

        ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

        ! One electron operator
         hmatel = get_one_e_int_real(i, a)

        ! Two electron operator
        ! < D | U | D_i^a > = 0 within the real space formulation of the
        ! Hubbard model.

        if (perm) hmatel = -hmatel

    end function slater_condon1_hub_real

    pure subroutine slater_condon1_hub_real_excit(f, connection, hmatel)

        ! Generate the matrix element between a determinant and a single
        ! excitation in the real space formulation of the Hubbard model.
        ! WARNING: this routine assumes that the excitation is allowed (i.e.
        ! the excitation is from an occupied orbital, i, to an unoccupied
        ! orbital, a, which is connected to i).  By skipping such checks, it is,
        ! however, faster.

        ! In:
        !    f: bit string representation of the Slater determinant, D.
        ! In/Out:
        !    connection: excit type describing the excitation between |D> and
        !    |D_i^a>.  On entry, only the from_orb and to_orb fields must be
        !    set.  On exit the perm field will also be set.
        ! Out:
        !    hmatel: < D | H | D_i^a >, the Hamiltonian matrix element between a 
        !    determinant and a single excitation of it in the real space
        !    formulation of the Hubbard model.

        use excitations, only: excit, find_excitation_permutation1

        integer(i0), intent(in) :: f(basis_length)
        type(excit), intent(inout) :: connection
        real(p), intent(out) :: hmatel

        ! a) Find out permutation required to line up determinants.
        call find_excitation_permutation1(f, connection)

        ! b) The matrix element connected |D> and |D_i^a> is <i|h|a> = -t.
        if (connection%perm) then
            hmatel = hubt
        else
            hmatel = -hubt
        end if

    end subroutine slater_condon1_hub_real_excit
    
    pure function diagonal_element_heisenberg(f) result(hmatel)

        ! In:
        !    f: bit string representation of the basis function.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Heisenberg Model.
        !
        ! Includes an uniform external field, if one is added.
        
        use basis, only: basis_length, basis_lookup
        use calc, only: ms_in
        use hubbard_real, only: connected_orbs
        use bit_utils, only: count_set_bits
        use system, only: nbonds, J_coupling, magnetic_field
        
        real(p) :: hmatel
        integer(i0), intent(in) :: f(basis_length)
        integer(i0) :: f_not(basis_length), g(basis_length)
        integer :: ipos, i, basis_find, counter
        
        counter = 0
        
        ! Count the number of 0-1 type bonds
        f_not = not(f)
        do i = 1, basis_length
            do ipos = 0, i0_end
                if (btest(f(i), ipos)) then
                    basis_find = basis_lookup(ipos, i)
                    g = iand(f_not, connected_orbs(:,basis_find))
                    counter = counter + sum(count_set_bits(g))
                end if
            end do
        end do
        
        ! Contribution to Hamiltonian from spin interactions:
        ! For a lattice the number of bonds is stored as nbonds.
        ! Bonds of type 0-0 or 1-1 will give a contribution of -J_coupling to the matrix
        ! element.  0-1 bonds will give +J_coupling contribution.
        ! The above loop counts the number of 0-1 bonds, so the remaining number
        ! of 0-0 or 1-1 bonds will be (nbonds-counter)
        ! So we have a contribution of -J_coupling*counter from the 0-1 bonds and 
        ! +J_coupling*(nbonds-counter) from the 1-1 and 0-0 bonds, so in total
        ! the matrix element is...
        hmatel = -J_coupling*(nbonds-2*counter)
        
        ! Contribution to Hamiltonian from external field:
        ! Each spin up (bit set) up gives a contirbution of -magnetic_field. Each spin down
        ! gives a contirbution of +magnetic_field. There are bits_set spins up, and 
        ! (nsites - bits_set) spins down, so the total contirbution is
        ! -magnetic_field*(2*bits_set-nsites) = -h_field*ms_in
        hmatel = hmatel - magnetic_field*ms_in

    end function diagonal_element_heisenberg
    
    pure function diagonal_element_heisenberg_staggered(f) result(hmatel)

        ! In:
        !    f: bit string representation of the basis function.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Heisenberg Model.
        !
        ! This function is for diagonal elements for the Hamiltonian which includes
        ! a staggered magnetization term.
        
        use basis, only: basis_length, basis_lookup
        use determinants, only: lattice_mask
        use hubbard_real, only: connected_orbs
        use bit_utils, only: count_set_bits
        use system, only: nbonds, nel, J_coupling, magnetic_field, staggered_magnetic_field
        
        real(p) :: hmatel
        integer(i0), intent(in) :: f(basis_length)
        integer(i0) :: f_not(basis_length), f_mask(basis_length), &
                       g(basis_length)
        integer :: ipos, i, basis_find, counter, sublattice1_up_spins
        
        counter = 0
        
        ! For non-frustrtated lattices where we want to add a staggered magnetization
        ! term, we need to calculate how many up spins are on each of the sublattices,
        ! so we first consider one particular sublattice and put 0's at all other sites:
        f_mask = iand(f, lattice_mask)
        sublattice1_up_spins = sum(count_set_bits(f_mask))
        
        ! Count the number of 0-1 type bonds
        f_not = not(f)
        do i = 1, basis_length
            do ipos = 0, i0_end
                if (btest(f(i), ipos)) then
                    basis_find = basis_lookup(ipos, i)
                    g = iand(f_not, connected_orbs(:,basis_find))
                    counter = counter + sum(count_set_bits(g))
                end if
            end do
        end do
        
        ! Contribution to Hamiltonian from spin interactions:
        ! For a lattice the number of bonds is stored as nbonds.
        ! Bonds of type 0-0 or 1-1 will give a contribution of -J_coupling to the matrix
        ! element.  0-1 bonds will give +J_coupling contribution.
        ! The above loop counts the number of 0-1 bonds, so the remaining number
        ! of 0-0 or 1-1 bonds will be (nbonds-counter)
        ! So we have a contribution of -J_coupling*counter from the 0-1 bonds and 
        ! +J_coupling*(nbonds-counter) from the 1-1 and 0-0 bonds, so in total
        ! the matrix element is...
        hmatel = -J_coupling*(nbonds-2*counter)
        
        ! Contibution to Hamiltonian from staggered field term.
        ! Split the lattice into a (+) sublattice and a (-) sublattice.
        ! For every up spin on a (+) site, add one to the hmatel. (add staggered_magnetic_field * 1)
        ! For every down spin on a (+) site, minus one.
        ! For every up spin on a (-) site, minus one.
        ! For every down spin on a (-) site, add one.
        ! sublattice1_up_spins gives the number of up spins on the (+) lattice, so there are
        ! (nel-sublattice1_up_spins) up spins on the (-) lattice, and a total of (nsites/2)
        ! sites on each sublattice. Putting all this together gives the following:
        hmatel = hmatel - staggered_magnetic_field*(4*sublattice1_up_spins - 2*nel)

    end function diagonal_element_heisenberg_staggered

    pure function offdiagonal_element_heisenberg(i, a) result(hmatel)

        ! In:
        !    i: index of the spin-orbital from which is changed from spin-up to
        !        spin-down in the reference basis function.
        !    a: index of the spin-orbital from which is changed from spin-down to
        !        spin-up in the reference basis function.
        ! Returns:
        !    < i | H | a >, the Hamiltonian matrix element between a basis
        !    function and a single excitation of it in the Heisenberg model.
        !
        ! This applies to all Heisenberg situations, including those with
        ! uniform and staggered external fields applied, since these additions
        ! only alter the diagonal elements.
        
        use basis, only: basis_length, bit_lookup
        use hubbard_real, only: connected_orbs
        use system, only: J_coupling
        
        real(p) :: hmatel
        integer, intent(in) :: i, a

        integer :: ipos, iel

        ipos = bit_lookup(1,a)
        iel = bit_lookup(2,a)

       if (btest(connected_orbs(iel,i), ipos)) then
!            ! If the two sites connected and of opposite spin, matrix element is -2J.
!            ! As get_excitation finds where a 'set bit' has moved from and
!            ! to, the latter condition is already met.
            hmatel = -2*J_coupling
        else
            hmatel = 0.0_p
        end if

    end function offdiagonal_element_heisenberg

end module hamiltonian
