module hamiltonian_heisenberg

! Module for evaluating Hamiltonian matrix elements for the Heisenberg model.

use const

implicit none

contains

    pure function get_hmatel_heisenberg(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f1, f2: bit string representation of the basis functions
        !        D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two basis functions,
        !    < D1 | H | D2 >, where the basis functions are formed from
        !    the list of spins on each site.

        ! Used in the Heisenberg model only.

        use excitations, only: excit, get_excitation
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(sys%basis%string_len), f2(sys%basis%string_len)
        type(excit) :: excitation

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! We assume Ms is conserved (ie has already been checked for).
        excitation = get_excitation(sys%nel, sys%basis, f1, f2)
        ! Connected determinants can differ by (at most) 1 spin orbitals.
        ! Space group symmetry not currently implemented.

        select case(excitation%nexcit)
        case(0)

            ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
            if (abs(sys%heisenberg%staggered_magnetic_field) > depsilon) then
                hmatel = diagonal_element_heisenberg_staggered(sys, f1)
            else
                hmatel = diagonal_element_heisenberg(sys, f1)
            end if

        case(1)

            hmatel = offdiagonal_element_heisenberg(sys, excitation%from_orb(1), excitation%to_orb(1))

        case default

            hmatel = 0.0_p

        end select

    end function get_hmatel_heisenberg

    pure function diagonal_element_heisenberg(sys, f) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of the basis function.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Heisenberg Model.
        !
        ! Includes an uniform external field, if one is added.

        use calc, only: ms_in
        use real_lattice, only: connected_orbs
        use bit_utils, only: count_set_bits
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer(i0) :: f_not(sys%basis%string_len), g(sys%basis%string_len)
        integer :: ipos, i, basis_find, counter

        counter = 0

        ! Count the number of 0-1 type bonds
        f_not = not(f)
        do i = 1, sys%basis%string_len
            do ipos = 0, i0_end
                if (btest(f(i), ipos)) then
                    basis_find = sys%basis%basis_lookup(ipos, i)
                    g = iand(f_not, connected_orbs(:,basis_find))
                    counter = counter + sum(count_set_bits(g))
                end if
            end do
        end do

        ! Contribution to Hamiltonian from spin interactions:
        ! For a lattice the number of bonds is stored as nbonds.
        ! Bonds of type 0-0 or 1-1 will give a contribution of -J/4 to the matrix
        ! element.  0-1 bonds will give +J/4 contribution.
        ! The above loop counts the number of 0-1 bonds, so the remaining number
        ! of 0-0 or 1-1 bonds will be (nbonds-counter)
        ! So we have a contribution of -(J/4)*counter from the 0-1 bonds and
        ! +(J/4)*(nbonds-counter) from the 1-1 and 0-0 bonds, so in total
        ! the matrix element is...
        hmatel = -(sys%heisenberg%J/4)*(sys%heisenberg%nbonds-2*counter)

        ! Contribution to Hamiltonian from external field:
        ! Each spin up (bit set) up gives a contribution of -magnetic_field/2. Each spin down
        ! gives a contribution of +magnetic_field/2. There are bits_set spins up, and
        ! (nsites - bits_set) spins down, so the total contribution is
        ! -(magnetic_field/2)*(2*bits_set-nsites) = -h_field*ms_in/2
        hmatel = hmatel - sys%heisenberg%magnetic_field*ms_in/2

    end function diagonal_element_heisenberg

    pure function diagonal_element_heisenberg_staggered(sys, f) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of the basis function.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Heisenberg Model.
        !
        ! This function is for diagonal elements for the Hamiltonian which includes
        ! a staggered magnetization term.

        use determinants, only: lattice_mask
        use real_lattice, only: connected_orbs
        use bit_utils, only: count_set_bits
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer(i0) :: f_not(sys%basis%string_len), f_mask(sys%basis%string_len), &
                       g(sys%basis%string_len)
        integer :: ipos, i, basis_find, counter, sublattice1_up_spins

        counter = 0

        ! For non-frustrtated lattices where we want to add a staggered magnetization
        ! term, we need to calculate how many up spins are on each of the sublattices,
        ! so we first consider one particular sublattice and put 0's at all other sites:
        f_mask = iand(f, lattice_mask)
        sublattice1_up_spins = sum(count_set_bits(f_mask))

        ! Count the number of 0-1 type bonds
        f_not = not(f)
        do i = 1, sys%basis%string_len
            do ipos = 0, i0_end
                if (btest(f(i), ipos)) then
                    basis_find = sys%basis%basis_lookup(ipos, i)
                    g = iand(f_not, connected_orbs(:,basis_find))
                    counter = counter + sum(count_set_bits(g))
                end if
            end do
        end do

        ! Contribution to Hamiltonian from spin interactions:
        ! For a lattice the number of bonds is stored as nbonds.
        ! Bonds of type 0-0 or 1-1 will give a contribution of -J/4 to the matrix
        ! element.  0-1 bonds will give +J/4 contribution.
        ! The above loop counts the number of 0-1 bonds, so the remaining number
        ! of 0-0 or 1-1 bonds will be (nbonds-counter)
        ! So we have a contribution of -(J/4)*counter from the 0-1 bonds and
        ! +(J/4)*(nbonds-counter) from the 1-1 and 0-0 bonds, so in total
        ! the matrix element is...
        hmatel = -(sys%heisenberg%J/4)*(sys%heisenberg%nbonds-2*counter)

        ! Contibution to Hamiltonian from staggered field term.
        ! Split the lattice into a (+) sublattice and a (-) sublattice.
        ! For every up spin on a (+) site, add one half to the hmatel. (add staggered_magnetic_field * 0.5)
        ! For every down spin on a (+) site, minus one half.
        ! For every up spin on a (-) site, minus one half.
        ! For every down spin on a (-) site, add one half.
        ! sublattice1_up_spins gives the number of up spins on the (+) lattice, so there are
        ! (nel-sublattice1_up_spins) up spins on the (-) lattice, and a total of (nsites/2)
        ! sites on each sublattice. Putting all this together gives the following:
        hmatel = hmatel - sys%heisenberg%staggered_magnetic_field*(4*sublattice1_up_spins - 2*sys%nel)/2

    end function diagonal_element_heisenberg_staggered

    pure function offdiagonal_element_heisenberg(sys, i, a) result(hmatel)

        ! In:
        !    sys: system to be studied.
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

        use real_lattice, only: connected_orbs
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, a

        integer :: ipos, iel

        ipos = sys%basis%bit_lookup(1,a)
        iel = sys%basis%bit_lookup(2,a)

        if (btest(connected_orbs(iel,i), ipos)) then
            ! If the two sites connected and of opposite spin, matrix element is -J/2.
            ! As get_excitation finds where a 'set bit' has moved from and
            ! to, the latter condition is already met.
            hmatel = -sys%heisenberg%J/2
        else
            hmatel = 0.0_p
        end if

    end function offdiagonal_element_heisenberg

end module hamiltonian_heisenberg
