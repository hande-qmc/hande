module excitations

! Module for generating random excitations.

use const

implicit none

contains

    pure function calc_pgen_hub_k(ij_sym, f, unocc_alpha, unocc_beta) result(pgen)

        ! Calculate the generation probability of a given excitation for the
        ! Hubbard model in momentum space.  The Hubbard model is a special case
        ! as it is a two-band system and so the generation probability is
        ! independent of the virtual spin-orbitals into which electrons are
        ! excited and depends only upon the spin-orbitals from which we excite.
        !
        ! Note that all the information required for input should be available
        ! during the FCIQMC algorithm and should not be needed to be calculated.
        !
        ! Further, we assume only allowed excitations are generated.
        !
        ! In:
        !    ij_sym: symmetry spanned by the (i,j) combination of occupied
        !        spin-orbitals from which electrons are excited.
        !    f: bit string representation of the determinant we're exciting
        !        from.
        !    unocc_alpha, unocc_beta: integer list of the unoccupied alpha and
        !        beta (respectively) spin-orbitals.
        ! Returns:
        !    pgen: the generation probability of the excitation.  See notes in
        !        spawning.

        use basis, only: basis_length, bit_lookup, nbasis
        use system, only: nvirt_alpha, nvirt_beta, nalpha, nbeta, nel, nsites
        use symmetry, only: inv_sym, sym_table

        real(dp) :: pgen
        integer, intent(in) :: ij_sym
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: unocc_alpha(nvirt_alpha), unocc_beta(nvirt_beta)

        integer :: forbidden_excitations, a, b, ab_sym, b_pos, b_el, mallow

        forbidden_excitations = 0

        ! The excitation i,j -> a,b is only allowed if k_i + k_j - k_a - k_b = 0
        ! (up to a reciprocal lattice vector).  We store k_i + k_j as ij_sym, so
        ! k_a + k_b must be the inverse of this.
        ab_sym = inv_sym(ij_sym)

        ! Fixing a means b is fixed in order to meet the symmetry requirements.
        ! We count the number of a orbitals which cannot be excited into due to
        ! the corresponding b orbital not being available.
        ! The Hubbard model is a 2-band system, which makes this pleasingly
        ! easy. :-)

        ! pgen then follows immediately.

        ! exciting from alpha, beta orbitals.
        do a = 1, nvirt_alpha
            b = 2*sym_table(ab_sym, unocc_beta(a)/2) ! b is a beta orbital.
            b_pos = bit_lookup(1,b)
            b_el = bit_lookup(2,b)
            if (btest(f(b_el), b_pos)) forbidden_excitations = forbidden_excitations + 1
        end do
        do a = 1, nvirt_beta
            b = 2*sym_table(ab_sym, (unocc_beta(a)+1)/2) + 1 ! b is an alpha orbital.
            b_pos = bit_lookup(1,b)
            b_el = bit_lookup(2,b)
            if (btest(f(b_el), b_pos)) forbidden_excitations = forbidden_excitations + 1
        end do
        mallow = nsites - nbeta

        ! N^2/4 is the number of ways of choosing i,j.
        ! An additional factor of two arises from the symmetry in the order of
        ! choosing a and b.
        pgen = 8.0_dp/(nel*nel*(mallow - forbidden_excitations))

    end function calc_pgen_hub_k

end module excitations
