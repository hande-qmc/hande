module excitations

! Module for generating random excitations.

use const

implicit none

contains

    pure function calc_pgen_hub_k(ij_sym, ispin, f, unocc_alpha, unocc_beta) result(pgen)

        use basis, only: basis_length, bit_lookup, nbasis
        use system, only: nvirt_alpha, nvirt_beta, nalpha, nbeta, nel, nsites
        use symmetry, only: inv_sym, sym_table

        real(dp) :: pgen
        integer, intent(in) :: ij_sym, ispin
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

        select case(ispin)
        case(-1)
            ! exciting from beta, beta orbitals.
            do a = 1, nvirt_beta
                b = 2*sym_table(ij_sym, unocc_beta(a)/2)
                if (a == b) then
                    ! Can't have the a double excitation into the *same* spin-orbital.
                    forbidden_excitations = forbidden_excitations + 1
                else
                    b_pos = bit_lookup(1,b)
                    b_el = bit_lookup(2,b)
                    if (btest(f(b_el), b_pos)) forbidden_excitations = forbidden_excitations + 1
                end if
            end do
            mallow = nsites - nalpha
        case(0)
            ! exciting from alpha, beta orbitals.
            do a = 1, nvirt_alpha
                b = 2*sym_table(ij_sym, unocc_beta(a)/2) ! b is a beta orbital.
                b_pos = bit_lookup(1,b)
                b_el = bit_lookup(2,b)
                if (btest(f(b_el), b_pos)) forbidden_excitations = forbidden_excitations + 1
            end do
            do a = 1, nvirt_beta
                b = 2*sym_table(ij_sym, (unocc_beta(a)+1)/2) + 1 ! b is a alpha orbital.
                b_pos = bit_lookup(1,b)
                b_el = bit_lookup(2,b)
                if (btest(f(b_el), b_pos)) forbidden_excitations = forbidden_excitations + 1
            end do
            mallow = nsites - nbeta
        case(1)
            ! exciting from alpha, alpha orbitals.
            do a = 1, nvirt_alpha
                b = 2*sym_table(ij_sym, (unocc_alpha(a)+1)/2) + 1
                if (a == b) then
                    ! Can't have the a double excitation into the *same* spin-orbital.
                    forbidden_excitations = forbidden_excitations + 1
                else
                    b_pos = bit_lookup(1,b)
                    b_el = bit_lookup(2,b)
                    if (btest(f(b_el), b_pos)) forbidden_excitations = forbidden_excitations + 1
                end if
            end do
            mallow = nbasis - nel
        end select

        ! Common factor.
        ! N(N-1)/2 is the number of ways of choosing i,j.
        ! An additional factor of two arises from the symmetry in the order of
        ! choosing a and b.
        pgen = 4.0_dp/(nel*(nel-1)*(mallow - forbidden_excitations))

    end function calc_pgen_hub_k

    pure subroutine choose_ij(occ_list, i ,j, isym)

        use symmetry, only: sym_table
        use system, only: nel

        integer, intent(in) :: occ_list(nel)
        integer, intent(out) :: i,j, isym
        integer :: ind

        ind = 5 ! replace with int(random number*nel*(nel-1)/2)+1.

        ! We use a triangular indexing scheme to compress 2 electron indices
        ! into 1.
        ! For i/=j and (for an arbitrary choice of i>j), a 1D index of 
        ! a strictly lower triangular array is:
        !   p = (i-1)(i-2)/2 + j,   where 1<=j<i and 1<=p<=n(n-1)/2
        ! This maps the indexing scheme as:
        !    .                  .
        !   2,1  .              1  .
        !   3,1 3,2  .      to  2  3  .
        !   4,1 4,2 4,3  .      3  4  5  .
        ! We want to do the reverse process in order to pick 2 electron labels
        ! from one random number.
        ! Consider the case where j=1.  i can trivially be determined from the
        ! quadratic equation:
        !   i = 3/2 + \sqrt(2p-1.75)
        ! As j<i and, for a fixed i, p increases monotonically with j, the
        ! integer part of i given by the above equation can never exceed the
        ! correct value.  Hence i can be found for arbitrary j by taking the
        ! floor of the preceeding equation.  j follows trivially.
        !
        ! See (for lower triangular arrays rather than strictly lower):
        ! Decoding the sequential indices of a (lower) triangular array
        ! SIS-75-1783,  S Rifkin, CERN report (CERN-DD-75-7).

        ! This might seem odd, but it enables us to pick the (i,j) pair to
        ! excite with half the calls to the random number generator, which
        ! represents a substantial saving. :-)

        i = int(1.50_dp + sqrt(2*ind-1.750_dp))
        j = ind - ((i-1)*(i-2))/2

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.

        i = occ_list(i)
        j = occ_list(j)

        ! Symmetry info is a simple lookup...

        isym = sym_table((i+1)/2,(j+1)/2)

    end subroutine choose_ij

end module excitations
