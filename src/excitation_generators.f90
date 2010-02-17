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
        ! 2/N(N+1) is the number of ways of choosing i,j.
        ! An additional factor of two arises from the symmetry in the order of
        ! choosing a and b.
        pgen = 4.0_dp/(nel*(nel-1)*(mallow - forbidden_excitations))

    end function calc_pgen_hub_k

end module excitations
