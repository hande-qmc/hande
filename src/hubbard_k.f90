module hubbard_k

! Momentum space formulation of the Hubbard model.

use const

implicit none

contains

    pure function get_one_e_int_k(phi1, phi2) result(one_e_int)

        ! In:
        !    phi1: index of a momentum-space basis function.
        !    phi2: index of a momentum-space basis function.
        ! Returns:
        !    <phi1 | T | phi2> where T is the kinetic energy operator.

        use basis

        real(p) :: one_e_int
        integer, intent(in) :: phi1, phi2

        ! T is diagonal in the basis of momentum-space functions.
        if (phi1 == phi2) then
            one_e_int = basis_fns(phi1)%sp_eigv
        else
            one_e_int = 0.0_p
        end if

    end function get_one_e_int_k

    elemental function get_two_e_int_k(phi1, phi2, phi3, phi4) result(two_e_int)

        ! In:
        !    phi1: index of a momentum-space basis function.
        !    phi2: index of a momentum-space basis function.
        !    phi3: index of a momentum-space basis function.
        !    phi4: index of a momentum-space basis function.
        ! Returns:
        !    The anti-symmetrized integral <phi1 phi2 || phi3 phi4>.

        use basis
        use system

        real(p) :: two_e_int
        integer, intent(in) :: phi1, phi2, phi3, phi4

        logical :: s12, s34

        ! The integral < k_1 k_2 | U | k_3 k_4 > = U/N \delta_{k_1+k2,k_3+k_4}
        ! where the delta function requires crystal momentum is conserved up to
        ! a reciprocal lattice vector.

        s12 = spin_symmetry(phi1, phi2)
        s34 = spin_symmetry(phi3, phi4)

        ! If phi1, phi2, phi3, phi4 all of the same spin, then < phi1 phi2 || phi3 phi4 > is zero
        ! as the Coulomb and exchange integrals will exactly cancel.

        if (.not.s12 .and. s12.eqv.s34 .and.  momentum_conserved(phi1, phi2, phi3, phi4)) then
            ! phi1, phi2 are of opposite spins and so are phi3 and phi4 and
            ! crystal momentum is conserved.
            ! Either the Coulomb integral or the exchange integral is
            ! non-zero.
            if (spin_symmetry(phi1, phi3)) then
                two_e_int = sys_global%hubbard%coulomb_k
            else
                two_e_int = -sys_global%hubbard%coulomb_k
            end if
        else
            two_e_int = 0.0_p
        end if

    end function get_two_e_int_k

    elemental function momentum_conserved(i, j, k, l) result(conserved)

        use basis

        ! In:
        !    i: index of a momentum-space basis function.
        !    j: index of a momentum-space basis function.
        !    k: index of a momentum-space basis function.
        !    l: index of a momentum-space basis function.
        ! Returns:
        !    True if crystal momentum is conserved in the integral <k_i k_j | U | k_k k_l>
        !    i.e. if k_i + k_j - k_k -k_l = 0 up to a reciprocal lattice vector.

        use momentum_symmetry

        logical :: conserved
        integer, intent(in) :: i, j, k, l
        integer :: delta_k

        delta_k = sym_table((i+1)/2,(j+1)/2)
        delta_k = sym_table(delta_k,inv_sym((k+1)/2))

        conserved = delta_k == (l+1)/2

    end function momentum_conserved

end module hubbard_k
