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
            one_e_int = basis_fns(phi1)%kinetic
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

        ! <phi1 phi2 || phi3 phi4>
        two_e_int = 0.0_p

        ! The integral < k_1 k_2 | U | k_3 k_4 > = U/N \delta_{k_1+k2,k_3+k_4}
        ! where the delta function requires crystal momentum is conserved up to
        ! a reciprocal lattice vector.

        ! <phi1 phi2 | phi3 phi4>
        if (spin_symmetry(phi1, phi3) .and. spin_symmetry(phi2, phi4)) then
            if (momentum_conserved(phi1, phi2, phi3, phi4)) then
                two_e_int = hub_k_coulomb
            end if
        end if

        ! <phi1 phi2 | phi4 phi3>
        if (spin_symmetry(phi1, phi4) .and. spin_symmetry(phi2, phi3)) then
            if (momentum_conserved(phi1, phi2, phi4, phi3)) then
                two_e_int = two_e_int - hub_k_coulomb
            end if
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

        use symmetry

        logical :: conserved
        integer, intent(in) :: i, j, k, l
        integer :: delta_k

        delta_k = sym_table((i+1)/2,(j+1)/2)
        delta_k = sym_table(delta_k,inv_sym((k+1)/2))
        delta_k = sym_table(delta_k,inv_sym((l+1)/2))

        conserved = delta_k == 1

    end function momentum_conserved

end module hubbard_k
