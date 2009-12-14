module hubbard_k

! Momentum space formulation of the Hubbard model.

implicit none

contains

    pure function get_one_e_int_k(phi1, phi2) result(one_e_int)

        ! In:
        !    phi1: index of a momentum-space basis function.
        !    phi2: index of a momentum-space basis function.
        ! Returns:
        !    <phi1 | T | phi2> where T is the kinetic energy operator.

        use basis

        real(dp) :: one_e_int
        integer, intent(in) :: phi1, phi2

        ! T is diagonal in the basis of momentum-space functions.
        if (phi1 == phi2) then
            one_e_int = basis_fns(phi1)%kinetic
        else
            one_e_int = 0.0_dp
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

        real(dp) :: two_e_int
        integer, intent(in) :: phi1, phi2, phi3, phi4

        ! <phi1 phi2 || phi3 phi4>
        two_e_int = 0.0_dp

        ! The integral < k_1 k_2 | U | k_3 k_4 > = U/N \delta_{k_1+k2,k_3+k_4}
        ! where the delta function requires crystal momentum is conserved up to
        ! a reciprocal lattice vector.

        ! <phi1 phi2 | phi3 phi4>
        if (spin_symmetry(phi1, phi3) .and. spin_symmetry(phi2, phi4)) then
            if (momentum_conserved(phi1, phi2, phi3, phi4)) then
                two_e_int = hubu/nsites
            end if
        end if

        ! <phi1 phi2 | phi4 phi3>
        if (spin_symmetry(phi1, phi4) .and. spin_symmetry(phi2, phi3)) then
            if (momentum_conserved(phi1, phi2, phi4, phi3)) then
                two_e_int = two_e_int - hubu/nsites
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

        logical :: conserved
        integer, intent(in) :: i, j, k, l
        integer :: delta_k(ndim)

        ! k_i + k_j - k_k -k_l in units of the reciprocal lattice vectors of the
        ! crystal cell.
        delta_k = basis_fns(i)%l + basis_fns(j)%l - basis_fns(k)%l - basis_fns(l)%l

        conserved = is_reciprocal_lattice_vector(delta_k)

    end function momentum_conserved

end module hubbard_k
