module hubbard_real

! Real space formulation of the Hubbard model.

use const

implicit none

contains

    pure function get_one_e_int_real(phi1, phi2) result(one_e_int)

        ! In:
        !    phi1: index of a real-space basis function.
        !    phi2: index of a real-space basis function.
        ! Returns:
        !    <phi1 | T | phi2> where T is the kinetic energy operator.

        use basis

        real(dp) :: one_e_int
        integer, intent(in) :: phi1, phi2

        one_e_int = 0.0_dp

    end function get_one_e_int_real

    elemental function get_two_e_int_real(phi1, phi2, phi3, phi4) result(two_e_int)

        ! In:
        !    phi1: index of a real-space basis function.
        !    phi2: index of a real-space basis function.
        !    phi3: index of a real-space basis function.
        !    phi4: index of a real-space basis function.
        ! Returns:
        !    The anti-symmetrized integral <phi1 phi2 || phi3 phi4>.

        use basis
        use system

        real(dp) :: two_e_int
        integer, intent(in) :: phi1, phi2, phi3, phi4

        ! <phi1 phi2 || phi3 phi4>
        two_e_int = 0.0_dp

    end function get_two_e_int_real

end module hubbard_real
