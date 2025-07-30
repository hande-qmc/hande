module hamiltonian

use const
use hamiltonian_data

implicit none

contains

    pure function get_hmatel(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system being studied.
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

        use hamiltonian_chung_landau, only: get_hmatel_chung_landau
        use hamiltonian_heisenberg, only: get_hmatel_heisenberg
        use hamiltonian_hub_k, only: get_hmatel_hub_k
        use hamiltonian_hub_real, only: get_hmatel_hub_real
        use hamiltonian_molecular, only: get_hmatel_mol
        use hamiltonian_ueg, only: get_hmatel_ueg
        use hamiltonian_ringium, only: get_hmatel_ringium
        use hamiltonian_trotter, only: get_hmatel_trotter
        use system

        type(hmatel_t) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(:), f2(:)

        select case(sys%system)
        case(chung_landau)
            hmatel = get_hmatel_chung_landau(sys, f1, f2)
        case(hub_k)
            hmatel = get_hmatel_hub_k(sys, f1, f2)
        case(hub_real)
            hmatel = get_hmatel_hub_real(sys, f1, f2)
        case(heisenberg)
            hmatel = get_hmatel_heisenberg(sys, f1, f2)
        case(read_in)
            if (sys%read_in%trotter) then
                hmatel = get_hmatel_trotter(sys, f1, f2)
            else
                hmatel = get_hmatel_mol(sys, f1, f2)
            end if
        case (ueg)
            hmatel = get_hmatel_ueg(sys, f1, f2)
        case(ringium)
            hmatel = get_hmatel_ringium(sys, f1, f2)
        end select

    end function get_hmatel

    pure function get_hmatel_complex(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system being studied.
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

        ! Use separate function for obtaining complex hmatel as:
        !  -Avoids having either dummy argument for imaginary component
        !   to real function (in form of extra result or complex component
        !   of existing). This would break all existing references to function
        !   in existing real codes due to either number of results or complex
        !   vs real types.
        !  -Limits number of bool checks to get hmatel.
        !  -Could make it easier to use pointers later (if complex just point
        !   to this one).
        !  -Alternatively have to define second hmatel factor in all real
        !   routines to allow return of two values from all pure values.
        !  -If use array of allocatable length, redefine existing for minimally
        !   useful changes (haven't got complex behaviour/arithmetic).

        use hamiltonian_periodic_complex, only: get_hmatel_periodic_complex
        use system

        type(hmatel_t) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(:), f2(:)

        select case(sys%system)
        case(read_in)
            hmatel = get_hmatel_periodic_complex(sys, f1, f2)
        end select

    end function get_hmatel_complex

end module hamiltonian
