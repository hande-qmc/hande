module hamiltonian

use const

implicit none

contains

    pure function get_hmatel_dets(sys, d1, d2) result(hmatel)

        ! In:
        !    sys: system being studied.
        !    d1, d2: integer labels of two determinants, as stored in the
        !            dets array.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >.

        ! This is just a wrapper around get_hmatel (which is itself a wrapper
        ! around system-specific functions) but is handy for computing matrix
        ! elements (slowly!) when we have the entire Hilbert space of
        ! determinants stored in dets_list.

        use determinant_enumeration, only: dets_list
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: d1, d2

        hmatel = get_hmatel(sys, dets_list(:,d1), dets_list(:,d2))

    end function get_hmatel_dets

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

        use determinants, only: basis_length
        use hamiltonian_chung_landau, only: get_hmatel_chung_landau
        use hamiltonian_heisenberg, only: get_hmatel_heisenberg
        use hamiltonian_hub_k, only: get_hmatel_hub_k
        use hamiltonian_hub_real, only: get_hmatel_hub_real
        use hamiltonian_molecular, only: get_hmatel_mol
        use hamiltonian_ueg, only: get_hmatel_ueg
        use system

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)

        select case(sys%system)
        case(chung_landau)
            hmatel = get_hmatel_chung_landau(f1, f2)
        case(hub_k)
            hmatel = get_hmatel_hub_k(f1, f2)
        case(hub_real)
            hmatel = get_hmatel_hub_real(f1, f2)
        case(heisenberg)
            hmatel = get_hmatel_heisenberg(f1, f2)
        case(read_in)
            hmatel = get_hmatel_mol(f1, f2)
        case (ueg)
            hmatel = get_hmatel_ueg(f1, f2)
        end select

    end function get_hmatel

end module hamiltonian
