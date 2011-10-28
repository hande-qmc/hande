module symmetry

! Module for symmetry routines common to all symmetries.

! This depends upon system-specific symmetry modules, so take care not to
! introduce circular dependencies by USEing it in the system-specific symmetry
! modules.

implicit none

contains

    pure function symmetry_orb_list(orb_list) result(isym)

        ! In:
        !    orb_list: list of orbitals (e.g. determinant).
        ! Returns:
        !    symmetry index of list (i.e. direct product of the representations
        !    of all the orbitals in the list).

        use momentum_symmetry, only: symmetry_orb_list_k
        use point_group_symmetry, only: symmetry_orb_list_mol
        use system, only: system_type, hub_k, read_in

        integer :: isym
        integer, intent(in) :: orb_list(:)

        select case(system_type)
        case(hub_k)
            isym = symmetry_orb_list_k(orb_list)
        case(read_in)
            isym = symmetry_orb_list_mol(orb_list)
        case default
            ! symmetry not implemented
            isym = 1
        end select

    end function symmetry_orb_list

end module symmetry
