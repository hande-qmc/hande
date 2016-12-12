module calc_system_init

! Functions required for final initialisation of system at start of calculations.

use system

implicit none

contains

    subroutine set_spin_polarisation(nbasis, sys)

        ! Set the spin polarisation information stored in components of sys.
        !    nalpha, nbeta: number of alpha, beta electrons.
        !    nvirt_alpha, nvirt_beta: number of alpha, beta virtual spin-orbitals.

        ! In:
        !    nbasis: number of single-particle basis functions.
        ! In/Out:
        !    sys: initialised system object describing system. On output
        !       components related to spin-polarisation are set.

        use errors, only: stop_all

        integer, intent(in) :: nbasis
        type(sys_t), intent(inout) :: sys

        character(len=*), parameter :: proc_name = 'set_spin_polarisation'
        character(len=*), parameter :: err_fmt = '("Required Ms not possible: ",i11, ".")'
        character(len=40) :: err

        select case(sys%system)

        case(heisenberg)

            ! Spin polarization is different (see comments in system) as the
            ! Heisenberg model is a collection of spins rather than electrons.
            ! See comments in init_system and at module-level.
            if (abs(sys%Ms) > sys%lattice%nsites) then
                write (err, err_fmt) sys%Ms
                call stop_all(proc_name, err)
            end if
            sys%nel = (sys%lattice%nsites + sys%Ms)/2
            sys%nvirt = (sys%lattice%nsites - sys%Ms)/2
            ! The Heisenberg model doesn't use these values, but they need to be
            ! initialized to something sensible as we allocate memory using them in
            ! alloc_det_info_t.
            sys%nalpha = 0
            sys%nbeta = 0
            sys%nvirt_alpha = 0
            sys%nvirt_beta = 0
            sys%max_number_excitations = min(sys%nel, (sys%lattice%nsites-sys%nel))

        case(chung_landau)

            ! Spinless system.  Similarly for the Heisenberg model but treat all fermions as alpha electrons.
            sys%nalpha = sys%nel
            sys%nvirt_alpha = sys%lattice%nsites - sys%nalpha
            sys%nbeta = 0
            sys%nvirt_beta = 0
            sys%max_number_excitations = min(sys%nel, (sys%lattice%nsites-sys%nel))

        case (ringium)

            if (sys%ms /= sys%nel) call stop_all(proc_name, "Ringium must be completely spin polarised.")

            sys%nalpha = sys%nel
            sys%nbeta = 0

            sys%nvirt_alpha = nbasis/2 - sys%nel
            sys%nvirt_beta = 0
            sys%max_number_excitations = min(sys%nalpha, sys%nvirt_alpha)

        case default

            ! Test whether the required spin is possible given the numbe of electrons.
            if (abs(mod(sys%Ms,2)) /= mod(sys%nel,2) .or. abs(sys%Ms) > sys%nel) then
                write (err, err_fmt) sys%Ms
                call stop_all(proc_name, err)
            end if

            sys%nbeta = (sys%nel - sys%Ms)/2
            sys%nalpha = (sys%nel + sys%Ms)/2

            sys%nvirt_alpha = nbasis/2 - sys%nalpha
            sys%nvirt_beta = nbasis/2 - sys%nbeta
            sys%max_number_excitations = min(sys%nel, sys%nvirt)

        end select

        call set_fermi_energy(sys)

    end subroutine set_spin_polarisation

    subroutine set_symmetry_aufbau(sys, io_unit)

        ! If initialising symmetry using Aufbau principle-chosen determinant
        ! find this determinant and determine its symmetry.

        ! In/Out:
        !   sys: object containing information about system under consideration.
        !       Must have initialised basis and symmetry.
        ! In:
        !   io_unit: io unit to write any additional info to.

        use reference_determinant, only: set_reference_det
        use symmetry, only: symmetry_orb_list
        use system, only: sys_t

        type(sys_t), intent(inout) :: sys
        integer, intent(in) :: io_unit
        integer, allocatable :: occ_list(:)

        ! Find the approximate lowest energy determinant.
        call set_reference_det(sys, occ_list, .false., sys%symmetry, io_unit)

        ! Set symmetry sector equal to that of lowest energy determinant.
        sys%symmetry= symmetry_orb_list(sys, occ_list)

    end subroutine set_symmetry_aufbau

end module calc_system_init
