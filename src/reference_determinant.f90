module reference_determinant

! Utility module for selecting and manipulating reference determinants.

! Please take care with adding module dependencies (directly or indirectly) via
! use statements.  This is intended to be a 'utility' module accessible from
! a variety of calculations (not just QMC) and hence care must be taken to avoid
! circular dependencies.

implicit none

contains

!--- Attempt to find a reasonable reference determinant ---

    subroutine set_reference_det(occ_list, override_input)

        ! Set the list of occupied orbitals in the reference determinant to be
        ! the spin-orbitals with the lowest kinetic energy which satisfy the
        ! spin polarisation.

        ! Note: this is for testing only!  The symmetry input is currently
        ! ignored.

        ! This should be used as a last resort if the user doesn't specify
        ! a reference determinant.  It attempts to generate a sensible (but by
        ! no means guaranteed best) suitable reference determinant.  The
        ! definition of 'best' is system dependent (see code below).

        ! In/Out:
        !   occ_list: allocatable integer array.  Contains the reference
        !   determinant (nel elements) on output.  Unchanged if already
        !   allocated unless override_input is set.
        ! In:
        !   override_input: if true, overwrite occ_list with the best guess of
        !   a reference determinant even if occ_list is allocated on input.

        use checking, only: check_allocate

        use errors, only: stop_all
        use system, only: nalpha, nbeta, nel, system_type, hub_k, hub_real, read_in, ueg, nsites, &
                          heisenberg, J_coupling
        use basis, only: bit_lookup
        use hubbard_real, only: connected_orbs

        integer, intent(inout), allocatable :: occ_list(:)
        logical, intent(in) :: override_input

        integer :: i, j, ierr, spins_set, connections
        integer :: bit_element, bit_pos
        logical :: set

        ! Leave the reference determinant unchanged if it's already been
        ! allocated (and presumably set).

        if (allocated(occ_list)) then
            ! Already set.  Just check that it's not totally insane.
            if (size(occ_list) /= nel) then
                select case(system_type)
                case(heisenberg)
                    call stop_all('set_reference_det', &
                        'Reference determinant supplied does not contain the &
                        &specified number of up spins.')
                case default
                    call stop_all('set_reference_det', &
                        'Reference determinant supplied does not contain the &
                        &specified number of electrons.')
                end select
            end if
            set = .true.
        else
            ! Allocate memory if required.
            allocate(occ_list(nel), stat=ierr)
            call check_allocate('occ_list',nel,ierr)
            set = .false.
        end if

        if (.not.set .or. override_input) then
            ! Attempt to find a dumb 'best guess' for choice of a reference
            ! determinant.
            select case(system_type)
            case(hub_k,read_in,ueg)
                ! Orbitals are ordered by their single-particle eigenvalues.
                ! Occupy the Fermi sphere/HF det.
                forall (i=1:nalpha) occ_list(i) = 2*i-1
                forall (i=1:nbeta) occ_list(i+nalpha) = 2*i
            case(hub_real)
                ! Attempt to keep electrons on different sites where possible.
                ! Sites 1, 3, 5, ... (occupy every other alpha orbital first, ie
                ! place a max of nsites/2 electrons.  (nsites+1)/2 accounts for
                ! the possibility that we have an odd number of sites.)
                forall (i=1:min(nalpha,(nsites+1)/2)) occ_list(i) = 4*i-3
                ! now occupy the alternate alpha orbitals
                forall (i=1:nalpha-min(nalpha,(nsites+1)/2)) &
                    occ_list(i+min(nalpha,(nsites+1)/2)) = 4*i-1
                ! Similarly for beta, but now occupying orbitals sites 2, 4,
                ! ..., preferentially.
                forall (i=1:min(nbeta,nsites/2)) occ_list(i+nalpha) = 4*i
                forall (i=1:nbeta-min(nbeta,nsites/2)) &
                    occ_list(i+nalpha+min(nbeta,nsites/2)) = 4*i-2
            case(heisenberg)
                ! Ferromagnetic case is easy: group identical spins together!
                if (J_coupling >= 0) then
                    forall (i=1:nel) occ_list(i) = i
                ! For the antiferromagnetic case, below. This is messy but should
                ! give a reasonable reference determinant for general cases, even
                ! for bizarre lattices. For bipartite lattices (eg 4x4, 6x6...)
                ! it will give the best possible reference determinant.
                else if (J_coupling < 0) then
                    ! Always set the first spin up
                    occ_list(1) = 1
                    spins_set = 1
                    ! Loop over other sites to find orbitals which are not connected to
                    ! the other sites previously chosen.
                    do i=2,nsites
                        bit_pos = bit_lookup(1,i)
                        bit_element = bit_lookup(2,i)
                        connections = 0
                        ! Loop over all chosen sites to see if they neighbour this site.
                        do j=1,spins_set
                            if (btest(connected_orbs(bit_element, occ_list(j)), bit_pos)) then
                                  connections = connections + 1
                            end if
                        end do
                        ! If this site has no neighbours which have been previously added
                        ! to the reference determinant, then we include it.
                        if (connections == 0) then
                            spins_set = spins_set + 1
                            occ_list(spins_set) = i
                        end if
                    end do
                    ! If, after finding all the sites which are not connected, we still haven't
                    ! chosen enough sites, we accept that we must have some neigbouring sites
                    ! included in the reference determinant and start choosing the remaining sites.
                    if (spins_set /= nel) then
                        ! Loop over all sites looking for extra spins to include in the
                        ! reference detereminant.
                        fill_sites: do i=2,nsites
                            connections = 0
                            ! Check if this site is already included.
                            do j=1,spins_set
                                if (occ_list(j) == i) connections = connections + 1
                            end do
                            ! If connection = 0, this site is not currently included in the
                            ! reference determinant, so add it.
                            if (connections == 0) then
                                spins_set = spins_set + 1
                                occ_list(spins_set) = i
                            end if
                            ! When the correct number of spins have been chosen to be up,
                            ! we are finished.
                            if (spins_set == nel) exit fill_sites
                        end do fill_sites
                    end if
                end if
            end select
        end if

    end subroutine set_reference_det


end module reference_determinant
