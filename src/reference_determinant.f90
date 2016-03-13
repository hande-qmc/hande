module reference_determinant

! Utility module for selecting and manipulating reference determinants.

! Please take care with adding module dependencies (directly or indirectly) via
! use statements.  This is intended to be a 'utility' module accessible from
! a variety of calculations (not just QMC) and hence care must be taken to avoid
! circular dependencies.

use const, only: i0, p

implicit none

type reference_t
    ! Bit string of reference determinant.
    integer(i0), allocatable :: f0(:)
    ! List of occupied orbitals in reference determinant.
    integer, allocatable :: occ_list0(:)
    ! Bit string of reference determinant used to generate the Hilbert space.
    ! This is usually identical to f0, but not necessarily (e.g. if we're doing
    ! a spin-flip calculation).
    integer(i0), allocatable :: hs_f0(:)
    ! hs_f0:hs_occ_list0 as f0:occ_list0.
    integer, allocatable :: hs_occ_list0(:)
    ! CCMC/CIQMC: max number of excitations from the reference to include in
    ! the calculation.
    ! DMQMC: permit density matrix elements to be non-zero only if the two
    ! determinants differ by at most ex_level excitations.
    ! Set to the number of electrons in the system to use the full space.
    integer :: ex_level = -1
    ! Energy of reference determinant.
    real(p) :: H00
    ! Value of <D0|O|D0>, where O is the operator we are sampling.
    ! (Applicable/set only if Hellmann--Feynman sampling is in operation.)
    real(p) :: O00
    ! Energy shift of the reference determinant i.e.
    ! <D0|H_new|D0>-<D0|H_old|D0>, where H_old and H_new are two different
    ! Hamiltonians. Used in ip-dmqmc when reweighing the initial density matrix.
    real(p) :: energy_shift = 0.0_p
end type reference_t

contains

!--- Attempt to find a reasonable reference determinant ---

    subroutine set_reference_det(sys, occ_list, override_input, ref_sym)

        ! Set the list of occupied orbitals in the reference determinant to be
        ! the spin-orbitals with the lowest kinetic energy which satisfy the
        ! spin polarisation.

        ! Note: this is for testing only!
        ! This should be used as a last resort if the user doesn't specify
        ! a reference determinant.  It attempts to generate a sensible (but by
        ! no means guaranteed best) suitable reference determinant.  The
        ! definition of 'best' is system dependent (see code below).

        ! In/Out:
        !   occ_list: allocatable integer array.  Contains a 'best guess' at
        !       a suitable reference determinant (nel elements; list of occupied
        !       orbitals in determinant) on output.  Note that this is not ordered.
        !       Unchanged if already allocated unless override_input is set.
        ! In:
        !   sys: system being studied.
        !   override_input: if true, overwrite occ_list with the best guess of
        !       a reference determinant even if occ_list is allocated on input.
        !   ref_sym (optional): if supplied, attempt to find the reference
        !       determinant with the lowest sum of single-particle energies with
        !       this symmetry index.  Ignored if less than sym0 or greater than
        !       sym_max.

        use const, only: i0, p, depsilon
        use checking, only: check_allocate
        use errors, only: stop_all

        use determinants, only: encode_det
        use symmetry, only: symmetry_orb_list
        use system

        type(sys_t), intent(in) :: sys
        integer, intent(inout), allocatable :: occ_list(:)
        logical, intent(in) :: override_input
        integer, intent(in), optional :: ref_sym

        integer :: i, j, ierr, spins_set, connections, iel, icore, jcore, ivirt, jvirt
        integer :: bit_element, bit_pos, tmp_occ_list(sys%nel), curr_occ_list(sys%nel), sym
        integer(i0) :: f(sys%basis%string_len)
        real(p) :: eigv_sum, sp_eigv_sum
        logical :: set

        ! Leave the reference determinant unchanged if it's already been
        ! allocated (and presumably set).

        if (allocated(occ_list)) then
            ! Already set.  Just check that it's not totally insane.
            if (size(occ_list) /= sys%nel) then
                select case(sys%system)
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
            allocate(occ_list(sys%nel), stat=ierr)
            call check_allocate('occ_list',sys%nel,ierr)
            set = .false.
        end if

        if (.not.set .or. override_input) then
            ! Attempt to find a dumb 'best guess' for choice of a reference
            ! determinant.
            select case(sys%system)
            case(hub_k,read_in,ueg,ringium)
                ! Orbitals are ordered by their single-particle eigenvalues.
                ! Occupy the Fermi sphere/HF det.
                forall (i=1:sys%nalpha) occ_list(i) = 2*i-1
                forall (i=1:sys%nbeta) occ_list(i+sys%nalpha) = 2*i

                ! Symmetry only implemented for these systems.  Do we need
                ! to find a determinant of different symmetry?
                ! This needs only be called for initialisation, so don't attempt
                ! to be clever and efficient...
                if (present(ref_sym)) then
                    if (ref_sym >= sys%sym0 .and. ref_sym <= sys%sym_max) then
                        call encode_det(sys%basis, occ_list, f)
                        ! If occ_list is already of the correct symmetry, then
                        ! have nothing to do.
                        sym = symmetry_orb_list(sys, occ_list)
                        if (sym /= ref_sym) then

                            ! Consider single excitations of our current
                            ! reference determinant, conserving only spin.
                            eigv_sum = huge(0.0_p)
                            do icore = 1, sys%nel
                                i = occ_list(icore)
                                do ivirt = 1, sys%basis%nbasis
                                    ! Ensure ivirt is not already in the
                                    ! determinant.
                                    if (.not.btest(f(sys%basis%bit_lookup(2,ivirt)), sys%basis%bit_lookup(1,ivirt)) .and. &
                                            sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(ivirt)%Ms) then
                                        tmp_occ_list = occ_list
                                        tmp_occ_list(icore) = ivirt
                                        if (symmetry_orb_list(sys, tmp_occ_list) == ref_sym) then
                                            sp_eigv_sum = 0.0_p
                                            do iel = 1, sys%nel
                                                sp_eigv_sum = sp_eigv_sum + sys%basis%basis_fns(tmp_occ_list(iel))%sp_eigv
                                            end do
                                            if (sp_eigv_sum+depsilon < eigv_sum) then
                                                curr_occ_list = tmp_occ_list
                                                eigv_sum = sp_eigv_sum
                                            end if
                                        end if
                                    end if
                                end do
                            end do

                            ! Consider double excitations of our current
                            ! reference determinant, conserving only spin.
                            associate(bit_lookup=>sys%basis%bit_lookup, nbasis=>sys%basis%nbasis, &
                                    basis_fns=>sys%basis%basis_fns)
                                do icore = 1, sys%nel
                                    i = occ_list(icore)
                                    do jcore = icore+1, sys%nel
                                        j = occ_list(jcore)
                                        do ivirt = 1, nbasis
                                            if (.not.btest(f(bit_lookup(2,ivirt)), bit_lookup(1,ivirt))) then
                                                do jvirt = ivirt+1, nbasis
                                                    if (.not.btest(f(bit_lookup(2,jvirt)), bit_lookup(1,jvirt)) .and. &
                                                            (basis_fns(i)%Ms + basis_fns(j)%Ms) == &
                                                            (basis_fns(ivirt)%Ms + basis_fns(jvirt)%Ms) ) then
                                                        tmp_occ_list = occ_list
                                                        tmp_occ_list(icore) = ivirt
                                                        tmp_occ_list(jcore) = jvirt
                                                        if (symmetry_orb_list(sys, tmp_occ_list) == ref_sym) then
                                                            sp_eigv_sum = 0.0_p
                                                            do iel = 1, sys%nel
                                                                sp_eigv_sum = sp_eigv_sum + &
                                                                    basis_fns(tmp_occ_list(iel))%sp_eigv
                                                            end do
                                                            if (sp_eigv_sum+depsilon < eigv_sum) then
                                                                curr_occ_list = tmp_occ_list
                                                                eigv_sum = sp_eigv_sum
                                                            end if
                                                        end if
                                                    end if
                                                end do
                                            end if
                                        end do
                                    end do
                                end do
                            end associate

                            occ_list = curr_occ_list
                            if (.not. eigv_sum < huge(0.0_p)) then
                                call stop_all('set_reference_det', &
                                    'Could not find determinant of required symmetry.')
                            end if

                        end if
                    end if
                end if
            case(hub_real)
                ! Attempt to keep electrons on different sites where possible.
                ! Sites 1, 3, 5, ... (occupy every other alpha orbital first, ie
                ! place a max of nsites/2 electrons.  (nsites+1)/2 accounts for
                ! the possibility that we have an odd number of sites.)
                forall (i=1:min(sys%nalpha,(sys%lattice%nsites+1)/2)) occ_list(i) = 4*i-3
                ! now occupy the alternate alpha orbitals
                forall (i=1:sys%nalpha-min(sys%nalpha,(sys%lattice%nsites+1)/2)) &
                    occ_list(i+min(sys%nalpha,(sys%lattice%nsites+1)/2)) = 4*i-1
                ! Similarly for beta, but now occupying orbitals sites 2, 4,
                ! ..., preferentially.
                forall (i=1:min(sys%nbeta,sys%lattice%nsites/2)) occ_list(i+sys%nalpha) = 4*i
                forall (i=1:sys%nbeta-min(sys%nbeta,sys%lattice%nsites/2)) &
                    occ_list(i+sys%nalpha+min(sys%nbeta,sys%lattice%nsites/2)) = 4*i-2
            case(chung_landau)
                ! As with the hub_real, attempt to keep fermions not on
                ! neighbouring sites.
                forall (i=1:sys%nel) occ_list(i) = 2*i-1
            case(heisenberg)
                ! Ferromagnetic case is easy: group identical spins together!
                if (sys%heisenberg%J >= 0) then
                    forall (i=1:sys%nel) occ_list(i) = i
                ! For the antiferromagnetic case, below. This is messy but should
                ! give a reasonable reference determinant for general cases, even
                ! for bizarre lattices. For bipartite lattices (eg 4x4, 6x6...)
                ! it will give the best possible reference determinant.
                else if (sys%heisenberg%J < 0) then
                    ! Always set the first spin up
                    occ_list(1) = 1
                    spins_set = 1
                    ! Loop over other sites to find orbitals which are not connected to
                    ! the other sites previously chosen.
                    do i=2,sys%lattice%nsites
                        bit_pos = sys%basis%bit_lookup(1,i)
                        bit_element = sys%basis%bit_lookup(2,i)
                        connections = 0
                        ! Loop over all chosen sites to see if they neighbour this site.
                        do j=1,spins_set
                            if (btest(sys%real_lattice%connected_orbs(bit_element, occ_list(j)), bit_pos)) then
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
                    if (spins_set /= sys%nel) then
                        ! Loop over all sites looking for extra spins to include in the
                        ! reference detereminant.
                        fill_sites: do i=2,sys%lattice%nsites
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
                            if (spins_set == sys%nel) exit fill_sites
                        end do fill_sites
                    end if
                end if
            end select
        end if

    end subroutine set_reference_det

!--- reference_t utils ---

    subroutine reference_t_json(js, ref, sys, dmqmc_energy_shift, terminal)

        ! Serialise a reference_t object in JSON format.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   reference: reference_t object containing the information about reference state (including any defaults set).
        !   sys (optional): system to which reference belongs.  If present, output the spin and symmetry information of the reference.
        !   dmqmc_energy_shift (optional): if true, print out the 'energy shift' for IP-DMQMC.  Default: false.
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out

        use determinants, only: spin_orb_list
        use system, only: sys_t, heisenberg, chung_landau
        use symmetry, only: symmetry_orb_list

        type(json_out_t), intent(inout) :: js
        type(reference_t), intent(in) :: ref
        type(sys_t), intent(in), optional :: sys
        logical, intent(in), optional :: dmqmc_energy_shift, terminal

        logical :: print_energy_shift

        print_energy_shift = .false.
        if (present(dmqmc_energy_shift)) print_energy_shift = dmqmc_energy_shift

        call json_object_init(js, 'reference')
        if (allocated(ref%occ_list0)) then
            call json_write_key(js, 'det', ref%occ_list0)
            if (present(sys)) then
                if (sys%system /= heisenberg .and. sys%system /= chung_landau) &
                    call json_write_key(js, 'det_ms', spin_orb_list(sys%basis%basis_fns, ref%occ_list0))
                call json_write_key(js, 'det_symmetry', symmetry_orb_list(sys, ref%occ_list0))
            end if
            call json_write_key(js, 'H00', ref%H00)
        end if
        if (allocated(ref%hs_occ_list0)) then
            call json_write_key(js, 'hilbert_space_det', ref%hs_occ_list0)
            if (present(sys)) then
                if (sys%system /= heisenberg .and. sys%system /= chung_landau) &
                    call json_write_key(js, 'hilbert_space_det_ms', spin_orb_list(sys%basis%basis_fns, ref%hs_occ_list0))
                call json_write_key(js, 'hilbert_space_det_symmetry', symmetry_orb_list(sys, ref%hs_occ_list0))
            end if
        end if
        call json_write_key(js, 'ex_level', ref%ex_level, .not.print_energy_shift)
        if (print_energy_shift) call json_write_key(js, 'shift', ref%energy_shift, .true.)
        call json_object_end(js, terminal)

    end subroutine reference_t_json

    subroutine dealloc_reference_t(ref)

        ! In/Out:
        !    ref: reference_t object.  Deallocated on output.

        use checking, only: check_deallocate

        type(reference_t), intent(inout) :: ref

        integer :: ierr

        if (allocated(ref%f0)) then
            deallocate(ref%f0, stat=ierr)
            call check_deallocate('ref%f0', ierr)
        end if
        if (allocated(ref%occ_list0)) then
            deallocate(ref%occ_list0, stat=ierr)
            call check_deallocate('ref%occ_list0', ierr)
        end if
        if (allocated(ref%hs_f0)) then
            deallocate(ref%hs_f0, stat=ierr)
            call check_deallocate('ref%hs_f0', ierr)
        end if
        if (allocated(ref%hs_occ_list0)) then
            deallocate(ref%hs_occ_list0, stat=ierr)
            call check_deallocate('ref%hs_occ_list0', ierr)
        end if

    end subroutine dealloc_reference_t

    subroutine copy_reference_t(ref1, ref2)

        ! Copy information about a reference determinant from one reference_t
        ! object to another.

        ! In:
        !    ref1: reference_t object to be copied.
        ! In/Out:
        !    ref2: destination reference_t object.  Reallocated to be the same
        !       size as ref1 before being set equal.

        use checking, only: check_allocate

        type(reference_t), intent(in) :: ref1
        type(reference_t), intent(inout) :: ref2

        integer :: ierr

        call dealloc_reference_t(ref2)

        if (allocated(ref1%f0)) then
            allocate(ref2%f0(size(ref1%f0)), stat=ierr)
            call check_allocate('ref2%f0', size(ref1%f0), ierr)
        end if
        if (allocated(ref1%occ_list0)) then
            allocate(ref2%occ_list0(size(ref1%occ_list0)), stat=ierr)
            call check_allocate('ref2%occ_list0', size(ref1%occ_list0), ierr)
        end if
        if (allocated(ref1%hs_f0)) then
            allocate(ref2%hs_f0(size(ref1%hs_f0)), stat=ierr)
            call check_allocate('ref2%hs_f0', size(ref1%hs_f0), ierr)
        end if
        if (allocated(ref1%hs_occ_list0)) then
            allocate(ref2%hs_occ_list0(size(ref1%hs_occ_list0)), stat=ierr)
            call check_allocate('ref2%hs_occ_list0', size(ref1%hs_occ_list0), ierr)
        end if

        ref2 = ref1

    end subroutine copy_reference_t

end module reference_determinant
