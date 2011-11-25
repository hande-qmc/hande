module hubbard

! Procedures relating to the Hubbard model (real- and momentum-space
! formulations).

use basis

implicit none

contains

    subroutine init_model_basis_fns()

        ! Produce the basis functions for model Hamiltonian systems.
        !
        ! The number of wavevectors is equal to the number of sites in the
        ! crystal cell (ie the number of k-points used to sample the FBZ of the
        ! primitive cell).  From the cell parameters and the "tilt" used (if
        ! any) generate the list of wavevectors and hence the kinetic energy
        ! associated with each basis function (two per wavevector to account for
        ! spin).

        use checking, only: check_allocate, check_deallocate
        use system
        use m_mrgref, only: mrgref
        use errors, only: stop_all
        use parallel, only: parent

        integer :: limits(3,3), nmax(3), kp(3) ! Support a maximum of 3 dimensions.
        integer :: i, j, k, ibasis, ierr
        type(basis_fn), allocatable, target :: tmp_basis_fns(:)
        type(basis_fn), pointer :: basis_fn_p
        integer, allocatable :: basis_fns_ranking(:)

        nbasis = 2*nsites
        nvirt = nbasis - nel

        ! Find basis functions.

        ! We use a minimal basis: the hubbard model consisting of two
        ! spin-orbitals per lattice site.

        ! In the momentum space formulation the basis functions consist of a
        ! set of wavevectors/k-points that lie within the first Brillouin zone.

        ! In the real space formulation the basis functions used are those
        ! residing at the lattice points: we just need to find which lattice
        ! points fall within the crystal cell.

        ! Momentum space:

        ! Fold the crystal cell into the FBZ.
        ! The k-points must be integer multiples of the reciprocal lattice
        ! vectors of the crystal cell (so that the wavefunction is periodic in
        ! the crystal cell) and fall within the first Brillouin zone of the
        ! primitive unit cell (so that a unique set of k-points are chosen).
        ! The volume of the FBZ is inversely proportional to the volume of the
        ! cell, and so the number of sites in the crystal cell is equal to the
        ! number of reciprocal crystal cells in the FBZ of the unit cell and 
        ! hence this gives the required number of wavevectors.

        ! Real space:

        ! The same procedure applies as for the momentum space: we find which
        ! lattice sites lie within the Wigner--Seitz cell.  In fact, due to the
        ! relationship between reciprocal space and real space (and due to how
        ! we store the wavevectors in terms of the reciprocal lattice vectors of
        ! the crystal cell), *exactly* the same approach is needed, so we're
        ! just going to abuse the same code. Shocking, I know.

        ! Maximum limits...
        ! [Does it show that I've been writing a lot of python recently?]
        nmax = 0 ! Set nmax(i) to be 0 for unused higher dimensions.
        limits = 0
        ! forall is a poor substitute for list comprehension. ;-)
        forall (i=1:ndim)
            forall (j=1:ndim, lattice(i,j) /= 0) 
                limits(i,j) = abs(nint(box_length(i)**2/(2*lattice(i,j))))
            end forall
            nmax(i) = maxval(limits(:,i))
        end forall

        allocate(basis_fns(nbasis), stat=ierr)
        call check_allocate('basis_fns',nbasis,ierr)
        allocate(tmp_basis_fns(nbasis/2), stat=ierr)
        call check_allocate('tmp_basis_fns',nbasis/2,ierr)
        allocate(basis_fns_ranking(nbasis/2), stat=ierr)
        call check_allocate('basis_fns_ranking',nbasis/2,ierr)

        ! Find all alpha spin orbitals.
        ibasis = 0
        do k = -nmax(3), nmax(3)
            do j = -nmax(2), nmax(2)
                do i = -nmax(1), nmax(1)
                    ! kp is the Wavevector in terms of the reciprocal lattice vectors of
                    ! the crystal cell.
                    kp = (/ i, j, k /)
                    if (in_FBZ(kp(1:ndim))) then
                        if (ibasis==nbasis) then
                            call stop_all('init_basis_fns','Too many basis functions found.')
                        else
                            ! Have found an allowed wavevector/site.
                            ! Add 2 spin orbitals to the set of the basis functions.
                            ibasis = ibasis + 1
                            call init_basis_fn(tmp_basis_fns(ibasis), l=kp(1:ndim), ms=1)
                        end if
                    end if
                end do
            end do
        end do

        if (ibasis /= nbasis/2) call stop_all('init_basis_fns','Not enough basis functions found.')

        ! Rank by kinetic energy (applies to momentum space formulation only).
        select case(system_type)
        case(hub_k)
            call mrgref(tmp_basis_fns(:)%sp_eigv, basis_fns_ranking)
        case(hub_real)
            forall (i=1:nbasis/2) basis_fns_ranking(i) = i
        end select

        ! Form the list of sorted basis functions with both alpha and beta
        ! spins.
        do i = 1, nbasis/2
            ! Can't set a kpoint equal to another kpoint as then the k pointers
            ! can be assigned whereas we want to *copy* the values.
            basis_fn_p => tmp_basis_fns(basis_fns_ranking(i))
            call init_basis_fn(basis_fns(2*i-1), l=basis_fn_p%l, ms=basis_fn_p%ms)
            call init_basis_fn(basis_fns(2*i), l=basis_fn_p%l, ms=-basis_fn_p%ms)
            deallocate(tmp_basis_fns(basis_fns_ranking(i))%l, stat=ierr)
            call check_deallocate('tmp_basis_fns(basis_fns_ranking(i',ierr)
        end do
        deallocate(tmp_basis_fns, stat=ierr)
        call check_deallocate('tmp_basis_fns',ierr)
        deallocate(basis_fns_ranking, stat=ierr)
        call check_deallocate('basis_fns_ranking',ierr)

        if (parent) then
            call write_basis_fn_header()
            do i = 1, nbasis
                call write_basis_fn(basis_fns(i), ind=i, new_line=.true.)
            end do
            write (6,'()')
        end if

    end subroutine init_model_basis_fns

end module hubbard
