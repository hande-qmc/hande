module hubbard

use system
use kpoints

implicit none

type(kpoint), allocatable :: basis_fns(:)
integer :: nbasis

contains

    subroutine init_basis_fns()

        ! Produce the basis functions.  The number of wavevectors is
        ! equal to the number of sites in the crystal cell (ie the number
        ! of k-points used to sample the FBZ of the primitive cell).
        ! From the cell parameters and the "tilt" used (if any) generate
        ! the list of wavevectors and hence the kinetic energy associated
        ! with each basis function (two per wavevector to account for spin).

        use errors, only: stop_all
        use parallel, only: iproc, parent

        integer :: nmax(3), kp(3) ! Support a maximum of 3 dimensions.
        integer :: i, j, k, ibasis, ierr
        type(kpoint) :: test_basis_fn
        logical :: t

        nbasis = 2*nsites

        ! Fold the crystal cell into the FBZ.
        ! The k-points must be integer multiples of the reciprocal lattice
        ! vectors of the crystal cell (so that the wavefunction is periodic in
        ! the crystal cell) and fall within the first Brillouin zone of the
        ! primitive unit cell (so that a unique set of k-points are chosen).
        ! The volume of the FBZ is inversely proportional to the volume of the
        ! cell, and so the number of sites in the crystal cell is equal to the
        ! number of reciprocal crystal cells in the FBZ of the unit cell and 
        ! hence this gives the required number of wavevectors.


        ! Maximum limits...
        ! [Does it show that I've been writing a lot of python recently?]
        nmax = 0 ! Set nmax(i) to be 0 for unused higher dimensions.
        forall (i=1:ndim) nmax(i) = maxval(abs(nint(box_length(i)**2/(2*lattice(:,i)))))

        allocate(basis_fns(2*nsites), stat=ierr)

        ibasis = 0
        ! The following can't be done as a forall due to compilers being crap...
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
                            ! Have found an allowed wavevector.
                            ! Add 2 spin orbitals to the set of the basis functions.
                            ibasis = ibasis + 1
                            call init_kpoint(basis_fns(ibasis), kp(1:ndim), 1)
                            ibasis = ibasis + 1
                            call init_kpoint(basis_fns(ibasis), kp(1:ndim), -1)
                        end if
                    end if
                end do
            end do
        end do

        if (ibasis /= nbasis) call stop_all('init_basis_fns','Not enough basis functions found.')

        if (iproc == parent) then
            write (6,'(/,1X,a15,/,1X,15("-"),/)') 'Basis functions'
            write (6,'(1X,a7)', advance='no') 'k-point'
            do i = 1,ndim
                write (6,'(3X)', advance='no')
            end do
            write (6,'(a,4X,a7)') 'ms','kinetic'
            do i = 1,nbasis
                call write_kpoint(basis_fns(i))
            end do
        end if

    end subroutine init_basis_fns

end module hubbard