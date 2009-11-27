module hubbard

use system
use kpoints

implicit none

type(kpoint), allocatable :: basis_fns(:)

contains

    subroutine init_basis_fns()

        ! Produce the basis functions.  The number of basis functions is
        ! equal to the number of sites in the crystal cell (ie the number
        ! of k-points used to sample the FBZ of the primitive cell).
        ! From the cell parameters and the "tilt" used (if any) generate
        ! the list of wavevectors and hence the kinetic energy associated
        ! with each basis function (one per wavevector).

        use errors, only: stop_all
        use parallel, only: iproc, parent

        integer :: nmax(3), kp(3) ! Support a maximum of 3 dimensions.
        integer :: i, j, k, nbasis, ierr
        type(kpoint) :: test_basis_fn
        logical :: t

        ! [Does it show that I've been writing a lot of python recently?]
        nmax = 0
        forall (i=1:ndim) nmax(i) = maxval(abs(nint(box_length(i)**2/(2*lattice(:,i)))))

        allocate(basis_fns(nsites), stat=ierr)

        nbasis = 0
        ! The following can't be done as a forall due to compilers being crap...
        do k = -nmax(3), nmax(3)
            do j = -nmax(2), nmax(2)
                do i = -nmax(1), nmax(1)
                    kp = (/ i, j, k /)
                    if (in_FBZ(kp(1:ndim))) then
                        nbasis = nbasis + 1
                        if (nbasis>nsites) then
                            call stop_all('init_basis_fns','Too many basis functions found.')
                        else
                            call init_kpoint(basis_fns(nbasis), kp(1:ndim))
                        end if
                    end if
                end do
            end do
        end do

        if (nbasis/=nsites) call stop_all('init_basis_fns','Not enough basis functions found.')

        if (iproc==parent) then
            write (6,'(/,1X,a15,/,1X,15("-"),/)') 'Basis functions'
            write (6,'(1X,a7)', advance='no') 'k-point'
            do i = 1,ndim
                write (6,'(3X)', advance='no')
            end do
            write (6,'(a7)') 'kinetic'
            do i = 1,nbasis
                call write_kpoint(basis_fns(i))
            end do
        end if

    end subroutine init_basis_fns

end module hubbard
