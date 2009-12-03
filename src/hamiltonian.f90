module hamiltonian

use const
use basis
use determinants
use hubbard
use parallel

implicit none

! Hamiltonian matrix.  Clearly the scaling of the memory demands with system
! size is horrendous.
real(dp), allocatable :: hamil(:,:) ! (ndets, ndets)

! If true, then the eigenvectors are found during exact diagonalisation as well
! as the eigenvalues.  Doing this is substantially more expensive.
logical :: find_eigenvectors = .false.

! If true then the non-zero elements of the Hamiltonian matrix are written to hamiltonian_file.
logical :: write_hamiltonian = .false.
character(255) :: hamiltonian_file = 'HAMIL'

contains

    subroutine generate_hamil()

        ! Generate the Hamiltonian matrix.
        ! Only generate the upper diagonal for use with Lapack routines.

        use utils, only: get_free_unit

        integer :: ierr, i, j, iunit
        type(excit) :: excitation

        allocate(hamil(ndets,ndets), stat=ierr)

        do i=1, ndets
            do j=i, ndets
                if (dets(i)%Ms == dets(j) % Ms) then
                    excitation = get_excitation(dets(i)%f, dets(j)%f)
                    if (excitation%nexcit <= 2) then
                        hamil(i,j) = get_hmatel(dets(i)%f, excitation)
                    else
                        hamil(i,j) = 0.0_dp
                    end if
                else
                    hamil(i,j) = 0.0_dp
                end if
            end do
        end do

        ! Fill in rest of Hamiltonian matrix.  In this case the Hamiltonian is
        ! symmetric (rather than "just" Hermitian).
!        do i=1,ndets
!            do j=1,i-1
!                hamil(i,j) = hamil(j,i)
!            end do
!        end do

        if (write_hamiltonian) then
            iunit = get_free_unit()
            open(iunit, file=hamiltonian_file, status='unknown')
            do i=1,ndets
                write (iunit,*) i,i,hamil(i,i)
                do j=i+1, ndets
                    if (abs(hamil(i,j)) > depsilon) write (iunit,*) i,j,hamil(i,j)
                end do
            end do
            close(iunit, status='keep')
        end if

    end subroutine generate_hamil

    subroutine exact_diagonalisation()
    
        ! Perform an exact diagonalisation of the Hamiltonian matrix.
        ! Note that this destroys the Hamiltonian matrix stored in hamil.

        real(dp), allocatable :: eigv(:), work(:)
        integer :: info, ierr, lwork
        integer :: i

        if (iproc == parent) then
            write (6,'(1X,a21,/,1X,21("-"))') 'Exact diagonalisation'
            write (6,'(/,1X,a35,/)') 'Performing exact diagonalisation...'
        end if

        lwork = 3*ndets - 1
        allocate(work(lwork), stat=ierr)
        allocate(eigv(ndets), stat=ierr)

        if (find_eigenvectors) then
            call dsyev('V', 'U', ndets, hamil, ndets, eigv, work, lwork, info)
        else
            call dsyev('N', 'U', ndets, hamil, ndets, eigv, work, lwork, info)
        end if

        deallocate(work)

        if (iproc == parent) then
            write (6,'(1X,a8,3X,a12)') 'State','Total energy'
            do i = 1, ndets
                write (6,'(1X,i8,f18.12)') i, eigv(i)
            end do
            write (6,'(/,1X,a13,f18.12,/)') 'Ground state:', eigv(1)
        end if

    end subroutine exact_diagonalisation

    subroutine lanczos_diagonalisation()

       use trl_info
       use trl_interface
       
       integer, parameter :: lohi = -1, ned = 5, maxlan = 40, mev = 10
       real(dp) :: eval(mev)
       real(dp), allocatable :: evec(:,:) ! (ndets, mev)
       type(trl_info_t) :: info
       integer :: i

       call trl_init_info(info, ndets, maxlan, lohi, ned)
       call trlan(hamil_vector, info, ndets, mev, eval, evec, ndets)
       call trl_print_info(info, ndets*2)

       contains

           subroutine hamil_vector(nrow, ncol, xin, ldx, yout, ldy)

               implicit None
               integer, intent(in) :: nrow, ncol, ldx, ldy
               real(dp), intent(in) :: xin(ldx,ncol)
               real(dp), intent(out) :: yout(ldy,ncol)
               ! local variables
               integer :: i, j, ioff, joff

               do j = 1, ncol
                   ioff = (j-1)*ldx
                   joff = (j-1)*ldy
                   do i = 1, nrow
                       yout(joff,i) = (i*i)*xin(ioff,i)
                   end do
               end do

           end subroutine hamil_vector

    end subroutine lanczos_diagonalisation
    
    pure function get_hmatel(root_det_f, excitation) result(hmatel)

        ! In:
        !    root_det_f(basis_length): bit string representation of the "root" Slater
        !        determinant (i.e. the determinant the excitation is from).
        !    excitation: excitation information linking another determinant to
        !        the root determinant.
        ! Returns:
        !    Hamiltonian matrix element between the root determinant and the
        !        excited determinant.

        real(dp) :: hmatel
        type(excit), intent(in) :: excitation
        integer(i0), intent(in) :: root_det_f(basis_length)
        integer :: root_det(nel)
        integer :: i, j

        hmatel = 0.0_dp

        select case(excitation%nexcit)
        case(0)

            root_det = decode_det(root_det_f)

            ! One electron operator
            do i = 1, nel
                hmatel = hmatel + get_one_e_int(root_det(i), root_det(i))
            end do

            ! Two electron operator
            do i = 1, nel
                do j = i+1, nel
                    hmatel = hmatel + get_two_e_int(root_det(i), root_det(j), root_det(i), root_det(j))
                end do
            end do

        case(1)

            root_det = decode_det(root_det_f)

            ! One electron operator
            hmatel = hmatel + get_one_e_int(excitation%from_orb(1), excitation%to_orb(1)) 

            ! Two electron operator
            do i = 1, nel
                hmatel = hmatel + get_two_e_int(root_det(i), excitation%from_orb(1), root_det(i), excitation%to_orb(1))
            end do

            if (excitation%perm) hmatel = -hmatel

        case(2)

            ! Two electron operator
            hmatel = get_two_e_int(excitation%from_orb(1), excitation%from_orb(2), excitation%to_orb(1), excitation%to_orb(2))

            if (excitation%perm) hmatel = -hmatel

        end select

    end function get_hmatel

end module hamiltonian
