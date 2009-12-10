module hamiltonian

use const
use basis
use determinants
use hubbard
use parallel

implicit none

! Flags for doing exact and/or Lanczos diagonalisation.
logical :: t_exact = .false., t_lanczos = .false.

! Hamiltonian matrix.  Clearly the scaling of the memory demands with system
! size is horrendous.
real(dp), allocatable :: hamil(:,:) ! (ndets, ndets)

! If true, then the eigenvectors are found during exact diagonalisation as well
! as the eigenvalues.  Doing this is substantially more expensive.
logical :: find_eigenvectors = .false.

! If true then the non-zero elements of the Hamiltonian matrix are written to hamiltonian_file.
logical :: write_hamiltonian = .false.
character(255) :: hamiltonian_file = 'HAMIL'

! Number of Lanczos eigenpairs to find.
integer :: nlanczos_eigv = 5
! Size of Lanczos basis.
integer :: lanczos_basis_length = 40

private :: hamil_vector

contains

    subroutine generate_hamil()

        ! Generate the Hamiltonian matrix.
        ! Only generate the upper diagonal for use with Lapack routines.

        use utils, only: get_free_unit

        integer :: ierr, i, j, iunit
        type(excit) :: excitation
        logical :: hamil_element_set

        allocate(hamil(ndets,ndets), stat=ierr)

        do i=1, ndets
            do j=i, ndets
                hamil_element_set = .false.
                if (dets(i)%Ms == dets(j) % Ms) then
                    excitation = get_excitation(dets(i)%f, dets(j)%f)
                    if (excitation%nexcit <= 2) then
                        if (is_reciprocal_lattice_vector(dets(i)%k-dets(j)%k)) then
                            hamil(i,j) = get_hmatel(dets(i)%f, excitation)
                            hamil_element_set = .true.
                        end if
                    end if
                end if
                if (.not.hamil_element_set) hamil(i,j) = 0.0_dp
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

    subroutine end_hamil()

        ! Clean up hamiltonian module.

        integer :: ierr

        deallocate(hamil, stat=ierr)

    end subroutine end_hamil

    subroutine exact_diagonalisation()
    
        ! Perform an exact diagonalisation of the Hamiltonian matrix.
        ! Note that this destroys the Hamiltonian matrix stored in hamil.

        real(dp), allocatable :: eigv(:), work(:)
        integer :: info, ierr, lwork
        integer :: i

        if (parent) then
            write (6,'(1X,a21,/,1X,21("-"))') 'Exact diagonalisation'
            write (6,'(/,1X,a35,/)') 'Performing exact diagonalisation...'
        end if

        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        call dsyev('N', 'U', ndets, hamil, ndets, eigv, work, -1, info)
        lwork = work(1)
        deallocate(work)

        ! Now perform the diagonalisation.
        allocate(work(lwork), stat=ierr)
        allocate(eigv(ndets), stat=ierr)

        if (find_eigenvectors) then
            call dsyev('V', 'U', ndets, hamil, ndets, eigv, work, lwork, info)
        else
            call dsyev('N', 'U', ndets, hamil, ndets, eigv, work, lwork, info)
        end if

        deallocate(work)

        if (parent) then
            write (6,'(1X,a8,3X,a12)') 'State','Total energy'
            do i = 1, ndets
                write (6,'(1X,i8,f18.12)') i, eigv(i)
            end do
            write (6,'(/,1X,a19,f18.12,/)') 'Exact ground state:', eigv(1)
        end if

    end subroutine exact_diagonalisation

    subroutine lanczos_diagonalisation()

        ! Perform a Lanczos diagonalisation of the Hamiltonian matrix.

        use trl_info
        use trl_interface
        
        integer, parameter :: lohi = -1
        integer :: mev
        real(dp), allocatable :: eval(:) ! (mev)
        real(dp), allocatable :: evec(:,:) ! (ndets, mev)
        type(trl_info_t) :: info
        integer :: i, ierr

        ! mev: number of eigenpairs that can be stored in eval and evec.
        ! twice the number of eigenvalues to be found is a reasonable default.
        mev = max(2*nlanczos_eigv, ndets)
       
        if (parent) then
            write (6,'(1X,a23,/,1X,23("-"))') 'Lanczos diagonalisation'
            write (6,'(/,1X,a37,/)') 'Performing lanczos diagonalisation...'
        end if
       
        ! Initialise trlan.
        ! info: type(trl_info_t).  Used by trl to store calculation info.
        ! ndets: number of rows of matrix on processor.
        ! lanczos_basis_length: maximum Lanczos basis size.
        ! lohi: -1 means calculate the smallest eigenvalues first (1 to calculate
        !       the largest).
        ! nlanczos_eigv: number of eigenvalues to compute.
        call trl_init_info(info, ndets, lanczos_basis_length, lohi, nlanczos_eigv)
       
        allocate(eval(mev), stat=ierr)
        allocate(evec(ndets,mev), stat=ierr)
       
        ! Call Lanczos diagonalizer.
        ! hamil_vector: matrix-vector multiplication routine.
        ! info: created in trl_init_info.
        ! ndets: number of rows of matrix on processor.
        ! mev: number of eigenpairs that can be stored in eval and evec.
        ! eval: array to store eigenvalue
        ! evec: array to store the eigenvectors
        ! lde: the leading dimension of evec (in serial case: ndets).
        call trlan(hamil_vector, info, ndets, mev, eval, evec, ndets)
       
        ! Get info...
        if (parent) then
            write (6,'(1X,a8,3X,a12)') 'State','Total energy'
            do i = 1, nlanczos_eigv
                write (6,'(1X,i8,f18.12)') i, eval(i)
            end do
            write (6,'(/,1X,a21,f18.12,/)') 'Lanczos ground state:', eval(1)
       
            write (6,'(1X,a27,/,1X,27("-"),/)') 'TRLan (Lanczos) information'
            call trl_print_info(info, ndets*2)
            write (6,'()')
        end if

        deallocate(eval, stat=ierr)
        deallocate(evec, stat=ierr)

    end subroutine lanczos_diagonalisation
       
    subroutine hamil_vector(nrow, ncol, xin, ldx, yout, ldy)
 
        ! Matrix-vector multiplication procedure for use with trlan.
        ! In:
        !    nrow: the number of rows on this processor if the problem is distributed 
        !        using MPI, otherwise the number of total rows in a Lanczos vector. 
        !    ncol: the number of vectors (columns in xin and yout) to be multiplied. 
        !    xin: the array to store the input vectors to be multiplied.
        !    ldx: the leading dimension of the array xin when it is declared as 
        !       two-dimensional array.
        !    ldy: the leading dimension of the array yout when it is declared as 
        !       two-dimensional array.
        ! Out:
        !    yout: the array to store results of the multiplication.
 
        implicit None
        integer, intent(in) :: nrow, ncol, ldx, ldy
        real(dp), intent(in) :: xin(ldx,ncol)
        real(dp), intent(out) :: yout(ldy,ncol)
        ! local variables
        integer :: i
 
        do i = 1, ncol
            ! y = H x,
            ! where H is the Hamiltonian matrix, x is the input Lanczos vector
            ! and y the output Lanczos vector.
            call dsymv('U', nrow, 1.0_dp, hamil, nrow, xin(:,i), 1, 0.0_dp, yout(:,i), 1)
        end do
 
    end subroutine hamil_vector
    
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
        ! Apply Slater--Condon rules.
        case(0)

            root_det = decode_det(root_det_f)

            ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >

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

            ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

            ! One electron operator
            hmatel = hmatel + get_one_e_int(excitation%from_orb(1), excitation%to_orb(1)) 

            ! Two electron operator
            do i = 1, nel
                hmatel = hmatel + get_two_e_int(root_det(i), excitation%from_orb(1), root_det(i), excitation%to_orb(1))
            end do

            if (excitation%perm) hmatel = -hmatel

        case(2)

            ! < D | H | D_{ij}^{ab} > = < ij || ab >

            ! Two electron operator
            hmatel = get_two_e_int(excitation%from_orb(1), excitation%from_orb(2), excitation%to_orb(1), excitation%to_orb(2))

            if (excitation%perm) hmatel = -hmatel

        end select

    end function get_hmatel

end module hamiltonian
