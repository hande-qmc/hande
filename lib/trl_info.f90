! $Id: trl_info.f90,v 1.6 1999/09/23 16:20:21 kewu Exp $
!!!
! a module that defines TRL_INFO_T type,
! TRL_INFO_T stors information for TRLAN which includes the definition
! of eigenvalue problem, the size of the problem, convergence
! tolerance, time used by various parts of the Lanczos routine, whether
! to print debug information and where to print them, etc.
!!!
Module trl_info
  Type TRL_INFO_T
     Integer :: stat    ! status  (error code) of TRLAN

     ! specification of eigenpairs wanted
     Integer :: lohi    ! which end of spectrum to compute
                        ! lohi < 0 --> the smallest eigenvalues
                        ! lohi = 0 --> whichever converge first
                        ! lohi > 0 --> the largest eigenvalues
     Integer :: ned     ! number of eigenpairs wanted
     Integer :: nec     ! number of eigenpairs converged
                        ! if the user has nec correct eigenvectors, then
                        ! they are expected to be stored at the beginning
                        ! of the eigenvector array
     Real(8) :: tol     ! an eigenpair is declared converged if its
                        ! residual norm is less than tol*||OP||

     ! specification of resource allowed to use by TRLAN
     Integer :: mpicom  ! the MPI communicator
     Integer :: maxlan  ! maximum basis size to be used
     Integer :: klan    ! the actual basis size, this value may be smaller
                        ! than maxlan.  It is set when restarting.
     Integer :: maxmv   ! maximum number of MATVEC allowed
                        ! one MATVEC == one application of the operator on
                        ! one vector

     Integer :: restart ! index of restarting schemes
     Integer :: locked  ! number of eigenvalue locked
     Integer :: guess   ! initial guess
     ! <= 0, user did not provide initial guess, use [1,1,..,1]
     ! =  1, user has supplied initial guess, will only use the first one
     ! >  1, restart with previous check-point file

     ! some information about the progress and resouce comsumption
     Integer :: matvec  ! number of MATVEC used by TRLAN
     Integer :: nloop   ! number of restart of the Lanczos iterations
     Integer :: north   ! number of full orthogonalization invoked
     Integer :: nrand   ! number of times random elements are introduced.
                        ! Random elements are introduced when an invariant
                        ! subspace is found but not all wanted eigenvalues
                        ! are computed.
     Integer :: flop    ! Floating-point operations count (EXCLUDING MATVEC)
     Integer :: flop_h  ! FLOPS used in re-orthogonalization
     Integer :: flop_r  ! FLOPS used in restarting
     Real(8) :: rflp
     Real(8) :: rflp_h
     Real(8) :: rflp_r

     ! variables to store timing results
     Integer :: clk_rate! system clock rate (SYSTEM_CLOCK)
     Integer :: clk_max ! maximum counter value
     Integer :: clk_tot ! total time spent in TRLAN (in clock ticks)
     Integer :: clk_op  ! time in applying the operator (MATVEC)
     Integer :: clk_orth! time in re-orthogonalization
     Integer :: clk_res ! time in restarting the Lanczos iterations
     Real(8) :: tick_t  ! the sum of clk_tot and tick_t is the actual time
     Real(8) :: tick_o
     Real(8) :: tick_h
     Real(8) :: tick_r
     Integer :: clk_in  ! time spent in reading input data file
     Integer :: wrds_in ! number of real(8) words read
     Integer :: clk_out ! time spent in writing output data file
     Integer :: wrds_out! number of real(8) words written to file

     Real(8) :: anrm    ! norm of the operator used

     Integer :: my_pe   ! the PE number of current processor (start with 0)
     Integer :: npes    ! number of PEs in the group
     Integer :: nloc    ! local problem size
     Integer :: ntot    ! global problem size

     ! how much inforation to output during the execution of TRLAN
     Integer :: verbose ! default only print information related to
                        ! fatal errors
                        ! if verbose > 0, more diagnostic messages
                        ! are printed to the following files
     Integer :: log_io  ! Fortran I/O unit number to be used for
                        ! debug output. Used if verbose > 0.
     Character(128) :: log_file
     ! base of the file names used by TRLAN to store debug information if
     ! verbose > 0.  The filenames are computed by appending 'PE#' to
     ! this basis.

     ! check-pointing parameters
     ! when cpflag is greater than 0, the basis vectors will be written
     ! out roughly 'cpflag' times.  For simplicitly, each PE writes its
     ! own portion of the basis vectors to a file with cpfile followed by
     ! the processor number.  The file is written as unformatted fortran
     ! files with the following content:
     ! nrow, kb(basis size)
     ! alpha(1:kb)
     ! beta(1:kb)
     ! 1st basis vector, 2nd basis vector, ..., last basis vector
     ! the residual vector
     Integer :: cpflag, cpio
     Character(128) :: cpfile, oldcpf

     ! variables needed to measure convergence factor (crat)
     ! convergence rate of the restarted Lanczos algorithm is measure by
     ! the reduction of residual norm per MATVEC.  The residual norm of
     ! the target is used.
     Real(8) :: crat
     Real(8) :: trgt    ! the Ritz value that might convege next
     Real(8) :: tres    ! residual norm of the target
     Integer :: tmv     ! MATVEC used when target and tres were recorded

     Integer :: dummy   ! a dummy variable to fill space
  End Type TRL_INFO_T
End Module trl_info
!!!
! the interface module -- stores the external user interfaces defined
! in the package
!!!
Module trl_interface
  Interface
!! initialization routine for trl_info_t
!
! Arguments:
! info   -- the information package to be used by TRLAN
! nrow   -- local dimension of the problem
! mxlan  -- maximum number of basis vectors to be used
! lohi   -- which end of the spectrum to compute
!           <0 : lower end, the smallest eigenvalues
!           >0 : high end, the largest eigenvalues
!            0 : either lower and or high end, whoever converges
!                first 
! ned    -- number of wanted eigenvalues and eigenvectors
! mxmv   -- (optional) maximum number of matrix-vector multiplications
!           allowed
!           default: info%ntot*info%ned
! trestart-- (optional) thick-restarting scheme, 1-4,
!           default: 1
! tol    -- (optional) tolerance on residual norm,
!           default: sqrt(epsilon)
! mpicom -- (optional) the MPI communicator,
!           default: a duplicate of MPI_COMM_WORLD
!           in sequential case, this is set to 0.
     Subroutine trl_init_info(info, nrow, mxlan, lohi, ned, tol,&
          & trestart, maxmv, mpicom)
       Use trl_info
       Implicit None 
       Integer, Intent(in) :: lohi, mxlan, ned, nrow
       Integer, Intent(in), Optional :: maxmv, mpicom, trestart
       Real(8), Intent(in), Optional :: tol
       Type(TRL_INFO_T), Intent(out) :: info
     End Subroutine trl_init_info
! the main program of TRLAN package
!
!\Arguments:
! op   -- the operator
! nrow -- number of rows that is one this processor (local problem size)
! mev  -- number of elements in eval and number of columns in evec
! eval -- the eigenvalues
! evec -- the eigenvectors
! lde  -- the leading dimension of evec
! info -- data structure to store the information about the eigenvalue
!         problem and the progress of TRLAN
! wrk  -- (optional) workspace, if it is provided and there is enough
!         space, the residual norm of the converged eigenpairs will be
!         stored at wrk(1:info%nec).
! lwrk -- (optional) the size of the optional workspace
     Subroutine trlan(op, info, nrow, mev, eval, evec, lde, wrk, lwrk)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T) :: info
       Integer, Intent(in) :: nrow, mev, lde
       Real(8) :: eval(mev), evec(lde,mev)
       Integer, Intent(in), Optional :: lwrk
       Real(8), Target, Optional :: wrk(*)
       External op
       !  Interface
       !     Subroutine OP(nrow, ncol, xin, ldx, yout, ldy)
       !       Integer, Intent(in) :: nrow, ncol, ldx, ldy
       !       Real(8), Dimension(ldx*ncol), Intent(in) :: xin
       !       Real(8), Dimension(ldy*ncol), Intent(out) :: yout
       !     End Subroutine OP
       !  End Interface
     End Subroutine trlan
! function to print relavent entries of info to Fortran I/O unit iou
! Arguments
!
! ** NOTE *** MUST be called on all PEs ***
!
! info  -- the TRL_INFO_T variable to be printed
! mvflop-- the number of floating-operations per MATVEC.  This
!          information has to be supplied by user, otherwise related
!          entries are left blank in the print out.
     Subroutine trl_print_info(info, mvflop)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in), Optional :: mvflop
     End Subroutine trl_print_info
! a terse version of trl_print_info
! this is a local routine, indivadual PE can call it without regard of
! whether other PEs do the same
     Subroutine trl_terse_info(info, iou)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: iou
     End Subroutine trl_terse_info
!!!
! trl_check_ritz
! performs a standard check on the computed Ritz pairs
!
! Arguments:
! op    -- the operator, aka the matrix-vector multiplication routine
! nrow  -- the problem size
! rvec  -- the array storing the Ritz vectors
! alpha -- The ritz values
! beta  -- The residual norms returned from a Lanczos routine (optional)
! eval  -- The actual eigenvalues (optional)
! wrk   -- Optional workspace
!!!
     Subroutine trl_check_ritz(op, info, nrow, rvec, alpha, beta, eval, wrk)
       Use trl_info
       Implicit None
       External op
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nrow
       Real(8), Dimension(:,:), Intent(in) :: rvec
       Real(8), Dimension(:), Intent(in) :: alpha
       Real(8), Dimension(:), Intent(in), Optional :: beta, eval
       Real(8), Dimension(:), Optional, Target :: wrk
     End Subroutine trl_check_ritz
!!!
! set information related to debugging, trl_init_info does not allow
! any debug information to be printed.
!
! MUST CALL TRL_INIT_INFO BEFORE CALLING THIS SUBROUTINE
!!!
! arguments:
! info   -- the info to be modified
! msglvl -- the new message level
!           < 0 : nothing printed
!           1:10 -- the high the level, the more debug information is
!                printed
! iou    -- Fortran I/O unit to be used to write out the debug
!           information
! file   -- leading part of the debug file name.
!           Each processor will generate its own debug file with the
!           file name formed by appending the processor number to
!           string file
     Subroutine trl_set_debug(info, msglvl, iou, file)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T) :: info
       Integer, Intent(in) :: msglvl, iou
       Character(*), Optional :: file
     End Subroutine trl_set_debug
!!!
! set check-pointing parameters
! it also checks the existence of check-point files if they are to be
! read.
!
! MUST CALL TRL_INIT_INFO BEFORE CALLING THIS SUBROUTINE
!!!
! arguments:
! info   -- the info to be modified
! cpflag -- check-pointing flag (how many times to write)
! cpio   -- the IO unit number to use
! file   -- the prefix of the check-point files
     Subroutine trl_set_checkpoint(info, cpflag, cpio, file)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T) :: info
       Integer, Intent(in) :: cpflag, cpio
       Character(*), Optional :: file
     End Subroutine trl_set_checkpoint
!!!
! set up parameters related to initial guesses of the Lanczos iterations
!
! It sets parameters and also tries to make sure that the checkpoint
! files exist when iguess>1 and oldcpf is specified.
!
! It uses info%cpio to open the checkpoint files.  If the default value
! of info%cpio is in use, make sure that the routine is called after
! trl_set_checkpoint function is called.
!
! MUST CALL TRL_INIT_INFO BEFORE CALLING THIS SUBROUTINE
!!!
     Subroutine trl_set_iguess(info, nec, iguess, oldcpf)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T) :: info
       Integer, Intent(in) :: iguess, nec
       Character(*), Intent(in), Optional :: oldcpf
     End Subroutine trl_set_iguess
!!!
! a subroutine to compute Rayleigh quotient
! given a set of Ritz vectors and Ritz values, normalize the Ritz
! vectors and compute their Rayleigh quotients to replace the Ritz
! values.
!
! Arguments:
! op   -- the matrix-vector multiplication routine
! info -- the data structure that stores infomation for TRLAN
! evec -- the array to store the portion of eigenvectors on this PE
! base -- the workspace used to store results of MATVEC
! eres -- store new Ritz values and new Rresidual norms
!         if there are NEV Ritz pairs, eres(1:NEV) stores the new
!         Rayleigh quotient and eres(nev+1:nev+nev) stores the
!         new residual norms.
!!!
     Subroutine trl_rayleigh_quotients(op, info, evec, eres, base)
       Use trl_info
       Implicit None
       External op
       Type(TRL_INFO_T) :: info
       Real(8), Dimension(:) :: eres
       Real(8), Dimension(:,:) :: evec
       Real(8), Dimension(:), Target, Optional :: base
     End Subroutine trl_rayleigh_quotients
!!!
! A separated Rayleigh-Ritz projection routine
! Given a set of approximately orthonormal vectors (V), this routine
! performs the following operations
! (1) V'*V ==> G
! (2) R'*R :=  G
! (3) V'*A*V => H1, inv(R')*H1*inv(R) => H
! (4) Y*D*Y' := H
! (5) V*inv(R)*Y => V, diag(D) => lambda,
!     r(i) = ||A*V(:,i)-lambda(i)*V(:,i)||
!
! Arguments:
! op   -- the operator (matrix-vector multiplication routine)
! info -- the structure for storing information about TRLAN
! evec -- the Ritz vectors to be manipulated
! eres -- the array to store new Ritz values and residual norms
! base -- the (optional) workspace to store the result of MATVEC
! wrk  -- the (optional) workspace to store projection matrix, etc..
!!!
     Subroutine trl_ritz_projection(op, info, evec, eres, wrk, base)
       Use trl_info
       External op
       Type(TRL_INFO_T) :: info
       Real(8), Dimension(:), Optional, Target :: wrk, base
       Real(8), Dimension(:, :) :: evec
       Real(8), Dimension(:) :: eres
     End Subroutine trl_ritz_projection
  End Interface
End Module trl_interface
