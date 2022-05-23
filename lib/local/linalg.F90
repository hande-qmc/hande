module linalg

! Set of helper wrappers around LAPACK calls.

implicit none

private
public :: syev, heev, geev, gemm, gemv, geqrf, orgqr, psyev, pheev, pgemm, pgeqrf, porgqr, plaprnt, syev_wrapper, &
            heev_wrapper, geev_wrapper, qr_wrapper, psyev_wrapper, pheev_wrapper, pqr_wrapper

interface syev
    module procedure ssyev_f90
    module procedure dsyev_f90
end interface

interface heev
    module procedure cheev_f90
    module procedure zheev_f90
end interface

interface geev
    module procedure sgeev_f90
    module procedure dgeev_f90
end interface geev

interface gemm
    module procedure sgemm_f90
    module procedure dgemm_f90
end interface gemm

interface gemv
    module procedure sgemv_f90
    module procedure dgemv_f90
end interface gemv

interface geqrf
    module procedure sgeqrf_f90
    module procedure dgeqrf_f90
end interface geqrf

interface orgqr
    module procedure sorgqr_f90
    module procedure dorgqr_f90
end interface orgqr

interface psyev
    module procedure pssyev_f90
    module procedure pdsyev_f90
end interface psyev

interface pheev
    module procedure pcheev_f90
    module procedure pzheev_f90
end interface pheev

interface pgemm
    module procedure psgemm_f90
    module procedure pdgemm_f90
end interface pgemm

interface pgeqrf
    module procedure psgeqrf_f90
    module procedure pdgeqrf_f90
end interface pgeqrf

interface porgqr
    module procedure psorgqr_f90
    module procedure pdorgqr_f90
end interface porgqr

interface plaprnt
    module procedure pslaprnt_f90
    module procedure pdlaprnt_f90
    module procedure pclaprnt_f90
    module procedure pzlaprnt_f90
end interface plaprnt

contains

    subroutine ssyev_f90(job, uplo, N, A, lda, W, work, lwork, info)

        ! See LAPACK ssyev procedure for details.

        use const, only: sp

        character, intent(in) :: job, uplo
        integer, intent(in) :: N, lda
        real(sp), intent(inout) :: A(lda,*)
        real(sp), intent(out) :: W(:), work(:)
        integer, intent(in) :: lwork
        integer, intent(out) :: info

        call ssyev(job, uplo, N, A, lda, W, work, lwork, info)

    end subroutine ssyev_f90

    subroutine dsyev_f90(job, uplo, N, A, lda, W, work, lwork, info)

        ! See LAPACK dsyev procedure for details.

        use const, only: dp

        character, intent(in) :: job, uplo
        integer, intent(in) :: N, lda
        real(dp), intent(inout) :: A(lda,*)
        real(dp), intent(out) :: W(:), work(:)
        integer, intent(in) :: lwork
        integer, intent(out) :: info

        call dsyev(job, uplo, N, A, lda, W, work, lwork, info)

    end subroutine dsyev_f90

    subroutine syev_wrapper(job, uplo, N, A, lda, W, info)

        ! Wrapper around ssyev/dsyev which automates using optimal workspace.

        use const, only: p
        use checking, only: check_allocate, check_deallocate

        character, intent(in) :: job, uplo
        integer, intent(in) :: N, lda
        real(p), intent(inout) :: A(lda,*)
        real(p), intent(out) :: W(:)
        integer, intent(out) :: info

        real(p), allocatable :: work(:)
        integer :: lwork, i, ierr

        lwork = -1
        do i = 1, 2
            allocate(work(abs(lwork)), stat=ierr)
            call check_allocate('work', abs(lwork), ierr)
            call syev(job, uplo, N, A, lda, W, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work)
            call check_deallocate('work',ierr)
        end do

    end subroutine syev_wrapper

    subroutine cheev_f90(job, uplo, N, A, lda, W, work, lwork, rwork, info)

        ! See LAPACK cheev procedure for details.

        use const, only: sp

        character, intent(in) ::  job, uplo
        integer, intent(in) :: N, lda
        complex(sp), intent(inout) :: A(lda,*)
        real(sp), intent(out) :: W(N)
        complex(sp), intent(out) :: work(:)
        integer, intent(in) :: lwork
        real(sp), intent(out) :: rwork(:)
        integer, intent(out) :: info

        call cheev(job, uplo, N, A, lda, W, work, lwork, rwork, info)

    end subroutine cheev_f90

    subroutine zheev_f90(job, uplo, N, A, lda, W, work, lwork, rwork, info)

        ! See LAPACK zheev procedure for details.

        use const, only: dp

        character, intent(in) ::  job, uplo
        integer, intent(in) :: N, lda
        complex(dp), intent(inout) :: A(lda,*)
        real(dp), intent(out) :: W(N)
        complex(dp), intent(out) :: work(:)
        integer, intent(in) :: lwork
        real(dp), intent(out) :: rwork(:)
        integer, intent(out) :: info

        call zheev(job, uplo, N, A, lda, W, work, lwork, rwork, info)

    end subroutine zheev_f90

    subroutine heev_wrapper(job, uplo, N, A, lda, W, rwork, info)

        ! Wrapper around cheev/zheev which automates using optimal workspace.

        use checking, only: check_allocate, check_deallocate
        use const, only: p

        character, intent(in) ::  job, uplo
        integer, intent(in) :: N, lda
        complex(p), intent(inout) :: A(lda,*)
        real(p), intent(out) :: W(N)
        real(p), intent(out) :: rwork(:)
        integer, intent(out) :: info

        complex(p), allocatable :: work(:)
        integer :: lwork, i, ierr

        lwork = -1
        do i = 1, 2
            allocate(work(abs(lwork)), stat=ierr)
            call check_allocate('work',abs(lwork),ierr)
            call heev(job, uplo, N, A, lda, W, work, lwork, rwork, info)
            lwork = nint(real(work(1)))
            deallocate(work)
            call check_deallocate('work',ierr)
        end do

    end subroutine heev_wrapper

    subroutine sgeev_f90(jobvl, jobvr, N, A, lda, WR, WI, VL, ldvl, VR, ldvr, work, lwork, info)

        ! See LAPACK sgeev procedure for details.

        use const, only: sp

        character, intent(in) :: jobvl, jobvr
        integer, intent(in) :: N, lda
        real(sp), intent(inout) :: A(lda,*)
        real(sp), intent(out) :: WR(:), WI(:)
        integer, intent(in) :: ldvl, ldvr, lwork
        real(sp), intent(out) :: VL(ldvl,*), VR(ldvr,*), WORK(:)
        integer, intent(out) :: info

        call sgeev(jobvl, jobvr, N, A, lda, WR, WI, VL, ldvl, VR, ldvr, work, lwork, info)

    end subroutine sgeev_f90

    subroutine dgeev_f90(jobvl, jobvr, N, A, lda, WR, WI, VL, ldvl, VR, ldvr, work, lwork, info)

        ! See LAPACK dgeev procedure for details.

        use const, only: dp

        character, intent(in) :: jobvl, jobvr
        integer, intent(in) :: N, lda
        real(dp), intent(inout) :: A(lda,*)
        real(dp), intent(out) :: WR(:), WI(:)
        integer, intent(in) :: ldvl, ldvr, lwork
        real(dp), intent(out) :: VL(ldvl,*), VR(ldvr,*), WORK(:)
        integer, intent(out) :: info

        call dgeev(jobvl, jobvr, N, A, lda, WR, WI, VL, ldvl, VR, ldvr, work, lwork, info)

    end subroutine dgeev_f90

    subroutine geev_wrapper(jobvl, jobvr, N, A, lda, WR, WI, VL, ldvl, VR, ldvr, info)

        ! Wrapper around sgeev/dgeev which automates using optimal workspace.

        use checking, only: check_allocate, check_deallocate
        use const, only: p

        character, intent(in) :: jobvl, jobvr
        integer, intent(in) :: N, lda
        real(p), intent(inout) :: A(lda,*)
        real(p), intent(out) :: WR(:), WI(:)
        integer, intent(in) :: ldvl, ldvr
        real(p), intent(out) :: VL(ldvl,*), VR(ldvr,*)
        integer, intent(out) :: info

        real(p), allocatable :: work(:)
        integer :: lwork, i, ierr

        lwork = -1
        do i = 1, 2
            allocate(work(abs(lwork)), stat=ierr)
            call check_allocate('work', abs(lwork), ierr)
            call geev(jobvl, jobvr, N, A, lda, WR, WI, VL, ldvl, VR, ldvr, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work)
            call check_deallocate('work',ierr)
        end do

    end subroutine geev_wrapper

    subroutine sgemm_f90(transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

        ! See LAPACK sgemm procedure for details

        use const, only: sp

        character, intent(in) :: transA, transB
        integer, intent(in) :: M, N, K, lda, ldb, ldc
        real(sp), intent(in) :: alpha, beta
        real(sp), intent(in) :: A(lda,*), B(ldb,*)
        real(sp), intent(inout) :: C(ldc,*)

        call sgemm(transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

    end subroutine sgemm_f90

    subroutine dgemm_f90(transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

        ! See LAPACK dgemm procedure for details

        use const, only: dp

        character, intent(in) :: transA, transB
        integer, intent(in) :: M, N, K, lda, ldb, ldc
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: A(lda,*), B(ldb,*)
        real(dp), intent(inout) :: C(ldc,*)

        call dgemm(transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)
        
    end subroutine dgemm_f90
    
    subroutine sgemv_f90(trans, M, N, alpha, A, lda, x, incx, beta, y, incy)

        ! See LAPACK sgemv procedure for details

        use const, only: sp

        character, intent(in) :: trans
        integer, intent(in) :: M, N, incx, incy, lda
        real(sp), intent(in) :: alpha, beta
        real(sp), intent(in) :: A(lda,*), x(*)
        real(sp), intent(inout) :: y(*)

        call sgemv(trans, M, N, alpha, A, lda, x, incx, beta, y, incy)

    end subroutine sgemv_f90

    subroutine dgemv_f90(trans, M, N, alpha, A, lda, x, incx, beta, y, incy)

        ! See LAPACK dgemv procedure for details

        use const, only: dp

        character, intent(in) :: trans
        integer, intent(in) :: M, N, incx, incy, lda
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: A(lda,*), x(*)
        real(dp), intent(inout) :: y(*)

        call dgemv(trans, M, N, alpha, A, lda, x, incx, beta, y, incy)

    end subroutine dgemv_f90

    subroutine sgeqrf_f90(M, N, A, lda, tau, work, lwork, info)

        ! See LAPACK sgeqrf procedure for details

        use const, only: sp

        integer, intent(in) :: M, N, lda, lwork
        real(sp), intent(inout) :: A(lda,*)
        real(sp), intent(out) :: tau(*)
        real(sp), intent(in) :: work(*)
        integer, intent(out) :: info

        call sgeqrf(M, N, A, lda, tau, work, lwork, info)

    end subroutine sgeqrf_f90

    subroutine dgeqrf_f90(M, N, A, lda, tau, work, lwork, info)

        ! See LAPACK dgeqrf procedure for details

        use const, only: dp

        integer, intent(in) :: M, N, lda, lwork
        real(dp), intent(inout) :: A(lda,*)
        real(dp), intent(out) :: tau(*)
        real(dp), intent(in) :: work(*)
        integer, intent(out) :: info

        call dgeqrf(M, N, A, lda, tau, work, lwork, info)

    end subroutine dgeqrf_f90

    subroutine sorgqr_f90(M, N, K, A, lda, tau, work, lwork, info)

        ! See LAPACK sorgqr procedure for details

        use const, only: sp

        integer, intent(in) :: M, N, K, lda, lwork
        real(sp), intent(inout) :: A(lda,*)
        real(sp), intent(out) :: tau(*)
        real(sp), intent(in) :: work(*)
        integer, intent(out) :: info

        call sorgqr(M, N, K, A, lda, tau, work, lwork, info)

    end subroutine sorgqr_f90

    subroutine dorgqr_f90(M, N, K, A, lda, tau, work, lwork, info)

        ! See LAPACK dorgqr procedure for details

        use const, only: dp

        integer, intent(in) :: M, N, K, lda, lwork
        real(dp), intent(inout) :: A(lda,*)
        real(dp), intent(out) :: tau(*)
        real(dp), intent(in) :: work(*)
        integer, intent(out) :: info

        call dorgqr(M, N, K, A, lda, tau, work, lwork, info)

    end subroutine dorgqr_f90

    subroutine qr_wrapper(M, N, A, lda, info)

        ! Wrapper around the geqrf (QR decomposition) and orgqr (extracting the Q matrix from the output of geqrf)
        ! subroutines, which automates using optimal workspace.

        use checking, only: check_allocate, check_deallocate
        use const, only: p

        integer, intent(in) :: M, N, lda
        real(p), intent(inout) :: A(lda,*)
        integer, intent(out) :: info

        integer :: lwork, ierr, i
        real(p), allocatable :: tau(:), work(:)

        allocate(tau(min(M,N)), stat=ierr)
        call check_allocate('tau', min(M,N), ierr)

        lwork = -1
        do i = 1, 2
            allocate(work(abs(lwork)), stat=ierr)
            call check_allocate('work', abs(lwork), ierr)
            call geqrf(M, N, A, lda, tau, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work, stat=ierr)
            call check_deallocate('work', ierr)
        end do

        lwork = -1
        do i = 1, 2
            allocate(work(abs(lwork)), stat=ierr)
            call check_allocate('work', abs(lwork), ierr)
            call orgqr(M, N, min(M,N), A, lda, tau, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work)
            call check_deallocate('work', ierr)
        end do

    end subroutine qr_wrapper

    subroutine pssyev_f90(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, info)

        ! See ScaLAPACK pdsyev procedure for details.

        use const, only: sp
        use errors, only: stop_all

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        real(sp), intent(inout) :: A(:,:)
        integer, intent(in) :: IA, JA, desca(:)
        real(sp), intent(out) :: W(:), Z(:,:)
        integer, intent(in) :: IZ, JZ, descz(:)
        real(sp), intent(out) :: work(:)
        integer, intent(in) :: lwork
        integer, intent(out) :: info

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pssyev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, info)
#else
        call stop_all('pssyev_f90', "Scalapack disabled at compile-time")
#endif

    end subroutine pssyev_f90

    subroutine pdsyev_f90(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, info)

        ! See ScaLAPACK pdsyev procedure for details.

        use const, only: dp
        use errors, only: stop_all

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        real(dp), intent(inout) :: A(:,:)
        integer, intent(in) :: IA, JA, desca(:)
        real(dp), intent(out) :: W(:), Z(:,:)
        integer, intent(in) :: IZ, JZ, descz(:)
        real(dp), intent(out) :: work(:)
        integer, intent(in) :: lwork
        integer, intent(out) :: info

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pdsyev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, info)
#else
        call stop_all('pdsyev_f90', "Scalapack disabled at compile-time")
#endif

    end subroutine pdsyev_f90

    subroutine psyev_wrapper(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, info)

        ! Wrapper around pssyev/pdsyev which automates using optimal workspace.

        use checking, only: check_allocate, check_deallocate
        use const, only: p
        use errors, only: stop_all

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        real(p), intent(inout) :: A(:,:)
        integer, intent(in) :: IA, JA, desca(:)
        real(p), intent(out) :: W(:), Z(:,:)
        integer, intent(in) :: IZ, JZ, descz(:)
        integer, intent(out) :: info

        real(p), allocatable :: work(:)
        integer :: lwork, i, ierr

        lwork = -1
        do i = 1, 2
            allocate(work(abs(lwork)), stat=ierr)
            call check_allocate('work', abs(lwork), ierr)
            call psyev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work, stat=ierr)
            call check_deallocate('work',ierr)
        end do

    end subroutine psyev_wrapper

    subroutine pcheev_f90(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, rwork, lrwork, info)

        ! See ScaLAPACK pcheev procedure for details.

        use const, only: sp
        use errors, only: stop_all

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        complex(sp), intent(inout) :: A(:,:)
        integer, intent(in) :: IA, JA, desca(:)
        real(sp), intent(out) :: W(:)
        complex(sp), intent(out) :: Z(:,:)
        integer, intent(in) :: IZ, JZ, descz(:)
        complex(sp), intent(out) :: work(:)
        integer, intent(in) :: lwork
        real(sp), intent(out) :: rwork(:)
        integer, intent(in) :: lrwork
        integer, intent(out) :: info

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pcheev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, rwork, lrwork, info)
#else
        call stop_all('pcheev_f90', "Scalapack disabled at compile-time")
#endif

    end subroutine pcheev_f90

    subroutine pzheev_f90(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, rwork, lrwork, info)

        ! See ScaLAPACK pzheev procedure for details.

        use const, only: dp
        use errors, only: stop_all

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        complex(dp), intent(inout) :: A(:,:)
        integer, intent(in) :: IA, JA, desca(:)
        real(dp), intent(out) :: W(:)
        complex(dp), intent(out) :: Z(:,:)
        integer, intent(in) :: IZ, JZ, descz(:)
        complex(dp), intent(out) :: work(:)
        integer, intent(in) :: lwork
        real(dp), intent(out) :: rwork(:)
        integer, intent(in) :: lrwork
        integer, intent(out) :: info

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pzheev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, rwork, lrwork, info)
#else
        call stop_all('pzheev_f90', "Scalapack disabled at compile-time")
#endif

    end subroutine pzheev_f90

    subroutine pheev_wrapper(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, rwork, lrwork, info)

        ! Wrapper around pcheev/pzheev which automates using optimal workspace.

        use checking, only: check_allocate, check_deallocate
        use const, only: p

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        complex(p), intent(inout) :: A(:,:)
        integer, intent(in) :: IA, JA, desca(:)
        real(p), intent(out) :: W(:)
        complex(p), intent(out) :: Z(:,:)
        integer, intent(in) :: IZ, JZ, descz(:)
        real(p), intent(out) :: rwork(:)
        integer, intent(in) :: lrwork
        integer, intent(out) :: info

        complex(p), allocatable :: work(:)
        integer :: lwork, i, ierr

        lwork = -1
        do i = 1,2
            allocate(work(abs(lwork)), stat=ierr)
            call check_allocate('work', abs(lwork), ierr)
            call pheev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, rwork, lrwork, info)
            lwork = nint(real(work(1)))
            deallocate(work, stat=ierr)
            call check_deallocate('work',ierr)
        end do

    end subroutine pheev_wrapper

    subroutine psgemm_f90(transA, transB, M, N, K, alpha, A, ia, ja, desca, B, ib, jb, descb, beta, C, ic, jc, descc)

        ! See ScaLAPACK psgemm procedure for details

        use errors, only: stop_all
        use const, only: sp

        character, intent(in) :: transA, transB
        integer, intent(in) :: M, N, K, ia, ja, ib, jb, ic, jc, desca(:), descb(:), descc(:)
        real(sp), intent(in) :: alpha, beta
        real(sp), intent(in) :: A(:,:), B(:,:)
        real(sp), intent(inout) :: C(:,:)

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call psgemm(transA, transB, M, N, K, alpha, A, ia, ja, desca, B, ib, jb, descb, beta, C, ic, jc, descc)
#else
        call stop_all('psgemm_f90', 'Scalapack disabled at compile-time')
#endif

    end subroutine psgemm_f90

    subroutine pdgemm_f90(transA, transB, M, N, K, alpha, A, ia, ja, desca, B, ib, jb, descb, beta, C, ic, jc, descc)

        ! See ScaLAPACK pdgemm procedure for details

        use errors, only: stop_all
        use const, only: dp

        character, intent(in) :: transA, transB
        integer, intent(in) :: M, N, K, ia, ja, ib, jb, ic, jc, desca(:), descb(:), descc(:)
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: A(:,:), B(:,:)
        real(dp), intent(inout) :: C(:,:)

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pdgemm(transA, transB, M, N, K, alpha, A, ia, ja, desca, B, ib, jb, descb, beta, C, ic, jc, descc)
#else
        call stop_all('pdgemm_f90', 'Scalapack disabled at compile-time')
#endif
        
    end subroutine pdgemm_f90

    subroutine psgeqrf_f90(M, N, A, ia, ja, desca, tau, work, lwork, info)

        ! See ScaLAPACK psgeqrf procedure for details

        use errors, only: stop_all
        use const, only: sp

        integer, intent(in) :: M, N, ia, ja, desca(:), lwork
        real(sp), intent(inout) :: A(:,:)
        real(sp), intent(out) :: tau(:)
        real(sp), intent(in) :: work(:)
        integer, intent(out) :: info

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call psgeqrf(M, N, A, ia, ja, desca, tau, work, lwork, info)
#else
        call stop_all('psgeqrf_f90', 'Scalapack disabled at compile-time')
#endif

    end subroutine psgeqrf_f90

    subroutine pdgeqrf_f90(M, N, A, ia, ja, desca, tau, work, lwork, info)

        ! See ScaLAPACK pdgeqrf procedure for details

        use errors, only: stop_all
        use const, only: dp

        integer, intent(in) :: M, N, ia, ja, desca(:), lwork
        real(dp), intent(inout) :: A(:,:)
        real(dp), intent(out) :: tau(:)
        real(dp), intent(in) :: work(:)
        integer, intent(out) :: info

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pdgeqrf(M, N, A, ia, ja, desca, tau, work, lwork, info)
#else
        call stop_all('pdgeqrf_f90', 'Scalapack disabled at compile-time')
#endif

    end subroutine pdgeqrf_f90

    subroutine psorgqr_f90(M, N, K, A, ia, ja, desca, tau, work, lwork, info)

        ! See ScaLAPACK psorgqr procedure for details

        use errors, only: stop_all
        use const, only: sp

        integer, intent(in) :: M, N, K, ia, ja, desca(:), lwork
        real(sp), intent(inout) :: A(:,:)
        real(sp), intent(out) :: tau(:)
        real(sp), intent(in) :: work(:)
        integer, intent(out) :: info

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call psorgqr(M, N, K, A, ia, ja, desca, tau, work, lwork, info)
#else
        call stop_all('psorgqr_f90', 'Scalapack disabled at compile-time')
#endif

    end subroutine psorgqr_f90

    subroutine pdorgqr_f90(M, N, K, A, ia, ja, desca, tau, work, lwork, info)

        ! See ScaLAPACK pdorgqr procedure for details

        use errors, only: stop_all
        use const, only: dp

        integer, intent(in) :: M, N, K, ia, ja, desca(:), lwork
        real(dp), intent(inout) :: A(:,:)
        real(dp), intent(out) :: tau(:)
        real(dp), intent(in) :: work(:)
        integer, intent(out) :: info

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pdorgqr(M, N, K, A, ia, ja, desca, tau, work, lwork, info)
#else
        call stop_all('pdorgqr_f90', 'Scalapack disabled at compile-time')
#endif

    end subroutine pdorgqr_f90

    subroutine pqr_wrapper(M, N, A, ia, ja, desca, info)

        ! Wrapper around the pgeqrf (QR decomposition) and porgqr (extracting the Q matrix from the output of pgeqrf)
        ! subroutines, which automates using optimal workspace.

        use checking, only: check_allocate, check_deallocate
        use const, only: p

        integer, intent(in) :: M, N, ia, ja, desca(:)
        real(p), intent(inout) :: A(:,:)
        integer, intent(out) :: info

        integer :: lwork, ierr, i
        real(p), allocatable :: tau(:), work(:)

        allocate(tau(min(M,N)), stat=ierr)
        call check_allocate('tau', min(M,N), ierr)

        lwork = -1
        do i = 1, 2
            allocate(work(abs(lwork)), stat=ierr)
            call check_allocate('work', abs(lwork), ierr)
            call pgeqrf(M, N, A, ia, ja, desca, tau, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work, stat=ierr)
            call check_deallocate('work', ierr)
        end do

        lwork = -1
        do i = 1, 2
            allocate(work(abs(lwork)), stat=ierr)
            call check_allocate('work', abs(lwork), ierr)
            call porgqr(M, N, min(M,N), A, ia, ja, desca, tau, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work)
            call check_deallocate('work', ierr)
        end do

    end subroutine pqr_wrapper

    subroutine pslaprnt_f90(M, N, A, IA, JA, desca, irprnt, icprnt, cmatnm, nout, work)

        ! See ScaLAPACK's pslaprnt for details.

        use const, only: sp
        use errors, only: stop_all

        integer, intent(in) :: M, N, IA, JA, desca(:), irprnt, icprnt, nout
        character(*), intent(in) :: cmatnm
        real(sp), intent(in) :: A(:,:)
        real(sp), intent(out) :: work(:)

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pslaprnt(M, N, A, IA, JA, desca, irprnt, icprnt, cmatnm, nout, work)
#else
        call stop_all('pslaprnt_f90', "Scalapack disabled at compile-time")
#endif

    end subroutine pslaprnt_f90

    subroutine pdlaprnt_f90(M, N, A, IA, JA, desca, irprnt, icprnt, cmatnm, nout, work)

        ! See ScaLAPACK's pdlaprnt for details.

        use const, only: dp
        use errors, only: stop_all

        integer, intent(in) :: M, N, IA, JA, desca(:), irprnt, icprnt, nout
        character(*), intent(in) :: cmatnm
        real(dp), intent(in) :: A(:,:)
        real(dp), intent(out) :: work(:)

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pdlaprnt(M, N, A, IA, JA, desca, irprnt, icprnt, cmatnm, nout, work)
#else
        call stop_all('pdlaprnt_f90', "Scalapack disabled at compile-time")
#endif

    end subroutine pdlaprnt_f90

    subroutine pclaprnt_f90(M, N, A, IA, JA, desca, irprnt, icprnt, cmatnm, nout, work)

        ! See ScaLAPACK's pclaprnt for details.

        use const, only: sp
        use errors, only: stop_all

        integer, intent(in) :: M, N, IA, JA, desca(:), irprnt, icprnt, nout
        character(*), intent(in) :: cmatnm
        complex(sp), intent(in) :: A(:,:)
        complex(sp), intent(out) :: work(:)

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pclaprnt(M, N, A, IA, JA, desca, irprnt, icprnt, cmatnm, nout, work)
#else
        call stop_all('pclaprnt_f90', "Scalapack disabled at compile-time")
#endif

    end subroutine pclaprnt_f90

    subroutine pzlaprnt_f90(M, N, A, IA, JA, desca, irprnt, icprnt, cmatnm, nout, work)

        ! See ScaLAPACK's pzlaprnt for details.

        use const, only: dp
        use errors, only: stop_all

        integer, intent(in) :: M, N, IA, JA, desca(:), irprnt, icprnt, nout
        character(*), intent(in) :: cmatnm
        complex(dp), intent(in) :: A(:,:)
        complex(dp), intent(out) :: work(:)

#if defined(PARALLEL) && ! defined(DISABLE_SCALAPACK)
        call pzlaprnt(M, N, A, IA, JA, desca, irprnt, icprnt, cmatnm, nout, work)
#else
        call stop_all('pzlaprnt_f90', "Scalapack disabled at compile-time")
#endif

    end subroutine pzlaprnt_f90

end module linalg
