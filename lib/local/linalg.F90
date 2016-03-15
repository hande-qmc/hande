module linalg

! Set of helper wrappers around LPAACK calls.

implicit none

private
public :: syev, heev, geev, psyev, pheev, plaprnt, syev_wrapper, &
            heev_wrapper, geev_wrapper, psyev_wrapper, pheev_wrapper

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

interface psyev
    module procedure pssyev_f90
    module procedure pdsyev_f90
end interface psyev

interface pheev
    module procedure pcheev_f90
    module procedure pzheev_f90
end interface pheev

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
