module linalg

! Set of helper wrappers around LPAACK calls.

implicit none

private
public :: syev, heev, geev, psyev, pheev

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

contains

    subroutine ssyev_f90(job, uplo, N, A, lda, W, work, lwork, info)

        ! See LAPACK ssyev procedure for details.

        use const, only: sp

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        real(sp), intent(inout) :: A(lda,*)
        integer, intent(in) :: lda
        real(sp), intent(out) :: W(:), work(:)
        integer, intent(in) :: lwork
        integer, intent(out) :: info

        call ssyev(job, uplo, N, A, lda, W, work, lwork, info)

    end subroutine ssyev_f90

    subroutine dsyev_f90(job, uplo, N, A, lda, W, work, lwork, info)

        ! See LAPACK dsyev procedure for details.

        use const, only: dp

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        real(dp), intent(inout) :: A(lda,*)
        integer, intent(in) :: lda
        real(dp), intent(out) :: W(:), work(:)
        integer, intent(in) :: lwork
        integer, intent(out) :: info

        call dsyev(job, uplo, N, A, lda, W, work, lwork, info)

    end subroutine dsyev_f90

    subroutine cheev_f90(job, uplo, N, A, lda, W, work, lwork, rwork, info)

        ! See LAPACK cheev procedure for details.

        use const, only: sp

        character, intent(in) ::  job, uplo
        integer, intent(in) :: N
        complex(sp), intent(inout) :: A(lda,*)
        integer, intent(in) :: lda
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
        integer, intent(in) :: N
        complex(dp), intent(inout) :: A(lda,*)
        integer, intent(in) :: lda
        real(dp), intent(out) :: W(N)
        complex(dp), intent(out) :: work(:)
        integer, intent(in) :: lwork
        real(dp), intent(out) :: rwork(:)
        integer, intent(out) :: info

        call zheev(job, uplo, N, A, lda, W, work, lwork, rwork, info)

    end subroutine zheev_f90

    subroutine sgeev_f90(jobvl, jobvr, N, A, lda, WR, WI, VL, ldvl, VR, ldvr, work, lwork, info)

        ! See LAPACK sgeev procedure for details.

        use const, only: sp

        character, intent(in) :: jobvl, jobvr
        integer, intent(in) :: N
        real(sp), intent(inout) :: A(lda,*)
        integer, intent(in) :: LDA
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
        integer, intent(in) :: N
        real(dp), intent(inout) :: A(lda,*)
        integer, intent(in) :: LDA
        real(dp), intent(out) :: WR(:), WI(:)
        integer, intent(in) :: ldvl, ldvr, lwork
        real(dp), intent(out) :: VL(ldvl,*), VR(ldvr,*), WORK(:)
        integer, intent(out) :: info

        call dgeev(jobvl, jobvr, N, A, lda, WR, WI, VL, ldvl, VR, ldvr, work, lwork, info)

    end subroutine dgeev_f90

    subroutine pssyev_f90(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, info)

        ! See ScaLAPACK pdsyev procedure for details.

        use const, only: sp

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        real(sp), intent(inout) :: A(:,:)
        integer, intent(in) :: IA, JA, desca(:)
        real(sp), intent(out) :: W(:), Z(:,:)
        integer, intent(in) :: IZ, JZ, descz(:)
        real(sp), intent(out) :: work(:)
        integer, intent(in) :: lwork
        integer, intent(out) :: info

#ifdef PARALLEL
        call pssyev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, info)
#endif

    end subroutine pssyev_f90

    subroutine pdsyev_f90(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, info)

        ! See ScaLAPACK pdsyev procedure for details.

        use const, only: dp

        character, intent(in) :: job, uplo
        integer, intent(in) :: N
        real(dp), intent(inout) :: A(:,:)
        integer, intent(in) :: IA, JA, desca(:)
        real(dp), intent(out) :: W(:), Z(:,:)
        integer, intent(in) :: IZ, JZ, descz(:)
        real(dp), intent(out) :: work(:)
        integer, intent(in) :: lwork
        integer, intent(out) :: info

#ifdef PARALLEL
        call pdsyev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, info)
#endif

    end subroutine pdsyev_f90

    subroutine pcheev_f90(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, rwork, lrwork, info)

        ! See ScaLAPACK pcheev procedure for details.

        use const, only: sp

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

#ifdef PARALLEL
        call pcheev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, rwork, lrwork, info)
#endif

    end subroutine pcheev_f90

    subroutine pzheev_f90(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, rwork, lrwork, info)

        ! See ScaLAPACK pzheev procedure for details.

        use const, only: dp

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

#ifdef PARALLEL
        call pzheev(job, uplo, N, A, IA, JA, desca, W, Z, IZ, JZ, descz, work, lwork, rwork, lrwork, info)
#endif

    end subroutine pzheev_f90

end module linalg
