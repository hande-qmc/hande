module utils

! Various utilities and tools...

implicit none

contains

    elemental function binom_i(m, n) result(binom)

        ! ACM Algorithm 160 translated to Fortran.
        ! Returns the binomial coefficient ^mC_n: the number of
        ! combinations of m things taken n at a time.

        integer :: binom
        integer, intent(in) :: m, n
        integer :: p,i, n1

        n1 = n
        p = m - n1
        if (n1 < p) then
            p = n1
            n1 = m - p
        end if
        binom = n1 + 1
        if (p == 0) binom = 1
        do i = 2, p
            binom = (binom*(n1+i))/i
        end do

    end function binom_i

    elemental function binom_r(m, n) result(binom)

        ! ACM Algorithm 160 translated to Fortran.
        ! Returns the binomial coefficient ^mC_n: the number of
        ! combinations of m things taken n at a time.

        use const, only: dp

        real(dp) :: binom
        integer, intent(in) :: m, n
        integer :: p,i, n1

        n1 = n
        p = m - n1
        if (n1 < p) then
            p = n1
            n1 = m - p
        end if
        binom = n1 + 1
        if (p == 0) binom = 1
        do i = 2, p
            binom = (binom*(n1+i))/i
        end do

    end function binom_r

    function get_free_unit() result(free_unit)

        ! Returns:
        !    The first free file unit above 10 and less than or equal to
        !    the paramater max_unit (currently set to 200).

        use errors, only: stop_all

        integer, parameter :: max_unit = 100
        integer :: free_unit
        integer :: i
        logical :: t_open, t_exist

        do i = 10, max_unit
            inquire(unit=i, opened=t_open, exist=t_exist)
            if (.not.t_open .and. t_exist) then
                free_unit = i
                exit
            end if
        end do
        if (i == max_unit+1) call stop_all('get_free_unit','Cannot find a free unit below max_unit.')

    end function get_free_unit

    elemental function int_fmt(i, padding) result(fmt1)
    
        ! In:
        !    i: an integer
        !    padding (optional): amount of padding to add to format statement.
        !        The default amount is 2.  The padding is used to include the
        !        sign if i is negative.
        ! Returns:
        !    fmt1: a format statement for an integer field which will hold
        !        i perfectly plus an amount of padding.
        
        ! This does take i/o formatting to a slightly OCD level addmittedly...

        character(4) :: fmt1
        integer, intent(in) :: i
        integer, intent(in), optional :: padding
        integer :: p
        real :: r

        if (present(padding)) then
            p = padding
        else
            p  = 2
        end if

        if (i == 0 .or. i==1) then
            r = 1.0
        else
            r = log10(real(abs(i)+1))
        end if
        p = ceiling(r+p)

        if (p < 10) then
            write (fmt1,'("i",i1)') p
        else if (p < 100) then
            write (fmt1,'("i",i2)') p
        else
            ! By this point we'll have hit integer overflow anyway...
            write (fmt1,'("i",i3)') p
        end if

    end function int_fmt

    subroutine append_ext(stem, n, s)

        ! Returns stem.n in s.

        character(*), intent(in) :: stem
        integer, intent(in) :: n
        character(*), intent(out) :: s
        character(10) :: ext

        write (ext,'('//int_fmt(n,0)//')') n
        s = stem//'.'//ext

    end subroutine append_ext

   subroutine get_unique_filename(stem, tnext, istart, filename)

        ! Find a filename which is either the "newest" or the next to be used.
        ! The filename is assumed to be stem.x, where x is an integer.

        ! In:
        !    stem: stem of the filename.
        !    tnext: the next unused filename is found if true, else the
        !        filename is set to be stem.x where stem.x exists and stem.x+1
        !        doesn't and x is greater than istart.
        !    istart: the integer of the first x value to check.
        !        If istart is negative, then the filename is set to be stem.x,
        !        where x = |istart+1|.  This overrides everything else.
        ! Out:
        !    filename.

        character(*), intent(in) :: stem
        logical, intent(in) :: tnext
        integer, intent(in) :: istart
        character(*), intent(out) :: filename

        integer :: i
        logical :: exists

        i = istart
        exists = .true.
        do while (exists)
            call append_ext(stem, i, filename)
            inquire(file=filename,exist=exists)
            i = i + 1
        end do

        if (.not.tnext) then
            ! actually want the last file which existed.
            ! this will return stem.istart if stem.istart doesn't exist.
            i = max(istart,i - 2)
            call append_ext(stem, i, filename)
        end if

        ! Have been asked for a specific file.
        if (istart < 0) then
            call append_ext(stem, abs(istart+1), filename)
        end if

    end subroutine get_unique_filename

    elemental function tri_ind(i,j) result(indx)

        ! Find the index corresponding to the (i,j)-th element of a lower
        ! triangular array.  This maps:
        !
        !   1,1                   1
        !   2,1 2,2               2  3
        !   3,1 3,2 3,3       to  4  5  6
        !   4,1 4,2 4,3 4,4       7  8  9 10
        !
        ! WARNING:
        ! We assume that i >= j.  It is the programmer's responsibility to check
        ! this and re-order i and j if required.
        !
        ! In:
        !    i: (1-indexed) row index
        !    j: (1-indexed) column index
        ! Returns:
        !    A combined (1-indexed) index for the corresponding element in
        !    a lower triangular array.

        integer :: indx
        integer, intent(in) :: i, j

        indx = ((i-1)*i)/2 + j

    end function tri_ind

end module utils
