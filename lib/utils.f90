module utils

! Various utilities and tools...

implicit none

contains

    function get_free_unit() result(free_unit)

        ! Returns:
        !    The first free file unit above 10 and less than or equal to
        !    the paramater max_unit (currently set to 200).

        use errors, only: stop_all

        integer, parameter :: max_unit = 200
        integer :: free_unit
        integer :: i
        logical :: t_open

        do i = 10, max_unit
            inquire(unit=i, opened=t_open)
            if (.not.t_open) then
                free_unit = i
                exit
            end if
        end do
        if (i == max_unit+1) call stop_all('get_free_unit','Cannot find a free unit below 200')

    end function get_free_unit

    elemental function int_fmt(i, padding) result(fmt1)
    
        ! In:
        !    i: an integer
        !    padding (optional): amount of padding to add to format statement.
        !        The default amount is 2.
        ! Returns:
        !    fmt1: a format statement for an integer field which will hold
        !        i perfectly plus an amount of padding.
        
        ! This does take i/o formatting to a slightly OCD level addmittedly...

        character(2) :: fmt1
        integer, intent(in) :: i
        integer, intent(in), optional :: padding
        integer :: p
        real :: r

        if (present(padding)) then
            p = padding
        else
            p  = 2
        end if

        r = log10(real(i))
        if (r < 10) then
            write (fmt1,'("i",i1)') nint(r+p)
        else if (r < 100) then
            write (fmt1,'("i",i2)') nint(r+p)
        else
            ! By this point we'll have hit integer overflow anyway...
            write (fmt1,'("i",i3)') nint(r+p)
        end if

    end function int_fmt

end module utils
