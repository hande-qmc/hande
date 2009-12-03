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
            if (t_open) then
                free_unit = i
                exit
            end if
        end do
        if (i == max_unit+1) call stop_all('get_free_unit','Cannot find a free unit below 200')

    end function get_free_unit

end module utils
