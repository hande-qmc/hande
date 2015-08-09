module json_out

! Helper procedures for creating json-structured output.

! Quick JSON summary:

! {
!   key_1: val,
!   key_2: [x_1, x_2, x_3],
!   key_3: {
!     subkey_1: v_1,
!     subkey_2: v_2
!   }
! }
!
! val can be an integer, float, string, array ([....]), or another object (created by { ... }).
! All strings must be quoted.
! Trailing commas are not allowed in the last record of each object.

! Usage:

! 1. Set json_out_t%io if other than STDOUT is desired.
! 2. call json_object_init.
! 3. call json_write_key for each key/value pair (possibly containing nested objects).
! 4. For the final record in each object, be sure to set terminal=.true. in the json_write_key call.
! 5. Call json_object_end.  (Note: this must be done when closing each object.)

implicit none

private
public :: json_out_t, json_object_init, json_object_end, json_write_key

type json_out_t
    ! Unit to output to.
    integer :: io = 6
    ! Indent level to use.  Indent (automatically) one level for each nested object.
    integer, private :: level = 0
end type json_out_t

interface json_write_key
    module procedure :: write_int_32
    module procedure :: write_int_64
    module procedure :: write_real_32
    module procedure :: write_real_64
    module procedure :: write_string
    module procedure :: write_bool
    module procedure :: write_int_32_arr
    module procedure :: write_int_64_arr
    module procedure :: write_real_32_arr
    module procedure :: write_real_64_arr
end interface json_write_key

contains

    !--- Internal helper procedures ---

    pure function indent_level(level) result(xc)

        ! Return nX, where n=4*level+1, to use to indent the JSON output.

        ! In:
        !    level: indent level.

        character(3) :: xc
        integer, intent(in), optional :: level

        if (level <= 2) then
            write (xc,'(i1,"X")') 4*level+1
        else
            write (xc,'(i2,"X")') 4*level+1
        end if

    end function indent_level

    function record_delim(terminal, fmt_smt) result(rec_delim)

        ! Create string to append to the end of each record.

        ! In, optional:
        !     terminal: if present and true, return a blank
        !     fmt_smt: if present and true, create a string to be included at the end of
        !       the format statement, otherwise create the string to be printed out in
        !       the write statement.

        character(4) :: rec_delim
        logical, intent(in), optional :: terminal, fmt_smt

        rec_delim = ','
        if (present(fmt_smt)) then
            if (fmt_smt) rec_delim = ',","'
        end if
        if (present(terminal)) then
            if (terminal) rec_delim = ''
        end if

    end function record_delim

    subroutine write_key(js, key)

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        character(3) :: xc
        integer :: io

        xc = indent_level(js%level)
        write (js%io,'('//xc//'a)', advance='no') '"'//trim(key)//'": '

    end subroutine write_key

    !--- Object initialisation and termination ---

    subroutine json_object_init(js, key, tag)

        ! Start a new JSON object.

        ! In:
        !    js: json_out_t object controlling the output unit to write to and (automatically)
        !        the indentation level.
        !    key (optional): key of JSON object.  Should probably be set unless creating a new JSON
        !        document/starting a new object from scratch.
        !   tag (optional): if true (default: false) print out a tag before the opening bracket
        !        to allow for easy identication of JSON blocks in a mixed format output file.

        type(json_out_t), intent(inout) :: js
        character(*), intent(in), optional :: key
        logical, intent(in), optional :: tag
        integer :: io, l
        character(3) :: xc

        if (present(tag)) then
            if (tag) write (js%io,'(1X,"-- Start JSON block --")')
        end if
        if (present(key)) then
            call write_key(js, key)
        else
            xc = indent_level(js%level)
            write (js%io,'('//trim(xc)//')', advance='no')
        end if
        write (js%io,'("{")')
        js%level = js%level + 1

    end subroutine json_object_init

    subroutine json_object_end(js, terminal, tag)

        ! Close a JSON object.

        ! In:
        !    js: json_out_t object controlling the output unit to write to and (automatically)
        !        the indentation level.
        !    terminal (optional): if true (default: false) create the last record of the current
        !        object.  Note: this should be true when closing the last nested object in a JSON
        !        document *and* when closing the containing JSON document/outermost JSON object
        !        itself.
        !   tag (optional): if true (default: false) print out a tag after the closing bracket
        !        to allow for easy identication of JSON blocks in a mixed format output file.

        type(json_out_t), intent(inout) :: js
        logical, intent(in), optional :: terminal, tag
        integer :: io
        character(3) :: xc
        character(4) :: term

        js%level = js%level - 1
        xc = indent_level(js%level)
        term = record_delim(terminal, .true.)
        write (js%io,'('//xc//',"}"'//term//')')
        if (present(tag)) then
            if (tag) write (js%io,'(1X,"-- End JSON block --")')
        end if

    end subroutine json_object_end

    !--- Output of key/value pairs ---

    ! All procedures have the following interface:

    ! In:
    !    js: json_out_t object controlling the output unit to write to and (automatically)
    !        the indentation level.
    !    key: key to print out.
    !    value: value to print out.
    !    terminal: if present and true, create the last record of the current object.

    subroutine write_int_32(js, key, val, terminal)

        ! Write out a key/pair for an integer(int_32) value.

        use const, only: int_32

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        integer(int_32), intent(in) :: val
        logical, intent(in), optional :: terminal

        call write_key(js, key)
        write (js%io,'(i0,a)') val, record_delim(terminal)

    end subroutine write_int_32

    subroutine write_int_64(js, key, val, terminal)

        ! Write out a key/pair for an integer(int_64) value.

        use const, only: int_64

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        integer(int_64), intent(in) :: val
        logical, intent(in), optional :: terminal

        call write_key(js, key)
        write (js%io,'(i0,a)') val, record_delim(terminal)

    end subroutine write_int_64

    subroutine write_real_32(js, key, val, terminal)

        ! Write out a key/pair for an real(sp) value.

        use const, only: sp

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        real(sp), intent(in) :: val
        logical, intent(in), optional :: terminal

        call write_key(js, key)
        write (js%io,'(f0.8,a)') val, record_delim(terminal)

    end subroutine write_real_32

    subroutine write_real_64(js, key, val, terminal)

        ! Write out a key/pair for an real(dp) value.

        use const, only: dp

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        real(dp), intent(in) :: val
        logical, intent(in), optional :: terminal

        call write_key(js, key)
        write (js%io,'(f0.8,a)') val, record_delim(terminal)

    end subroutine write_real_64

    subroutine write_string(js, key, val, terminal)

        ! Write out a key/pair for a character string.

        use const, only: dp

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        character(*), intent(in) :: val
        logical, intent(in), optional :: terminal

        call write_key(js, key)
        write (js%io,*) '"'//trim(val)//'"', record_delim(terminal)

    end subroutine write_string

    subroutine write_bool(js, key, val, terminal)

        ! Write out a key/pair for a logical.

        use const, only: dp

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        logical, intent(in) :: val
        logical, intent(in), optional :: terminal

        call write_key(js, key)
        if (val) then
            write (js%io,*) 'true'//record_delim(terminal)
        else
            write (js%io,*) 'false'//record_delim(terminal)
        end if

    end subroutine write_bool

    subroutine write_int_32_arr(js, key, val, terminal)

        ! Write out a key/pair for a integer(int_32) 1D array.

        use const, only: int_32
        use utils, only: int_fmt

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        integer(int_32), intent(in) :: val(:)
        logical, intent(in), optional :: terminal
        integer :: i

        call write_key(js, key)
        write (js%io, '("[")', advance='no')
        do i = 1, size(val)
            write (js%io,'('//int_fmt(val(i),1)//')', advance='no') val(i)
            if (i/=size(val)) write (js%io,'(",")', advance='no')
        end do
        write (js%io,'("]"'//record_delim(terminal, .true.)//')')

    end subroutine write_int_32_arr

    subroutine write_int_64_arr(js, key, val, terminal)

        ! Write out a key/pair for a integer(int_64) 1D array.

        use const, only: int_64
        use utils, only: int_fmt

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        integer(int_64), intent(in) :: val(:)
        logical, intent(in), optional :: terminal
        integer :: i

        call write_key(js, key)
        write (js%io, '("[")', advance='no')
        do i = 1, size(val)
            write (js%io,'('//int_fmt(val(i),1)//')', advance='no') val(i)
            if (i/=size(val)) write (js%io,'(",")', advance='no')
        end do
        write (js%io,'("]"'//record_delim(terminal, .true.)//')')

    end subroutine write_int_64_arr

    subroutine write_real_32_arr(js, key, val, terminal)

        ! Write out a key/pair for a real(sp) 1D array.

        use const, only: sp

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        real(sp), intent(in) :: val(:)
        logical, intent(in), optional :: terminal
        integer :: i

        call write_key(js, key)
        write (js%io, '("[")', advance='no')
        do i = 1, size(val)
            write (js%io,'(f0.8)', advance='no') val(i)
            if (i/=size(val)) write (js%io,'(",")', advance='no')
        end do
        write (js%io,'("]"'//record_delim(terminal, .true.)//')')

    end subroutine write_real_32_arr

    subroutine write_real_64_arr(js, key, val, terminal)

        ! Write out a key/pair for a real(dp) 1D array.

        use const, only: dp

        type(json_out_t), intent(in) :: js
        character(*), intent(in) :: key
        real(dp), intent(in) :: val(:)
        logical, intent(in), optional :: terminal
        integer :: i

        call write_key(js, key)
        write (js%io, '("[")', advance='no')
        do i = 1, size(val)
            write (js%io,'(f0.8)', advance='no') val(i)
            if (i/=size(val)) write (js%io,'(",")', advance='no')
        end do
        write (js%io,'("]"'//record_delim(terminal, .true.)//')')

    end subroutine write_real_64_arr

end module json_out
