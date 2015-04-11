module utils

! Various utilities and tools...

implicit none

interface int_fmt
    module procedure int_fmt_int_32
    module procedure int_fmt_int_64
end interface int_fmt

contains

! --- Combinatorics ---

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

    elemental function factorial(n) result(fac)

        ! In:
        !   n: integer
        ! Returns:
        !   factorial of n, n!

        ! WARNING:
        ! This does *not* safeguard against integer overflow.  As such, it is
        ! only suitably for 0 <= n <= 12.  Investigating using log(n!) or the
        ! Gamma function if required for larger values.

        integer :: fac
        integer, intent(in) :: n

        integer :: i

        fac = 1
        do i = 2, n
            fac = fac*i
        end do

    end function factorial

    elemental function factorial_combination_1(m,n) result (combination)

        ! Given m and n input, this function returns
        ! combination = n!m!/(n+m+1)!
        ! Required to calculate amplitudes of the Neel singlet
        ! trial function for the Heisenberg model

        use const, only: dp

        real(dp) :: combination
        integer, intent(in) :: m, n
        integer :: m1, n1, i

        ! Choose m1 to be the larger of m and n
        if (m >= n) then
            m1 = m
            n1 = n
        else
            m1 = n
            n1 = m
        end if
        combination = 1
        do i = 1, n1
            combination = combination * i/(i+m1)
        end do
        combination = combination/(m1+n1+1)

    end function factorial_combination_1

    pure subroutine next_comb(n, k, comb, ierr)

        ! Create the next combination of k objects from a set of size n.

        ! In:
        !    n: size of set.
        !    k: size of subset/combination.
        ! In/Out:
        !    comb: contains the previous combination on input and the next
        !        combination on output.  Use (0,1,2,...,k-1) for the first
        !        combination.  The returned combination is 0-indexed.
        ! Out:
        !    ierr: 0 if a combination is found and 1 if there are no more
        !        combinations.

        ! Translated from the C implementation at:
        ! https://compprog.wordpress.com/2007/10/17/generating-combinations-1/.
        ! The author, Alexandru Scvortov (scvalex@gmail.com), gave permission
        ! via email that we are 'free to use/modify/redistribute that code
        ! however [we] like'.  His code is provided as-is.

        integer, intent(in) :: n, k
        integer, intent(inout) :: comb(0:k-1) ! 0-indexed for easy translation from original C.
        integer, intent(out) :: ierr

        integer :: i, j

        ! The final element of the previous combination must, at very least, be
        ! incremented.
        i = k - 1
        comb(i) = comb(i) + 1

        do while (comb(i) >= n - k + 1 + i)
            ! The i-th element is greater than the maximum allowed value (note
            ! that combinations are generated as ascending sequences).  Hence
            ! there are no more combinations which begin comb(0:i-1).  Increment
            ! comb(i-1) and use that as the start of the next combination.
            i = i - 1
            if (i < 0) exit
            comb(i) = comb(i) + 1
        end do

        if (comb(0) > n - k) then
            ! combination (n-k, n-k+1, ..., n) reached.
            ! no more combinations can be geerated.
            ierr = 1
        else
            ierr = 0
            ! comb now looks like (..., x, n_1, n_2, n_3, ..., n_n),
            ! where n_i ! >= max value of comb.
            ! Turn it into (..., x, x+1, x+1, ...), ie the first combination
            ! which starts with ( ..., x).  Subsequent calls to next_comb will
            ! then iterate over 'higher' combinations with the same starting
            ! sequence.
            ! The i-th position at the end of the previous loop holds the
            ! x entry.
            do j = i+1, k-1
                comb(j) = comb(j-1) + 1
            end do
        end if

    end subroutine next_comb

!--- format statement formatting ---

    elemental function int_fmt_int_32(i, padding) result(fmt1)

        ! In:
        !    i: an integer
        !    padding (optional): amount of padding to add to format statement.
        !        Default: 2.
        ! Returns:
        !    fmt1: a format statement for an integer field which will hold
        !        i perfectly plus an amount of padding.

        ! This does take i/o formatting to a slightly OCD level addmittedly...

        use const, only: dp

        character(4) :: fmt1
        integer, intent(in) :: i
        integer, intent(in), optional :: padding
        real(dp) :: logi

        if (i == 0 .or. i==1) then
            logi = 1.0
        else
            logi = log10(real(abs(i)+1,dp))
        end if
        if (i < 0) logi = logi + 1

        fmt1 = int_fmt_helper(logi, padding)

    end function int_fmt_int_32

    elemental function int_fmt_int_64(i, padding) result(fmt1)

        ! In:
        !    i: a long integer
        !    padding (optional): amount of padding to add to format statement.
        !        Default: 2.
        ! Returns:
        !    fmt1: a format statement for an integer field which will hold
        !        i perfectly plus an amount of padding.

        ! This does take i/o formatting to a slightly OCD level addmittedly...

        use const, only: int_64, dp

        character(4) :: fmt1
        integer(int_64), intent(in) :: i
        integer, intent(in), optional :: padding
        real(dp) :: logi

        if (i == 0 .or. i==1) then
            logi = 1.0
        else
            logi = log10(real(abs(i)+1,dp))
        end if
        if (i < 0) logi = logi + 1

        fmt1 = int_fmt_helper(logi, padding)

    end function int_fmt_int_64

    elemental function int_fmt_helper(logi, padding) result(fmt1)

        ! In:
        !    logi: log10 of an integer.
        !    padding (optional): amount of padding to add to format statement.
        !        The default amount is 2.  The padding is used to include the
        !        sign if i is negative.
        ! Returns:
        !    fmt1: a format statement for an real field using the G format
        !       statement which will hold i perfectly plus an amount of padding.

        use const, only: dp

        character(4) :: fmt1
        real(dp), intent(in) :: logi
        integer, intent(in), optional :: padding
        integer :: p

        if (present(padding)) then
            p = padding
        else
            p = 2
        end if

        p = ceiling(logi+p)

        if (p < 10) then
            write (fmt1,'("i",i1)') p
        else if (p < 100) then
            write (fmt1,'("i",i2)') p
        else
            ! By this point we'll have hit integer overflow (for 32- and 64-bit
            ! integers) anyway...
            write (fmt1,'("i",i3)') p
        end if

    end function int_fmt_helper

! --- File names and file handling ---

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

    elemental subroutine append_ext(stem, n, s)

        ! Returns stem.n in s.

        character(*), intent(in) :: stem
        integer, intent(in) :: n
        character(*), intent(out) :: s
        character(10) :: ext

        write (ext,'('//int_fmt(n,0)//')') n
        s = stem//'.'//ext

    end subroutine append_ext

   subroutine get_unique_filename(stem, suffix, tnext, istart, filename)

        ! Find a filename which is either the "newest" or the next to be used.
        ! The filename is assumed to be stem.xsuffix, where x is an integer.

        ! In:
        !    stem: stem of the filename.
        !    suffix: suffix of filename
        !    tnext: the next unused filename is found if true, else the
        !        filename is set to be stem.x where stem.x exists and stem.x+1
        !        doesn't and x is greater than istart.
        !    istart: the integer of the first x value to check.
        !        If istart is negative, then the filename is set to be stem.x,
        !        where x = |istart+1|.  This overrides everything else.
        ! Out:
        !    filename.

        character(*), intent(in) :: stem, suffix
        logical, intent(in) :: tnext
        integer, intent(in) :: istart
        character(*), intent(out) :: filename

        integer :: i
        logical :: exists

        if (istart >= 0) then

            i = istart
            exists = .true.
            do while (exists)
                call append_ext(stem, i, filename)
                filename = trim(filename)//suffix
                inquire(file=filename,exist=exists)
                i = i + 1
            end do

            if (.not.tnext) then
                ! actually want the last file which existed.
                ! this will return stem.istart if stem.istart doesn't exist.
                i = max(istart,i - 2)
                call append_ext(stem, i, filename)
                filename = trim(filename)//suffix
            end if

        else
            ! Have been asked for a specific file.
            call append_ext(stem, abs(istart+1), filename)
            filename = trim(filename)//suffix
        end if

    end subroutine get_unique_filename

    subroutine read_file_to_buffer(buffer, fname, in_unit)

        ! Read the entire contents of a file into a character buffer.

        ! In (optional):
        !    fname: name of file to be read.  If supplied, the file is opened, read
        !        and then closed.
        !    in_unit: unit of already opened file.  Reads the rest of the file and
        !        returns with the file at the *final* record.  Overrides fname.
        ! Out:
        !    buffer: set to the contents in the file.  Enlarged to hold the entire
        !        contents of the file.

        ! NOTE: for GCC 4.7 and below, we take a fixed buffer and assume that
        ! it is large enough to hold the entire file and lines are (at most)
        ! 1024 characters long.  For all other compilers we take an allocatalbe
        ! character string and dynamically resize it.

        use, intrinsic :: iso_fortran_env

        use errors, only: stop_all

        character(*), optional, intent(in) :: fname
        integer, optional, intent(in) :: in_unit
        ! GCC 4.7 and below don't have sufficient F2003 support for
        ! allocatable strings.
#if __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ > 7))
        character(:), allocatable, intent(out) :: buffer

        character(:), allocatable :: line
#else
        character(*), intent(out) :: buffer

        character(1024) :: line
#endif
        character(256) :: chunk
        integer :: chunk_size, stat, iunit

        if (present(in_unit)) then
            iunit = in_unit
        else if (present(fname)) then
            iunit = get_free_unit()
            open(iunit, file=fname, status='old', form='formatted')
        else
            call stop_all('read_file_to_buffer', 'Neither file nor file handle supplied to read_file_to_buffer.')
        end if

        ! NOTE: we use Fortran 2003's left-hand-side reallocation feature for both line and buffer...
        stat = 0
#if (__GNUC__ == 4 && (__GNUC_MINOR__ < 8))
        buffer = ''
#endif
        do
            line = ''
            do
                read(iunit, '(A)', advance='no', iostat=stat, size=chunk_size) chunk
                if (stat /= 0 .and. .not. (stat == iostat_eor .or. stat == iostat_end)) then
                    if (present(fname)) call stop_all('read_file_to_buffer', 'Problem reading file: '//trim(fname))
                    if (present(in_unit)) call stop_all('read_file_to_buffer', 'Problem reading file from given unit.')
                end if
                line = trim(line) // chunk(:chunk_size)
                if (stat == iostat_eor .or. stat == iostat_end) exit
            end do
            if (stat == iostat_end) exit
#if __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ > 7))
            if (allocated(buffer)) then
                buffer = buffer // new_line(line) // line
            else
                buffer = line
            end if
#else
            if (trim(buffer) /= '') then
                buffer = trim(buffer) // new_line(line) // line
            else
                buffer = line
            end if
#endif
        end do

        if (present(fname)) close(iunit, status='keep')

    end subroutine read_file_to_buffer

! --- Array indexing ---

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

    elemental function tri_ind_reorder(i,j) result(indx)

        ! Find the index corresponding to the (i,j)-th element of a lower
        ! triangular array.
        !
        ! We assume that i >= j.  If this is not the case (i.e. (i,j) refers to
        ! an element in the upper triangular array) then the index of the
        ! transpose element (i.e. (j,i)) is returned.
        !
        ! This maps:
        !
        !   1,1 1,2 1,3 1,4        1  2  4  7
        !   2,1 2,2 2,3 2,4        2  3  5  8
        !   3,1 3,2 3,3 3,4   to   4  5  6  9
        !   4,1 4,2 4,3 4,4        7  8  9 10
        !
        ! In:
        !    i: (1-indexed) index
        !    j: (1-indexed) index
        ! Returns:
        !    A combined (1-indexed) index for the corresponding element in
        !    a lower triangular array.

        integer :: indx
        integer, intent(in) :: i, j

        if (i>=j) then
            indx = tri_ind(i,j)
        else
            indx = tri_ind(j,i)
        end if

    end function tri_ind_reorder

! --- String conversion ---

    pure function fstring_to_carray(string_f03) result(array_c)

        ! Convert a Fortran string into a C string.  The result can be passed to
        ! C procedures which require a *char object, which corresponds to an
        ! array in Fortran.

        ! In:
        !    string_f03: Fortran character string.
        ! Returns:
        !    The equivalent (null-terminated) C string in a character array.

        use, intrinsic :: iso_c_binding, only: c_char, c_null_char

        character(len=*), intent(in) :: string_f03
        character(kind=c_char) :: array_c(len(string_f03)+1)

        integer :: i

        do i = 1, len(string_f03)
            array_c(i) = string_f03(i:i)
        end do
        array_c(i) = c_null_char

    end function fstring_to_carray

    pure function carray_to_fstring(array_c) result(string_f03)

        ! Convert a C string into a Fortran string.  The input can be from a
        ! C procedure which uses a *char object (which corresponds to an array
        ! in Fortran) and the output is the equivalent standard Fortran string.

        ! In:
        !    array_c: null-terminated C string in a character array.
        ! Returns:
        !    The equivalent Fortran character string (without null termination).

        use, intrinsic :: iso_c_binding, only: c_char, c_null_char

        character(kind=c_char), intent(in) :: array_c(:)
        character(len=size(array_c)-1) :: string_f03

        integer :: i

        ! Don't copy any (presumably garbage) from beyond the null character.
        string_f03 = ''
        do i = 1, size(array_c)-1
            if (array_c(i) == c_null_char) exit
            string_f03(i:i) = array_c(i)
        end do

    end function carray_to_fstring

!--- Informative output ---

    subroutine rng_init_info(seed)

        integer, intent(in) :: seed

        write (6,'(1X,a3,/,1X,3("^"),/)') 'RNG'
        write (6,'(1X,a51,'//int_fmt(seed,1)//',a1,/)') 'Initialised random number generator with a seed of:', seed, '.'

    end subroutine rng_init_info

    subroutine print_matrix(matrix)

        ! Print out a given real matrix in a neat format.

        ! In:
        !    matrix: The matrix which is to be output to the screen. 

        use const, only: p

        real(p), intent(in) :: matrix(:,:)
        integer :: i,j

        do i = 1, ubound(matrix,1)
            do j = 1, ubound(matrix,2)
                write(6, '(es17.10,2X)', advance = 'no') matrix(i,j)
            end do
            write(6,'(1X)', advance = 'yes')
        end do

    end subroutine print_matrix

end module utils
