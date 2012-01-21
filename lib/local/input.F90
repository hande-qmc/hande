MODULE input

!  Fortran90 input parsing module
!
!  Copyright (C) 2005 Anthony J. Stone
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to
!  the Free Software Foundation, Inc., 51 Franklin Street,
!  Fifth Floor, Boston, MA 02110-1301, USA, or see
!  http://www.gnu.org/copyleft/gpl.html


#ifdef NAGF95
use f90_unix_env, ONLY: getarg
#endif

IMPLICIT NONE

CHARACTER(LEN=800), SAVE :: char=''
LOGICAL, SAVE :: skipbl=.false., clear=.true., echo=.false.,           &
    debug=.false., more
INTEGER, SAVE :: item=0, nitems=0, loc(0:80)=0, end(80)=0,               &
    line(0:10)=0, level=0, nerror=0, ir=5, last=0, unit(0:10)

CHARACTER(LEN=26), PARAMETER ::                                        &
    upper_case="ABCDEFGHIJKLMNOPQRSTUVWXYZ",                           &
    lower_case="abcdefghijklmnopqrstuvwxyz"
CHARACTER, PARAMETER :: space = " ", bra = "(", ket = ")",             &
    comma = ",", squote = "'", dquote = '"', tab=achar(9),             &
    plus="+", minus="-", dot="."

CHARACTER(LEN=8), SAVE :: concat = "+++"
CHARACTER(LEN=40) :: file(10)=""

INTEGER, SAVE :: lc=3

INTEGER, PARAMETER :: sp=kind(1.0),dp=kind(1.d0)!, qp=selected_real_kind(30)
INTEGER, PARAMETER :: l=selected_int_kind(15)

INTERFACE readf
  MODULE PROCEDURE read_single, read_double!, read_quad
END INTERFACE

PRIVATE
PUBLIC :: item, nitems, read_line, stream, reada, readu, readl,        &
    readf, readi, readli, getf, geta, geti, getli, reread, input_options,             &
    upcase, locase, report, die, assert, find_io, read_colour,         &
    getargs, parse, char, ir
!  AJWT - added , ir to above
!  Free-format input routines

!     CALL READ_LINE(eof[,inunit])
!  Read next input record from unit IR into buffer, and analyse for data
!  items, null items (enclosed by commas), and comment (enclosed in
!  parentheses, which may be nested, and occurring where a new item may
!  occur). Data items are terminated by space or comma, unless enclosed
!  in single or double quotes.
!  If the optional argument inunit is provided, the record is read from
!  there instead.
!  If a line ends with the concatenation sequence (default "+++")
!  possibly followed by spaces, the sequence and any following spaces
!  are removed and the next line concatenated in its place.
!  This is repeated if necessary until a line is read that does not end
!  with the concatenation sequence.
!  The logical variable eof is set TRUE if end of file is encountered,
!  otherwise FALSE.
!  The public module variable ITEM is set to zero and NITEMS to the number
!  of items in the record.

!     CALL PARSE(string)
!  Parse the string provided in the same way as the read_line routine
!  (which uses this routine itself) and leave the details in the buffer
!  as for a line read from the input. Input directives are not
!  interpreted.

!     CALL READx(V)
!  Read an item of type x from the buffer into variable V:
!     CALL READF   single or double precision, depending on the type of V
!     CALL READI   integer
!     CALL READLI  long integer
!     CALL READA   character string
!     CALL READU   character string, uppercased
!     CALL READL   character string, lowercased
!  All of these return null or zero if there are no more items on the
!  current line. Otherwise ITEM is incremented so that it gives the
!  number of the item that has just been read.
!     CALL READF(V,factor)
!  If the optional second argument is supplied to READF, the variable V
!  is divided by it after it has been read. This is convenient for converting
!  from external to internal units.
!     CALL REREAD(K)
!  K>0  Prepare to read item K
!  K<0  Go back |K| items
!  K=0  Same as K=-1, i.e. reread last item.

!     CALL GETx
!  Same as the corresponding READx, but a new record is read if there are
!  no more items in the current one.

!     CALL READ_COLOUR(fmt,col,clamp)
!  Read a colour definition, in a form specified by FMT:
!  FMT="GREY": read a single number between 0 and 1
!  FMT="RGB": read 3 numbers, which are RGB colour values between 0 and 1
!  FMT="RGB255": read 3 numbers, which are RGB colour values between 0 and 255
!  FMT="RGBX": read a single 6-character string giving the RGB colour
!        values in hexadecimal.
!  In the last two cases the colour values are scaled to the range from 0 to 1.
!  CLAMP is an optional logical argument. If present and true, the colour
!        values are clamped (after scaling, if appropriate) to the range 0 to 1.

!  ITEM    is the number of the last item read from the buffer

!  NITEMS  is the number of items in the buffer

!  char(I:I) is the Ith character in the buffer

!  LOC(I)  is the starting position of the Ith item
!  END(I)  is the position of the last character of the Ith item

!  LINE    is the number of lines (records) read

!  If SKIPBL is set to TRUE, lines containing no items other than
!  comment are skipped automatically, so that INPUT will always return
!  a line with at least one item (possibly null).  If SKIPBL is FALSE
!  (default) no skipping occurs and the next data line is returned
!  regardless of what it contains.

!  If CLEAR is set to TRUE (default) then null items will be returned
!  as zero or blank. If an attempt is made to read more than
!  NITEMS items from a line, the items are treated as null. If
!  CLEAR is FALSE, a variable into which a null item is read is
!  left unaltered.

!  NERROR specifies the treatment of errors found while reading numbers:
!    0  Hard error - print message and stop (default).
!    1  Soft error - print message and return zero.
!    2  Soft error - no message, return zero, set NERROR to -1.
!  If NERROR is set to 2 it is important to test and reset it after reading
!  a number.

!  IR is the input stream from which the data are to be read (default 5).

!  If ECHO is TRUE, the input line will be reflected to standard output.

!  LAST gives the position of the last non-blank character in the line.



!     CALL REPORT(message[,REFLECT])
!  Error-message routine. Prints the error message, followed by the
!  current input buffer if REFLECT is present and TRUE, and stops.

CONTAINS

SUBROUTINE read_line(eof,inunit)

!  Read next input record from unit IR into buffer, and analyse for data
!  items, null items (enclosed by commas), and comment (enclosed in
!  parentheses, which may be nested, and occurring where spaces may
!  occur).
!  If the optional argument inunit is specified, read a line from that
!  unit instead of IR.

!  Stream-switching commands may occur in the data:
!      #include file-name
!         switches input to be from the file specified;
!      #revert
!         (or end-of-file) reverts to the previous file.
!  However it does not make sense for stream-switching commands to
!  be used when the input unit is specified explicitly. If they are given
!  in this case, they apply to the default input stream.
!    Also
!      #concat "string"
!         sets the concatenation string; e.g.
!         #concat "\"
!      #width <n>
!         Input line width limited to n characters. Default is 128; may
!         need to be reduced to 80 for some IBM computers. Maximum is 255.


!  CONCAT is the line concatenation string: if it is found as the last
!  non-blank character string on a line, the following line is read and
!  appended to the current line, overwriting the concatenation string.
!  This procedure is repeated if necessary until a line is found that does
!  not end with the concatenation string.

!  The concatenation string may be modified by changing its declaration
!  above. If it is null, no concatentation occurs. The concatenation string
!  may also be changed by the #concat directive in the data file.

LOGICAL, INTENT(OUT) :: eof
INTEGER, INTENT(IN), OPTIONAL :: inunit

CHARACTER(LEN=255) :: w, f
CHARACTER :: term

INTEGER, SAVE :: lrecl = 128
INTEGER :: in, fail, i, k, l, m

eof=.false.
if (present(inunit)) then
  in=inunit
else
  in=ir
end if

char=""
lines: do
  more=.true.
  m=1
  do while (more)
    last=m+lrecl-1
    line(level)=line(level)+1
    read (in,"(a)",end=900) char(m:last)
    go to 10
!  End of file
900 if (more .and. m .gt. 1) then
      print "(a)", "Apparently concatenating at end-of-file"
      call report("Unexpected end of data file",.true.)
    endif
    if (level .gt. 0) then
      !  Revert to previous input
      close(in)
      level=level-1
      ir=unit(level)
      in=ir
      cycle lines
    else
      !  End of input
      char(1:last)=""
      item=0
      nitems=-1
      eof=.true.
      return
    endif

!  Find last non-blank character
10  last=verify(char,space//tab,back=.true.)
    if (echo) print "(a)", char(m:last)
!  Look for concatenation string
    if (lc .gt. 0 .and. last .ge. lc) then
      more=(char(last-lc+1:last) .eq. concat)
      if (more) then
        m=last-lc+1
      endif
    else
      more=.false.
    endif
  end do  ! while (more)

!  Replace tab by single space
  do while (index(char,tab) .gt. 0)
    L=index(char,tab)
    char(L:L)=space
  end do

!  Logical line assembled. First look for input directives
  L=1
  do while (char(L:L) .eq. space .and. L .lt. last)
    L=L+1
  end do
  if (char(L:L) .eq. "#") then
    !       M=index(char(L:),space)+L-1
    M=L
    do while(char(M:M) .ne. space .and. M .le. last)
      M=M+1
    end do
    w=char(L:M-1)
    call upcase(w)
    if (M .gt. last) then
      f=" "
    else
      do while (char(M:M) .eq. space)
        M=M+1
      end do
      if (char(M:M) .eq. squote .or. char(M:M) .eq. dquote) then
        term=char(M:M)
        M=M+1
      else
        term=space
      endif
      !         L=index(char(M:),term)+M-1
      L=M
      do while(char(M:M) .ne. term .and. M .le. last)
        M=M+1
      end do
      f=char(L:M-1)
    endif
    select case(w)
    case("#")          ! Comment -- ignore
    case("#INCLUDE")   ! Take input from specified file
      if (f .eq. " ") call report                                      &
          ("No filename given in #include directive",.true.)
      if (level .eq. 0) unit(0)=ir
      level=level+1
      line(level)=0
      ir=find_io(91)
      unit(level)=ir
      open(unit=ir,file=f,status="old",iostat=fail)
      if (fail .ne. 0) then
        call report(trim(f)//" could not be opened",.true.)
      end if
      in=ir
      file(level)=f(:40)
    case("#CONCAT")
      concat=f(:8)
      lc=len(trim(concat))
    case("#REVERT")
      close(in)
      file(level)=""
      level=level-1
      ir=unit(level)
      in=ir
    case("#WIDTH")
      read (unit=f,fmt="(i6)") lrecl
    case("#ECHO")
      echo=.true.
      if (f .eq. "OFF") echo=.false.
    case("#DEBUG")
      debug=.true.
      if (f .eq. "OFF") debug=.false.
    case default
      call report("Unrecognized directive "//trim(w)//"in input",.true.)
    end select
    cycle lines
  endif

  call parse

!  Blank except for comment?
  if (nitems .eq. 0 .and. skipbl) then
    cycle lines   !  Read another line
  else
    exit lines    !  Finished
  end if

end do lines

if (debug) then
  !       print "(8(I6,I4))", (loc(i), end(i), i=1,nitems)
  if (echo .and. nitems .gt. 0) then
    print "(100a1)", (" ", i=1,loc(1)-1),                              &
        (("+", i=loc(k), end(k)), (" ", i=end(k)+1, loc(k+1)-1), k=1,nitems-1), &
        ("+", i=loc(nitems), end(nitems))
  endif
  print "(I2,A)", nitems, " items"
endif

END SUBROUTINE read_line

!-----------------------------------------------------------------------

SUBROUTINE parse(string)

CHARACTER(LEN=*), OPTIONAL :: string

INTEGER :: L, state, nest
LOGICAL :: tcomma
CHARACTER :: term, c

if (present(string)) then
  char=string
  last=len_trim(string)
end if

!  Analyse input
!  State numbers:
!  0  Looking for item
!  1  Reading quoted string
!  2  Reading unquoted item
!  3  Reading comment
!  4  Expecting space or comma after quoted string

state=0
item=0
nitems=0
L=0            ! Position in input buffer
tcomma=.true.  !  True if last item was terminated by comma
               !  and also at start of buffer


chars: do

!  Read next character
  L=L+1
  if (L .gt. last) then
!  End of line
    if (nitems .gt. 0) then
      select case(state)
      case(1)
        call report("Closing quote missing",.true.)
      case(2)
        end(nitems)=L-1
      end select
    endif
    return
  endif

  c=char(L:L)
!   if (debug) print "(A,I3,A,A,A,I1)",                                &
!       "L = ", L, "  Character ", c, "  state ", state
  select case (state)
  case(0)                ! Looking for next item
    select case(c)
    case(space,tab)      ! Keep on looking
      !           cycle chars
    case(bra)            ! Start of comment
      nest=1
      state=3
      !           cycle chars
    case(squote,dquote)  ! Start of quoted string
      nitems=nitems+1
      loc(nitems)=L
      term=c
      state=1
    case(comma)
      if (tcomma) then   ! Null item between commas
        nitems=nitems+1
        loc(nitems)=0
      endif
      tcomma=.true.
    case default         ! Start of unquoted item
      nitems=nitems+1
      loc(nitems)=L
      state=2
    end select

  case(1)                ! Reading through quoted string
    if (c .eq. term) then ! Closing quote found
      end(nitems)=L
      state=4
      !         else
      !           cycle chars
    endif

  case(2)                ! Reading through unquoted item
    select case(c)
    case(space,tab)      ! Terminator
      end(nitems)=L-1
      state=0
      tcomma=.false.
!   case(bra)            ! Start of comment -- treat as space
!     end(nitems)=L-1    ! This code allows for parenthesised comments
!     tcomma=.false.     ! to be embedded in unquoted strings. This is
!     state=3            ! not a good idea. Such comments now have to occur
!     nest=1             ! only where a new item might begin.
    case(comma)          ! Comma-terminated
      end(nitems)=L-1
      state=0
      tcomma=.true.
      !         case default
      !           cycle chars
    end select

  case(3)                ! Reading through comment
    select case(c)
    case(bra)            ! Nested parenthesis
      nest=nest+1
    case(ket)
      nest=nest-1
      if (nest .eq. 0) then  ! End of comment
        state=0          ! Space or comma not required after comment
      endif
      !         case default
      !           cycle chars
    end select

  case(4)                ! Expecting space or comma
    select case(c)
    case(space,tab)
      tcomma=.false.
      state=0
    case(comma)
      tcomma=.true.
      state=0
    case(bra)            ! Start of comment -- treat as space
      tcomma=.false.
      state=3
      nest=1
    case default
      call report("Space or comma needed after quoted string",.true.)
    end select

  end select

end do chars

END SUBROUTINE parse

!-----------------------------------------------------------------------

SUBROUTINE getargs

CHARACTER(LEN=80) :: word

nitems=0
last=-1

do
  nitems=nitems+1
  call getarg(nitems,word)
  if (word .eq. "") then
    nitems=nitems-1
    exit
  else
    loc(nitems)=last+2
    char(last+2:)=word
    last=len_trim(char)
    end(nitems)=last
  endif
end do

END SUBROUTINE getargs

!-----------------------------------------------------------------------

SUBROUTINE input_options(default,clear_if_null,skip_blank_lines,       &
    echo_lines, error_flag, concat_string)

IMPLICIT NONE
LOGICAL, OPTIONAL :: default, clear_if_null, skip_blank_lines,         &
    echo_lines
INTEGER, OPTIONAL :: error_flag
CHARACTER(LEN=*), optional :: concat_string

if (present(default)) then
  if (default) then
    clear=.true.
    skipbl=.false.
    echo=.false.
    nerror=0
    concat="+++"
  endif
endif
if (present(clear_if_null)) then
  clear=clear_if_null
endif
if (present(skip_blank_lines)) then
  skipbl=skip_blank_lines
endif
if (present(echo_lines)) then
  echo=echo_lines
endif
if (present(error_flag)) then
  nerror=error_flag
endif
if (present(concat_string)) then
  if (len(trim(concat_string)) .gt. 8) call report                     &
      ("Concatenation string must be 8 characters or fewer",.false.)
  concat=concat_string
  lc=len(trim(concat_string))
endif

END SUBROUTINE input_options

!-----------------------------------------------------------------------

SUBROUTINE stream(n)

INTEGER, INTENT(IN) :: n
!  Set the input stream for subsequent data to be N.

ir=n

END SUBROUTINE stream

!-----------------------------------------------------------------------

SUBROUTINE reada(m)

!  Copy characters from the next item into the character variable M.
!  If the first character is a single or double quote, the string is
!  terminated by a matching quote and the quotes are removed.

CHARACTER(LEN=*), INTENT(INOUT) :: m
INTEGER :: l

if (clear) m=""
!  If there are no more items on the line, M is unchanged
if (item .ge. nitems) return

item=item+1
!  Null item?
if (loc(item) .eq. 0) return

l=loc(item)
if (char(l:l) .eq. squote .or. char(l:l) .eq. dquote) then
  m=char(l+1:end(item)-1)
else
  m=char(l:end(item))
endif

END SUBROUTINE reada

!-----------------------------------------------------------------------

! SUBROUTINE read_quad(A,factor)
!
! !  Read the next item from the buffer as a real (quadruple precision) number.
! !  If the optional argument factor is present, the value read should be
! !  divided by it. (External value = factor*internal value)
!
! REAL(KIND=qp), INTENT(INOUT) :: a
! REAL(KIND=qp), INTENT(IN), OPTIONAL :: factor
!
! CHARACTER(LEN=50) :: string
!
! if (clear) a=0.0_qp
!
! !  If there are no more items on the line, I is unchanged
! if (item .ge. nitems) return
!
! string=""
! call reada(string)
! !  If the item is null, I is unchanged
! if (string == "") return
! read (unit=string,fmt=*,err=99) a
! if (present(factor)) then
!   a=a/factor
! endif
! return
!
! 99 a=0.0_qp
! select case(nerror)
! case(-1,0)
!   call report("Error while reading real number",.true.)
! case(1)
!   print "(2a)", "Error while reading real number. Input is ", trim(string)
! case(2)
!   nerror=-1
! end select
!
! END SUBROUTINE read_quad

!-----------------------------------------------------------------------

SUBROUTINE read_double(A,factor)

!  Read the next item from the buffer as a real (double precision) number.
!  If the optional argument factor is present, the value read should be
!  divided by it. (External value = factor*internal value)

DOUBLE PRECISION, INTENT(INOUT) :: a
DOUBLE PRECISION, INTENT(IN), OPTIONAL :: factor

CHARACTER(LEN=50) :: string

if (clear) a=0d0

!  If there are no more items on the line, I is unchanged
if (item .ge. nitems) return

string=""
call reada(string)
!  If the item is null, I is unchanged
if (string == "") return
read (unit=string,fmt=*,err=99) a
if (present(factor)) then
  a=a/factor
endif
return

99 a=0d0
select case(nerror)
case(-1,0)
  call report("Error while reading real number",.true.)
case(1)
  print "(2a)", "Error while reading real number. Input is ", trim(string)
case(2)
  nerror=-1
end select

END SUBROUTINE read_double

!-----------------------------------------------------------------------

SUBROUTINE read_single(A,factor)

!  Read the next item from the buffer as a real (double precision) number.
!  If the optional argument factor is present, the value read should be
!  divided by it. (External value = factor*internal value)

REAL(kind=sp), INTENT(INOUT) :: a
REAL(kind=sp), INTENT(IN), OPTIONAL :: factor

DOUBLE PRECISION :: aa

if (present(factor)) then
  call read_double(aa,real(factor,dp))
else
  call read_double(aa)
endif
a=real(aa,sp)

END SUBROUTINE read_single

!-----------------------------------------------------------------------

SUBROUTINE readi(I)
!  Read an integer from the current record

INTEGER, INTENT(INOUT) :: i

CHARACTER(LEN=50) :: string

if (clear) i=0

!  If there are no more items on the line, I is unchanged
if (item .ge. nitems) return

string=""
call reada(string)
!  If the item is null, I is unchanged
if (string == "") return
read (unit=string,fmt=*,err=99) i
return

99 i=0
select case(nerror)
case(-1,0)
  call report("Error while reading integer",.true.)
case(1)
  print "(2a)", "Error while reading integer. Input is ", trim(string)
case(2)
  nerror=-1
end select

END SUBROUTINE readi

!-----------------------------------------------------------------------

SUBROUTINE readli(I)
!  Read an integer from the current record

INTEGER(l), INTENT(INOUT) :: i

CHARACTER(LEN=50) :: string

if (clear) i=0

!  If there are no more items on the line, I is unchanged
if (item .ge. nitems) return

string=""
call reada(string)
!  If the item is null, I is unchanged
if (string == "") return
read (unit=string,fmt=*,err=99) i
return

99 i=0
select case(nerror)
case(-1,0)
  call report("Error while reading integer",.true.)
case(1)
  print "(2a)", "Error while reading integer. Input is ", trim(string)
case(2)
  nerror=-1
end select

END SUBROUTINE readli

!-----------------------------------------------------------------------

SUBROUTINE readu(m)
CHARACTER(LEN=*) m

call reada(m)
call upcase(m)

END SUBROUTINE readu

!-----------------------------------------------------------------------

SUBROUTINE readl(m)
CHARACTER(LEN=*) m

call reada(m)
call locase(m)

END SUBROUTINE readl

!-----------------------------------------------------------------------

SUBROUTINE getf(A,factor)
!  Read the next item as a double-precision number, reading new data
!  records if necessary.
!  If the optional argument factor is present, the value read should be
!  divided by it. (External value = factor*internal value)

DOUBLE PRECISION, INTENT(INOUT) :: A
DOUBLE PRECISION, INTENT(IN), OPTIONAL :: factor

LOGICAL :: eof

do
  if (item .lt. nitems) then
    call readf(a,factor)
    exit
  else
    call read_line(eof)
    if (eof) then
      print "(A)", "End of file while attempting to read a number"
      stop
    endif
  endif
end do

END SUBROUTINE getf

!-----------------------------------------------------------------------

SUBROUTINE geti(I)
!  Get an integer, reading new data records if necessary.
INTEGER, INTENT(INOUT) :: i
LOGICAL :: eof

do
  if (item .lt. nitems) then
    call readi(i)
    exit
  else
    call read_line(eof)
    if (eof) then
      print "(A)", "End of file while attempting to read a number"
      stop
    endif
  endif
end do

END SUBROUTINE geti

!-----------------------------------------------------------------------

SUBROUTINE getli(I)
!  Get a long integer, reading new data records if necessary.
INTEGER(l), INTENT(INOUT) :: i
LOGICAL :: eof

do
  if (item .lt. nitems) then
    call readli(i)
    exit
  else
    call read_line(eof)
    if (eof) then
      print "(A)", "End of file while attempting to read a number"
      stop
    endif
  endif
end do

END SUBROUTINE getli

!-----------------------------------------------------------------------

SUBROUTINE geta(m)
!  Get a character string
CHARACTER(LEN=*) m

LOGICAL eof

do
  if (item .lt. nitems) then
    call reada(m)
    exit
  else
    call read_line(eof)
    if (eof) then
      print "(A)", "End of file while attempting to read a character string"
      stop
    endif
  endif
end do

END SUBROUTINE geta

!-----------------------------------------------------------------------

SUBROUTINE reread(k)

INTEGER, INTENT(IN) :: k
!  k>0  Reread from item k
!  k<0  Go back |k| items
!  k=0  Same as k=-1, i.e. reread last item.

if (k .lt. 0) then
  item=item+k
else if (k .eq. 0) then
  item=item-1
else
  item=k-1
endif
if (item .lt. 0) item=0

END SUBROUTINE reread

!-----------------------------------------------------------------------

SUBROUTINE upcase(word)
CHARACTER(LEN=*), INTENT(INOUT) :: word
INTEGER :: i,k

do i=1,len(word)
  k=index(lower_case,word(i:i))
  if (k .ne. 0) word(i:i)=upper_case(k:k)
end do

END SUBROUTINE upcase

!-----------------------------------------------------------------------

SUBROUTINE locase(word)
CHARACTER(LEN=*), INTENT(INOUT) :: word
INTEGER :: i,k

do i=1,len(word)
  k=index(upper_case,word(i:i))
  if (k .ne. 0) word(i:i)=lower_case(k:k)
end do

END SUBROUTINE locase

!-----------------------------------------------------------------------

SUBROUTINE die(c,reflect)

CHARACTER(LEN=*), INTENT(IN) :: c
LOGICAL, INTENT(IN), OPTIONAL :: reflect

call report(c,reflect)

END SUBROUTINE die

!-----------------------------------------------------------------------

SUBROUTINE report(c,reflect)

CHARACTER(LEN=*), INTENT(IN) :: c
LOGICAL, INTENT(IN), OPTIONAL :: reflect
INTEGER :: i, i1, i2, l

CHARACTER(LEN=3) s1, s2

print "(a)", c
if (present(reflect)) then
  if (reflect) then
    l=loc(item)
    i2=min(last,l+20)
    i1=max(1,i2-70)
    s1=" "
    if (i1 .gt. 1) s1="..."
    s2=" "
    if (i2 .lt. last) s2="..."
    if (level .gt. 0) then
      print "(a, I5, a,a)", "Input line ", line(level),                &
          " in file ", trim(file(level))
    else
      print "(a, I5)", "Input line ", line(level)
    endif
    print "(a3,1x,a,1x,a3)", s1, char(i1:i2), s2
    print "(3x,80a1)", (" ", i=i1,l), "*"
  end if
end if
stop

END SUBROUTINE report

!----------------------------------------------------------------

SUBROUTINE assert(test,string,reflect)

LOGICAL, INTENT(IN) :: test
CHARACTER(LEN=*), INTENT(IN) :: string
LOGICAL, INTENT(IN), OPTIONAL :: reflect

if (.not. test) call report(string,reflect)

END SUBROUTINE assert

!----------------------------------------------------------------

INTEGER FUNCTION find_io(start)

!  Find an unused unit number for input or output. Unit n=start is used
!  if available; otherwise n is incremented until an unused unit is found.
!  Unit numbers are limited to the range 1-100; if n reaches 100 the
!  search starts again at 1.

IMPLICIT NONE
INTEGER, INTENT(IN) :: start
LOGICAL :: in_use, exists
CHARACTER(LEN=40) :: string
INTEGER :: n, n0
INTEGER, PARAMETER :: max_unit=99

n0=start
if (n0 .le. 1 .or. n0 .gt. max_unit) n0=1
n=n0
in_use=.true.
do while (in_use)
  inquire(n,opened=in_use,exist=exists)
  if (exists) then
    if (.not. in_use) exit
  else
    write (unit=string,fmt="(a,i3,a)") "Unit number", n, " out of range"
    call report (string)
  endif
  n=n+1
  if (n > max_unit) n=1
  if (n == n0) then
    call report ("No i/o unit available")
  end if
end do
find_io=n

END FUNCTION find_io

!----------------------------------------------------------------

SUBROUTINE read_colour(fmt, colour, clamp)

CHARACTER(LEN=*), INTENT(IN) :: fmt
REAL(kind=sp), INTENT(OUT) :: colour(3)
LOGICAL, INTENT(IN), OPTIONAL :: clamp
CHARACTER(LEN=6) :: x
INTEGER :: i, r, g, b
DOUBLE PRECISION :: c

select case(fmt)
case("GREY","GRAY")
  call readf(c)
  colour(:)=real(c,sp)
case("RGB")
  call readf(colour(1))
  call readf(colour(2))
  call readf(colour(3))
case("RGB255")
  call readf(colour(1),255.0_sp)
  call readf(colour(2),255.0_sp)
  call readf(colour(3),255.0_sp)
case("RGBX")
  call readu(x)
  read (x(1:2),"(z2)") r
  colour(1)=real(r/255d0,sp)
  read (x(3:4),"(z2)") g
  colour(2)=real(g/255d0,sp)
  read (x(5:6),"(z2)") b
  colour(3)=real(b/255d0,sp)
case default
  call die('COLOUR keyword not recognised',.true.)
end select

if (present(clamp)) then
  if (clamp) then
    do i=1,3
      if (colour(i)>1d0) colour(i)=1d0
      if (colour(i)<0d0) colour(i)=0d0
    end do
  end if
end if

END SUBROUTINE read_colour

END MODULE input
