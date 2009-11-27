module comb_m

! Module for calculating combination sets.

! For further details see:
!    Algorithm 515: Generation of a vector from the Lexicographical index [G6]
!    B. P. Beckles and M. Lybanon.
!    ACM Transactions on Mathematical Software, 3 180-182 (1977).

! Translated into Fortran 90 by James Spencer.

implicit none

contains

    pure function comb(n, p, l) result(c)

        ! ACM Algorithm 515
        ! Find the combination set of n things taken p at a time for a given
        ! lexicographical index, l.
        ! In:
        !    n: number of items in set.
        !    p: number of items in each combination. 1 <= p <= n.
        !    l: lexicographical index. 1 <= l <= ^nC_p.
        ! Returns:
        !    c: output array of length p containing the combination set.


        integer, intent(in) :: n, p, l
        integer :: c(p)
        integer :: k, r, p1,i

        ! Initialise lower bound at zero.
        k = 0
        ! Loop to select elements in ascending order.
        p1 = p-1
        do i = 1, p1
            ! Set lower bound for next element value.
            c(i) = 0
            if (i /= 1) c(i) = c(i-1)
            ! Loop to check validity of each element value.
            do
                c(i) = c(i) + 1
                r = binom(n-c(i), p-i)
                k = k + r
                if (k >= l) exit
            end do
            k = k - r
        end do
        c(p) = c(p1) + l - k

    end function comb

    pure function binom(m, n)
        ! ACM Algorithm 160 translated to Fortran.
        ! Returns the binomial coefficient ^mC_n: the number of
        ! combinations of m things taken n at a time.

        integer :: binom
        integer, intent(in) :: m, n
        integer :: p,i, n1, r

        n1 = n
        p = m - n1
        if (n1 < p) then
            p = n1
            n1 = m - p
        end if
        r = n1 + 1
        if (p == 0) r = 1
        do i = 2, p
            r = (r*(n1+i))/i
        end do
        binom = r

    end function binom

end module comb_m

program test_comb

    use comb_m

    integer :: i
    integer :: c(2)

    write (6,*) binom(10,4)

    do i = 1,6
        write (6,*) comb(4,2,i)
    end do

end program test_comb
