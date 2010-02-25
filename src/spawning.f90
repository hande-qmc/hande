module spawning

! Module for procedures involved in the spawning step of the FCIQMC algorithm.

use const
implicit none

contains

    subroutine choose_ij(occ_list, i ,j, isym)

        ! Randomly choose a pair of spin-orbitals.
        ! In:
        !    occ_list: Integer list of occupied spin-orbitals.
        ! Out:
        !    i, j: randomly selected spin-orbitals.
        !    isym: symmetry label of the (i,j) combination.

        use symmetry, only: sym_table
        use system, only: nel
        use dSFMT_interface, only: genrand_real2

        integer, intent(in) :: occ_list(nel)
        integer, intent(out) :: i,j, isym
        integer :: ind
        real(dp) :: r

        ! We use a triangular indexing scheme to compress 2 electron indices
        ! into 1.
        ! For i/=j and (for an arbitrary choice of i>j), a 1D index of 
        ! a strictly lower triangular array is:
        !   p = (i-1)(i-2)/2 + j,   where 1<=j<i and 1<=p<=n(n-1)/2
        ! This maps the indexing scheme as:
        !    .                  .
        !   2,1  .              1  .
        !   3,1 3,2  .      to  2  3  .
        !   4,1 4,2 4,3  .      3  4  5  .
        ! We want to do the reverse process in order to pick 2 electron labels
        ! from one random number.
        ! Consider the case where j=1.  i can trivially be determined from the
        ! quadratic equation:
        !   i = 3/2 + \sqrt(2p-1.75)
        ! As j<i and, for a fixed i, p increases monotonically with j, the
        ! integer part of i given by the above equation can never exceed the
        ! correct value.  Hence i can be found for arbitrary j by taking the
        ! floor of the preceeding equation.  j follows trivially.
        !
        ! See (for lower triangular arrays rather than strictly lower):
        ! Decoding the sequential indices of a (lower) triangular array
        ! SIS-75-1783,  S Rifkin, CERN report (CERN-DD-75-7).

        ! This might seem odd, but it enables us to pick the (i,j) pair to
        ! excite with half the calls to the random number generator, which
        ! represents a substantial saving. :-)

        r = genrand_real2()
        ind = int(r*nel*(nel-1)/2)+1.

        i = int(1.50_dp + sqrt(2*ind-1.750_dp))
        j = ind - ((i-1)*(i-2))/2

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.

        i = occ_list(i)
        j = occ_list(j)

        ! Symmetry info is a simple lookup...

        isym = sym_table((i+1)/2,(j+1)/2)

    end subroutine choose_ij

end module spawning
