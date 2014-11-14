module molecular_integral_types

! Module for derived types used to contain molecular/generic integrals.

! We only define the types here; all procedures are in molecular_integrals.

use const, only: p
use base_types, only: alloc_rp1d

implicit none

! NOTE: interaction with integral stores is best done using the store_* and get_*
! procedures provided in molecular_integrals rather than directly accessing them.

! Store for one-body integrals, <i|o|j>, where i,j are spin basis functions and
! o is a one-electron operator.
type one_body
    ! integrals(ispin, isym)%v(indx) corresponds to the <i|o|j> integral (assuming
    ! i,j conserve spin and spatial symmetry), where ispin and isym index the spin
    ! and spatial symmetry of i and j and indx is the combined (triangular) index of
    ! i and j within that spin and symmetry block.
    ! See access procedures for this in practice.
    ! This data structure makes it possible and relative easy to only store the
    ! integrals which are non-zero by symmetry (ie a small fraction of the possible
    ! integrals).
    ! Note that only one spin channel is needed (and stored) in RHF calculations.
    type(alloc_rp1d), allocatable :: integrals(:,:)
    ! bit string representations of irreducible representations
    integer :: op_sym
    ! From a UHF calculation?
    logical :: uhf
end type one_body

! Store for two-body integrals, <ij|o|ab>, where i,j,a,b are spin basis functions and
! o is a two-electron operator.
type two_body
    ! integrals(ispin)%v(indx) gives the integral <ij|o_2|ab>, where ispin depends upon
    ! the spin combination (ie all alpha, all beta, and haf alpha, half beta) and
    ! indx is related to i,j,a,b.  As we deal with real orbitals only (except
    ! see below), we can use permutation symmetry to reduce the number of
    ! integrals by a factor of 8.  See access procedures for this in action.
    ! Note that only one spin channel is needed (and stored) in RHF calculations.

    ! L_z symmetry:
    ! The use of L_z symmetry complicates matters as, whilst the integrals remain
    ! real, the orbitals are complex and so we lose 8-fold permutation symmetry
    ! and are left with 4-fold permutation symmetry.  However, consider the
    ! integrals < i j | k l > and < k j | i l >, where the orbitals are labelled
    ! by their L_z value.  In order to conserve angular momentum, i+j=k+l and
    ! k+j=i+l => i-j=k-l.  There are two options: either (at most) one integral
    ! is non-zero, in which case we can simply test for whether the integral is
    ! non-zero by symmetry and then look it up using the same scheme as above
    ! (using 8-fold permuation symmetry) or both integrals are non-zero, in
    ! which case i=k and j=l.  This means we have an integral like < i j | i' j' >,
    ! where i and i' are different orbitals with the same L_z.  Due to the
    ! integral transformation, the complex phase factor comes solely from L_z.
    ! Since i and i' have the same L_z, the product (i^* i') is real and hence
    ! (i^* i') = (i'^* i).
    ! In summary < i j | k l> = < k j | i l > if i.Lz=k.Lz and j.Lz=l.Lz,
    ! otherwise only one of < i j | k l> and < k j | i l > is non-zero by
    ! symmetry.  Hence by testing whether an integral is non-zero by symmetry
    ! first, we can continue to use the same lookup scheme for integrals
    ! involving orbitals which have L_z symmetry as we do for the standard
    ! (real) orbitals.

    ! TODO:
    ! * can compress coulomb integral store by ensuring integrand is totally
    !   symmetric, as is done for the one-body integrals.
    type(alloc_rp1d), allocatable :: integrals(:)
    ! bit string representations of irreducible representations
    integer :: op_sym
    ! From a UHF calculation?
    logical :: uhf
end type

end module molecular_integral_types
