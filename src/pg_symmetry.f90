module point_group_symmetry

! Module for handling point group symmetry, as read in from FCIDUMP files.

! This was made much easier thanks to conversations with Alex Thom...

! NOTE:
! It seems physicists are far less obsessed with point group symmetries than my
! undergraduate chemistry course.  I recommend a thorough reading of the
! relevant sections of the classic book 'Group Theory and Quantum Mechanics' by
! Tinkham.

! Point group symmetry
! --------------------
!
! The quantum chemistry packages we use to generate FCIDUMP files only implement
! D2h symmetry (and subgroups thereof).  Whilst this means some symmetries are
! not considered, the advantage for us is that all point groups we will consider
! are real and Abelian.  Thus:
!
! * all irreducible representations are 1D.
! * all operations are their own inverse
! * \Gamma_i \cross \Gamma_i = \Gamma_1, where \Gamma_i is an arbitrary
!   irreducible representation, \cross indicates direct product and \Gamma_1 is
!   the totally-symmetric representation.
! * all irreducible representations can be represented by (at most) three
!   generators and the behaviour of a function under those generators is
!   sufficient to completely determine the symmetry of that function.
!
! Furthermore, we assume that each basis function spans one (and only one)
! irreducible representation---this can always be done when working in D2h (or
! a subgroup thereof).
!
! An irreducible representation is labelled by its behaviour under the
! generators of the point group using a bit string where the i-th bit
! corresponds to the i-th generator.  A set bit indicates that the
! representation is *odd* with respect to that generator and an unset bit
! indicates that the representation is *even* with respect to that generator.
! Thus the totally symmetric representation is always labelled by the 0 bit
! string.
!
! Direct products are easily evaluated: even \cross even = even,
! odd \cross odd = even and odd \cross even = odd.  Hence the direct product can
! be found by taking XOR of the representations involved.
!
! As we consider (at most) D2h (3 generators) we can just use a standard integer
! to represent the irreducible representations.  The higher bits are wasted, but
! the memory used is minimal.

use const

implicit none

contains

end module point_group_symmetry
