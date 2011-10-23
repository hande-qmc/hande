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
! Direct products are easily evaluated by considering the characters of the
! generator operations in the product: even \cross even = even, odd \cross odd = even 
! and odd \cross even = odd.  Hence the direct product can be found by
! taking XOR of the representations involved.
!
! As we consider (at most) D2h (3 generators) we can just use a standard integer
! to represent the irreducible representations.  The higher bits are wasted, but
! the memory used is minimal.
!
! For a point group containing n generators, there are 2^n irreducible
! representations.  Due to the bit representation described above, these
! representations are labelled by the set of integers {0,1,...2^n-1}.

use const

implicit none

! Following the above discussion, the totally symmetric representation is given
! by the null bit string.
integer, parameter :: gamma_sym = 0

! nbasis_sym(i) gives the number of (spin) basis functions in the i-th symmetry,
! where i is the bit string describing the irreducible representation.
integer, allocatable :: nbasis_sym(:)

! nbasis_sym_spin(1,i) gives the number of spin-down basis functions in the i-th
! symmetry where i is the bit string describing the irreducible representation.
! Similarly, j=2 gives the analagous quantity for spin-up basis functions.
! For RHF calculations nbasis_sym_spin(:,i) = nbasis_sym(i)/2.
integer, allocatable :: nbasis_sym_spin(:,:)

contains

    subroutine init_pg_symmetry()

        ! Initialise point group symmetry information.
        ! *Must* be called after basis functions are initialised and have their
        ! symmetries set from the FCIDUMP file.

        use checking, only: check_allocate

        use basis, only: basis_fns, nbasis
        use symmetry, only: nsym

        integer :: i, ierr

        nsym = maxval(basis_fns(:)%sym) + 1
        
        allocate(nbasis_sym(0:nsym-1), stat=ierr)
        call check_allocate('nbasis_sym', nsym, ierr)
        allocate(nbasis_sym_spin(2,0:nsym-1), stat=ierr)
        call check_allocate('nbasis_sym_spin', 2*nsym, ierr)

        nbasis_sym = 0
        nbasis_sym_spin = 0

        do i = 1, nbasis
            nbasis_sym(basis_fns(i)%sym) = nbasis_sym(basis_fns(i)%sym) + 1
            basis_fns(i)%sym_index = nbasis_sym(basis_fns(i)%sym)
            if (basis_fns(i)%ms == -1) then
                nbasis_sym_spin(1,basis_fns(i)%sym) = nbasis_sym_spin(1,basis_fns(i)%sym) + 1
                basis_fns(i)%sym_spin_index = nbasis_sym_spin(1,basis_fns(i)%sym)
            else
                nbasis_sym_spin(2,basis_fns(i)%sym) = nbasis_sym_spin(2,basis_fns(i)%sym) + 1
                basis_fns(i)%sym_spin_index = nbasis_sym_spin(2,basis_fns(i)%sym)
            end if
        end do

    end subroutine init_pg_symmetry

    subroutine end_pg_symmetry()

        ! Deallocate arrays containing point group symmetry information.

        use checking, only: check_deallocate

        integer :: ierr

        if (allocated(nbasis_sym)) then
            deallocate(nbasis_sym, stat=ierr)
            call check_deallocate('nbasis_sym', ierr)
        end if

    end subroutine end_pg_symmetry

end module point_group_symmetry
