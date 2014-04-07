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
! by the null bit string (with some Lz symmetry which is added later)
integer :: gamma_sym = 0

! nbasis_sym(i) gives the number of (spin) basis functions in the i-th symmetry,
! where i is the bit string describing the irreducible representation.
integer, allocatable :: nbasis_sym(:) ! (sym0:sym_max)

! nbasis_sym_spin(1,i) gives the number of spin-down basis functions in the i-th
! symmetry where i is the bit string describing the irreducible representation.
! Similarly, j=2 gives the analagous quantity for spin-up basis functions.
! For RHF calculations nbasis_sym_spin(:,i) = nbasis_sym(i)/2.
integer, allocatable :: nbasis_sym_spin(:,:) ! (2,sym0:sym_max)

! sym_spin_basis_fns(:,ims,isym) gives the list of spin functions (ims=1 for
! down, ims=2 for up) with symmetry isym.  We merrily waste some memory (not all
! symmetries will have the same number of basis functions), so a 0 entry
! indicates no more basis functions with the given spin and spatial symmetries.
integer, allocatable :: sym_spin_basis_fns(:,:,:) ! (max(nbasis_sym_spin),2,sym0:sym_max)

integer :: pg_mask= 7  !The 3 symmetry generators.

!Lz symmetry is stored in the higher bits of the symmetry.
integer :: Lz_max=10  !This is set to I functions
integer :: Lz_mask=63*8 !We allow -32...31 as our range for Lz in the bits above the Abelian symmetry
! Used to offset the Lz values when encoding into a symmetry
integer :: Lz_offset
contains

    subroutine init_pg_symmetry(sys)

        ! Initialise point group symmetry information.
        ! *Must* be called after basis functions are initialised and have their
        ! symmetries set from the FCIDUMP file.

        ! In/Out:
        !    sys: system being studied.  On output the symmetry fields are set.

        use checking, only: check_allocate

        use basis, only: basis_fns, nbasis
        use system, only: sys_t
        use calc, only: sym_in

        type(sys_t), intent(inout) :: sys

        integer :: i, ierr, ims, ind

        integer :: maxsym, maxLz, maxpg
        ! molecular systems use symmetry indices starting from 0.
        sys%sym0 = 0

        ! Given n generators, there must be 2^n irreducible representations.
        ! in the working space.  (Note that the wavefunctions might only span
        ! a subspace of the point group considered by the generating quantum
        ! chemistry package.)
        ! Set sys%nsym to be 2^n so that we always work in the smallest group
        ! spanned by the wavefunctions.
        ! +1 comes from the fact that the basis_fn%sym gives the bit string
        ! representation.

        maxsym = 2**ceiling(log(real(maxval(basis_fns(:)%sym)+1))/log(2.0))
        write(6,*) "Maximum symmetry found", maxsym
        maxLz = maxval(basis_fns(:)%lz)
        write(6,*) "Maximum Lz found", maxLz 
        pg_mask=maxsym-1
        write(6,*) "Symmetry Mask:", pg_mask

        !The maximum encountered value of Lz will be maxLz*3 (and the minimum the negative of this)
        !We allocate bits for this above the pointgroup symmetry bits.
        Lz_mask=(2**ceiling(log(real(6*maxLz+1))/log(2.0))-1)*maxsym
        write(6,*) "Lz Mask:", Lz_mask
        !In Lz world, 0 means Lz=-3*maxLz*maxsym, so Lz_offset means Lz=0
        Lz_offset=3*maxLz*maxsym
        write(6,*) "Lz offset (corresponds to Lz=0):", Lz_offset
        gamma_sym=Lz_offset*maxsym
        write(6,*) "Totally symmetric symmetry: ", gamma_sym
        if(sym_in /= huge(1)) then
            !Need to modify to include Lz:
            sym_in=sym_in+Lz_offset
        endif
        sys%nsym = (6*maxLz+1)*maxsym
        sys%sym_max = sys%nsym-1

        allocate(nbasis_sym(0:sys%nsym-1), stat=ierr)
        call check_allocate('nbasis_sym', sys%nsym, ierr)
        allocate(nbasis_sym_spin(2,0:sys%nsym-1), stat=ierr)
        call check_allocate('nbasis_sym_spin', 2*sys%nsym, ierr)

        nbasis_sym = 0
        nbasis_sym_spin = 0

        do i = 1, nbasis
            !Encode the Lz into the symmetry.
            basis_fns(i)%sym = basis_fns(i)%sym + (basis_fns(i)%lz*maxsym+Lz_offset)
            nbasis_sym(basis_fns(i)%sym) = nbasis_sym(basis_fns(i)%sym) + 1
            basis_fns(i)%sym_index = nbasis_sym(basis_fns(i)%sym)

            ims = (basis_fns(i)%Ms+3)/2 ! Ms=-1,1 -> ims=1,2

            nbasis_sym_spin(ims,basis_fns(i)%sym) = nbasis_sym_spin(ims,basis_fns(i)%sym) + 1
            basis_fns(i)%sym_spin_index = nbasis_sym_spin(ims,basis_fns(i)%sym)

        end do

        allocate(sym_spin_basis_fns(maxval(nbasis_sym_spin),2,0:sys%nsym-1), stat=ierr)
        call check_allocate('sym_spin_basis_fns', size(sym_spin_basis_fns), ierr)
        sym_spin_basis_fns = 0

        do i = 1, nbasis
            ims = (basis_fns(i)%Ms+3)/2 ! Ms=-1,1 -> ims=1,2
            ind = minloc(sym_spin_basis_fns(:,ims,basis_fns(i)%sym), dim=1) ! first non-zero element
            sym_spin_basis_fns(ind, ims, basis_fns(i)%sym) = i
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
        if (allocated(nbasis_sym_spin)) then
            deallocate(nbasis_sym_spin, stat=ierr)
            call check_deallocate('nbasis_sym_spin', ierr)
        end if
        if (allocated(sym_spin_basis_fns)) then
            deallocate(sym_spin_basis_fns, stat=ierr)
            call check_deallocate('sym_spin_basis_fns', ierr)
        end if

    end subroutine end_pg_symmetry

    subroutine print_pg_symmetry_info(sys)

        ! Write out point group symmetry information.

        ! In:
        !    sys: system being studied.

        use system, only: sys_t
        use parallel, only: parent

        type(sys_t), intent(in) :: sys

        integer :: i, j

        if (parent) then
            write (6,'(1X,a20,/,1X,20("-"),/)') "Symmetry information"

            write (6,'(1X,a78,/)') 'The matrix below gives the direct products of the irreducible representations.'
            ! Note that we never actually store this.
            do i = 0, sys%nsym-1
                do j = 0, sys%nsym-1
                    write (6,'(1X,i2)',advance='no') cross_product_pg_sym(i,j)
                end do
                write (6,'()')
            end do

            write (6,'(/,1X,a93,/)') 'The table below gives the number of basis functions spanning each irreducible representation.'

            write (6,'(1X,"irrep  nbasis  nbasis_up  nbasis_down")')
            do i = 0, sys%nsym-1
                write (6,'(1X,i3,4X,i5,3X,i5,6X,i5)') i, nbasis_sym(i), nbasis_sym_spin(:,i)
            end do

            write (6,'()')

        end if

    end subroutine print_pg_symmetry_info

    elemental function cross_product_pg_basis(i, j) result(sym_ij)

        ! In:
        !    i,j: (indices of) spin-orbitals.
        ! Returns:
        !    The bit string representation of the irreducible representation
        !    formed from the direct product i%sym \cross j%sym, where i%sym is
        !    the irreducible representation spanned by the i-th spin-orbital.

        use basis, only: basis_fns

        integer :: sym_ij
        integer, intent(in) :: i, j

        sym_ij = cross_product_pg_sym(basis_fns(i)%sym, basis_fns(j)%sym)

    end function cross_product_pg_basis

    elemental function cross_product_pg_sym(sym_i, sym_j) result(sym_ij)

        ! In:
        !    sym_i,sym_j: bit string representations of irreducible
        !    representations of a point group.
        ! Returns:
        !    The bit string representation of the irreducible representation
        !    formed from the direct product sym_i \cross sym_j.

        integer :: sym_ij
        integer, intent(in) :: sym_i, sym_j

        sym_ij = ior(iand(ieor(sym_i, sym_j),pg_mask), &
                iand(sym_i,Lz_mask)+iand(sym_j,Lz_mask)-Lz_offset)
    end function cross_product_pg_sym


    elemental function pg_sym_conj(sym) result(rsym)
        integer, intent(in) :: sym
        integer :: rsym

        !Take the symmetry conjugate.  The point group part is the same.
        !The Lz needs to become -Lz, but
        rsym  = ior(iand(sym,pg_mask), &
                iand(2*Lz_offset-iand(sym,Lz_mask),Lz_mask))

    end function pg_sym_conj
    elemental function is_gamma_irrep_pg_sym(sym) result(is_gamma)

        ! In:
        !    sym: bit string representation of an irreducible representation of
        !    a point group.
        ! Returns:
        !    True if sym represents the Gamma_1 (totally symmetric) irreducible
        !    representation.

        logical :: is_gamma
        integer, intent(in) :: sym

        is_gamma = gamma_sym == sym

    end function is_gamma_irrep_pg_sym

    pure function symmetry_orb_list_mol(orb_list) result(isym)

        ! In:
        !    orb_list: list of orbitals (e.g. determinant).
        ! Returns:
        !    symmetry index of list (i.e. direct product of the representations
        !    of all the orbitals in the list).

        use basis, only: basis_fns

        integer :: isym
        integer, intent(in) :: orb_list(:)

        integer :: i

        isym = gamma_sym
        do i = lbound(orb_list, dim=1), ubound(orb_list, dim=1)
            isym = cross_product_pg_sym(isym, basis_fns(orb_list(i))%sym)
        end do

    end function symmetry_orb_list_mol

end module point_group_symmetry
