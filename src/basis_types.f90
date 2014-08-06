module basis_types

    use const, only: p, i0

    ! Structs for holding basis information.
    implicit none

    ! Define a spin orbital.
    type basis_fn_t
        ! Set of quantum numbers describing the basis function.
        ! l is used in two different contexts depending upon whether the orbitals
        ! are defined in momentum space or in real space.  Applies only to model
        ! Hamiltonians (e.g. Hubbard model).
        ! Momentum space:
        !     l is the wavevector in terms of the reciprocal lattice vectors of the crystal cell.
        ! Real space:
        !     l is the position of the basis function within the crystal cell in
        !     units of the lattice vectors of the primitive unit cell.
        ! Obviously we should not convert between descriptions within one
        ! calculation! ;-)
        ! For molecular systems l is not used but is allocated with size 0.  This is
        ! solely a legacy 'feature' and could be safely changed, if someone is so
        ! inclined...
        integer, allocatable :: l(:)
        integer :: spatial_index
        ! Index of the irreducible representation spanned by the orbital.  Used only
        ! in systems where point group symmetry is used (e.g.  molecules).  See
        ! notes in pg_symmetry.
        integer :: sym = 0
        ! Index of basis function within the symmetry block.  sym_index = n
        ! indicates the basis function is the fifth in basis_fns array to have the
        ! symmetry given by sym.
        ! Used only with point_group symmetry.
        integer :: sym_index = 0
        ! Index of basis function within the symmetry block.  sym_spin_index = n
        ! indicates the basis function is the fifth in basis_fns array to have the
        ! symmetry given by sym *and* with the spin given by ms.
        ! Used only with point_group symmetry.
        integer :: sym_spin_index = 0
        ! Spin of the electron (1 or -1).
        integer :: ms
        ! single-particle energy of basis function.
        ! model Hamiltonians in momentum space:
        !     sp_eigv is the kinetic energy of the basis function.
        ! model Hamiltonians in real space:
        !     sp_eigv is not set/used.
        ! molecular systems:
        !     sp_eigv is the single-particle energy read in from the FCIDUMP file
        !     (e.g. Hartree--Fock or Kohn--Sham eigenvalue).
        real(p) :: sp_eigv
        !     Store the Lz of the an orbital when reading in from an FCIDUMP
        !     This is later encoded in the symmetry information.
        integer :: lz=0
    end type basis_fn_t

    ! Struct for holding information about the entire basis set of a system.
    type basis_t
        ! Store of information about the (spin) basis functions of the system.
        ! The *odd* indices contain the alpha (spin up) functions.  This is in
        ! contrast to the bit strings used to refer to determinants where the *even*
        ! bits refer to alpha (spin up) functions.  This difference arises because
        ! fortran numbers bits from 0...
        type(basis_fn_t), allocatable :: basis_fns(:) ! (nbasis)

        ! Number of basis functions.
        ! For the Hubbard model this is equal to twice the number of sites as there are
        ! 2 spin orbitals per site, for the Heisenberg model to the number of sites,
        ! for UEG equal to twice the number of k-points within the energy cutoff and for
        ! read in (e.g. molecular) systems the number of single-particle states read in.
        integer :: nbasis

        ! The determinants are stored as a bit string.  Each element of an array is
        ! an integer of kind i0 (containing i0_length bits).
        ! (The bit type has just been deleted from the forthcoming F2008 standard, so we
        ! won't hold our breath until we can use bits directly......)
        ! basis_length is the length of the byte array necessary to contain a bit for
        ! each basis function, i.e. ceiling(nbasis/i0_length).
        ! If separate_strings is true, then we actually store the alpha and beta
        ! strings separately, and so basis_length is 2*ceiling(nbasis/(2*i0_length)).
        integer :: basis_length

        ! DMQMC uses two determinants for each psip to refer to the two indices
        ! of the relevant matrix element. Hence, the byte array used in DMQMC has
        ! 2*basis_length components. There are some procedures which require basis_length
        ! when used for standard FCIQMC but 2*basis_length when used for DMQMC. It is
        ! therefore useful to have a quantity which equal to 2*basis_length for DMQMC and
        ! equal to basis_length for other methods. Then a procedure can use this quantity
        ! and will work for both methods, making it general. This quantity is total_basis_length.
        ! total_basis_length can then be used when we want to refer to *both* determinants
        ! in DMQMC, and hence the entire bitstring.
        ! [review] - AJWT: This looks ok for now but does make me feel a little uneasy.
        ! [review] - AJWT:   It feels like either it should be the whole basis for DMQMC
        ! [review] - AJWT:   or perhaps we need to abstract the concept of the space we're
        ! [review] - AJWT:   working in.  Either way it's something to keep an eye on.
        integer :: total_basis_length

        ! A determinant is stored in the array f(nbasis).  A basis function is occupied
        ! in the determinant if the relevant bit is set.  The relevant bit is given by
        ! bit_element, the element of the array which contains the bit corresponding to
        ! the basis function, and bit_position, which contains the position of the bit
        ! within the given element.  bit_lookup(:,i) gives the (/ bit_position,
        ! bit_element /) of the i-th basis function.
        ! Note fortran numbers bits starting from 0.
        integer, allocatable :: bit_lookup(:,:) ! (2, nbasis)

        ! The reverse lookup to bit_lookup.
        ! basis_lookup(i,j) gives the basis function corresponding to
        ! the i-th bit in the j-th element of a determinant array.
        integer, allocatable :: basis_lookup(:,:) ! (0:i0_end, basis_length)
    end type basis_t

    contains

        subroutine init_basis_strings(b, separate_strings)

            ! Initialise the string information in a basis_t object for
            ! converting a bit-string representation of a list of orbitials to/from the
            ! integer list.

            ! In/Out:
            !   b: basis_t object to be set.  On input b%nbasis must be set.  On output
            !   the bit string look up tables
            ! In:
            !   separate_strings: true if the alpha and beta orbitals are stored in
            !       separate strings rather than in alternating within one string.

            use calc, only: doing_calc, dmqmc_calc
            use const, only: i0_end, i0_length
            use checking, only: check_allocate

            type(basis_t), intent(inout) :: b
            logical, intent(in) :: separate_strings

            integer :: i, bit_element, bit_pos, ierr

            if (separate_strings) then
                b%basis_length = 2*ceiling(real(b%nbasis)/(2*i0_length))
            else
                b%basis_length = ceiling(real(b%nbasis)/i0_length)
            end if

            if(doing_calc(dmqmc_calc)) then
                b%total_basis_length = 2*b%basis_length
            else
                b%total_basis_length = b%basis_length
            end if

            ! Lookup arrays.
            allocate(b%bit_lookup(2,b%nbasis), stat=ierr)
            call check_allocate('b%bit_lookup',2*b%nbasis,ierr)
            allocate(b%basis_lookup(0:i0_end,b%basis_length), stat=ierr)
            call check_allocate('b%basis_lookup',i0_length*b%basis_length,ierr)
            b%basis_lookup = 0

            if (separate_strings) then
                do i = 1, b%nbasis-1, 2
                    ! find position of alpha orbital
                    bit_pos = mod((i+1)/2, i0_length) - 1
                    if (bit_pos == -1) bit_pos = i0_end
                    bit_element = ((i+1)/2+i0_end)/i0_length
                    b%bit_lookup(:,i) = (/ bit_pos, bit_element /)
                    b%basis_lookup(bit_pos, bit_element) = i
                    ! corresponding beta orbital is in the same position in the
                    ! second half of the string.
                    bit_element = bit_element + b%basis_length/2
                    b%bit_lookup(:,i+1) = (/ bit_pos, bit_element /)
                    b%basis_lookup(bit_pos, bit_element) = i+1
                end do
            else
                do i = 1, b%nbasis
                    bit_pos = mod(i, i0_length) - 1
                    if (bit_pos == -1) bit_pos = i0_end
                    bit_element = (i+i0_end)/i0_length
                    b%bit_lookup(:,i) = (/ bit_pos, bit_element /)
                    b%basis_lookup(bit_pos, bit_element) = i
                end do
            end if

        end subroutine init_basis_strings

        subroutine print_basis_metadata(b, nel, heisenberg_system)

            ! Print out metadata regarding the basis (but not the basis itself).

            ! In:
            !    b: basis_t to be (partly) printed.
            !    nel: number of electrons in the system
            !    heisenberg_system: true if the system of interest is a Heisenberg model.

            use utils, only: int_fmt
            use parallel, only: parent
            use const, only: i0_length

            type(basis_t), intent(in) :: b
            integer, intent(in) :: nel
            logical, intent(in) :: heisenberg_system

            character(4) :: fmt1(4)

            if (parent) then
                fmt1 = int_fmt((/nel, b%nbasis, i0_length, b%basis_length/), padding=1)
                if (heisenberg_system) then
                    write (6,'(1X,a22,'//fmt1(1)//')') 'Number of alpha spins:', nel
                else
                    write (6,'(1X,a20,'//fmt1(1)//')') 'Number of electrons:', nel
                end if
                write (6,'(1X,a26,'//fmt1(2)//')') 'Number of basis functions:', b%nbasis
                write (6,'(/,1X,a61,'//fmt1(3)//')') 'Bit-length of integers used to store determinant bit-strings:', i0_length
                write (6,'(1X,a57,'//fmt1(4)//',/)') &
                    'Number of integers used to store determinant bit-strings:', b%basis_length
            end if

        end subroutine print_basis_metadata

        subroutine copy_basis_t(b1, b2)

            ! Copy a basis_t object into a new basis_t object.

            ! In:
            !   b1: basis_t to be copied.
            ! Out:
            !   b2: destination basis_t.

            use const, only: i0_end

            type(basis_t), intent(in) :: b1
            type(basis_t), intent(out) :: b2

            integer :: i, n

            b2%basis_length = b1%basis_length
            b2%total_basis_length = b1%total_basis_length
            b2%nbasis = b1%nbasis

            if (allocated(b1%bit_lookup)) then
                allocate(b2%bit_lookup(2,b2%nbasis), source=b1%bit_lookup)
            end if
            if (allocated(b1%basis_lookup)) then
                allocate(b2%basis_lookup(0:i0_end,b2%basis_length), source=b1%basis_lookup)
            end if

            if (allocated(b1%basis_fns)) then
                allocate(b2%basis_fns(b2%nbasis))
                if (allocated(b1%basis_fns(1)%l)) then
                    n = size(b1%basis_fns(1)%l)
                    do i = 1, b2%nbasis
                        allocate(b2%basis_fns(i)%l(n))
                        b2%basis_fns(i) = b1%basis_fns(i)
                    end do
                end if
            end if

        end subroutine copy_basis_t

        subroutine dealloc_basis_t(b)

            ! In/Out:
            !    b: basis_t object to be deallocated.

            use checking, only: check_deallocate

            type(basis_t), intent(inout) :: b
            integer :: ierr

            if (allocated(b%bit_lookup)) then
                deallocate(b%bit_lookup, stat=ierr)
                call check_deallocate('b%bit_lookup', ierr)
            end if
            if (allocated(b%basis_lookup)) then
                deallocate(b%basis_lookup, stat=ierr)
                call check_deallocate('b%basis_lookup', ierr)
            end if
            call dealloc_basis_fn_t_array(b%basis_fns)

        end subroutine dealloc_basis_t

        subroutine dealloc_basis_fn_t_array(basis_fns)

            ! In/Out:
            !    basis_fns: array of basis_fn_t to be deallocated.

            use checking, only: check_deallocate

            type(basis_fn_t), intent(inout), allocatable :: basis_fns(:)
            integer :: i, ierr

            if (allocated(basis_fns)) then
                do i = lbound(basis_fns, dim=1), ubound(basis_fns, dim=1)
                    if (allocated(basis_fns(i)%l)) then
                        deallocate(basis_fns(i)%l, stat=ierr)
                        call check_deallocate('basis_fns(i)%l', ierr)
                    end if
                end do
                deallocate(basis_fns, stat=ierr)
                call check_deallocate('basis_fns', ierr)
            end if

        end subroutine dealloc_basis_fn_t_array

end module basis_types
