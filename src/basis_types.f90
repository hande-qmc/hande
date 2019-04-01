module basis_types

    use const, only: p, i0

    ! Structs for holding basis information.
    implicit none

    ! Define a spin orbital.
    type basis_fn_t
        ! Set of quantum numbers describing the basis function.
        ! l is used in two different contexts depending upon whether the orbitals
        ! are defined in momentum space or in real space.  Applies to model
        ! Hamiltonians (e.g. Hubbard model) and momentum space read_in
        ! Hamiltonians (ie. periodic systems).
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
        ! contrast to the bit strings (see below) used to refer to determinants
        ! where the *even* bits refer to alpha (spin up) functions.  This
        ! difference arises because fortran numbers bits from 0 but arrays start
        ! (by default) from 1...
        type(basis_fn_t), allocatable :: basis_fns(:) ! (nbasis)

        ! Number of basis functions.
        ! For the Hubbard model this is equal to twice the number of sites as there are
        ! 2 spin orbitals per site, for the Heisenberg model to the number of sites,
        ! for UEG equal to twice the number of k-points within the energy cutoff and for
        ! read in (e.g. molecular) systems the number of single-particle states read in.
        integer :: nbasis

        ! We commonly store the many-particle basis functions (e.g. spin
        ! products for the Heisenberg model, determinants for fermions) as a bit
        ! string.  The bit string is stored in an array of i0 integers.
        ! tot_string_len gives the size of this array.
        ! The alpha and beta orbitals are interleaved and bit_string_len is
        ! ceiling(nbasis/i0_length).
        !
        ! Note it's much more efficient to do operations on 32-bit or 64-bit
        ! integers than individual bits (or, indeed smaller integers) on modern
        ! architectures, so there's no need to cry about the bit type being
        ! removed from proposed F2008 standard or the wasted space at the end of
        ! the bit array.
        integer :: bit_string_len
        ! We can store additional information in the bit string to enable different
        ! sorting orders. This takes the form of an additional info_string_len type-i0
        ! integers at the end of the bit string.
        ! NB This information should be removed when manipulating bit strings, eg. in
        ! excitation generators. As it is currently only used within CCMC this simply
        ! involves removing this information when selecting a cluster, but for any
        ! other use it may be more involved.
        integer :: info_string_len
        ! The overall length of the bit string in type-i0 integers is the sum of
        ! bit_string_len and info_string_len
        integer :: tot_string_len

        ! The QMC algorithms implemented essentially sample a tensor of arbitrary rank.
        ! For example, FCIQMC samples a vector (rank 1) and DMQMC a matrix (rank 2).
        ! It is sometimes convenient to store/do operations on a bit string
        ! formed from concatenating the bit strings of the individual tensor
        ! labels.  Then length of this array is given by tensor_label_len.
        ! tensor_label_len = (rank of tensor) * tot_string_len.
        ! NOTE: this must be set before running a QMC algorithm where the rank is not 1.
        integer :: tensor_label_len

        ! Bit masks to reveal the list of alpha basis functions and beta functions occupied
        ! in a Slater determinant.
        integer(i0) :: alpha_mask, beta_mask

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
        integer, allocatable :: basis_lookup(:,:) ! (0:i0_end, bit_string_len)

        ! excit_mask(:,i) is a bit field with bits corresponding to all orbitals with
        ! a higher index than i set.
        integer(i0), allocatable :: excit_mask(:,:) ! (tot_string_len, nbasis)

    end type basis_t

    contains

        subroutine init_basis_strings(b)

            ! Initialise the string information in a basis_t object for
            ! converting a bit-string representation of a list of orbitials to/from the
            ! integer list.

            ! In/Out:
            !   b: basis_t object to be set.  On input b%nbasis must be set.
            !      On output the bit string look-up tables are also set.

            use const, only: i0_end, i0_length
            use checking, only: check_allocate

            type(basis_t), intent(inout) :: b

            integer :: i, bit_element, bit_pos, ierr

            ! Assume not storing any additional information within bit string initially.
            b%bit_string_len = ceiling(real(b%nbasis) / i0_length)
            b%info_string_len = 0
            b%tot_string_len = b%bit_string_len + b%info_string_len
            b%tensor_label_len = b%tot_string_len

            ! Lookup arrays.
            allocate(b%bit_lookup(2,b%nbasis), stat=ierr)
            call check_allocate('b%bit_lookup',2*b%nbasis,ierr)
            allocate(b%basis_lookup(0:i0_end,b%bit_string_len), stat=ierr)
            call check_allocate('b%basis_lookup',i0_length*b%bit_string_len,ierr)
            b%basis_lookup = 0

            do i = 1, b%nbasis
                bit_pos = mod(i, i0_length) - 1
                if (bit_pos == -1) bit_pos = i0_end
                bit_element = (i+i0_end)/i0_length
                b%bit_lookup(:,i) = (/ bit_pos, bit_element /)
                b%basis_lookup(bit_pos, bit_element) = i
            end do

            ! Bit masks...
            ! Alpha basis functions are in the even bits.  alpha_mask = 01010101...
            ! Beta basis functions are in the odd bits.    beta_mask  = 10101010...
            b%alpha_mask = 0_i0
            b%beta_mask = 0_i0
            do i = 0, i0_end
                if (mod(i,2)==0) then
                    b%alpha_mask = ibset(b%alpha_mask,i)
                else
                    b%beta_mask = ibset(b%beta_mask,i)
                end if
            end do

        end subroutine init_basis_strings

        pure subroutine reset_extra_info_bit_string(basis, f)

            type(basis_t), intent(in) :: basis
            integer(i0), intent(inout) :: f(:)

            ! Just clear all information stored within information section
            ! of bit string.
            if (basis%info_string_len/=0) then
                f(basis%bit_string_len+1:) = 0_i0
            end if

        end subroutine reset_extra_info_bit_string

        subroutine print_basis_metadata(b, nel, heisenberg_system, io_unit)

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
            integer, intent(in), optional :: io_unit

            character(4) :: fmt1(4)
            integer :: iunit

            iunit = 6
            if (present(io_unit)) iunit = io_unit

            if (parent) then
                fmt1 = int_fmt((/nel, b%nbasis, i0_length, b%tot_string_len/), padding=1)
                if (heisenberg_system) then
                    write (iunit,'(1X,a22,'//fmt1(1)//')') 'Number of alpha spins:', nel
                else
                    write (iunit,'(1X,a20,'//fmt1(1)//')') 'Number of electrons:', nel
                end if
                write (iunit,'(1X,a26,'//fmt1(2)//')') 'Number of basis functions:', b%nbasis
                write (iunit,'(/,1X,a61,'//fmt1(3)//')') 'Bit-length of integers used to store determinant bit-strings:', i0_length
                write (iunit,'(1X,a57,'//fmt1(4)//',/)') &
                    'Number of integers used to store determinant bit-strings:', b%tot_string_len
            end if

        end subroutine print_basis_metadata

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
