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
        integer, pointer :: l(:) => NULL()
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

        ! number of basis functions.
        ! For the Hubbard model is equal to twice the number of sites as there are
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

        ! DMQMC uses two determinants for each psip to refer to the two components
        ! of the relevant matrix element. Hence, the bitstring which is stored in DMQMC has
        ! 2*basis_length components. There are some procedures which required basis_length
        ! when used for stahndard FCIQMC but 2*basis_length when used for DMQMC. It is
        ! therefore useful to have a quantity which equal to 2*basis_length for DMQMC and
        ! equal to basis_length for other methods. Then a procedure can use this quantity
        ! and will work for both methods, making it general. This quantity is total_basis_length.
        ! total_basis_length can then be used when we want to refer to *both* determinants
        ! in DMQMC, and hence the entire bitstring.
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

end module basis_types
