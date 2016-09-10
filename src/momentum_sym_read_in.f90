module momentum_sym_read_in

! Module for handing crystal momentum symmetry routines unique to real, periodic systems.
!
! We provide here a very brief overview of the theory underlying the general area
! of solid state electronic structure, as well as a more thorough description of
! the approach taken to impement this within HANDE.
!
! Momentum symmetry for non-model periodic systems
! ------------------------------------------------
!
! When performing calculations within a periodic system we use a basis constructed
! of Bloch functions with periodicity dictated by the underlying Bravais lattice.
! These functions satisfy the Bloch condition given by
!
!           Psi_k (x + t_n) = e^{i k.t_n} \Psi_k(x)     (1)
!
! where t_n is an integer multiple of the underlying lattice periodicity and k is
! a vector describing the behaviour of the function under translation. This gives
! a general functional form of
!
!           Psi(r) =  e^{i k.r} u(r)               (2)
!
! where u(r) is function with the periodicity of the lattice.
!
! Defining a reciprocal (kspace) lattice in the usual way, we have certain
! allowed values of k and can relate all point within the space to those within
! the first Brillouin zone.
!
! Utilising expression (2) and the properties of periodic functions and their
! integrals we obtain a neccessary but not sufficient condition that an integral
! between such functions be nonzero- that the sum of all kvectors be 0 (modulo
! the periodicity of the lattice). Bearing in mind the effect of conjugation and
! lattice periodicity, this condition is:
!
!           <i|a> = 0           if                  k_a - k_i /= 0
!           <ij|ab> = 0         if      k_a + k_b - k_i - k_j /= 0
!
! We can use this knowledge of integral values to reduce integral storage costs and
! improve our excitation generation.
!
! Implementation within HANDE
! ---------------------------
!
! To implement this in a way compatible with the existing symmetry framework within
! hande requires a mapping of all possible kpoint values onto a single (preferably
! contiguous) index. This minimises code duplication between momentum and pg sym.
! This mapping is defined within get_kpoint_index, with the reverse mapping defined
! in get_kpoint_vector. This index is what will be utilised to within our
! calculations, enabling us to avoid actually manipulating kpoint vectors.
!
! With our indexing defined we can then store our inverses and cross products within
! lookup tables (sys%read_in%mom_sym%inv_sym and sys%read_in%mom_sym%sym_table
! respectively) initialised within init_momentum_symmetry within momentum_symmetry.f90
! and obtain the relevant values from a lookup rather than calculation. The size of
! these lookup tables scales as N_sym and N_sym^2 respectively, so for a hypothetical
! 10x10x10 kpoint grid we would only require storing 1 million integer values within
! sym_table. At this point we would have much bigger issues with our two-body integral
! storage and application of our actual method to the system.

use system

implicit none

contains

    subroutine init_basis_momentum_symmetry_info(sys)

        ! Initialises all required information for use of basis kpoint symmetry.

        ! Specifically, allocates and sets pg_sym%nbasis_sym,
        ! pg_sym%nbasis_sym_spin and pg_sym%sym_spin_basis_fns appropriately.

        ! In/Out:
        !   sys: system being studied. On output all information about basis function
        !       symmetry in read_in%mom_sym set as required.
! [review] - AJWT: It is not immediately obious what information this routine expects
! [review] - AJWT: to have already been set in sys (or where it is set).
        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all

        type(sys_t), intent(inout) :: sys
        integer, allocatable :: current_index(:,:)
        integer :: i, ierr, ims

        allocate(sys%read_in%pg_sym%nbasis_sym_spin(2,1:sys%nsym), stat=ierr)
        call check_allocate('nbasis_sym_spin',sys%nsym,ierr)

        allocate(sys%read_in%pg_sym%nbasis_sym(1:sys%nsym), stat=ierr)
        call check_allocate('nbasis_sym',sys%nsym,ierr)

        sys%read_in%pg_sym%nbasis_sym = 0
        sys%read_in%pg_sym%nbasis_sym_spin = 0

        associate(basis_fns=>sys%basis%basis_fns, &
                    nbasis_sym_spin=>sys%read_in%pg_sym%nbasis_sym_spin, &
                    nbasis_sym=>sys%read_in%pg_sym%nbasis_sym)
            do i = 1, sys%basis%nbasis
                ims = (basis_fns(i)%Ms+3)/2
                nbasis_sym_spin(ims, basis_fns(i)%sym) = &
                            nbasis_sym_spin(ims, basis_fns(i)%sym) + 1
                nbasis_sym(basis_fns(i)%sym) = nbasis_sym(basis_fns(i)%sym) + 1
            end do
        end associate

        allocate(sys%read_in%pg_sym%sym_spin_basis_fns(maxval(sys%read_in%pg_sym%nbasis_sym_spin), &
                        2, sys%sym0_tot:sys%sym_max_tot), stat=ierr)
        call check_allocate('nbasis_sym_spin', &
                sys%nsym*maxval(sys%read_in%pg_sym%nbasis_sym_spin)*2, ierr)


        allocate(current_index(sys%sym0_tot:sys%sym_max_tot, 2), stat=ierr)
        call check_allocate('current_index',sys%nsym,ierr)
        current_index = 1

        associate(sym_spin_basis_fns=>sys%read_in%pg_sym%sym_spin_basis_fns, basis_fns=>sys%basis%basis_fns)
            do i = 1, sys%basis%nbasis
                ims = (basis_fns(i)%Ms+3)/2
                sym_spin_basis_fns(current_index(basis_fns(i)%sym, ims), ims, basis_fns(i)%sym) = i
                basis_fns(i)%sym_spin_index = current_index(basis_fns(i)%sym, ims)
                current_index(basis_fns(i)%sym, ims) = current_index(basis_fns(i)%sym, ims) + 1
            end do
        end associate

        deallocate(current_index, stat=ierr)
        call check_deallocate('current_index',ierr)

        sys%read_in%cross_product_sym_ptr => cross_product_periodic_read_in
        sys%read_in%sym_conj_ptr => mom_sym_conj

    end subroutine init_basis_momentum_symmetry_info

! [review] - AJWT: Below, sym is a 3-vector for the kpoint, but in mom_sym_conj it is a symmetry index
! [review] - AJWT: Some consistent naming conventions to distinguish the two would be helpful.
! [review] - AJWT: [Later] Of course neither to be confused with isym, the FCIDUMP symmetry index.
    pure function is_gamma_sym_periodic_read_in(mom_sym, sym) result(is_gamma_sym)

        ! Checks if symmetry given is the gamma point symmetry.
        ! In:
        !   mom_sym: basis function symmetry information.
        !   sym: kpoint to compare, expressed via 3 integers.
        ! Returns:
        !   true: if symmetry provided is gamma sym.
        !   false: otherwise.

        ! For momentum symmetry in real (read in from an FCIDUMP), periodic
        ! systems.

        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        integer, intent(in) :: sym(3)
        logical :: is_gamma_sym

        is_gamma_sym = all(modulo(sym, mom_sym%nprop) == mom_sym%gamma_point)

    end function is_gamma_sym_periodic_read_in

    pure function mom_sym_conj(read_in, sym) result(rsym)

        ! Returns symmetry index of complex conjugate of provided
        ! sym index.
        ! Since using complex plane waves, e^(ik.r), this will in
        ! general be e^(-ik.r) so just return inverse symmetry.

        ! In:
        !   sys: information about system being studied. We use the
        !       momentum symmetry information.
        !   sym: symmetry to compare.
        ! Returns:
        !   Symmetry of complex conjugate of sym.

        ! For momentum symmetry in real (read in from an FCIDUMP), periodic
        ! systems.

        use system, only: sys_read_in_t

        type(sys_read_in_t), intent(in) :: read_in
        integer, intent(in) :: sym
        integer :: rsym

        rsym = read_in%mom_sym%inv_sym(sym)

    end function mom_sym_conj

!--- Cross products ---

    pure function cross_product_periodic_read_in(read_in, a1, a2) result(prod)

        ! Return the index of the symmetry resulting from the cross
        ! product of the symmetries of indexes a1 and a2

        ! In:
        !   mom_sym: basis function symmetry information.
        !   a1, a2: symmetry indexes to return cross product of.
        ! Returns:
        !   Index of direct product of symmetries a1 & a2.

        use system, only: sys_read_in_t

        type(sys_read_in_t), intent(in) :: read_in
        integer, intent(in) :: a1, a2
        integer :: prod

        prod = read_in%mom_sym%sym_table(a1, a2)

    end function cross_product_periodic_read_in

!--- Indexing conversion routines: ---

    pure subroutine decompose_trans_sym(isym, propbitlen, abel_sym)

        ! Takes symmetry index for translationally symmetric wavefunction and
        ! returns representation of three "quantum numbers". In accordance
        ! with approach used in NECI, values stored according to:
        !   isym = 1 + \sum_i sym(i) * 2 ** (propbitlen * (i-1))

        ! For an example of this in practice see symmetry_types.f90

        ! In:
        !   isym: symmetry index provided in FCIDUMP file.
        !   propbitlen: length of bit representation of each "quantum number"
        !       within isym.
        ! Out:
        !   abel_sym: array of 3 quantum numbers stored in isym.

        use const, only: int_32, int_64

        integer(int_64), intent(in) :: isym
        integer, intent(in) :: propbitlen
        integer, intent(out) :: abel_sym(3)

        ! Use Iand and mask to select only bits in first propbitlen bits of
        ! isym.
        abel_sym(1) = int(Iand(isym, 2_int_64 ** propbitlen - 1), int_32)
        ! Bit shift to access correct bits of isym.
        abel_sym(2) = int(Iand(Ishft(isym, -propbitlen), &
                            2_int_64 ** (propbitlen) - 1), int_32)
        abel_sym(3) = int(Ishft(isym, -(propbitlen * 2)), int_32)

    end subroutine decompose_trans_sym

    pure function get_kpoint_index(a, nprop) result(ind)

        ! Converts from kpoint symmetry vector into unique index.
        ! If we know size of unit cell, can calculate unique index by tiling
        ! first along axis 1, then 2, then 3 in 3 dimensions.

        ! By defining a unique mapping between the kpoint grid and a contiguous
        ! indexes we can refer to symmetry via only this index and arrays
        ! specifying the resultant products and conjugates/inverses. This
        ! enables easy compatibility with much of the pg_sym functionality
        ! once these arrays are constructed.

        ! As such this functonality is mainly for use when initialising
        ! symmetry information within read_in%mom_sym.

        ! In:
        !   a: array containing representation of kpoint wavevectors.
        !   nprop: condition of periodic bounary conditions used, ie the
        !       supercell dimension.
        ! Returns:
        !   Index of given kpoint within indexing scheme.

        integer, intent(in) :: a(3), nprop(3)
        integer :: ind

        ! Want to start from index 1 at gamma point (0,0,0)
        ind = 1 + a(1) + nprop(1) * a(2) + nprop(1) * nprop(2) * a(3)

    end function get_kpoint_index

    pure subroutine get_kpoint_vector(ind, nprop, a)

        ! Get kpoint index, in 3D array, from index defined in get_kpoint_index.

        ! this functonality is mainly for use when initialising symmetry
        ! information within read_in%mom_sym.

        ! In:
        !   ind: index to decode.
        !   nprop: condition of periodic boundary conditions used, ie the
        !       supercell dimension.
        ! Out:
        !   a: array containing kpoint identifier in terms of 3 "quantum numbers".

        integer, intent(in) :: ind, nprop(3)
        integer, intent(out) :: a(3)
        integer :: scratch

! [review] - AJWT: Will the real arithmetic ever end up with a number which rounds the wrong way?
! [review] - AJWT: eg. for large nprop, scratch below could end up being 1.999999999 which rownds
! [review] - AJWT: down to 1 rather than getting the correct value 2.
        scratch = real(ind-1)/real(nprop(1)*nprop(2))
        a(3) = int(scratch)
        scratch = ind - 1 - a(3) * nprop(1) * nprop(2)
        a(2) = int(real(scratch)/real(nprop(1)))
        a(1) = scratch - a(2) * nprop(1)

    end subroutine get_kpoint_vector

end module momentum_sym_read_in
