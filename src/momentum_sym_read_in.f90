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
! in get_kpoint_vector. This index is what will be utilised within our calculations,
! enabling us to avoid actually manipulating kpoint vectors.
!
! With our indexing defined we can then store our inverses and cross products within
! lookup tables (sys%read_in%mom_sym%inv_sym and sys%read_in%mom_sym%sym_table
! respectively) initialised within init_read_in_momentum_symmetry in this module
! and obtain the relevant values from a lookup rather than calculation. The size of
! these lookup tables scales as N_sym and N_sym^2 respectively, so for a hypothetical
! 10x10x10 kpoint grid we would only require storing 1 million integer values within
! sym_table. At this point we would have much bigger issues with our two-body integral
! storage and application of our actual methods to the system.

use system

implicit none

contains

    subroutine init_read_in_momentum_symmetry(sys)

        ! Construct the symmetry tables.

        ! In/Out:
        !    sys: system to be studied.  On output the symmetry components are set.

        use system, only: sys_t
        use checking, only: check_allocate
        use errors, only: stop_all

        type(sys_t), intent(inout) :: sys
        ! [VAN]: Why is this a(3) but ksum(sys%lattice%ndim)?
        integer :: i, j, k, ierr, a(3)
        integer :: ksum(sys%lattice%ndim)

        ! Use 1-index in common with model periodic systems.
        sys%sym0 = 1
        sys%nsym = product(sys%read_in%mom_sym%nprop)
        sys%sym_max = sys%nsym
        sys%sym0_tot = sys%sym0
        sys%nsym_tot = sys%nsym
        sys%sym_max_tot = sys%sym_max

        ! Multiple wavevectors in each irrep. Number of wavevectors
        ! depends on kpoint grid used, but absolute maximum less than
        ! 1000 (10x10x10) due to hard coded limits in read_in.

        ! Feasible to calculate and store product and inverse tables,
        ! so we do so to avoid repeated calculation and enable abstraction
        ! of interfaces for various functions. This gives easier sharing of
        ! logic with pg_sym functionality in eg. excit_gen_mol.f90 or
        ! molecular _integrals.F90.
        allocate(sys%read_in%mom_sym%sym_table(sys%nsym, sys%nsym), stat=ierr)
        call check_allocate('sym_table',sys%nsym*sys%nsym,ierr)
        allocate(sys%read_in%mom_sym%inv_sym(sys%nsym), stat=ierr)
        call check_allocate('inv_sym',sys%nsym,ierr)

        call init_basis_momentum_symmetry_info(sys)

        sys%read_in%mom_sym%gamma_sym = 0
        do i = 1, sys%nsym
            call get_kpoint_vector(i, sys%read_in%mom_sym%nprop, a)

            if (all(a == 0)) sys%read_in%mom_sym%gamma_sym = i

        end do
        sys%read_in%pg_sym%gamma_sym = sys%read_in%mom_sym%gamma_sym
        if (sys%tot_sym) sys%symmetry = sys%read_in%mom_sym%gamma_sym

        if (sys%read_in%mom_sym%gamma_sym == 0) call stop_all('init_momentum_symmetry', 'Gamma-point symmetry not found. ')

        do i = sys%sym0, sys%nsym
            do j = i, sys%nsym
                call get_kpoint_vector(i, sys%read_in%mom_sym%nprop, a)
                call get_kpoint_vector(j, sys%read_in%mom_sym%nprop, ksum)
                ksum = modulo(ksum + a, sys%read_in%mom_sym%nprop)
                do k = 1, sys%nsym
                    call get_kpoint_vector(k, sys%read_in%mom_sym%nprop, a)
                    if (is_gamma_sym_periodic_read_in(sys%read_in%mom_sym, ksum - a)) then
                        sys%read_in%mom_sym%sym_table(i,j) = k
                        sys%read_in%mom_sym%sym_table(j,i) = k
                        if (k == sys%read_in%mom_sym%gamma_sym) then
                            sys%read_in%mom_sym%inv_sym(i) = j
                            sys%read_in%mom_sym%inv_sym(j) = i
                        end if
                        exit
                    end if
                end do
            end do
        end do

    end subroutine init_read_in_momentum_symmetry

    subroutine init_basis_momentum_symmetry_info(sys)

        ! Initialises all required information for use of basis kpoint symmetry.

        ! Specifically, allocates and sets pg_sym%nbasis_sym,
        ! pg_sym%nbasis_sym_spin and pg_sym%sym_spin_basis_fns appropriately.

        ! Requires that the sys%nsym, sys%sym0, sys%sym_max etc be set, as well
        ! as symmetry basis function symmetries be set. This is done within
        ! init_momentum_symmetry in momentum_symmetry.f90 and when reading in
        ! system information respectively.

        ! In/Out:
        !   sys: system being studied. On output all information about basis function
        !       symmetry in read_in%mom_sym set as required.

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

    subroutine print_mom_sym_info(sys, io_unit)

        ! Write out momentum symmetry info (for non-model periodic systems).
        ! In:
        !    sys: system to be studied, with all symmetry components set.
        !    io_unit (optional): io unit to print info to.

        use system, only: sys_t
        use parallel, only: parent

        type(sys_t), intent(in) :: sys
        integer, intent(in), optional :: io_unit
        integer :: i, j, iunit, k_vector(sys%lattice%ndim)

        iunit = 6
        if (present(io_unit)) iunit = io_unit

        if (parent) then
            write (iunit,'(1X,a20,/,1X,20("-"),/)') "Symmetry information"
            write (iunit,'(1X,a63,/)') 'The table below gives the label and inverse of each wavevector.'
            write (iunit,'(1X,a5,4X,a7)', advance='no') 'Index','k-point'
            do i = 1, sys%lattice%ndim
                write (iunit,'(3X)', advance='no')
            end do
            write (iunit,'(a7)') 'Inverse'
            do i = 1, sys%nsym
                write (iunit,'(i4,5X)', advance='no') i
                call get_kpoint_vector(i, sys%read_in%mom_sym%nprop, k_vector)
                write (iunit,'(1X,"(")', advance='no')
                write (iunit,'(i3)',advance='no') k_vector(1)
                do j = 2,sys%lattice%ndim
                    write (iunit,'(",",i3)',advance='no') k_vector(j)
                end do
                write (iunit,'(")")', advance='no')
                write (iunit,'(5X,i4)') sys%read_in%mom_sym%inv_sym(i)
            end do
            write (iunit,'()')
            write (iunit,'(1X,a83,/)') &
                "The matrix below gives the result of k_i+k_j to within a reciprocal lattice vector."
            do i = 1, sys%nsym
                do j = 1, sys%nsym
                    write (iunit,'(2X,i0)', advance='no') sys%read_in%mom_sym%sym_table(j,i)
                end do
                write (iunit,'()')
            end do
            write (iunit,'()')
        end if

    end subroutine print_mom_sym_info

    pure function is_gamma_sym_periodic_read_in(mom_sym, kpoint_vector) result(is_gamma_sym)

        ! Checks if symmetry given is the gamma point symmetry.
        ! In:
        !   mom_sym: basis function symmetry information.
        !   kpoint_vector: kpoint to compare, expressed via 3 integers.
        ! Returns:
        !   true: if symmetry provided is gamma sym.
        !   false: otherwise.

        ! For momentum symmetry in real (read in from an FCIDUMP), periodic
        ! systems.

        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        integer, intent(in) :: kpoint_vector(3)
        logical :: is_gamma_sym

        is_gamma_sym = all(modulo(kpoint_vector, mom_sym%nprop) == mom_sym%gamma_point)

    end function is_gamma_sym_periodic_read_in

    pure function mom_sym_conj(read_in, sym_index) result(conj_index)

        ! Returns symmetry index of complex conjugate of provided
        ! sym index.
        ! Since using complex plane waves, e^(ik.r), this will in
        ! general be e^(-ik.r) so just return inverse symmetry.

        ! In:
        !   sys: information about system being studied. We use the
        !       momentum symmetry information.
        !   sym_index: symmetry index to conjgate.
        ! Returns:
        !   Index of complex conjugate of symmetry sym_index.

        ! For momentum symmetry in real (read in from an FCIDUMP), periodic
        ! systems.

        use system, only: sys_read_in_t

        type(sys_read_in_t), intent(in) :: read_in
        integer, intent(in) :: sym_index
        integer :: conj_index

        conj_index = read_in%mom_sym%inv_sym(sym_index)

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

    pure subroutine decompose_trans_sym(orbsym, propbitlen, kpoint_vector)

        ! Takes orbsym value for translationally symmetric wavefunction and
        ! returns representation of three "quantum numbers" in a vector. In
        ! accordance with widely used approach, values stored according to:
        !
        !   orbsym = 1 + \sum_i sym(i) * 2 ** (propbitlen * (i-1))

        ! For an example of this in practice see symmetry_types.f90

        ! In:
        !   sym_index: symmetry index provided in FCIDUMP file.
        !   propbitlen: length of bit representation of each "quantum number"
        !       within isym.
        ! Out:
        !   kpoint_vector: array of 3 quantum numbers stored in isym.

        use const, only: int_32, int_64

        integer(int_64), intent(in) :: orbsym
        integer, intent(in) :: propbitlen
        integer, intent(out) :: kpoint_vector(3)

        ! Use Iand and mask to select only bits in first propbitlen bits of
        ! isym.
        kpoint_vector(1) = int(Iand(orbsym, 2_int_64 ** propbitlen - 1), int_32)
        ! Bit shift to access correct bits of isym.
        kpoint_vector(2) = int(Iand(Ishft(orbsym, -propbitlen), &
                            2_int_64 ** (propbitlen) - 1), int_32)
        kpoint_vector(3) = int(Ishft(orbsym, -(propbitlen * 2)), int_32)

    end subroutine decompose_trans_sym

    pure function get_kpoint_index(kpoint_vector, nprop) result(sym_index)

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
        !   kpoint_vector: array containing representation of kpoint wavevectors.
        !   nprop: condition of periodic bounary conditions used, ie the
        !       supercell dimension.
        ! Returns:
        !   Index of given kpoint within indexing scheme.

        integer, intent(in) :: kpoint_vector(3), nprop(3)
        integer :: sym_index

        ! Want to start from index 1 at gamma point (0,0,0)
        sym_index = 1 + kpoint_vector(1) + nprop(1) * kpoint_vector(2) + &
                    nprop(1) * nprop(2) * kpoint_vector(3)

    end function get_kpoint_index

    pure subroutine get_kpoint_vector(ind, nprop, a)

        ! Get kpoint vector, in 3D, from index defined in get_kpoint_index.

        ! This functonality is mainly for use when initialising symmetry
        ! information within read_in%mom_sym.

        ! In:
        !   ind: index to decode.
        !   nprop: condition of periodic boundary conditions used, ie the
        !       supercell dimension.
        ! Out:
        !   a: array containing kpoint identifier in terms of 3 "quantum numbers".

        use const, only: dp

        integer, intent(in) :: ind, nprop(3)
        integer, intent(out) :: a(3)
        integer :: scratch

        scratch = floor(real(ind-1, kind=dp)/real(nprop(1)*nprop(2), kind=dp))
        a(3) = scratch
        scratch = ind - 1 - a(3) * nprop(1) * nprop(2)
        a(2) = floor(real(scratch, kind=dp)/real(nprop(1), kind=dp))
        a(1) = scratch - a(2) * nprop(1)

    end subroutine get_kpoint_vector

end module momentum_sym_read_in
