module molecular_integrals

! Module for storing and accessing integrals for molecular systems.
! These integrals are previously calculated using a quantum chemistry package
! (e.g. PSI4, MOLPRO or QChem).

use const, only: p
use molecular_integral_types
use const, only: int_64

implicit none

! Indexing type for two_body_t integral stores.
! i.e. an encoding of the 4-index quartet a,b,c,d into an index for the integral store
type, private :: int_indx
    integer :: spin_channel     !If alpha and beta spin-orbitals differ, we store
                                ! different combinations of these in different spin channels
                                ! 1=bbbb, 2=aaaa, 3=baba, 4=ababa
    integer(int_64) :: indx     ! The index within the spin channel
    logical :: conjugate        ! Indicates if we need to take the complex conjugate of the integral
end type int_indx

type, private :: int_ex_indx
    integer :: spin_channel
    integer :: repeat_ind
    integer(int_64) :: triind
    logical :: conjugate
end type int_ex_indx

interface get_one_body_int_mol
    module procedure get_one_body_int_mol_real
    module procedure get_one_body_int_mol_complex
end interface

interface get_two_body_int_mol
    module procedure get_two_body_int_mol_real
    module procedure get_two_body_int_mol_complex
end interface

contains

!--- Memory allocation and deallocation ---

    subroutine init_one_body_t(sys, op_sym, imag, store)

        ! Allocate memory required for the integrals involving a one-body
        ! operator.

        ! In:
        !    uhf: whether integral store is from a UHF calculation or RHF
        !       calculation.
        !    op_sym: bit string representations of irreducible representations
        !       of a point group.  See point_group_symmetry.
        !    imag: whether integral store contains the imaginary component of
        !       complex integrals.
        ! Out:
        !    store: one-body integral store with components allocated to hold
        !       interals.  Note that the integral store is *not* zeroed.

        use checking, only: check_allocate
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: op_sym
        logical, intent(in) :: imag
        type(one_body_t), intent(out) :: store

        integer :: ierr, i, s, ispin, nspin

        store%op_sym = op_sym
        store%uhf = sys%read_in%uhf
        store%imag = imag
        ! if rhf then need to store only integrals for spatial orbitals.
        ! ie < i,alpha j,beta | a,alpha b,beta > = < i,alpha j,alpha | a,alpha b,alpha >
        if (store%uhf) then
            nspin = 2
        else
            nspin = 1
        end if

        associate(nbasis_sym_spin=>sys%read_in%pg_sym%nbasis_sym_spin)
            ! Allocate general store for the one-electron integrals.
            allocate(store%integrals(nspin,lbound(nbasis_sym_spin, dim=2):ubound(nbasis_sym_spin, dim=2)), stat=ierr)
            call check_allocate('one_body_store', nspin*size(nbasis_sym_spin, dim=2), ierr)

            ! <i|o|j> is only non-zero if the integrand is totally symmetric, ie
            ! \Gamma_i \cross \Gamma_o \cross \Gamma_j = \Gamma_1.
            ! Currently all operators we consider are Hermitian.
            ! => Store only lower triangle of each symmetry block in o_{ij}.
            ! Within the block each state is labelled by its symmetry and spin
            ! index.  If \Gamma_o \= \Gamma_1, then i and j are of different
            ! symmetries.  We get around this by arranging index_i=<index_j.  In this
            ! case some memory is wasted (as we store the diagonal elements in both
            ! the i and j symmetry blocks) and if the number of states with the same
            ! symmetry as i is greater than those with with same symmetry as j, but
            ! this effect will be small.
            !
            ! o_{ij} is only non-zero if i and j are of the same spin.
            ! Furthermore, o_{i,alpha, j,alpha} = o_{i,beta, j, beta} in RHF
            ! calculations.
            ! => store spin blocks separately and only store both in UHF
            ! calculations.
            do ispin = 1, nspin
                do i = lbound(store%integrals, dim=2), ubound(store%integrals, dim=2)
                    s = (nbasis_sym_spin(ispin,i)*(nbasis_sym_spin(ispin,i)+1))/2
                    allocate(store%integrals(ispin,i)%v(s), stat=ierr)
                    call check_allocate('one_body_store_component', s, ierr)
                end do
            end do
        end associate

    end subroutine init_one_body_t

    subroutine end_one_body_t(store)

        ! Deallocate components of a store of integrals involving a one-body operator.

        ! In/Out:
        !    store: one-body integral store with components allocated to hold
        !    integrals which are deallocated upon exit.

        use checking, only: check_deallocate

        type(one_body_t), intent(inout) :: store
        integer :: i, ierr, ispin

        if (allocated(store%integrals)) then
            do ispin = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
                do i = lbound(store%integrals, dim=2), ubound(store%integrals, dim=2)
                    deallocate(store%integrals(ispin,i)%v, stat=ierr)
                    call check_deallocate('one_body_store_component', ierr)
                end do
            end do
            deallocate(store%integrals, stat=ierr)
            call check_deallocate('one_body_store', ierr)
        end if

    end subroutine end_one_body_t

    subroutine init_two_body_t(sys, op_sym, imag, store)

        ! Allocate memory required for the integrals involving a two-body
        ! operator.

        ! In:
        !    sys: sys_t object containing info on current system. We use:
        !           -sys%read_in%uhf
        !           -sys%basis%nbasis
        !           -sys%read_in%comp
        !    op_sym: bit string representations of irreducible representations
        !       of a point group.  See point_group_symmetry.
        !    imag: whether integral store contains imaginary component of complex
        !       integrals.
        ! Out:
        !    store: two-body integral store with components allocated to hold
        !       interals.  Note that the integral store is *not* zeroed.

        ! NB with shared memory, this will return the same chunk of memory for
        !    each given spin/imag combination, even if called a second time.
        !    i.e. you cannot have more than one type of 2-body integral store.

        use checking, only: check_allocate
        use const, only: int_64
        use parallel, only: parent
        use system, only: sys_t
        use shmem, only: allocate_shared

        logical, intent(in) :: imag
        integer, intent(in) :: op_sym
        type(sys_t), intent(in) :: sys
        type(two_body_t), intent(out) :: store

        integer :: ierr, ispin, nspin, mem_reqd, iunit
        integer(int_64):: npairs, nintgrls
        character(30) :: int_name

        iunit = 6

        store%op_sym = op_sym
        store%uhf = sys%read_in%uhf
        store%imag = imag
        store%comp = sys%read_in%comp
        ! if rhf then need to store only integrals for spatial orbitals.
        ! ie < i,alpha j,beta | a,alpha b,beta > = < i,alpha j,alpha | a,alpha b,alpha >
        if (store%uhf) then
            nspin = 4
        else
            nspin = 1
        end if

        ! Allocate general store for each spin-channel the two-electron integrals.
        allocate(store%integrals(nspin), stat=ierr)
        call check_allocate('two_body_store', nspin**2, ierr)

        ! Allocate component of store for each spin-channel.
        ! The spatial parts are identical in RHF, thus need store only one
        ! spin-channel.
        ! In UHF need to store <a a|a a>, <a b|a b>, <b a|b a> and <b b|b b>
        ! (where a==alpha spin-orbital and b==beta spin-orbital).
        ! For the integral <i j|a b>, where (i,j,a,b) are spatial-orbitals,
        ! there are M(M+1)/2=N_p (i,a) pairs (and similarly for (j,b) pairs).
        ! Due to permutation symmetry there are thus N_p(N_p+1)/2 integrals per
        ! spin-channel, where 2M is the number of spin-orbitals.
        ! NOTE:
        ! Compression due to spatial symmetry not yet implemented.
        npairs = ((sys%basis%nbasis/2)*(sys%basis%nbasis/2 + 1))/2
        ! If complex twice as many integrals so need twice the size of array
        ! as <ia|jb> != <ib|ja>
        if (store%comp) then
            nintgrls = (npairs*(npairs+1))
        else
            nintgrls = (npairs*(npairs+1))/2
        end if

        if (parent .and. (.not.store%comp .or. .not.store%imag)) then
#ifdef SINGLE_PRECISION
            mem_reqd = int((nintgrls*4*nspin)/10**6)
#else
            mem_reqd = int((nintgrls*8*nspin)/10**6)
#endif
            if (store%comp) mem_reqd = 2 * mem_reqd
            write(iunit,'(1X,a,i0)') 'Memory required for all two body integrals (MB) on each processor: ', &
                            mem_reqd
            write(iunit,'(1X, a,/)') 'It is left to the user to ensure that this does not exceed available resources.'
        end if

        do ispin = 1, nspin
            if (.not.imag) then
                write (int_name, '("two_body_store_component",i1)') ispin
            else
                write (int_name, '("two_body_store_component",i1,"imag")') ispin
            end if
            associate(int_store=>store%integrals(ispin))
                call allocate_shared(int_store%v, int_name, int_store%shmem_handle, nintgrls)
            end associate
        end do

    end subroutine init_two_body_t

    subroutine end_two_body_t(store)

        ! Deallocate comptwonts of a store of integrals involving a two-body operator.

        ! In/Out:
        !    store: two-body integral store with comptwonts allocated to hold
        !    integrals which are deallocated upon exit.

        use checking, only: check_deallocate
        use shmem, only: deallocate_shared

        type(two_body_t), intent(inout) :: store
        integer :: ierr, ispin

        if (allocated(store%integrals)) then
            do ispin = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
                associate(int_store=>store%integrals(ispin))
                    call deallocate_shared(int_store%v, int_store%shmem_handle)
                end associate
            end do
            deallocate(store%integrals, stat=ierr)
            call check_deallocate('two_body_store', ierr)
        end if

    end subroutine end_two_body_t

    subroutine init_two_body_exchange_t(sys, op_sym, store)

        ! Allocate memory required for the integrals involving a two-body
        ! operator.

        ! In:
        !    sys: sys_t object containing info on current system. We use:
        !           -sys%read_in%uhf
        !           -sys%basis%nbasis
        !    op_sym: bit string representations of irreducible representations
        !       of a point group.  See point_group_symmetry.
        ! Out:
        !    store: two-body integral store with components allocated to hold
        !       interals.  Note that the integral store is *not* zeroed.

        use checking, only: check_allocate
        use const, only: int_64
        use parallel, only: parent
        use system, only: sys_t

        integer, intent(in) :: op_sym
        type(sys_t), intent(in) :: sys
        type(two_body_exchange_t), intent(out) :: store

        integer :: ierr, ispin, nspin, mem_reqd, iunit
        integer(int_64):: npairs, nintgrls

        iunit = 6

        store%op_sym = op_sym
        store%uhf = sys%read_in%uhf
        store%comp = sys%read_in%comp
        ! if rhf then need to store only integrals for spatial orbitals.
        ! ie < i,alpha j,beta | a,alpha b,beta > = < i,alpha j,alpha | a,alpha b,alpha >
        ! Note that integral form in spinorbitals must be
        ! <ij|ai>. This requires that all spinorbitals have
        ! the same spin, so we have two spin channels in the
        ! UHF case (all up or all down) versus one in the RHF case.

        if (store%uhf) then
            nspin = 2
        else
            nspin = 1
        end if

        ! Allocate general store for each spin-channel the two-electron integrals.
        allocate(store%integrals(nspin), stat=ierr)
        call check_allocate('two_body_store', nspin**2, ierr)

        ! Allocate component of store for each spin-channel.
        ! The spatial parts are identical in RHF, thus need store only one
        ! spin-channel.
        ! In UHF need to store <a a|a a> and <b b|b b>
        ! (where a==alpha spin-orbital and b==beta spin-orbital).
        ! For the integral <i j|a i>, where (i,j,a,i) are spatial-orbitals,
        ! there are M(M+1)/2=N_p (j,a) pairs and each pair can go with a single
        ! value of i.
        ! NOTE:
        ! Compression due to symmetry not yet implemented.
        npairs = ((sys%basis%nbasis/2)*(sys%basis%nbasis/2 + 1))/2
        nintgrls = npairs * sys%basis%nbasis/2

        if (parent ) then
#ifdef SINGLE_PRECISION
            mem_reqd = int((nintgrls*4*nspin)/10**6)
#else
            mem_reqd = int((nintgrls*8*nspin)/10**6)
#endif
            write(iunit,'(1X,a,i0)') 'Memory required for all additional exchange integrals (MB) on each processor: ', &
                            mem_reqd
            write(iunit,'(1X, a,/)') 'It is left to the user to ensure that this does not exceed available resources.'
        end if

        do ispin = 1, nspin
            allocate(store%integrals(ispin)%v(sys%basis%nbasis/2,npairs), stat=ierr)
            call check_allocate('two_body_store_component', nintgrls, ierr)
        end do

    end subroutine init_two_body_exchange_t

    subroutine end_two_body_exchange_t(store)

        ! Deallocate comptwonts of a store of integrals involving a two-body operator.

        ! In/Out:
        !    store: two-body integral store with comptwonts allocated to hold
        !    integrals which are deallocated upon exit.

        use checking, only: check_deallocate

        type(two_body_exchange_t), intent(inout) :: store
        integer :: ierr, ispin

        if (allocated(store%integrals)) then
            do ispin = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
                deallocate(store%integrals(ispin)%v, stat=ierr)
                call check_deallocate('two_body_store_component', ierr)
            end do
            deallocate(store%integrals, stat=ierr)
            call check_deallocate('two_body_store', ierr)
        end if

    end subroutine end_two_body_exchange_t

!--- Zeroing ---

     pure subroutine zero_one_body_int_store(store)

        ! Zero a one-body integral store.

        ! In:
        !    store: one-body integral store with components allocated to hold
        !    interals.
        ! Out:
        !    store: one-body integral store with integral array now set to zero.

        type(one_body_t), intent(inout) :: store

        integer :: i, ispin

        do i = lbound(store%integrals, dim=2), ubound(store%integrals, dim=2)
            do ispin = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
                store%integrals(ispin,i)%v = 0.0_p
            end do
        end do

     end subroutine zero_one_body_int_store

     pure subroutine zero_two_body_int_store(store)

        ! Zero a two-body integral store.

        ! In:
        !    store: two-body integral store with components allocated to hold
        !    interals.
        ! Out:
        !    store: two-body integral store with integral array now set to zero.

        type(two_body_t), intent(inout) :: store

        integer :: ispin

        do ispin = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
            store%integrals(ispin)%v = 0.0_p
        end do

     end subroutine zero_two_body_int_store

     pure subroutine zero_two_body_exchange_int_store(store)

        ! Zero a two-body integral store.

        ! In:
        !    store: two-body integral store with components allocated to hold
        !    interals.
        ! Out:
        !    store: two-body integral store with integral array now set to zero.

        type(two_body_exchange_t), intent(inout) :: store

        integer :: ispin

        do ispin = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
            store%integrals(ispin)%v = 0.0_p
        end do

     end subroutine zero_two_body_exchange_int_store

!--- Integral access ---

! [todo] - fast and specific 'get' functions for UHF and RHF cases

! First, symmetry check functions:

    pure function check_one_body_sym(i, j, sys, op_sym) result(sym_allowed)

        ! Check if one body integral symmetry allowed, using appropriate symmetry
        ! (pg_sym vs mom_sym).

        ! In:
        !    i,j: (indices of) spin-orbitals.
        !    sys: sys_t object containing info on current system.
        !    op_sym: symmetry of operator.
        ! Out:
        !    sym_allowed: boolean value denoting whether integral can be
        !       nonzero. false if must be zero by symmetry.

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, op_sym
        logical :: sym_allowed

        ! If dealing with complex plane waves or Lz conj(sym_i) = inv(sym_i).
        ! If pg_sym conj(sym_i) = inv(sym_i) = sym_i
        ! So in general we need:
        !       conj(sym_i) = inv(op_sym x sym_j)
        !             sym_i = op_sym x sym_j
        sym_allowed = (sys%basis%basis_fns(i)%sym == &
                sys%read_in%cross_product_sym_ptr(sys%read_in, &
                op_sym, sys%basis%basis_fns(j)%sym))

    end function check_one_body_sym

    pure function check_two_body_sym(i, j, a, b, sys, op_sym) result (allowed)

        ! Function to check if a two-body integral is symmetry allowed.

        ! In:
        !    i,j,a,b: (indices of) spin-orbitals.
        !    op_sym:  The symmetry of the operator
        !    sys: information on system being studied.
        ! Out:
        !    allowed: boolean value denoting whether or not integral can
        !       be nonzero by symmetry.

        use system, only: sys_t
        use read_in_symmetry, only: cross_product_basis_read_in

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b, op_sym
        integer :: sym_ij, sym_ab
        logical :: allowed

        sym_ij = cross_product_basis_read_in(sys, i, j)
        sym_ab = cross_product_basis_read_in(sys, a, b)
        ! As dealing with complex plane waves conj(sym_i) = inv_sym(sym_i).
        ! So need:
        !       conj(sym_ij) = inv_sym(op_sym x sym_ab)
        !            sym_ij = op_sym x sym_ab
        allowed = (sym_ij == sys%read_in%cross_product_sym_ptr(sys%read_in, &
                                           sym_ab, op_sym))

    end function check_two_body_sym

! 1. < i | o_1 | j >

    subroutine store_one_body_int(i, j, intgrl, sys, suppress_err_msg, store, ierr)
        ! Store <i|o_1|j> in the appropriate slot in the one-body integral
        ! store.  The store does not have room for non-zero integrals, so it is
        ! assumed that <i|o_1|j> is non-zero by spin and spatial symmetry.
        !
        ! Checks point group symmetry before storage. For use with molecular systems.
        !
        ! In:
        !    i,j: (indices of) spin-orbitals.
        !    intgrl: <i|o_1|j>, where o_1 is a one-body operator.
        !    suppress_err_msg: if true, don't print out any error messages.
        !    sys: information of system being studied. Must have symmetry
        !       (pg or momentum as appropriate) and basis function initialised.
        ! In/out:
        !    store: one-body integral store.  On exit the <i|o_1|j> is also
        !       stored.
        ! Out:
        !    ierr: 0 if no error is encountered, 1 if integral should be non-zero
        !       by symmetry but is larger than depsilon.

        use system, only: sys_t
        use molecular_integral_types, only: one_body_t
        use errors, only: warning
        use const, only: depsilon

        integer, intent(in) :: i, j
        real(p), intent(in) :: intgrl
        type(sys_t), intent(in) :: sys
        logical, intent(in) :: suppress_err_msg
        type(one_body_t) :: store
        integer, intent(out) :: ierr
        character(255) :: error

        ierr = 0

        if (check_one_body_sym(i, j, sys, store%op_sym) .and. &
                        sys%basis%basis_fns(i)%ms == sys%basis%basis_fns(j)%ms) then
            call store_one_body_int_nonzero(i, j, intgrl, sys%basis%basis_fns, &
                                        suppress_err_msg, store)

        else if (abs(intgrl) > depsilon) then
            if (.not.suppress_err_msg) then
                write (error, '("<i|o|j> should be zero by symmetry: &
                                &<",i3,"|o|",i3,"> =",f16.10)') i, j, intgrl
                call warning('store_one_body_int', trim(error), -1)
            end if
            ierr = 1
        end if

    end subroutine store_one_body_int

    subroutine store_one_body_int_nonzero(i, j, intgrl, basis_fns, suppress_err_msg, store)

        ! Store <i|o_1|j> in the appropriate slot in the one-body integral
        ! store.  The store does not have room for non-zero integrals, so it is
        ! assumed that <i|o_1|j> is non-zero by spin and spatial symmetry.

        ! NB should be called once symmetry and spin have been checked.

        ! In:
        !    i,j: (indices of) spin-orbitals.
        !    intgrl: <i|o_1|j>, where o_1 is a one-body operator.
        !    suppress_err_msg: if true, don't print out any error messages.
        !    basis_fns: list of single-particle basis functions.
        ! In/out:
        !    store: one-body integral store.  On exit the <i|o_1|j> is also
        !       stored.
        use basis_types, only: basis_fn_t
        use errors, only: warning
        use utils, only: tri_ind
        type(basis_fn_t), intent(in) :: basis_fns(:)
        integer, intent(in) :: i, j
        real(p), intent(in) :: intgrl
        logical, intent(in) :: suppress_err_msg
        type(one_body_t), intent(inout) :: store

        integer :: ii, jj, spin
        ! Integral is (should be!) non-zero by symmetry.
        if (store%uhf) then
            spin = (basis_fns(i)%Ms+3)/2 ! Ms=-1,1 -> spin=1,2
        else
            spin = 1
        end if
        ii = basis_fns(i)%sym_spin_index
        jj = basis_fns(j)%sym_spin_index
        if (ii == jj) then
            ! See note about how operators which are no symmetric are
            ! handled in init_one_body_t.
            store%integrals(spin,basis_fns(i)%sym)%v(tri_ind(ii,jj)) = intgrl
            store%integrals(spin,basis_fns(j)%sym)%v(tri_ind(jj,ii)) = intgrl
        else if (ii > jj) then
            store%integrals(spin,basis_fns(i)%sym)%v(tri_ind(ii,jj)) = intgrl
        else if (store%imag) then
            ! j>i.  We store <i|o|j> with (i>j) so store the complex conjugate of
            ! <j|o|i> which requires a sign change only if this is the imaginary part
            ! of an integral.
            store%integrals(spin,basis_fns(j)%sym)%v(tri_ind(jj,ii)) = - intgrl
        else
            store%integrals(spin,basis_fns(j)%sym)%v(tri_ind(jj,ii)) = intgrl
        end if

    end subroutine store_one_body_int_nonzero

    pure function get_one_body_int_mol_real(store, i, j, sys) result(intgrl)

        ! In:
        !    store: one-body integral store.
        !    i,j: (indices of) spin-orbitals.
        !    sys: information on system under consideration.
        ! Returns:
        !    <i|o|j>, the corresponding one-body matrix element, where o is a
        !    one-body operator given by store.
        !
        ! NOTE:
        !    If <i|o|j> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_one_body_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use system, only: sys_t

        real(p) :: intgrl
        type(sys_t), intent(in) :: sys
        type(one_body_t), intent(in) :: store
        integer, intent(in) :: i, j

        if (check_one_body_sym(i, j, sys, store%op_sym) .and. &
                sys%basis%basis_fns(i)%ms == sys%basis%basis_fns(j)%ms) then
            intgrl = get_one_body_int_mol_nonzero(store, i, j, sys%basis%basis_fns)
        else
            intgrl = 0.0_p
        end if
    end function get_one_body_int_mol_real

    pure function get_one_body_int_mol_complex(store, im_store, i, j, sys) result(intgrl)

        ! In:
        !    store: one-body integral real component store.
        !    i,j: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        !    pg_sym: information on the symmetries of the basis functions.
        !    im_store: one-body integral imaginary component store.
        ! Returns:
        !    <i|o|j>, the corresponding one-body matrix element, where o is a
        !    one-body operator given by store.
        !
        ! NOTE:
        !    If <i|o|j> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_one_body_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use system, only: sys_t

        complex(p) :: intgrl
        real(p) :: re, im
        type(sys_t), intent(in) :: sys
        type(one_body_t), intent(in) :: store, im_store
        integer, intent(in) :: i, j

        ! As dealing with complex plane waves conj(sym_i) = inv(sym_i).
        ! So need:
        !       conj(sym_i) = inv(op_sym cross sym_j)
        !             sym_i = op_sym cross sym_j
        if (check_one_body_sym(i, j, sys, store%op_sym) .and. &
                sys%basis%basis_fns(i)%ms == sys%basis%basis_fns(j)%ms) then
            re = get_one_body_int_mol_nonzero(store, i, j, sys%basis%basis_fns)
            im = get_one_body_int_mol_nonzero(im_store, i, j, sys%basis%basis_fns)
            intgrl = cmplx(re, im, p)
        else
            intgrl = cmplx(0.0_p, 0.0_p, p)
        end if

    end function get_one_body_int_mol_complex

    pure function get_one_body_int_mol_nonzero(store, i, j, basis_fns) result(intgrl)

        ! In:
        !    store: one-body integral store.
        !    i,j: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        ! Returns:
        !    <i|o|j>, the corresponding one-body matrix element, where o is a
        !    one-body operator given by store.
        !
        ! NOTE:
        !    This assumes that <i|h|j> is known the be non-zero by spin and
        !    spatial symmetry.  If this is not true then this routine will return
        !    either an incorrect value or cause an array-bounds error.  If
        !    <i|h|j> might be zero by symmetry, get_one_body_int_mol must be called
        !    instead.
        !    It is faster to call RHF- or UHF-specific routines.

        use basis_types, only: basis_fn_t
        use utils, only: tri_ind

        real(p) :: intgrl
        type(basis_fn_t), intent(in) :: basis_fns(:)
        type(one_body_t), intent(in) :: store
        integer, intent(in) :: i, j

        integer :: ii, jj, spin

        if (store%uhf) then
            spin = (basis_fns(i)%Ms+3)/2 ! Ms=-1,1 -> spin=1,2
        else
            spin = 1
        end if
        ii = basis_fns(i)%sym_spin_index
        jj = basis_fns(j)%sym_spin_index

        if (ii >= jj) then
            intgrl = store%integrals(spin, basis_fns(i)%sym)%v(tri_ind(ii,jj))
        else if (store%imag) then
                ! j>i.  We store <i|o|j> with (i>j) so return its complex conjugate
                ! <j|o|i> which requires a sign change only if this is the imaginary part
                ! of an integral
            intgrl = - store%integrals(spin, basis_fns(j)%sym)%v(tri_ind(jj,ii))
        else
            intgrl = store%integrals(spin, basis_fns(j)%sym)%v(tri_ind(jj,ii))
        end if

    end function get_one_body_int_mol_nonzero

! 2. < i j | o_2 | a b >

    pure function two_body_int_indx(uhf, i, j, a, b, basis_fns) result(indx)

        ! In:
        !    uhf: whether integral store is from a UHF calculation (and hence is
        !        stored using spin orbital labels rather than spatial orbitals).
        !    i,j,a,b: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        ! Returns:
        !    indx: spin-channel and index of a two_body_t integral store which contains the
        !        <ij|o_2|ab> integral, assuming the integral is non-zero by spin
        !        and spatial symmetry.

        ! NOTE:
        !     This is not optimised for RHF systems, where the spin-channel is
        !     always 1.

        !   This should NOT be used for complex systems.

        use basis_types, only: basis_fn_t
        use utils, only: tri_ind
        use const, only: int_64

        type(int_indx) :: indx
        type(basis_fn_t), intent(in) :: basis_fns(:)
        logical, intent(in) :: uhf
        integer, intent(in) :: i, j, a, b

        integer(int_64) :: ia, jb
        integer :: ii, jj, aa, bb

        ! Use permutation symmetry to find unique indices corresponding to the
        ! desired integral.
        ! As orbitals and integrals are real, <ij|o_2|ab> = <ji|o_2|ba> = <ab|o_2|ij> = <ba|o_2|ji>
        !                                                 = <ib|o_2|aj> = <bi|o_2|ja> = <aj|o_2|ib> = <ja|o_2|bi>
        ! Obviously we wish to use this permutation symmetry to reduce the
        ! storage space required by a factor of 8.

        ! For UHF systems we must also keep track of the spin of the orbitals
        ! during the permutations so we know which spin channel the integral is
        ! in.

        ! Require i>=a and j>=b.
        if (i < a) then
            ii = a
            aa = i
        else
            ii = i
            aa = a
        end if
        if (j < b) then
            jj = b
            bb = j
        else
            jj = j
            bb = b
        end if
        ia = tri_ind(basis_fns(ii)%spatial_index, basis_fns(aa)%spatial_index)
        jb = tri_ind(basis_fns(jj)%spatial_index, basis_fns(bb)%spatial_index)

        ! Combine ia and jb in a unique way.
        ! This amounts to requiring (i,a) > (j,b), i.e. i>j || (i==j && a>b),
        ! for example.
        ! Hence find overall index after applying 3-fold permutation symmetry.
        ! NOTE: this test *only* looks at the spatial indices so it is not
        ! sufficient for detecting the case where (e.g.) i and j are different
        ! spin-orbitals with the same spatial index (see below).
        if (ia < jb) then
            indx%indx = tri_ind(jb, ia)
        else
            indx%indx = tri_ind(ia, jb)
        end if

        ! Find spin channel.
        if (uhf) then

            ! Due to overall index depending on spatial indices, we must check
            ! the spin indices to determine whether there's another flip in
            ! order to obtain the unique set of indices for this integral.
            ! If ia and jb are different, then the choice of ordering is trivial.
            ! If ia = jb (ie i,a and j,b both have one spin orbital from each of a pair
            ! of spatial orbitals) then we need to make an arbitrary choice as to which
            ! permutation to look up.
            if ( ia < jb .or. ( ia == jb .and. ii < jj) ) then
                aa = ii ! don't need aa and bb any more; use as scratch space
                ii = jj
                jj = aa
            end if

            if (basis_fns(ii)%ms == -1) then
                if (basis_fns(jj)%ms == -1) then
                    ! down down down down
                    indx%spin_channel = 1
                else
                    ! down up down up
                    indx%spin_channel = 3
                end if
            else
                if (basis_fns(jj)%ms == 1) then
                    ! up up up up
                    indx%spin_channel = 2
                else
                    ! up down up down
                    indx%spin_channel = 4
                end if
            end if
        else
            indx%spin_channel = 1
        end if

    end function two_body_int_indx

    pure function two_body_int_indx_complex(uhf, i, j, a, b, basis_fns) result(indx)

        ! In:
        !    uhf: whether integral store is from a UHF calculation (and hence is
        !        stored using spin orbital labels rather than spatial orbitals).
        !    i,j,a,b: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        ! Returns:
        !    indx: spin-channel and index of a two_body_t integral store which contains the
        !        <ij|o_2|ab> integral, assuming the integral is non-zero by spin
        !        and spatial symmetry.

        ! NOTE:
        !     This is not optimised for RHF systems, where the spin-channel is
        !     always 1.

        !   This MUST be used for complex integral stores.

        use basis_types, only: basis_fn_t
        use utils, only: tri_ind, tri_ind_reorder
        use const, only: int_64

        type(int_indx) :: indx
        type(basis_fn_t), intent(in) :: basis_fns(:)
        logical, intent(in) :: uhf
        integer, intent(in) :: i, j, a, b
        integer :: orbs(4), oldorbs(4)

        integer(int_64) :: ia, jb
        integer :: ii, jj, aa, bb
        logical :: conj, eswap

        integer :: maxv, scratch(2)

        ! Use permutation symmetry to find unique indices corresponding to the
        ! desired integral.
        ! As orbitals and integrals are complex,
        !       <ij|o_2|ab> = <ji|o_2|ba> = <ab|o_2|ij>* = <ba|o_2|ji>*
        !       =/= <ib|o_2|aj> = <bi|o_2|ja> = <aj|o_2|ib>* = <ja|o_2|bi>*
        ! Obviously we wish to use this permutation symmetry to reduce the
        ! storage space required by a factor of 4.

        ! For UHF systems we must also keep track of the spin of the orbitals
        ! during the permutations so we know which spin channel the integral is
        ! in.
        ! We're given integral <ij|o_2|ab>
        ! If doing complex, can't simply reorder spinorbital indexes naively; have to use spatial
        ! orbitals to ensure correct result.
        ! We require i>=a, but cannot always also guarantee j>=b. Instead can only specify i>=j,a,b.
        ! If i == a can have j >= b, but must ensure ia >= jb. This can lead to complications in
        ! cases of multiple values being equal to the maximum index. As such, we aim to ensure
        ! ia >= jb then j >= b:
        ! - For instance, in the case of <ij|ai>, where i is the largest spatial index.
        !   Conventionally we would seek i >= a (satified by definition), then seek to ensure the
        !   largest of j and a was in the second position from the left. If j > a, this would lead
        !   to ia < jb and giving us a nonsensical final index iajb. If instead we seek to maximise
        !   the value of the a position, we will obtain the correct answer.

        oldorbs = (/i, j, a, b/)

        ii = basis_fns(i)%spatial_index
        jj = basis_fns(j)%spatial_index
        aa = basis_fns(a)%spatial_index
        bb = basis_fns(b)%spatial_index

        orbs = (/ii, jj, aa, bb/)
        maxv = max(ii,jj,aa,bb)
        ! Use counting here to avoid degeneracy issues from complex notation.
        ! If only one maximum value, have unique choice of ordering. In degenerate cases, have to
        ! be more careful
        ! First figure out if we have to conjugate or swap electrons for best alignment
        if (count(orbs == maxv) == 1) then
            conj = (maxv == aa .or. maxv == bb)
            eswap = (maxv == jj .or. maxv == bb)
        else
            if (ii == jj .and. ii == maxv) then
                ! Have <ii|ab>, <ii|ai>, <ii|ib> or <ii|ii>
                conj = .false.
                eswap = (bb > aa)
            else if (aa == bb .and. aa == maxv) then
                ! Have <ij|aa>, <aj|aa> or <ia|aa>
                conj = .true.
                eswap = (jj > ii)
            else if (ii == aa .and. ii == maxv) then
                ! <ij|ib>
                conj = (bb > jj)
                eswap = .false.
            else if (ii == bb .and. ii == maxv) then
                ! <ij|ai>
                conj = (jj > aa)
                eswap = conj
            else if (jj == aa .and. jj == maxv) then
                ! <ij|jb>
                conj = (ii > bb)
                eswap = (.not.conj)
            else if (jj == bb .and. jj == maxv) then
                ! <ij|aj>
                conj = (ii < aa)
                eswap = .true.
            end if
        end if

        ! Perform operations we've decided we need.

        if (conj) then
            scratch = (/aa, ii/)
            ii = scratch(1)
            aa = scratch(2)
            scratch = (/bb, jj/)
            jj = scratch(1)
            bb = scratch(2)
        end if

        if (eswap) then
            scratch = (/ii, jj/)
            jj = scratch(1)
            ii = scratch(2)
            scratch = (/aa, bb/)
            bb = scratch(1)
            aa = scratch(2)

        end if

        ! Combine values uniquely, with reordering for jb as we can't guarantee j >= b
        ia = tri_ind(ii, aa)
        jb = tri_ind_reorder(jj, bb)

        ! Combine ia and jb in a unique way.
        ! This amounts to requiring (i,a) >= (j,b), i.e. i>=j || (i==j && a>=b),
        ! for example.
        ! Hence find overall index after applying 3-fold permutation symmetry.
        ! NOTE: this test *only* looks at the spatial indices so it is not
        ! sufficient for detecting the case where (e.g.) i and j are different
        ! spin-orbitals with the same spatial index (see below).
        ! As two possible permutations give same ia and jb values, need to
        ! ensure give different values.

        ! We observe that each previously unique index for a real integral store
        ! has a pair of values associated with it in the complex case, <ij|ab> and
        ! <ib|aj> for i >= j,a,b, j > b (see note on j == b case below).
        ! If we then double our inital index we can store the <ij|ab> value at the
        ! associated odd-value index and the <ib|aj> value at the even-value index.
        ! This gives a unique index for each integral in the complex system.

        ! In various cases of equality the assumed pair of integrals will be identical
        ! and so one of the indexes will be unused. Depending on system size this will
        ! contribute a degree of inefficiency.

        if (jj < bb) then
            indx%indx = 2 * tri_ind(ia, jb)
        else
            indx%indx = 2 * tri_ind(ia, jb) - 1
        end if

        indx%conjugate = conj

        ! Find spin channel.
        if (uhf) then

            ! Due to overall index depending on spatial indices, we must check
            ! the spin indices to determine whether there's another flip in
            ! order to obtain the unique set of indices for this integral.
            ! If ia and jb are different, then the choice of ordering is trivial.
            ! If ia == jb (ie i,a and j,b both have one spin orbital from each of a pair
            ! of spatial orbitals) then we need to make an arbitrary choice as to which
            ! permutation to look up.
            ! First apply operations we know we need to spinorbitals.
            if (conj) then
                scratch = (/oldorbs(1), oldorbs(3)/)
                oldorbs(1) = scratch(2)
                oldorbs(3) = scratch(1)
                scratch = (/ oldorbs(2), oldorbs(4) /)
                oldorbs(2) = scratch(2)
                oldorbs(4) = scratch(1)
            end if

            if (eswap) then
                scratch = (/oldorbs(1), oldorbs(2)/)
                oldorbs(1) = scratch(2)
                oldorbs(2) = scratch(1)
                scratch = (/ oldorbs(3), oldorbs(4) /)
                oldorbs(3) = scratch(2)
                oldorbs(4) = scratch(1)
            end if

            ! From previous operations, should have ensured ia >= jb & i >= j for
            ! spatial orbitals. As such, we now just have to figure out which spin
            ! channel the resultant arrangement is in. As we've already ensured a
            ! unique arrangement, we can ensure any rearrangement of same spinorbitals
            ! will give the same index and spin channel.

            ii = oldorbs(1)
            jj = oldorbs(2)


            if (basis_fns(ii)%ms == -1) then
                if (basis_fns(jj)%ms == -1) then
                    ! down down down down
                    indx%spin_channel = 1
                else
                    ! down up down up
                    indx%spin_channel = 3
                end if
            else
                if (basis_fns(jj)%ms == 1) then
                    ! up up up up
                    indx%spin_channel = 2
                else
                    ! up down up down
                    indx%spin_channel = 4
                end if
            end if
        else
            indx%spin_channel = 1
        end if

    end function two_body_int_indx_complex

    subroutine store_two_body_int(i, j, a, b, intgrl, sys, suppress_err_msg, store, ierr)

        ! Store <ij|o_2|ab> in the appropriate slot in the two-body integral store.
        ! The store does not have room for non-zero integrals, so it is assumed
        ! that <ij|o_2|ab> is non-zero by spin and spatial symmetry.
        !
        ! (Note that compression by spatial symmetry is currently not
        ! implemented.)

        ! In:
        !    i,j,a,b: (indices of) spin-orbitals.
        !    intgrl: <ij|o_2|ab>, where o_2 is a two-electron operator.  Note
        !       the integral is expressed in *PHYSICIST'S NOTATION*.
        !    sys: information on the system being studied. Requires symmetry
        !       information (pg or momentum as appropriate) and basis function
        !       info to be set.
        !    suppress_err_msg: if true, don't print out any error messages.
        ! In/out:
        !    store: two-body integral store.  On exit the <ij|o_2|ab> is also
        !       stored.
        ! Out:
        !    ierr: 0 if no error is encountered, 1 if integral should be non-zero
        !       by symmetry but is larger than depsilon.

        use system, only: sys_t
        use const, only: depsilon
        use errors, only: warning

        integer, intent(in) :: i, j, a, b
        real(p), intent(in) :: intgrl
        type(sys_t), intent(in) :: sys
        logical, intent(in) :: suppress_err_msg
        type(two_body_t), intent(inout) :: store
        integer, intent(out) :: ierr

        character(255) :: error

        ierr = 0

        ! Should integral be non-zero by symmetry?

        if (check_two_body_sym(i, j, a, b, sys, store%op_sym) .and. &
                sys%basis%basis_fns(i)%ms == sys%basis%basis_fns(a)%ms .and. &
                sys%basis%basis_fns(j)%ms == sys%basis%basis_fns(b)%ms) then
            call store_two_body_int_nonzero(i, j, a, b, intgrl, sys%basis%basis_fns, store, ierr)
        else if (abs(intgrl) > depsilon) then
            if (.not.suppress_err_msg) then
                write (error, '("<ij|o|ab> should be zero by symmetry: &
                                &<",2i3,"|",2i3,"> =",f16.10)') i, j, a, b, intgrl
                call warning('store_two_body_int_pg_sym', trim(error), -1)
            end if
            ierr = 1
        end if

    end subroutine store_two_body_int

    subroutine store_two_body_int_nonzero(i, j, a, b, intgrl, basis_fns, store, ierr)
        use basis_types, only: basis_fn_t
        integer, intent(in) :: i, j, a, b
        real(p), intent(in) :: intgrl
        type(basis_fn_t), intent(in) :: basis_fns(:)
        type(two_body_t), intent(inout) :: store
        integer, intent(out) :: ierr

        type(int_indx) :: indx

        ! Store integral
        if (store%comp) then
            indx = two_body_int_indx_complex(store%uhf, i, j, a, b, basis_fns)
            if (indx%conjugate.and.store%imag) then
                store%integrals(indx%spin_channel)%v(indx%indx) = - intgrl
            else
                store%integrals(indx%spin_channel)%v(indx%indx) = intgrl
            end if
        else
            indx = two_body_int_indx(store%uhf, i, j, a, b, basis_fns)
            store%integrals(indx%spin_channel)%v(indx%indx) = intgrl
        end if

    end subroutine store_two_body_int_nonzero

    pure function get_two_body_int_mol_real(store, i, j, a, b, sys) result(intgrl)
        ! In:
        !    store: two-body integral store.
        !    i,j,a,b: (indices of) spin-orbitals.
        !    sys: information on system under consideration.
        ! Returns:
        !    < i j | o_2 | a b >, the integral between the (i,a) co-density and
        !    the (j,b) co-density involving a two-body operator o_2 given by
        !    store.
        !
        ! NOTE:
        !    If <ij|ab> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_two_body_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use system, only: sys_t

        real(p) :: intgrl
        type(two_body_t), intent(in) :: store
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b

        if (check_two_body_sym(i, j, a, b, sys, store%op_sym) .and. &
                sys%basis%basis_fns(i)%ms == sys%basis%basis_fns(a)%ms .and. &
                sys%basis%basis_fns(j)%ms == sys%basis%basis_fns(b)%ms) then
            intgrl = get_two_body_int_mol_nonzero(store, i, j, a, b, sys%basis%basis_fns)
        else
            intgrl = 0.0_p
        end if

    end function get_two_body_int_mol_real

    pure function get_two_body_int_mol_complex(store, im_store, i, j, a, b, sys) result(intgrl)

        ! In:
        !    store: store for two-body integral real component.
        !    i,j,a,b: (indices of) spin-orbitals.
        !    im_store: store for two-body integral imaginary component.
        !    sys: information on the system. We use symmetry and basis function
        !       info.
        ! Returns:
        !    < i j | o_2 | a b >, the integral between the (i,a) co-density and
        !    the (j,b) co-density involving a two-body operator o_2 given by
        !    store.
        !
        ! NOTE:
        !    If <ij|ab> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_two_body_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use system, only: sys_t

        complex(p) :: intgrl
        real(p) :: re, im
        type(two_body_t), intent(in) :: store, im_store
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b

        if (check_two_body_sym(i, j, a, b, sys, store%op_sym) .and. &
                    sys%basis%basis_fns(j)%ms == sys%basis%basis_fns(b)%ms .and. &
                    sys%basis%basis_fns(i)%ms == sys%basis%basis_fns(a)%ms) then
            re = get_two_body_int_mol_nonzero(store, i, j, a, b, sys%basis%basis_fns)
            im = get_two_body_int_mol_nonzero(im_store, i, j, a, b, sys%basis%basis_fns)
            intgrl = cmplx(re, im, p)
        else
            intgrl = cmplx(0.0_p,0.0_p, p)
        end if

    end function get_two_body_int_mol_complex

    pure function get_two_body_int_mol_nonzero(store, i, j, a, b, basis_fns) result(intgrl)

        ! In:
        !    store: two-body integral store.
        !    i,j,a,b: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        ! Returns:
        !    < i j | o_2 | a b >, the integral between the (i,a) co-density and
        !    the (j,b) co-density involving a two-body operator o_2 given by
        !    store.
        !
        ! NOTE:
        !    This assumes that <ij|ab> is known the be non-zero by spin and
        !    spatial symmetry.  If this is not true then this routine will return
        !    either an incorrect value or cause an array-bounds error.  If
        !    <ij|ab> might be zero by symmetry, get_two_body_int_mol must be called
        !    instead.
        !    It is faster to call RHF- or UHF-specific routines.

        use basis_types, only: basis_fn_t

        real(p) :: intgrl
        type(two_body_t), intent(in) :: store
        integer, intent(in) :: i, j, a, b
        type(basis_fn_t), intent(in) :: basis_fns(:)

        type(int_indx) :: indx

        if (store%comp) then
            indx = two_body_int_indx_complex(store%uhf, i, j, a, b, basis_fns)
            if (indx%conjugate.and.store%imag) then
                intgrl = - store%integrals(indx%spin_channel)%v(indx%indx)
            else
                intgrl = store%integrals(indx%spin_channel)%v(indx%indx)
            end if
        else
            indx = two_body_int_indx(store%uhf, i, j, a, b, basis_fns)
            intgrl = store%integrals(indx%spin_channel)%v(indx%indx)
        end if

    end function get_two_body_int_mol_nonzero

    pure function get_two_body_exchange_pbc_int_real(store, i, j, a, b, sys) result(intgrl)

        ! In:
        !   store: store for two body pbc exchange integrals with modified electron
        !       interaction potential.
        !   i,j,a,b: (indicies of) spin-orbitals.
        !    sys: information on the system. We use symmetry and basis function
        !       info.
        ! Returns:
        !    < i j | o_2 | a b >_x, the integral between the (i,a) co-density and
        !    the (j,b) co-density involving a two-body operator o_2 given by
        !    store using a modified electron-electron interaction potential. This
        !    is stored separately to other integrals as the unmodified interaction
        !    is sometimes required.

        use system, only: sys_t

        type(two_body_exchange_t), intent(in) :: store
        integer, intent(in) :: i,j,a,b
        type(sys_t), intent(in) :: sys
        real(p) :: intgrl

        if (check_two_body_sym(i,j,a,b, sys, store%op_sym) .and. &
                    sys%basis%basis_fns(j)%ms == sys%basis%basis_fns(b)%ms .and. &
                    sys%basis%basis_fns(i)%ms == sys%basis%basis_fns(a)%ms) then
            intgrl = get_two_body_exchange_pbc_int_nonzero(store, i, j, a, b, sys%basis%basis_fns)
        else
            intgrl = 0.0_p
        end if

    end function get_two_body_exchange_pbc_int_real

    pure function get_two_body_exchange_pbc_int_complex(store, im_store, i, j, a, b, sys) result(intgrl)

        ! In:
        !   store: store for real component of the two body pbc exchange
        !       integrals with modified electron interaction potential.
        !   im_store: store for imag component of the two body pbc exchange
        !       integrals with modified electron interaction potential.
        !   i,j,a,b: (indicies of) spin-orbitals.
        !    sys: information on the system. We use symmetry and basis function
        !       info.
        ! Returns:
        !    < i j | o_2 | a b >_x, the integral between the (i,a) co-density and
        !    the (j,b) co-density involving a two-body operator o_2 given by
        !    store using a modified electron-electron interaction potential. This
        !    is stored separately to other integrals as the unmodified interaction
        !    is sometimes required.

        use system, only: sys_t

        type(two_body_exchange_t), intent(in) :: store, im_store
        integer, intent(in) :: i,j,a,b
        type(sys_t), intent(in) :: sys
        complex(p) :: intgrl
        real(p) :: re, im

        if (check_two_body_sym(i,j,a,b, sys, store%op_sym) .and. &
                    sys%basis%basis_fns(j)%ms == sys%basis%basis_fns(b)%ms .and. &
                    sys%basis%basis_fns(i)%ms == sys%basis%basis_fns(a)%ms) then
            re = get_two_body_exchange_pbc_int_nonzero(store, i, j, a, b, sys%basis%basis_fns)
            im = get_two_body_exchange_pbc_int_nonzero(im_store, i, j, a, b, sys%basis%basis_fns)
            intgrl = cmplx(re, im, p)
        else
            intgrl = cmplx(0.0_p, 0.0_p, p)
        end if

    end function get_two_body_exchange_pbc_int_complex

    pure function get_two_body_exchange_pbc_int_nonzero(store, i, j, a, b, basis_fns) result(intgrl)

        ! In:
        !    store: two-body integral store.
        !    i,j,a,b: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        ! Returns:
        !    < i j | o_2 | a b >, the integral between the (i,a) co-density and
        !    the (j,b) co-density involving a two-body operator o_2 given by
        !    store.
        !
        ! NOTE:
        !    This assumes that <ij|ab> is known the be non-zero by spin and
        !    spatial symmetry.  If this is not true then this routine will return
        !    either an incorrect value or cause an array-bounds error.  If
        !    <ij|ab> might be zero by symmetry, get_two_body_int_mol must be called
        !    instead.
        !    It is faster to call RHF- or UHF-specific routines.

        use basis_types, only: basis_fn_t

        type(two_body_exchange_t), intent(in) :: store
        integer, intent(in) :: i,j,a,b
        type(basis_fn_t), intent(in) :: basis_fns(:)
        real(p) :: intgrl
        type(int_ex_indx) :: indx

        indx = pbc_ex_int_indx(store%uhf, i, j, a, b, basis_fns)

        intgrl = store%integrals(indx%spin_channel)%v(indx%repeat_ind,indx%triind)

        if (store%imag .and. indx%conjugate) intgrl = -intgrl

    end function get_two_body_exchange_pbc_int_nonzero

! 3. additional < i j | o_2 | a i > for PBC.

    subroutine store_pbc_int_mol(i, j, a, b, intgrl, basis_fns, store, ierr)

        use basis_types, only: basis_fn_t
        use errors, only: warning
        use utils, only: tri_ind

        integer, intent(in) :: i,j,a,b
        integer, intent(out) :: ierr
        real(p), intent(in) :: intgrl
        real(p) :: intgrl_loc
        type(basis_fn_t), intent(in) :: basis_fns(:)
        type(two_body_exchange_t), intent(inout) :: store
        type(int_ex_indx) :: indx

        indx = pbc_ex_int_indx(store%uhf, i, j, a, b, basis_fns)


        if (indx%conjugate .and. store%imag) then
            intgrl_loc = -intgrl
        else
            intgrl_loc = intgrl
        end if

        store%integrals(indx%spin_channel)%v(indx%repeat_ind,indx%triind) = intgrl_loc

    end subroutine store_pbc_int_mol

    pure function pbc_ex_int_indx(uhf, i, j, a, b, basis_fns) result(indx)

        use basis_types, only: basis_fn_t
        use errors, only: warning
        use utils, only: tri_ind

        type(int_ex_indx) :: indx
        logical, intent(in) :: uhf
        integer, intent(in) :: i, j, a, b
        type(basis_fn_t), intent(in) :: basis_fns(:)

        logical :: conjugate, elec_swap
        integer :: ii,jj,aa,bb,scratch

        ! Index goes integrals(ispin)%v(repeated_basis, triind of remainder)
        conjugate = .false.
        elec_swap = .false.

        if (i == b) then
            if (j==a) then
                ! edge case- need to choose sensible representation
                ! as no guarentee both forms will appear in fcidump.
                ! Always choose to use largest pair for pair indexing.
                if (i < j) then
                    elec_swap = .true.
                end if
            else
                ! Conventional case; choose such that j>a.
                if (j < a) then
                    conjugate = .true.
                    elec_swap = .true.
                end if
            end if
        else if (j == a) then
            ! Can't have edge case. Just swap indexes and ensure satisfy
            ! conditions for tri ind.
            if (i < b) then
                conjugate = .true.
            else
                elec_swap = .true.
            end if
        end if

        ii = i
        jj = j
        aa = a
        bb = b

        if (conjugate) then
            scratch = ii  
            ii = aa
            aa = scratch
            scratch = jj
            jj = bb
            bb = scratch
        end if
        if (elec_swap) then
            scratch = ii
            ii = jj
            jj = ii
            scratch = aa
            aa = bb
            bb = scratch
        end if

        indx%repeat_ind = ii
        indx%triind = tri_ind(basis_fns(jj)%spatial_index, basis_fns(aa)%spatial_index)

        ! Note that integral form in spinorbitals must be
        ! <ij|ai>. This requires that all spinorbitals have
        ! the same spin, so we have two spin channels in the
        ! UHF case versus one in the RHF case.
        ! In this case spin channels are
        !  1 = down down down down
        !  2 = up up up up 

        if (uhf) then
            if (basis_fns(i)%ms == -1) then
                indx%spin_channel = 1
            else
                indx%spin_channel = 2
            end if
        else
            indx%spin_channel = 1
        end if

        indx%conjugate = conjugate

    end function pbc_ex_int_indx

!--- Parallel broadcasting ---

    subroutine broadcast_one_body_t(store, data_proc)

        ! Broadcast the integral store from data_proc to all processors.
        ! In/Out:
        !    store: one-body integral store.  On input the integrals are only
        !    stored on data_proc.  On output all processors have an identical copy of
        !    the integral store.
        ! In:
        !    data_proc: processor on which the integral store is already filled.

        use parallel

        type(one_body_t), intent(inout) :: store
        integer, intent(in) :: data_proc
#ifdef PARALLEL
        integer :: i, j, ierr

        ! Yes, I know I *could* use an MPI derived type, but coding this took 10
        ! minutes rather than several hours and the loss of elegance is minimal.
        call MPI_BCast(store%op_sym, 1, mpi_integer, data_proc, MPI_COMM_WORLD, ierr)
        do i = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
            do j = lbound(store%integrals, dim=2), ubound(store%integrals, dim=2)
                call MPI_BCast(store%integrals(i,j)%v, size(store%integrals(i,j)%v), mpi_preal, data_proc, MPI_COMM_WORLD, ierr)
            end do
        end do
#else
        integer :: i
        ! A null operation so I can use -Werror when compiling
        ! without this procedure throwing an error.
        i = data_proc
        i = size(store%integrals)
#endif

    end subroutine broadcast_one_body_t

    subroutine broadcast_two_body_t(store, data_proc, max_broadcast_chunk)

        ! Broadcast the integral store from data_proc to all processors.
        ! In/Out:
        !    store: two-body integral store.  On input the integrals are only
        !       stored on data_proc.  On output all processors have an identical
        !       copy of the integral store.
        ! In:
        !    data_proc: processor on which the integral store is already filled.
        !    max_broadcast_chunk: maximum number of integral values to broadcast together in
        !       contiguous mpi datatype.

        use parallel
        use const, only: p, dp, int_64
        use checking, only: check_allocate, check_deallocate
        use errors, only: warning

        type(two_body_t), intent(inout) :: store
        integer, intent(in) :: data_proc, max_broadcast_chunk
#ifdef PARALLEL
        integer :: iunit
        integer :: i, ierr, nblocks, nnext, mpi_preal_block, optimal_block_size, j
        integer(int_64) :: nmain, ncurr

        iunit = 6
        call MPI_BCast(store%op_sym, 1, mpi_integer, data_proc, MPI_COMM_WORLD, ierr)

        ! Broadcast integrals between **nodes** such that each processor can access the integrals.
        ! Note intra_node_comm is the MPI_COMM_WORLD internode communicator if MPI-3 shared memory is not in use
        ! (i.e. each processor is effectively its own node).  See comments in parallel for more details.
        do i = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
            ! Only use chunked broadcasting if have more elements than can broadcast in
            ! single MPI call.
            if (size(store%integrals(i)%v, kind=int_64) > max_broadcast_chunk) then
                ! Broadcasting more elements than mpi supports by default- can only take 32-bit
                ! size parameter in MPI_BCast.

                associate(ints=>store%integrals(i)%v)
                    ! Instead use custom contiguous type and calculate optimal number/size of
                    ! blocks to minimise size of remaining elements to be broadcast. Passing
                    ! an array slice into mpi_bcast can lead to formation of a temporary array
                    ! the size of the full array (depending on compiler), so we have to work
                    ! around this and minimise size of this call.
                    ! Using a constant contiguous type size block_size gives a maximum remainder size
                    ! of block_size - 1.
                    ! In comparison, calculating the optimal block size reduces this maximum to
                    ! nblocks - 1 (this is currently used).
                    call get_optimal_integral_block(size(ints, kind=int_64), max_broadcast_chunk, nblocks, &
                                                optimal_block_size, mpi_preal_block)
                    nmain = int(nblocks, kind=int_64) * int(optimal_block_size, kind=int_64)
                    nnext = int(size(ints, kind=int_64) - nmain)

                    if (parent) then
                        write(iunit,'(1X,"Integral Broadcasting",/,1X,21("-"),/)')
                        write(iunit,'(1X,"Integral array larger than max_broadcast_chunk ",i0,".",&
                                    &/,1X,"Using contiguous MPI types for broadcast.",/)') max_broadcast_chunk
                        write(iunit,'(1X,"Broadcasting coulomb integrals using ",i4," blocks of size ",&
                                    &es11.4E3,".")') nblocks, real(optimal_block_size)
                        write(iunit,'(1X,"This corresponds to ", es11.4E3," integrals in the main broadcast "&
                                    &,/,1X,"and ", es11.4E3," in the remainder.")') real(nmain), real(nnext)
                    end if

                    ncurr = 1
                    if (intra_node_comm /= MPI_COMM_NULL) then
                        do j = 1,nblocks
                           call MPI_BCast(ints(ncurr), optimal_block_size, mpi_preal, data_proc, intra_node_comm, ierr)
                           ncurr = ncurr + optimal_block_size
                        end do
                    end if

                    ! Finally broadcast the remaining values not included in previous block.
                    ! In some compilers (intel) passing array slices into mpi calls leads to
                    ! creation of a temporary array the size of the full integral list. To
                    ! avoid this we instead use assumed-size property of arrays in MPI.
                    ! If we move to the mpi_f08 binding in future assumed-shape arrays are used
                    ! so this approach will no longer work. However, the explicit interface
                    ! to MPI_BCast will enable array slicing without risk of array temporary
                    ! creation.

                    if (nnext > 0 .and. intra_node_comm /= MPI_COMM_NULL) &
                        call MPI_BCast(ints(nmain+1), nnext, mpi_preal, data_proc, intra_node_comm, ierr)
                    ! Finally tidy up mpi types.
                    call mpi_type_free(mpi_preal_block, ierr)
                    if (parent) then
                        write(iunit, '(/,1X,"Broadcasting completed.")')
                        write(iunit,'(/,1X,21("-"),/)')
                    end if
                end associate
            else if (intra_node_comm /= MPI_COMM_NULL) then
                call MPI_BCast(store%integrals(i)%v, size(store%integrals(i)%v), mpi_preal, data_proc, intra_node_comm, ierr)
            end if
        end do
        if (inter_node_comm /= MPI_COMM_NULL) call MPI_Barrier(inter_node_comm, ierr)
#else
        integer :: i
        ! A null operation so I can use -Werror when compiling
        ! without this procedure throwing an error.
        i = data_proc
        i = size(store%integrals)
#endif

    end subroutine broadcast_two_body_t

    subroutine broadcast_two_body_exchange_t(store, data_proc)

        use parallel
        use errors, only: warning
        use const, only: p, dp, int_64

        type(two_body_exchange_t), intent(inout) :: store
        integer, intent(in) :: data_proc

#ifdef PARALLEL
        integer :: i, j, ierr
        do i = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
            do j = lbound(store%integrals(i)%v,dim=1), ubound(store%integrals(i)%v,dim=1)

                call MPI_BCAST(store%integrals(i)%v(j,:), size(store%integrals(i)%v(j,:)), &
                                mpi_preal, data_proc, MPI_COMM_WORLD, ierr)
            end do
        end do
#else
        integer :: i
        ! A null operation so I can use -Werror when compiling
        ! without this procedure throwing an error.
        i = data_proc
        i = size(store%integrals)
#endif

    end subroutine broadcast_two_body_exchange_t

#ifdef PARALLEL
    subroutine get_optimal_integral_block(nints, max_broadcast_chunk, nblocks, optimal_block_size, &
                                        mpi_preal_block)

        ! For a given number of integrals and maximum block size calculate block number and
        ! size that (approximately) minimises the size of integral remainder to be broadcast.
        ! In:
        !   nints: total number of integrals to be broadcast.
        !   max_broadcast_chunk: maximum number of integral values to broadcast in a single
        !       mpi contiguous type.
        ! Out:
        !   nblocks: number of contiguous type blocks to broadcast in.
        !   optimal_block_size: optimal number of elements to broadcast in each block.
        !   mpi_preal_block: custom mpi type identifier for type of size specified in
        !       optimal_block_size.  Must be freed with MPI_Type_Free when no longer needed.

        use const, only: int_32, int_64
        use parallel
        integer(int_64), intent(in) :: nints
        integer, intent(in) :: max_broadcast_chunk
        integer, intent(out) :: nblocks, optimal_block_size, mpi_preal_block
        integer :: ierr

        ! First calculate the largest number of max_broadcast_chunk that don't exceed nints.
        nblocks = int(nints / max_broadcast_chunk)
        ! Now use more blocks than this.
        nblocks = nblocks + 1
        ! Finally calculate the block size that gives the smallest residual array size to be
        ! broadcast.
        optimal_block_size = int(nints / nblocks, kind=int_32)
        ! Initialise the mpi custom type of appropriate size for this block.
        call mpi_type_contiguous(optimal_block_size, mpi_preal, mpi_preal_block, ierr)
        call mpi_type_commit(mpi_preal_block, ierr)

    end subroutine get_optimal_integral_block
#endif

end module molecular_integrals
