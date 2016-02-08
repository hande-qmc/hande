module molecular_integrals

! Module for storing and accessing integrals for molecular systems.
! These integrals are previously calculated using a quantum chemistry package
! (e.g. PSI4, MOLPRO or QChem).

! [todo] - (compile-time option) allocate arrays using shmem.

use const, only: p
use molecular_integral_types

implicit none

! Indexing type for two_body_t integral stores.
! i.e. an encoding of the 4-index quartet a,b,c,d into an index for the integral store
type, private :: int_indx
    integer :: spin_channel     !If alpha and beta spin-orbitals differ, we store
                                ! different combinations of these in different spin channels
                                ! 1=bbbb, 2=aaaa, 3=baba, 4=ababa
    integer :: indx             ! The index within the spin channel
    logical :: conjugate        ! Indicates if we need to take the complex conjugate of the integral
end type int_indx

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

    subroutine init_one_body_t(uhf, op_sym, nbasis_sym_spin, imag, store)

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

        logical, intent(in) :: uhf
        integer, intent(in) :: op_sym
        integer, allocatable, intent(in) :: nbasis_sym_spin(:,:)
        logical, intent(in) :: imag
        type(one_body_t), intent(out) :: store

        integer :: ierr, i, s, ispin, nspin

        store%op_sym = op_sym
        store%uhf = uhf
        store%imag = imag
        ! if rhf then need to store only integrals for spatial orbitals.
        ! ie < i,alpha j,beta | a,alpha b,beta > = < i,alpha j,alpha | a,alpha b,alpha >
        if (uhf) then
            nspin = 2
        else
            nspin = 1
        end if

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

    subroutine init_two_body_t(uhf, nbasis, op_sym, comp, imag, store)

        ! Allocate memory required for the integrals involving a two-body
        ! operator.

        ! In:
        !    uhf: whether integral store is from a UHF calculation or RHF
        !       calculation.
        !    nbasis: number of single-particle spin functions.
        !    op_sym: bit string representations of irreducible representations
        !    of a point group.  See point_group_symmetry.
        !    comp: whether integral store is from a calculation using complex
        !       orbitals and integrals. If true, the store will contain either 
        !       the real or imaginary components of the complex two body integrals.
        !    imag: whether integral store contains imaginary component of complex
        !       integrals.
        ! Out:
        !    store: two-body integral store with components allocated to hold
        !    interals.  Note that the integral store is *not* zeroed.

        use checking, only: check_allocate

        logical, intent(in) :: uhf, comp, imag
        integer, intent(in) :: nbasis, op_sym
        type(two_body_t), intent(out) :: store

        integer :: ierr, ispin
        integer :: nspin, npairs, nintgrls

        store%op_sym = op_sym
        store%uhf = uhf
        store%imag = imag
        store%comp = comp
        ! if rhf then need to store only integrals for spatial orbitals.
        ! ie < i,alpha j,beta | a,alpha b,beta > = < i,alpha j,alpha | a,alpha b,alpha >
        if (uhf) then
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
        npairs = ((nbasis/2)*(nbasis/2 + 1))/2
        ! If complex twice as many integrals so need twice the size of array
        ! as <ia|jb> != <ib|ja>
        if (comp) then
            nintgrls = (npairs*(npairs+1))
        else
            nintgrls = (npairs*(npairs+1))/2
        end if
        do ispin = 1, nspin
            allocate(store%integrals(ispin)%v(nintgrls), stat=ierr)
            call check_allocate('two_body_store_component', nintgrls, ierr)
        end do

    end subroutine init_two_body_t

    subroutine end_two_body_t(store)

        ! Deallocate comptwonts of a store of integrals involving a two-body operator.

        ! In/Out:
        !    store: two-body integral store with comptwonts allocated to hold
        !    integrals which are deallocated upon exit.

        use checking, only: check_deallocate

        type(two_body_t), intent(inout) :: store
        integer :: ierr, ispin

        if (allocated(store%integrals)) then
            do ispin = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
                deallocate(store%integrals(ispin)%v, stat=ierr)
                call check_deallocate('two_body_store_component', ierr)
            end do
            deallocate(store%integrals, stat=ierr)
            call check_deallocate('two_body_store', ierr)
        end if

    end subroutine end_two_body_t

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

!--- Integral access ---

! [todo] - fast and specific 'get' functions for UHF and RHF cases

! 1. < i | o_1 | j >

    subroutine store_one_body_int_mol(i, j, intgrl, basis_fns, pg_sym, suppress_err_msg, store, ierr)

        ! Store <i|o_1|j> in the appropriate slot in the one-body integral
        ! store.  The store does not have room for non-zero integrals, so it is
        ! assumed that <i|o_1|j> is non-zero by spin and spatial symmetry.
        !
        ! In:
        !    i,j: (indices of) spin-orbitals.
        !    intgrl: <i|o_1|j>, where o_1 is a one-body operator.
        !    suppress_err_msg: if true, don't print out any error messages.
        !    basis_fns: list of single-particle basis functions.
        !    pg_sym: information on the symmetries of the basis functions.
        ! In/out:
        !    store: one-body integral store.  On exit the <i|o_1|j> is also
        !       stored.
        ! Out:
        !    ierr: 0 if no error is encountered, 1 if integral should be non-zero
        !       by symmetry but is larger than depsilon.

        use basis_types, only: basis_fn_t
        use point_group_symmetry, only: cross_product_pg_sym, is_gamma_irrep_pg_sym, pg_sym_conj
        use symmetry_types, only: pg_sym_t

        use const, only: depsilon
        use errors, only: warning
        use utils, only: tri_ind

        type(basis_fn_t), intent(in) :: basis_fns(:)
        type(pg_sym_t), intent(in) :: pg_sym
        integer, intent(in) :: i, j
        real(p), intent(in) :: intgrl
        logical, intent(in) :: suppress_err_msg
        type(one_body_t) :: store
        integer, intent(out) :: ierr

        integer :: ii, jj, spin
        integer :: sym
        character(255) :: error

        ierr = 0

        sym = cross_product_pg_sym(pg_sym, pg_sym_conj(pg_sym, basis_fns(i)%sym), basis_fns(j)%sym)
        sym = cross_product_pg_sym(pg_sym, sym, store%op_sym)

        ! Currently inefficient for complex as will check symmetry when storing real and 
        ! imaginary components of same integral.

        if (is_gamma_irrep_pg_sym(pg_sym, sym) .and. basis_fns(i)%ms == basis_fns(j)%ms) then

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
        else if (abs(intgrl) > depsilon) then
            if (.not.suppress_err_msg) then
                write (error, '("<i|o|j> should be zero by symmetry: &
                                &<",i3,"|o|",i3,"> =",f16.10)') i, j, intgrl
                call warning('store_one_body_int_mol', trim(error), -1)
            end if
            ierr = 1
        end if

    end subroutine store_one_body_int_mol

    pure function get_one_body_int_mol_real(store, i, j, basis_fns, pg_sym) result(intgrl)

        ! In:
        !    store: one-body integral store.
        !    i,j: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        !    pg_sym: information on the symmetries of the basis functions.
        ! Returns:
        !    <i|o|j>, the corresponding one-body matrix element, where o is a
        !    one-body operator given by store.
        !
        ! NOTE:
        !    If <i|o|j> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_one_body_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use basis_types, only: basis_fn_t
        use point_group_symmetry, only: pg_sym_conj, cross_product_pg_sym, is_gamma_irrep_pg_sym
        use symmetry_types, only: pg_sym_t

        real(p) :: intgrl
        type(basis_fn_t), intent(in) :: basis_fns(:)
        type(pg_sym_t), intent(in) :: pg_sym
        type(one_body_t), intent(in) :: store
        integer, intent(in) :: i, j

        integer :: sym

        sym = cross_product_pg_sym(pg_sym, pg_sym_conj(pg_sym, basis_fns(i)%sym),basis_fns(j)%sym)
        sym = cross_product_pg_sym(pg_sym, sym, store%op_sym)

        if (is_gamma_irrep_pg_sym(pg_sym, sym) .and. basis_fns(i)%ms == basis_fns(j)%ms) then
            intgrl = get_one_body_int_mol_nonzero(store, i, j, basis_fns)
        else
            intgrl = 0.0_p
        end if

    end function get_one_body_int_mol_real

    pure function get_one_body_int_mol_complex(store, i, j, basis_fns, pg_sym, im_store) result(intgrl)

        ! In:
        !    store: one-body integral real component store.
        !    i,j: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        !    pg_sym: information on the symmetries of the basis functions.
        !    im_store (optional): one-body integral imaginary component store.
        ! Returns:
        !    <i|o|j>, the corresponding one-body matrix element, where o is a
        !    one-body operator given by store.
        !
        ! NOTE:
        !    If <i|o|j> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_one_body_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use basis_types, only: basis_fn_t
        use point_group_symmetry, only: pg_sym_conj, cross_product_pg_sym, is_gamma_irrep_pg_sym
        use symmetry_types, only: pg_sym_t

        complex(p) :: intgrl
        real(p) :: re, im
        type(basis_fn_t), intent(in) :: basis_fns(:)
        type(pg_sym_t), intent(in) :: pg_sym
        type(one_body_t), intent(in) :: store, im_store
        integer, intent(in) :: i, j

        integer :: sym

        sym = cross_product_pg_sym(pg_sym, pg_sym_conj(pg_sym, basis_fns(i)%sym),basis_fns(j)%sym)
        sym = cross_product_pg_sym(pg_sym, sym, store%op_sym)

        if (is_gamma_irrep_pg_sym(pg_sym, sym) .and. basis_fns(i)%ms == basis_fns(j)%ms) then
            re = get_one_body_int_mol_nonzero(store, i, j, basis_fns)
            im = get_one_body_int_mol_nonzero(im_store, i, j, basis_fns)
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

        type(int_indx) :: indx
        type(basis_fn_t), intent(in) :: basis_fns(:)
        logical, intent(in) :: uhf
        integer, intent(in) :: i, j, a, b

        integer :: ia, jb, ii, jj, aa, bb

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

        ! Comine ia and jb in a unique way.
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

        type(int_indx) :: indx
        type(basis_fn_t), intent(in) :: basis_fns(:)
        logical, intent(in) :: uhf
        integer, intent(in) :: i, j, a, b
        integer :: orbs(4), oldorbs(4)

        integer :: ia, jb, ii, jj, aa, bb
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
        ! in. This is WIP currently- due to the use of spatial orbitals for 
        ! ordering in complex we will have to change implementation slightly.

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

    subroutine store_two_body_int_mol(i, j, a, b, intgrl, basis_fns, pg_sym, suppress_err_msg, store, ierr)

        ! Store <ij|o_2|ab> in the appropriate slot in the two-body integral store.
        ! The store does not have room for non-zero integrals, so it is assumed
        ! that <ij|o_2|ab> is non-zero by spin and spatial symmetry.
        !
        ! (Note that compression by spatial symmetry is currently not
        ! implemented.)
        !
        ! In:
        !    i,j,a,b: (indices of) spin-orbitals.
        !    intgrl: <ij|o_2|ab>, where o_2 is a two-electron operator.  Note
        !       the integral is expressed in *PHYSICIST'S NOTATION*.
        !    basis_fns: list of single-particle basis functions.
        !    pg_sym: information on the symmetries of the basis functions.
        !    suppress_err_msg: if true, don't print out any error messages.
        ! In/out:
        !    store: two-body integral store.  On exit the <ij|o_2|ab> is also
        !       stored.
        ! Out:
        !    ierr: 0 if no error is encountered, 1 if integral should be non-zero
        !       by symmetry but is larger than depsilon.

        use basis_types, only: basis_fn_t
        use point_group_symmetry, only: cross_product_pg_basis, cross_product_pg_sym, is_gamma_irrep_pg_sym, pg_sym_conj
        use symmetry_types, only: pg_sym_t

        use const, only: depsilon
        use errors, only: warning

        integer, intent(in) :: i, j, a, b
        real(p), intent(in) :: intgrl
        type(basis_fn_t), intent(in) :: basis_fns(:)
        type(pg_sym_t), intent(in) :: pg_sym
        logical, intent(in) :: suppress_err_msg
        type(two_body_t), intent(inout) :: store
        integer, intent(out) :: ierr

        integer :: sym_ij, sym_ab, sym
        type(int_indx) :: indx
        character(255) :: error

        ierr = 0

        ! Should integral be non-zero by symmetry?
        sym_ij = cross_product_pg_basis(pg_sym, i, j, basis_fns)
        sym_ab = cross_product_pg_basis(pg_sym, a, b, basis_fns)
        sym = cross_product_pg_sym(pg_sym, pg_sym_conj(pg_sym, sym_ij), sym_ab)
        sym = cross_product_pg_sym(pg_sym, sym, store%op_sym)

        if (is_gamma_irrep_pg_sym(pg_sym, sym) .and. basis_fns(i)%ms == basis_fns(a)%ms &
                                       .and. basis_fns(j)%ms == basis_fns(b)%ms) then
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
        else if (abs(intgrl) > depsilon) then
            if (.not.suppress_err_msg) then
                write (error, '("<ij|o|ab> should be zero by symmetry: &
                                &<",2i3,"|",2i3,"> =",f16.10)') i, j, a, b, intgrl
                call warning('store_two_body_int_mol', trim(error), -1)
            end if
            ierr = 1
        end if

    end subroutine store_two_body_int_mol

    pure function get_two_body_int_mol_real(store, i, j, a, b, basis_fns, pg_sym) result(intgrl)

        ! In:
        !    store: two-body integral store.
        !    i,j,a,b: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        !    pg_sym: information on the symmetries of the basis functions.
        ! Returns:
        !    < i j | o_2 | a b >, the integral between the (i,a) co-density and
        !    the (j,b) co-density involving a two-body operator o_2 given by
        !    store.
        !
        ! NOTE:
        !    If <ij|ab> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_two_body_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use basis_types, only: basis_fn_t
        use point_group_symmetry, only: cross_product_pg_basis, cross_product_pg_sym, is_gamma_irrep_pg_sym, pg_sym_conj
        use symmetry_types, only: pg_sym_t

        real(p) :: intgrl
        type(two_body_t), intent(in) :: store
        type(basis_fn_t), intent(in) :: basis_fns(:)
        type(pg_sym_t), intent(in) :: pg_sym
        integer, intent(in) :: i, j, a, b

        integer :: sym_ij, sym_ab, sym

        sym_ij = pg_sym_conj(pg_sym, cross_product_pg_basis(pg_sym, i, j, basis_fns))
        sym_ab = cross_product_pg_basis(pg_sym, a, b, basis_fns)
        sym = cross_product_pg_sym(pg_sym, sym_ij, sym_ab)
        sym = cross_product_pg_sym(pg_sym, sym, store%op_sym)
        if (is_gamma_irrep_pg_sym(pg_sym, sym) .and. basis_fns(i)%ms == basis_fns(a)%ms &
                                       .and. basis_fns(j)%ms == basis_fns(b)%ms) then
            intgrl = get_two_body_int_mol_nonzero(store, i, j, a, b, basis_fns)
        else
            intgrl = 0.0_p
        end if

    end function get_two_body_int_mol_real

    pure function get_two_body_int_mol_complex(store, i, j, a, b, basis_fns, pg_sym, im_store) result(intgrl)

        ! In:
        !    store: store for two-body integral real component.
        !    i,j,a,b: (indices of) spin-orbitals.
        !    basis_fns: list of single-particle basis functions.
        !    pg_sym: information on the symmetries of the basis functions.
        !    im_store: store for two-body integral imaginary component.
        ! Returns:
        !    < i j | o_2 | a b >, the integral between the (i,a) co-density and
        !    the (j,b) co-density involving a two-body operator o_2 given by
        !    store.
        !
        ! NOTE:
        !    If <ij|ab> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_two_body_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use basis_types, only: basis_fn_t
        use point_group_symmetry, only: cross_product_pg_basis, cross_product_pg_sym, is_gamma_irrep_pg_sym, pg_sym_conj
        use symmetry_types, only: pg_sym_t

        complex(p) :: intgrl
        real(p) :: re, im
        type(two_body_t), intent(in) :: store, im_store
        type(basis_fn_t), intent(in) :: basis_fns(:)
        type(pg_sym_t), intent(in) :: pg_sym
        integer, intent(in) :: i, j, a, b

        integer :: sym_ij, sym_ab, sym

        sym_ij = pg_sym_conj(pg_sym, cross_product_pg_basis(pg_sym, i, j, basis_fns))
        sym_ab = cross_product_pg_basis(pg_sym, a, b, basis_fns)
        sym = cross_product_pg_sym(pg_sym, sym_ij, sym_ab)
        sym = cross_product_pg_sym(pg_sym, sym, store%op_sym)
        if (is_gamma_irrep_pg_sym(pg_sym, sym) .and. basis_fns(i)%ms == basis_fns(a)%ms &
                                       .and. basis_fns(j)%ms == basis_fns(b)%ms) then
            re = get_two_body_int_mol_nonzero(store, i, j, a, b, basis_fns)
            im = get_two_body_int_mol_nonzero(im_store, i, j, a, b, basis_fns)
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

    subroutine broadcast_two_body_t(store, data_proc)

        ! Broadcast the integral store from data_proc to all processors.
        ! In/Out:
        !    store: two-body integral store.  On input the integrals are only
        !    stored on data_proc.  On output all processors have an identical copy of
        !    the integral store.
        ! In:
        !    data_proc: processor on which the integral store is already filled.

        use parallel

        type(two_body_t), intent(inout) :: store
        integer, intent(in) :: data_proc
#ifdef PARALLEL
        integer :: i, ierr
        ! Yes, I know I *could* use an MPI derived type, but coding this took 10
        ! minutes rather than several hours and the loss of elegance is minimal.
        call MPI_BCast(store%op_sym, 1, mpi_integer, data_proc, MPI_COMM_WORLD, ierr)
        do i = lbound(store%integrals, dim=1), ubound(store%integrals, dim=1)
            call MPI_BCast(store%integrals(i)%v, size(store%integrals(i)%v), mpi_preal, data_proc, MPI_COMM_WORLD, ierr)
        end do
#else
        integer :: i
        ! A null operation so I can use -Werror when compiling
        ! without this procedure throwing an error.
        i = data_proc
        i = size(store%integrals)
#endif

    end subroutine broadcast_two_body_t

end module molecular_integrals
