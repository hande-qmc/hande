module molecular_integrals

! Module for storing and accessing integrals for molecular systems.
! These integrals are previously calculated using a quantum chemistry package
! (e.g. MOLPRO or QChem).

! TODO:
! * (compile-time option) allocate arrays using shmem.

use const, only: p
use base_types, only: alloc_rp1d

implicit none

! Indexing type for two_e_int store.
type int_indx
    integer :: spin_channel, indx
end type

! Interaction with the integral stores is best done using the store_* and get_*
! procedures provided below rather than directly accessing them.

! Store for one-body integrals, <i|h|j>, where i,j are spin basis functions and
! h is the one-electron operator.
! one_e_int(ispin, isym)%v(indx) corresponds to the <i|h|j> integral (assuming
! i,j conserve spin and spatial symmetry), where ispin and isym index the spin
! and spatial symmetry of i and j and indx is the combined (triangular) index of
! i and j within that spin and symmetry block.
! See access procedures for this in practice.
! This data structure makes it possible and relative easy to only store the
! integrals which are non-zero by symmetry (ie a small fraction of the possible
! integrals).
! Note that only one spin channel is needed (and stored) in RHF calculations.
type(alloc_rp1d), allocatable :: one_e_int(:,:)

! Store for the two-body integrals, <ij|1/r_12|ab>, where i,j,a,b are spin basis
! functions and 1/r_12 is the Coulomb operator.
! two_e_int(ispin)%v(indx) gives the integral <ij|ab>, where ispin depends upon
! the spin combination (ie all alpha, all beta, and haf alpha, half beta) and
! indx is related to i,j,a,b.  As we deal with real orbitals only, we can use
! permutation symmetry to reduce the number of integrals by a factor of 8.
! See access procedures for this in action.
! Note that only one spin channel is needed (and stored) in RHF calculations.
! TODO:
! * can compress coulomb integral store by ensuring integrand is totally
!   symmetric, as is done for the one-body integrals.
type(alloc_rp1d), allocatable :: two_e_int(:)

contains

!--- Memory allocation and deallocation ---

    subroutine init_molecular_integrals()

        ! Initialise integral stores for molecular integrals (subsequently read
        ! in from an FCIDUMP file).

        ! *Must* be called after point group symmetry is initialised.

        use basis, only: nbasis
        use point_group_symmetry, only: nbasis_sym_spin
        use system, only: uhf

        use checking, only: check_allocate

        integer :: ierr, i, s, ispin
        integer :: nspin, npairs, nintgrls

        ! if rhf then need to store only integrals for spatial orbitals.
        ! ie < i,alpha j,beta | a,alpha b,beta > = < i,alpha j,alpha | a,alpha b,alpha >
        if (uhf) then
            nspin = 2
        else
            nspin = 1
        end if

        ! Allocate general store for the one-electron integrals.
        allocate(one_e_int(nspin,lbound(nbasis_sym_spin, dim=2):ubound(nbasis_sym_spin, dim=2)), stat=ierr)
        call check_allocate('one_e_int', size(one_e_int), ierr)

        ! <i|h|j> is only non-zero if i and j are the same symmetry.
        ! Furthermore, h_{ij} is Hermitian.
        ! => Store only lower triangle of each symmetry block in h_{ij}.
        ! h_{ij} is only non-zero if i and j are of the same spin.
        ! Furthermore, h_{i,alpha, j,alpha} = h_{i,beta, j, beta} in RHF
        ! calculations.
        ! => store spin blocks separately and only store both in UHF
        ! calculations.
        do ispin = 1, nspin
            do i = lbound(one_e_int, dim=2), ubound(one_e_int, dim=2)
                s = (nbasis_sym_spin(ispin,i)*(nbasis_sym_spin(ispin,i)+1))/2
                allocate(one_e_int(ispin,i)%v(s), stat=ierr)
                call check_allocate('one_e_int_component', s, ierr)
            end do
        end do

        ! Allocate general store for the two-electron integrals.
        allocate(two_e_int(nspin**2), stat=ierr)
        call check_allocate('two_e_int', nspin**2, ierr)

        ! Allocate component of two_e_int for each spin-channel.
        ! The spatial parts are identical in RHF, thus need store only one
        ! spin-channel.
        ! In UHF need to store <a a|a a>, <a b|a b>, <b a|b a> and <b b|b b>
        ! (where a==alpha spin-orbital and b==beta spin-orbital).
        ! For the integral <i j|a b>, where (i,j,a,b) are spatial-orbitals,
        ! there are M(M+1)/2=N_p (i,a) pairs (and similarly for (j,b) pairs).
        ! Due to permutation symmetry there are thus N_p(N_p+1)/2 integrals per
        ! spin-channel, where 2M is the number of spin-orbitals.
        npairs = ((nbasis/2)*(nbasis/2 + 1))/2
        nintgrls = (npairs*(npairs+1))/2
        do ispin = 1, nspin**2
            allocate(two_e_int(ispin)%v(nintgrls), stat=ierr)
            call check_allocate('two_e_int_component', nintgrls, ierr)
        end do

    end subroutine init_molecular_integrals

    subroutine end_molecular_integrals()

        ! Deallocate arrays containing molecular integrals.

        use checking, only: check_deallocate

        integer :: i, ierr, ispin

        if (allocated(one_e_int)) then
            do ispin = lbound(one_e_int, dim=1), ubound(one_e_int, dim=1)
                do i = lbound(one_e_int, dim=2), ubound(one_e_int, dim=2)
                    deallocate(one_e_int(ispin,i)%v, stat=ierr)
                        call check_deallocate('one_e_int_component', ierr)
                end do
            end do
            deallocate(one_e_int, stat=ierr)
            call check_deallocate('one_e_int', ierr)
        end if
        if (allocated(two_e_int)) then
            deallocate(two_e_int, stat=ierr)
            call check_deallocate('two_e_int', ierr)
        end if

    end subroutine end_molecular_integrals

!--- Integral access ---

! TODO:
! fast and specific 'get' functions for UHF and RHF cases

! 1. < i | h | j >

    subroutine store_one_e_int_mol(i, j, intgrl)

        ! Store <i|h|j> in the appropriate slot in one_e_int.
        ! one_e_int does not have room for non-zero integrals, so it is assumed
        ! that <i|h|j> is non-zero by spin and spatial symmetry.
        !
        ! In:
        !    i,j: (indices of) spin-orbitals.
        !    intgrl: <i|h|j>, where h is the one-body Hamiltonian operator.

        use basis, only: basis_fns
        use system, only: uhf

        use const, only: depsilon
        use errors, only: stop_all
        use utils, only: tri_ind

        integer, intent(in) :: i, j
        real(p), intent(in) :: intgrl

        integer :: ii, jj, spin
        character(255) :: error

        if (basis_fns(i)%sym == basis_fns(j)%sym .and. basis_fns(i)%ms == basis_fns(j)%ms) then
            ! Integral is (should be!) non-zero by symmetry.
            if (uhf) then
                if (basis_fns(i)%ms > 0) then
                    spin = 1
                else
                    spin = 2
                end if
            else
                spin = 1
            end if
            ii = basis_fns(i)%sym_spin_index
            jj = basis_fns(j)%sym_spin_index
            if (ii >= jj) then
                one_e_int(spin,basis_fns(i)%sym)%v(tri_ind(ii,jj)) = intgrl
            else
                one_e_int(spin,basis_fns(i)%sym)%v(tri_ind(jj,ii)) = intgrl
            end if
        else if (abs(intgrl) > depsilon) then
            write (error, '("<i|h|j> should be non-zero by symmetry: &
                            &<",i3,"|h|",i3,"> =",f16.10)') i, j, intgrl
            call stop_all('store_one_e_int_mol', error)
        end if

    end subroutine store_one_e_int_mol

    elemental function get_one_e_int_mol(i, j) result(intgrl)

        ! In:
        !    i,j: (indices of) spin-orbitals.
        ! Returns:
        !    <i|h|j>, the corresponding one-body matrix element, where h is the
        !    one-body Hamiltonian operator.
        !
        ! NOTE:
        !    If <i|h|j> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_one_e_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use basis, only: basis_fns

        real(p) :: intgrl
        integer, intent(in) :: i, j

        if (basis_fns(i)%sym == basis_fns(j)%sym .and. basis_fns(i)%ms == basis_fns(j)%ms) then
            intgrl = get_one_e_int_mol_nonzero(i, j)
        else
            intgrl = 0.0_p
        end if

    end function get_one_e_int_mol

    elemental function get_one_e_int_mol_nonzero(i, j) result(intgrl)

        ! In:
        !    i,j: (indices of) spin-orbitals.
        ! Returns:
        !    <i|h|j>, the corresponding one-body matrix element, where h is the
        !    one-body Hamiltonian operator.
        !
        ! NOTE:
        !    This assumes that <i|h|j> is known the be non-zero by spin and
        !    spatial symmetry.  If this is not true then this routine will return
        !    either an incorrect value or cause an array-bounds error.  If
        !    <i|h|j> might be zero by symmetry, get_one_e_int_mol must be called
        !    instead.
        !    It is faster to call RHF- or UHF-specific routines.

        use basis, only: basis_fns
        use system, only: uhf

        use utils, only: tri_ind

        real(p) :: intgrl
        integer, intent(in) :: i, j

        integer :: ii, jj, spin

        if (uhf) then
            if (basis_fns(i)%ms > 0) then
                spin = 1
            else
                spin = 2
            end if
        else
            spin = 1
        end if
        ii = basis_fns(i)%sym_spin_index
        jj = basis_fns(j)%sym_spin_index

        if (ii >= jj) then
            intgrl = one_e_int(spin, basis_fns(i)%sym)%v(tri_ind(ii,jj))
        else
            intgrl = one_e_int(spin, basis_fns(i)%sym)%v(tri_ind(jj,ii))
        end if

    end function get_one_e_int_mol_nonzero

! 2. < i j | a b >

    elemental function two_e_int_indx(i, j, a, b) result(indx)

        ! In:
        !    i,j,a,b: (indices of) spin-orbitals.
        ! Returns:
        !    indx: spin-channel and index of two_e_int store which contains the
        !    <ij|ab> integral.

        use basis, only: basis_fns
        use system, only: uhf

        use utils, only: tri_ind

        type(int_indx) :: indx
        integer, intent(in) :: i, j, a, b

        integer :: ii, jj, aa, bb, tmp

        ii = i
        jj = j
        aa = a
        bb = b

        ! Use permutation symmetry to find unique integral.
        ! Require i<a and j<b.
        if (ii > aa) then
            tmp = aa
            aa = ii
            ii = tmp
        end if
        if (jj > bb) then
            tmp = bb
            bb = jj
            jj = tmp
        end if
        ! Require (i,a) < (j,b), i.e. i<j || (i==j && a<b)
        if (ii > jj .or. (ii==jj .and. aa > bb) ) then
            tmp = ii
            ii = jj
            jj = tmp
            tmp = aa
            aa = bb
            bb = tmp
        end if

        ! Find spin channel.
        if (uhf) then
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

        ! Convert to spatial indices
        ii = basis_fns(ii)%spatial_index
        jj = basis_fns(jj)%spatial_index
        aa = basis_fns(aa)%spatial_index
        bb = basis_fns(bb)%spatial_index

        ! Find index.
        indx%indx = tri_ind(tri_ind(ii,aa),tri_ind(jj,bb))

    end function two_e_int_indx

    subroutine store_two_e_int_mol(i, j, a, b, intgrl)

        ! Store <ij|1/r_12|ab> in the appropriate slot in two_e_int.
        ! two_e_int does not have room for non-zero integrals, so it is assumed
        ! that <ij|1/r_12|ab> is non-zero by spin and spatial symmetry.
        !
        ! (Note that compression by spatial symmetry is currently not
        ! implemented.)
        !
        ! In:
        !    i,j,a,b: (indices of) spin-orbitals.
        !    intgrl: <ij|1/r_12|ab>, where 1/r_12 is the two-electron Coulomb
        !    operator.  Note the integral is expressed in *PHYSICIST'S
        !    NOTATION*.

        use basis, only: basis_fns
        use point_group_symmetry, only: cross_product_pg_basis, cross_product_pg_sym, is_gamma_irrep_pg_sym

        use const, only: depsilon
        use errors, only: stop_all

        integer, intent(in) :: i, j, a, b
        real(p), intent(in) :: intgrl

        integer :: sym_ij, sym_ab
        type(int_indx) :: indx
        character(255) :: error

        ! Should integral be non-zero by symmetry?
        sym_ij = cross_product_pg_basis(i, j)
        sym_ab = cross_product_pg_basis(a, b)
        if (is_gamma_irrep_pg_sym(cross_product_pg_sym(sym_ij, sym_ab)) &
            .and. basis_fns(i)%ms == basis_fns(a)%ms &
            .and. basis_fns(j)%ms == basis_fns(b)%ms) then
            ! Store integral
            indx = two_e_int_indx(i, j, a, b)
            two_e_int(indx%spin_channel)%v(indx%indx) = intgrl
        else if (abs(intgrl) > depsilon) then
            write (error, '("<ij|ab> should be non-zero by symmetry: &
                            &<",2i3,"|",2i3,"> =",f16.10)') i, j, a, b, intgrl
        end if

        indx = two_e_int_indx(i, j, a, b)

    end subroutine store_two_e_int_mol

    elemental function get_two_e_int_mol(i, j, a, b) result(intgrl)

        ! In:
        !    i,j,a,b: (indices of) spin-orbitals.
        ! Returns:
        !    intgrl: < i j | 1/r_12 | a b >, the Coulomb integral between the
        !    (i,a) co-density and the (j,b) co-density.
        !
        ! NOTE:
        !    If <ij|ab> is known the be non-zero by spin and spatial symmetry,
        !    then it is faster to call get_one_e_int_mol_nonzero.
        !    It is also faster to call RHF- or UHF-specific routines.

        use basis, only: basis_fns
        use point_group_symmetry, only: cross_product_pg_basis, cross_product_pg_sym, is_gamma_irrep_pg_sym

        real(p) :: intgrl
        integer, intent(in) :: i, j, a, b

        integer :: sym_ij, sym_ab

        sym_ij = cross_product_pg_basis(i, j)
        sym_ab = cross_product_pg_basis(a, b)
        if (is_gamma_irrep_pg_sym(cross_product_pg_sym(sym_ij, sym_ab)) &
            .and. basis_fns(i)%ms == basis_fns(a)%ms &
            .and. basis_fns(j)%ms == basis_fns(b)%ms) then
            intgrl = get_two_e_int_mol_nonzero(i, j, a, b)
        else
            intgrl = 0.0_p
        end if

    end function get_two_e_int_mol

    elemental function get_two_e_int_mol_nonzero(i, j, a, b) result(intgrl)

        ! In:
        !    i,j,a,b: (indices of) spin-orbitals.
        ! Returns:
        !    intgrl: < i j | 1/r_12 | a b >, the Coulomb integral between the
        !    (i,a) co-density and the (j,b) co-density.
        !
        ! NOTE:
        !    This assumes that <ij|ab> is known the be non-zero by spin and
        !    spatial symmetry.  If this is not true then this routine will return
        !    either an incorrect value or cause an array-bounds error.  If
        !    <ij|ab> might be zero by symmetry, get_two_e_int_mol must be called
        !    instead.
        !    It is faster to call RHF- or UHF-specific routines.

        use basis, only: basis_fns

        real(p) :: intgrl
        integer, intent(in) :: i, j, a, b

        type(int_indx) :: indx

        indx = two_e_int_indx(i, j, a, b)
        intgrl = two_e_int(indx%spin_channel)%v(indx%indx)

    end function get_two_e_int_mol_nonzero

end module molecular_integrals
