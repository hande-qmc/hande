module moulecular_integrals

!

use const, only: p
use base_types, only: alloc_rp1d

implicit none

type(alloc_rp1d), allocatable :: one_e_int(:,:)

real(p), allocatable :: two_e_int(:)

contains

!--- Memory allocation and deallocation ---

    subroutine init_molecular_integrals()

        ! Initialise integral stores for molecular integrals (subsequently read
        ! in from an FCIDUMP file).

        ! *Must* be called after point group symmetry is initialised.

        use point_group_symmetry, only: nbasis_sym_spin
        use read_in, only: uhf

        use checking, only: check_allocate

        integer :: ierr, i, s, ispin
        integer :: nspin

        ! TODO:
        ! * can compress coulomb integral store by ensuring integrand is totally
        !   symmetric.

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

        ! TODO:
        ! two-electron integral store.

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

    subroutine store_one_e_int_mol(i, j, intgrl)

        ! Store <i|h|j> in the appropriate slot in one_e_int.
        ! one_e_int does not have room for non-zero integrals, so it is assumed
        ! that <i|h|j> is non-zero by spin and spatial symmetry.
        !
        ! In:
        !    i,j: (indices of) spin-orbitals.
        !    intgrl: <i|h|j>, where h is the one-body Hamiltonian operator.

        use basis, only: basis_fns
        use read_in, only: uhf

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
        use read_in, only: uhf

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

end module moulecular_integrals
