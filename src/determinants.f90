module determinants

! Generation, inspection and manipulation of Slater determinants.

use ccmc_data, only: cluster_t
use const
use system
use basis
use parallel, only: parent

implicit none

! --- Slater determinants ---

! Bit string representation of the Slater determinant.
! f is used throughout to indicate a Slater determinant
! represented as a bit string.
! The *even* bits contain the alpha (spin up) functions.  This is in
! contrast to the list of basis functions, basis_fns, where the *odd*
! indices refer to alpha (spin up) functions.  This difference arises because
! fortran numbers bits from 0...
! If separate_strings is turned on, then the first basis_length/2 integers
! represent the alpha orbitals and the second half of the bit array the beta
! orbitals.

! --- Bit masks ---

! Bit masks to reveal the list of alpha basis functions and beta functions
! occupied in a Slater determinant, required for Hubbard model.
! If separate_strings is false, then:
!     Alpha basis functions are in the even bits.  alpha_mask = 01010101...
!     Beta basis functions are in the odd bits.    beta_mask  = 10101010...
! otherwise:
!     Alpha basis functions are stored in the first basis_length/2 integers
!     followed by the beta orbitals in the next basis_length/2 integers.
integer(i0) :: alpha_mask, beta_mask

! For the Heisenberg model, certain lattices can be split into two
! sublattices such that all the sites on one sublattice only have neighbors
! on the other sublattice. This is important when finding the staggered
! magnetisation:
!
! \hat{M} = \sum_{i}(-1)^{\zeta}\sigma_{i}^{z}
!
! Here zeta will be +1 for sites on one sublattice, and -1 for sites on the
! other sublattice. This is the standard measure of antiferromagnetism.
integer(i0), allocatable :: lattice_mask(:)

! If true the determinant bit string is formed from concatenating the strings
! for the alpha and beta orbitals rather than interleaving them.
! Note that this in general uses more memory due to padding at the end of the
! alpha and beta strings.
! Whilst some memory could be saved by having no padding between the end of the
! alpha string and the start of the beta string (ie not have the beta string
! start on a new integer) this would make certain operations (eg finding
! <D|U|D> in the real space Hubbard model) much harder and slower.
! WARNING: the vast majority of procedures assume this to be false.  It is the
! developer's responsibility to ensure required procedures can handle the case
! when it is true.
! Further, decoding the bit string does not produce an ordered list of
! orbitals but several procedures assume that the list *is* ordered.  It is
! therefore not sufficient to only check the procedures which operate directly
! upon bit strings.
logical :: separate_strings = .false.

!--- Info for FCI calculations ---

! Only used in FCI calculations, where we can be certain that we have fewer
! determinants than 2**31-1 (ie no overflow).
! Whilst it's set (and frequently overflows) in FCIQMC calculations, we never
! actually use it then.  See the estimate_hilbert_space option to obtain an
! estimate estimate (or, in real-space systems, exact to a certain precision)
! for the size of the Hilbert space for a given symmetry which avoids overflows.
integer :: tot_ndets

! --- FCIQMC info ---

! A handy type for containing a lot of information about a determinant.
! This is convenient for passing around different amounts of info when
! we need consistent interfaces.
! Not all compenents are necessarily allocated: only those needed at the time.
type det_info
    ! bit representation of determinant.
    integer(i0), pointer :: f(:)  => NULL()  ! (basis_length)
    ! List of occupied spin-orbitals.
    integer, pointer :: occ_list(:)  => NULL()  ! (nel)
    ! List of occupied alpha/beta spin-orbitals
    integer, pointer :: occ_list_alpha(:), occ_list_beta(:)
    ! List of unoccupied alpha/beta spin-orbitals
    integer, pointer :: unocc_list_alpha(:), unocc_list_beta(:)
    ! Number of unoccupied orbitals with each spin and symmetry.
    ! The first index maps to spin using (Ms+3)/2, where Ms=-1 is spin-down and
    ! Ms=1 is spin-up.
    integer, pointer :: symunocc(:,:) ! (2,sym0:sym_max)
    ! is the determinant an initiator determinant or not? (used only in
    ! i-FCIQMC).
    integer :: initiator_flag
    ! Pointer (never allocated) to corresponding elements in walker_data array.
    real(p), pointer :: data(:) => NULL()
    ! Pointer to an existing cluster_t variable.  Used *only* in CCMC and so
    ! should *not* be used in generic routines.  In particular, great care
    ! should be taken with excitation generators which are designed for both
    ! FCIQMC and CCMC.
    type(cluster_t), pointer :: cluster
end type det_info

interface operator(.detgt.)
    module procedure det_gt
end interface

contains

!--- Initialisation and finalisation of module-level variables ---

    subroutine init_determinants()

        ! Initialise determinant information: number of determinants, how
        ! they are stored as bit strings, lookup arrays for converting from
        ! integer list of orbitals to bit strings and vice versa.

        use checking, only: check_allocate
        use utils, only: binom_i
        use utils, only: get_free_unit, int_fmt
        use calc, only: doing_calc, exact_diag, lanczos_diag, doing_calc, &
                        dmqmc_calc, ras, ras1, ras3, ras1_min, truncation_level

        integer :: i, j, k, bit_pos, bit_element, ierr, site_index
        character(4) :: fmt1(5)

        tot_ndets = binom_i(nbasis, nel)

        ! See note in basis.
        if (separate_strings) then
            basis_length = 2*ceiling(real(nbasis)/(2*i0_length))
            last_basis_ind = nbasis/2 - i0_length*(basis_length/2-1) - 1
        else
            basis_length = ceiling(real(nbasis)/i0_length)
            last_basis_ind = nbasis - i0_length*(basis_length-1) - 1
        end if

        if(doing_calc(dmqmc_calc)) then
            total_basis_length = 2*basis_length
        else
            total_basis_length = basis_length
        end if

        if (parent) then
            fmt1 = int_fmt((/nel, nbasis, tot_ndets, i0_length, basis_length/), padding=1)
            if (system_type == heisenberg) then
                write (6,'(1X,a22,'//fmt1(1)//')') 'Number of alpha spins:', nel
            else
                write (6,'(1X,a20,'//fmt1(1)//')') 'Number of electrons:', nel
            end if
            write (6,'(1X,a26,'//fmt1(2)//')') 'Number of basis functions:', nbasis
            if (doing_calc(exact_diag+lanczos_diag)) &
                write (6,'(1X,a32,'//fmt1(3)//')') 'Total size of determinant space:', tot_ndets
            write (6,'(/,1X,a61,'//fmt1(4)//')') 'Bit-length of integers used to store determinant bit-strings:', i0_length
            write (6,'(1X,a57,'//fmt1(5)//',/)') 'Number of integers used to store determinant bit-strings:', basis_length
        end if

        ! Lookup arrays.
        allocate(bit_lookup(2,nbasis), stat=ierr)
        call check_allocate('bit_lookup',2*nbasis,ierr)
        allocate(basis_lookup(0:i0_end,basis_length), stat=ierr)
        call check_allocate('basis_lookup',i0_length*basis_length,ierr)
        basis_lookup = 0

        if (separate_strings) then
            do i = 1, nbasis-1, 2
                ! find position of alpha orbital
                bit_pos = mod((i+1)/2, i0_length) - 1
                if (bit_pos == -1) bit_pos = i0_end
                bit_element = ((i+1)/2+i0_end)/i0_length
                bit_lookup(:,i) = (/ bit_pos, bit_element /)
                basis_lookup(bit_pos, bit_element) = i
                ! corresponding beta orbital is in the same position in the
                ! second half of the string.
                bit_element = bit_element + basis_length/2
                bit_lookup(:,i+1) = (/ bit_pos, bit_element /)
                basis_lookup(bit_pos, bit_element) = i+1
            end do
        else
            do i = 1, nbasis
                bit_pos = mod(i, i0_length) - 1
                if (bit_pos == -1) bit_pos = i0_end
                bit_element = (i+i0_end)/i0_length
                bit_lookup(:,i) = (/ bit_pos, bit_element /)
                basis_lookup(bit_pos, bit_element) = i
            end do
        end if

        ! Alpha basis functions are in the even bits.  alpha_mask = 01010101...
        ! Beta basis functions are in the odd bits.    beta_mask  = 10101010...
        ! This is assumming separate_strings is off...
        alpha_mask = 0
        beta_mask = 0
        do i = 0, i0_end
            if (mod(i,2)==0) then
                alpha_mask = ibset(alpha_mask,i)
            else
                beta_mask = ibset(beta_mask,i)
            end if
        end do

        ! For Heisenberg systems, to include staggered fields and to calculate
        ! the staggered magnetisation, we require lattice_mask. Here we find
        ! lattice_mask for a gerenal bipartite lattice.
        if (system_type == heisenberg .and. bipartite_lattice) then
            allocate (lattice_mask(basis_length), stat=ierr)
            call check_allocate('lattice_mask',basis_length,ierr)
            lattice_mask = 0
            do k = 1,lattice_size(3)
                do j = 1,lattice_size(2)
                    do i = 1,lattice_size(1),2
                        site_index = (lattice_size(2)*lattice_size(1))*(k-1) + &
                                      lattice_size(1)*(j-1) + mod(j+k,2) + i
                        bit_pos = bit_lookup(1, site_index)
                        bit_element = bit_lookup(2, site_index)
                        lattice_mask(bit_element) = ibset(lattice_mask(bit_element), bit_pos)
                    end do
                end do
            end do
        end if

        if (all(ras > 0)) then
            allocate(ras1(basis_length), stat=ierr)
            call check_allocate('ras1', basis_length, ierr)
            allocate(ras3(basis_length), stat=ierr)
            call check_allocate('ras3', basis_length, ierr)
            ras1 = 0
            ras3 = 0
            ras1_min = nel - truncation_level
            do i = 1, 2*ras(1) ! RAS is in *spatial* orbitals.
                bit_pos = bit_lookup(1,i)
                bit_element = bit_lookup(2,i)
                ras1(bit_element) = ibset(ras1(bit_element), bit_pos)
            end do
            do i = nbasis-2*ras(2)+1, nbasis
                bit_pos = bit_lookup(1,i)
                bit_element = bit_lookup(2,i)
                ras3(bit_element) = ibset(ras3(bit_element), bit_pos)
            end do
        else
            ras = -1
        end if

    end subroutine init_determinants

    subroutine end_determinants()

        ! Clean up after determinants.

        use checking, only: check_deallocate

        use calc, only: ras1, ras3

        integer :: ierr

        deallocate(bit_lookup, stat=ierr)
        call check_deallocate('bit_lookup',ierr)
        deallocate(basis_lookup, stat=ierr)
        call check_deallocate('basis_lookup',ierr)
        if (allocated(lattice_mask)) then
            deallocate(lattice_mask, stat=ierr)
            call check_deallocate('lattice_mask',ierr)
        end if
        if (allocated(ras1)) then
            deallocate(ras1, stat=ierr)
            call check_deallocate('ras1',ierr)
            deallocate(ras3, stat=ierr)
            call check_deallocate('ras3',ierr)
        end if

    end subroutine end_determinants

!--- Initialisation and finalisation of det_info objects ---

    subroutine alloc_det_info(det_info_t)

        ! Allocate the components of a det_info variable.
        ! Out:
        !    det_info_t: det_info variable with components allocated to the
        !    appropriate sizes.

        use checking, only: check_allocate

        type(det_info), intent(inout) :: det_info_t
        integer :: ierr

        allocate(det_info_t%f(basis_length), stat=ierr)
        call check_allocate('det_info_t%f',basis_length,ierr)
        allocate(det_info_t%occ_list(nel), stat=ierr)
        call check_allocate('det_info_t%occ_list',nel,ierr)
        allocate(det_info_t%occ_list_alpha(nalpha), stat=ierr)
        call check_allocate('det_info_t%occ_list_alpha',nalpha,ierr)
        allocate(det_info_t%occ_list_beta(nbeta), stat=ierr)
        call check_allocate('det_info_t%occ_list_beta',nbeta,ierr)
        allocate(det_info_t%unocc_list_alpha(nvirt_alpha), stat=ierr)
        call check_allocate('det_info_t%unocc_list_alpha',nvirt_alpha,ierr)
        allocate(det_info_t%unocc_list_beta(nvirt_beta), stat=ierr)
        call check_allocate('det_info_t%unocc_list_beta',nvirt_beta,ierr)
        allocate(det_info_t%symunocc(2,sym0:sym_max), stat=ierr)
        call check_allocate('det_info_t%symunocc',size(det_info_t%symunocc),ierr)

    end subroutine alloc_det_info

    subroutine dealloc_det_info(det_info_t)

        ! Deallocate the components of a det_info variable.
        ! Out:
        !    det_info_t: det_info variable with all components deallocated.

        use checking, only: check_deallocate

        type(det_info), intent(inout) :: det_info_t
        integer :: ierr

        deallocate(det_info_t%f, stat=ierr)
        call check_deallocate('det_info_t%f',ierr)
        deallocate(det_info_t%occ_list, stat=ierr)
        call check_deallocate('det_info_t%occ_list',ierr)
        deallocate(det_info_t%occ_list_alpha, stat=ierr)
        call check_deallocate('det_info_t%occ_list_alpha',ierr)
        deallocate(det_info_t%occ_list_beta, stat=ierr)
        call check_deallocate('det_info_t%occ_list_beta',ierr)
        deallocate(det_info_t%unocc_list_alpha, stat=ierr)
        call check_deallocate('det_info_t%unocc_list_alpha',ierr)
        deallocate(det_info_t%unocc_list_beta, stat=ierr)
        call check_deallocate('det_info_t%unocc_list_beta',ierr)
        deallocate(det_info_t%symunocc, stat=ierr)
        call check_deallocate('det_info_t%symunocc',ierr)

    end subroutine dealloc_det_info

    subroutine set_spin_polarisation(Ms)

        ! Set the spin polarisation information stored in module-level
        ! variables:
        !    nalpha, nbeta: number of alpha, beta electrons.
        !    nvirt_alpha, nvirt_beta: number of alpha, beta virtual spin-orbitals.
        ! In:
        !    Ms: spin of determinants that are being considered.

        use errors, only: stop_all

        integer, intent(in) :: Ms

        select case(system_type)

        case(heisenberg)

            ! Spin polarization is different (see comments in system) as the
            ! Heisenberg model is a collection of spins rather than electrons.

        case default

            ! Find the number of determinants with the required spin.
            if (abs(mod(Ms,2)) /= mod(nel,2)) call stop_all('set_spin_polarisation','Required Ms not possible.')

            nbeta = (nel - Ms)/2
            nalpha = (nel + Ms)/2

            nvirt_alpha = nbasis/2 - nalpha
            nvirt_beta = nbasis/2 - nbeta

        end select

    end subroutine set_spin_polarisation

!--- Encode determinant bit strings ---

    pure subroutine encode_det(occ_list, bit_list)

        ! In:
        !    occ_list(nel): integer list of occupied orbitals in the Slater determinant.
        ! Out:
        !    bit_list(basis_length): a bit string representation of the occupied
        !        orbitals.   The first element contains the first i0_length basis
        !        functions, the second element the next i0_length and so on.  A basis
        !        function is occupied if the relevant bit is set.

        integer, intent(in) :: occ_list(nel)
        integer(i0), intent(out) :: bit_list(basis_length)
        integer :: i, orb, bit_pos, bit_element

        bit_list = 0
        do i = 1, nel
            orb = occ_list(i)
            bit_pos = bit_lookup(1,orb)
            bit_element = bit_lookup(2,orb)
            bit_list(bit_element) = ibset(bit_list(bit_element), bit_pos)
        end do

    end subroutine encode_det

!--- Decode determinant bit strings ---

    pure subroutine decode_det(f, occ_list)

        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    occ_list(nel): integer list of occupied orbitals in the Slater determinant.

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(out) :: occ_list(nel)
        integer :: i, j, iorb

        iorb = 1
        outer: do i = 1, basis_length
            do j = 0, i0_end
                if (btest(f(i), j)) then
                    occ_list(iorb) = basis_lookup(j, i)
                    if (iorb == nel) exit outer
                    iorb = iorb + 1
                end if
            end do
        end do outer

    end subroutine decode_det

    pure subroutine decode_det_occ(f, d)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info variable.  The following components are set:
        !        occ_list(nel): integer list of occupied spin-orbitals in the
        !            Slater determinant.

        integer(i0), intent(in) :: f(basis_length)
        type(det_info), intent(inout) :: d
        integer :: i, j, iocc, iunocc_a, iunocc_b

        iocc = 0
        iunocc_a = 0
        iunocc_b = 0

        do i = 1, basis_length
            do j = 0, i0_end
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    d%occ_list(iocc) = basis_lookup(j, i)
                end if
                if (iocc == nel) exit
            end do
        end do

    end subroutine decode_det_occ

    pure subroutine decode_det_occ_spinunocc(f, d)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.  The unoccupied alpha and beta
        ! orbitals are given separately, as this is convenient for FCIQMC.
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info variable.  The following components are set:
        !        occ_list(nel): integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        unocc_list_alpha(nvirt_alpha): integer list of unoccupied alpha
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_beta(nvirt_beta): integer list of unoccupied beta
        !            spin-orbitals in the Slater determinant.

        integer(i0), intent(in) :: f(basis_length)
        type(det_info), intent(inout) :: d
        integer :: i, j, iocc, iunocc_a, iunocc_b

        iocc = 0
        iunocc_a = 0
        iunocc_b = 0

        do i = 1, basis_length
            do j = 0, i0_end
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    d%occ_list(iocc) = basis_lookup(j, i)
                else
                    if (mod(j,2)==0) then
                        ! alpha state (even bit index, odd basis function index)
                        iunocc_a = iunocc_a + 1
                        d%unocc_list_alpha(iunocc_a) = basis_lookup(j, i)
                    else
                        ! beta state (odd bit index, even basis function index)
                        iunocc_b = iunocc_b + 1
                        d%unocc_list_beta(iunocc_b) = basis_lookup(j, i)
                    end if
                end if
                ! Have we covered all basis functions?
                ! This avoids examining any "padding" at the end of f.
                if (iocc+iunocc_a+iunocc_b==nbasis) exit
            end do
        end do

    end subroutine decode_det_occ_spinunocc

    pure subroutine decode_det_spinocc_spinunocc(f, d)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.
        !
        ! We return the lists for alpha and beta electrons separately.
        !
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info variable.  The following components are set:
        !        occ_list(nel): integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha(nalpha): integer list of occupied alpha
        !            spin-orbitals in the Slater determinant.
        !        occ_list_beta(nbeta): integer list of occupied beta
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_alpha(nvirt_alpha): integer list of unoccupied alpha
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_beta(nvirt_beta): integer list of unoccupied beta
        !            spin-orbitals in the Slater determinant.

        integer(i0), intent(in) :: f(basis_length)
        type(det_info), intent(inout) :: d
        integer :: i, j, iocc, iocc_a, iocc_b, iunocc_a, iunocc_b, orb

        iocc = 0
        iocc_a = 0
        iocc_b = 0
        iunocc_a = 0
        iunocc_b = 0
        orb = 0

        do i = 1, basis_length - 1
            ! Manual unrolling allows us to avoid 2 mod statements
            ! and some branching.
            do j = 0, i0_end, 2
                ! Test alpha orbital.
                orb = orb + 1
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    iocc_a = iocc_a + 1
                    d%occ_list(iocc) = orb
                    d%occ_list_alpha(iocc_a) = orb
                else
                    iunocc_a = iunocc_a + 1
                    d%unocc_list_alpha(iunocc_a) = orb
                end if
                ! Test beta orbital.
                orb = orb + 1
                if (btest(f(i), j+1)) then
                    iocc = iocc + 1
                    iocc_b = iocc_b + 1
                    d%occ_list(iocc) = orb
                    d%occ_list_beta(iocc_b) = orb
                else
                    iunocc_b = iunocc_b + 1
                    d%unocc_list_beta(iunocc_b) = orb
                end if
            end do
        end do

        ! Deal with the last element in the determinant bit array separately.
        ! Note that decoding a bit string is surprisingly slow (or, more
        ! importantly, adds up when doing billions of times).
        ! Treating the last element as a special case rather than having an if
        ! statement in the above loop results a speedup of the Hubbard k-space
        ! FCIQMC calculations of 1.5%.
        do j = 0, last_basis_ind, 2
            ! Test alpha orbital.
            orb = orb + 1
            if (btest(f(i), j)) then
                iocc = iocc + 1
                iocc_a = iocc_a + 1
                d%occ_list(iocc) = orb
                d%occ_list_alpha(iocc_a) = orb
            else
                iunocc_a = iunocc_a + 1
                d%unocc_list_alpha(iunocc_a) = orb
            end if
            ! Test beta orbital.
            orb = orb + 1
            if (btest(f(i), j+1)) then
                iocc = iocc + 1
                iocc_b = iocc_b + 1
                d%occ_list(iocc) = orb
                d%occ_list_beta(iocc_b) = orb
            else
                iunocc_b = iunocc_b + 1
                d%unocc_list_beta(iunocc_b) = orb
            end if
        end do

    end subroutine decode_det_spinocc_spinunocc

    pure subroutine decode_det_occ_symunocc(f, d)

        !0 Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info variable.  The following components are set:
        !        occ_list(nel): integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        symunocc(2, sym0:symmax): number of unoccupied orbitals of each
        !            spin/symmetry.  The same indexing scheme is used for
        !            nbasis_sym_spin.

        use point_group_symmetry, only: nbasis_sym_spin

        integer(i0), intent(in) :: f(basis_length)
        type(det_info), intent(inout) :: d
        integer :: i, j, iocc, iunocc_a, iunocc_b, orb, ims, isym

        iocc = 0
        iunocc_a = 0
        iunocc_b = 0

        d%symunocc = nbasis_sym_spin

        do i = 1, basis_length
            do j = 0, i0_end
                if (btest(f(i), j)) then
                    orb = basis_lookup(j, i)
                    ims = (basis_fns(orb)%ms+3)/2
                    isym = basis_fns(orb)%sym
                    iocc = iocc + 1
                    d%occ_list(iocc) = orb
                    d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
                end if
                if (iocc == nel) exit
            end do
        end do

    end subroutine decode_det_occ_symunocc

!--- Extract information from bit strings ---

    pure function det_momentum(occ_list) result(ksum)

        ! In:
        !    occ_list(nel): integer list of occupied orbitals in the Slater determinant.
        ! Returns:
        !    ksum: Sum of wavevectors of the occupied orbitals in the Slater determinant.

        integer :: ksum(ndim)
        integer, intent(in) :: occ_list(nel)
        integer :: i

        ksum = 0
        do i = 1, nel
            ksum = ksum + basis_fns(occ_list(i))%l
        end do

    end function det_momentum

    pure function det_spin(f) result(Ms)

        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Returns:
        !    Ms: total spin of the determinant in units of electron spin (1/2).

        use bit_utils, only: count_set_bits

        integer :: Ms
        integer(i0), intent(in) :: f(basis_length)
        integer(i0) :: a, b
        integer :: i

        Ms = 0
        do i = 1, basis_length
            ! Find bit string of all alpha orbitals.
            a = iand(f(i), alpha_mask)
            ! Find bit string of all beta orbitals.
            b = iand(f(i), beta_mask)
            Ms = Ms + count_set_bits(a) - count_set_bits(b)
        end do

    end function det_spin

    pure function spin_orb_list(orb_list) result(ms)

        ! In:
        !    orb_list: list of orbitals (e.g. determinant).
        ! Returns:
        !    Ms: total spin of the determinant in units of electron spin (1/2).

        use basis, only: basis_fns

        integer :: ms
        integer, intent(in) :: orb_list(:)

        integer :: i

        ms = 0
        do i = lbound(orb_list, dim=1), ubound(orb_list, dim=1)
            ms = ms + basis_fns(orb_list(i))%Ms
        end do

    end function spin_orb_list

!--- Manipulate/transform determinant bitstrings ---

    pure function det_invert_spin(f) result(f_inv)

        ! Applies the spin inversion operator to a determinant.
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Returns:
        !    f_inv(basis_length): bit string representation of the Slater
        !        determinant after the application of the spin inversion
        !        operator.

        integer(i0) :: f_inv(basis_length)
        integer(i0), intent(in) :: f(basis_length)

        integer(i0) :: a,b
        integer :: i

        do i = 1, basis_length
            ! Find bit string of all alpha orbitals.
            a = iand(f(i), alpha_mask)
            ! Find bit string of all beta orbitals.
            b = iand(f(i), beta_mask)
            f_inv(i) = ishft(a,1) + ishft(b,-1)
        end do

    end function det_invert_spin

!--- Comparison of determinants ---

    pure function det_gt(f1, f2) result(gt)

        ! In:
        !    f1(total_basis_length): bit string representation of the Slater
        !        determinant.
        !    f2(total_basis_length): bit string representation of the Slater
        !        determinant.
        !    (For DMQMC this bitstring contains information for both determinants)
        ! Returns:
        !    True if the first element of f1 which is not equal to the
        !    corresponding element of f2 is greater than the corresponding
        !    element in f2.

        logical :: gt
        integer(i0), intent(in) :: f1(total_basis_length), f2(total_basis_length)

        integer :: i

        gt = .false.
        do i = 1, total_basis_length
            if (f1(i) > f2(i)) then
                gt = .true.
                exit
            else if (f1(i) < f2(i)) then
                gt = .false.
                exit
            end if
        end do

    end function det_gt

    pure function det_compare(f1, f2) result(compare)

        ! In:
        !    f1(total_basis_length): bit string representation of the Slater
        !        determinant.
        !    f2(total_basis_length): bit string representation of the Slater
        !        determinant.
        !    (For DMQMC this bitstring contains information for both determinants)
        ! Returns:
        !    0 if f1 and f2 are identical;
        !    1 if the first non-identical element in f1 is smaller than the
        !    corresponding element in f2;
        !    -1 if the first non-identical element in f1 is greater than the
        !    corresponding element in f2;

        integer :: compare
        integer(i0), intent(in) :: f1(total_basis_length), f2(total_basis_length)

        integer :: i

        compare = 0
        do i = 1, total_basis_length
            if (f1(i) < f2(i)) then
                compare = 1
                exit
            else if (f1(i) > f2(i)) then
                compare = -1
                exit
            end if
        end do

    end function det_compare

!--- Output ---

    subroutine write_det(f, iunit, new_line)

        ! Write out a determinant as a list of occupied orbitals in the
        ! Slater determinant.
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        !    iunit (optional): io unit to which the output is written.
        !        Default: 6 (stdout).
        !    new_line (optional): if true, then a new line is written at
        !        the end of the list of occupied orbitals.  Default: no
        !        new line.

        use utils, only: int_fmt

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in), optional :: iunit
        logical, intent(in), optional :: new_line
        integer :: occ_list(nel), io, i
        character(4) :: fmt1

        if (present(iunit)) then
            io = iunit
        else
            io = 6
        end if

        call decode_det(f, occ_list)
        fmt1 = int_fmt(nbasis,1)

        write (io,'("|")', advance='no')
        do i = 1, nel
            write (io,'('//fmt1//')', advance='no') occ_list(i)
        end do
        write (io,'(1X,">")', advance='no')
        if (present(new_line)) then
            if (new_line) write (io,'()')
        end if

    end subroutine write_det

end module determinants
