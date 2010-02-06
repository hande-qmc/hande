module determinants

! Generation, inspection and manipulation of Slater determinants. 

use const
use system
use basis
use parallel

implicit none

type det
    ! Bit string representation of the Slater determinant.
    ! f is used throughout to indicate a Slater determinant
    ! represented as a bit string.
    ! The *even* bits contain the alpha (spin up) functions.  This is in
    ! contrast to the list of basis functions, basis_fns, where the *odd*
    ! indices refer to alpha (spin up) functions.  This difference arises because 
    ! fortran numbers bits from 0...
    integer(i0), pointer :: f(:)  => NULL()  ! (basis_length)
    ! Total spin of the determinant in units of electron spin (1/2).   
    integer(i0), pointer :: Ms => NULL()
    ! Sum of wavevectors of the occupied orbitals in the Slater determinant.
    integer, pointer :: k(:) => NULL()       ! (ndim)
end type det

! Store of determinant information.
! This will quickly become a memory issue, but for dealing with the FCI of small
! systems it is ok.  It is possible to store only determinants of the same spin
! at the same time, reducing memory requirements.
! For momentum-space systems the list of determinants is grouped by wavevector
! (up to a reciprocal lattice vector).

! Rather than creating an array of type(det), which leads to a serious
! memory overhead due to the need for pointers/allocatable arrays in the
! derived type, we instead create 3 separate arrays.  A variable of type det
! can be used to point to the appropriate information, if needed.

! Bit list of the Slater determinant.  See note for f in det type.
integer(i0), allocatable, target :: dets_list(:,:) ! (basis_length,ndets)

!  Total spin of each determinant in units of electron spin (1/2).
integer(i0), allocatable, target :: dets_Ms(:) ! (ndets)

! Sum of wavevectors of the occupied orbitals in the Slater determinant.
integer(i0), allocatable, target :: dets_ksum(:,:) ! (ndim,ndets)

! Number of determinants stored in dets.
! This is the number of determinants enumerated in enumerate_determinants with
! the desired spin.
integer :: ndets

! Total size of determinant space.
integer :: tot_ndets

! A handy type for containing the excitation information needed to connect one
! determinant to another.
type excit
    ! Excitation level.
    integer :: nexcit
    ! Orbitals which are excited from and to.
    ! Only used for single and double excitations.
    integer :: from_orb(2), to_orb(2)
    ! True if a total odd number of permutations is required to align
    ! the determinants.  Only used for single and double excitations.
    logical :: perm
end type excit

! If true then the determinant list is written to determinant_file.
logical :: write_determinants = .false.
character(255) :: determinant_file = 'DETS'
integer :: det_unit

! Bit masks to reveal the list of alpha basis functions and beta functions
! occupied in a Slater determinant.
! Alpha basis functions are in the even bits.  alpha_mask = 01010101...
! Beta basis functions are in the odd bits.    beta_mask  = 10101010...
integer(i0) :: alpha_mask, beta_mask

contains

    subroutine init_determinants()

        ! Initialise determinant information: number of determinants, how
        ! they are stored as bit strings, lookup arrays for converting from
        ! integer list of orbitals to bit strings and vice versa.

        use comb_m, only: binom
        use utils, only: get_free_unit, int_fmt

        integer :: i, bit_pos, bit_element, ierr
        character(2) :: fmt1(3)

        tot_ndets = binom(nbasis, nel)

        ! See note in basis.
        basis_length = nbasis/i0_length
        if (mod(nbasis,i0_length) /= 0) basis_length = basis_length + 1

        if (parent) then
            fmt1 = int_fmt((/nel, nbasis, tot_ndets/), padding=1)
            write (6,'(1X,a20,'//fmt1(1)//')') 'Number of electrons:', nel
            write (6,'(1X,a26,'//fmt1(2)//')') 'Number of basis functions:', nbasis
            write (6,'(1X,a26,'//fmt1(3)//',/)') 'Size of determinant space:', tot_ndets
        end if

        ! Lookup arrays.
        allocate(bit_lookup(2,nbasis), stat=ierr)
        allocate(basis_lookup(0:i0_end,basis_length), stat=ierr)
        basis_lookup = 0

        do i = 1, nbasis
            bit_pos = mod(i, i0_length) - 1
            if (bit_pos == -1) bit_pos = i0_end
            bit_element = (i+i0_end)/i0_length
            bit_lookup(:,i) = (/ bit_pos, bit_element /)
            basis_lookup(bit_pos, bit_element) = i
        end do

        ! Alpha basis functions are in the even bits.  alpha_mask = 01010101...
        ! Beta basis functions are in the odd bits.    beta_mask  = 10101010...
        alpha_mask = 0
        beta_mask = 0
        do i = 0, i0_end
            if (mod(i,2)==0) then
                alpha_mask = ibset(alpha_mask,i)
            else
                beta_mask = ibset(beta_mask,i)
            end if
        end do

        det_unit = get_free_unit()
        open(det_unit, file=determinant_file, status='unknown')

    end subroutine init_determinants

    subroutine end_determinants()

        ! Clean up after determinants.

        integer :: ierr, i

        deallocate(dets_list, stat=ierr)
        deallocate(dets_Ms, stat=ierr)
        if (allocated(dets_ksum)) deallocate(dets_ksum, stat=ierr)
        deallocate(bit_lookup, stat=ierr)
        deallocate(basis_lookup, stat=ierr)

        close(det_unit, status='keep')

    end subroutine end_determinants

    subroutine enumerate_determinants(Ms)
    
        ! Find the Slater determinants that can be formed from the
        ! basis functions.  The list of determinants is stored in the
        ! module level dets array.  For momentum-space systems the list
        ! of determinants is grouped by wavevector (up to a reciprocal lattice
        ! vector).
        ! In:
        !   Ms: spin of determinants to be found.  If not given then
        !       all determinants are enumerated.
        !       The determinant list is stored in the dets array.
        !       If Ms is not given then all possible determinants are 
        !       enumerated.

        use comb_m, only: binom, comb
        use errors, only: stop_all
        use m_mrgref, only: mrgref
        use utils, only: get_free_unit, int_fmt

        integer, intent(in), optional :: Ms

        integer :: i, j, idet, c(nel), ierr, iunit, ibasis
        integer :: nalpha,  nbeta, nalpha_combinations, nbeta_combinations
        character(2) :: fmt1
        integer(i0), allocatable :: dets_ksum_tmp(:,:)
        integer(i0), allocatable, target :: dets_list_tmp(:,:)
        integer(i0), pointer :: dets_p(:,:)
        integer, allocatable :: dets_sym(:), ranking(:)
        type(det) :: d

        if (allocated(dets_list)) deallocate(dets_list, stat=ierr)
        if (allocated(dets_Ms)) deallocate(dets_Ms, stat=ierr)
        if (allocated(dets_ksum)) deallocate(dets_ksum, stat=ierr)

        if (present(Ms)) then
            ! Find the number of determinants with the required spin.
            if (mod(Ms,2) /= mod(nel,2)) call stop_all('enumerate_dets','Required Ms not possible.')
            nbeta = (nel - Ms)/2
            nalpha = (nel + Ms)/2
            nbeta_combinations = binom(nbasis/2, nbeta)
            nalpha_combinations = binom(nbasis/2, nalpha)
            ndets = nalpha_combinations*nbeta_combinations
        else
            ! Total size of determinant space.  (This will be painful!)
            ndets = binom(nbasis, nel)
        end if

        allocate(dets_list(basis_length, ndets), stat=ierr)
        allocate(dets_Ms(ndets), stat=ierr)

        if (system_type == hub_real) then
            dets_p => dets_list
        else
            allocate(dets_list_tmp(basis_length, ndets), stat=ierr)
            allocate(dets_ksum_tmp(ndim, ndets), stat=ierr)
            allocate(dets_ksum(ndim, ndets), stat=ierr)
            allocate(dets_sym(ndets), stat=ierr)
            dets_p => dets_list_tmp
        end if

        if (present(Ms)) then
            do i = 1, nbeta_combinations
                ! comb(nbasis/2, nbeta, i) will give the sites occupied by
                ! electrons in the beta spin orbital.
                ! beta orbitals are defined to be the even numbered basis
                ! functions, hence the conversion.
                c = 2*comb(nbasis/2, nbeta, i)
                do j = 1, nalpha_combinations
                    ! alpha orbitals are defined to be the odd numbered basis functions, hence
                    ! the conversion.
                    c(nbeta+1:nel) = 2*comb(nbasis/2, nalpha, j) - 1
                    idet = (i-1)*nalpha_combinations + j
                    dets_p(:,idet) = encode_det(c)
                    dets_ksum_tmp(:,idet) = det_momentum(c)
                end do
            end do
        else
            ! For the entire list it's easy: we just loop over all determinants
            ! and comb returns each possible combination in lexicographical
            ! order.
            do i = 1, ndets
                c = comb(nbasis, nel, i)
                dets_p(:,idet) = encode_det(c)
            end do
        end if

        if (system_type /= hub_real) then
            ! Rank by wavevector.
            ! All determinants have an overall wavevector that is one of the
            ! wavevectors sampled in the Brillouin zone (up to a reciprocal
            ! lattice vector).
            do idet = 1, ndets
                do i = 1, nbasis, 2
                    if (is_reciprocal_lattice_vector(basis_fns(i)%l-dets_ksum_tmp(:,idet))) then
                        dets_sym(idet) = i
                        exit
                    end if
                end do
            end do

            allocate(ranking(ndets), stat=ierr)
            call mrgref(dets_sym, ranking)

            deallocate(dets_sym, stat=ierr)

            ! Store in dets in block format: group wavevectors together.
            do idet = 1, ndets
                dets_list(:,idet) = dets_p(:,ranking(idet))
                dets_ksum(:,idet) = dets_ksum_tmp(:,ranking(idet))
                dets_Ms(idet) = det_spin(dets_list(:,idet))
            end do

            deallocate(ranking, stat=ierr)
            deallocate(dets_list_tmp, stat=ierr)
            deallocate(dets_ksum_tmp, stat=ierr)
        else
            do idet = 1, ndets
                dets_Ms(idet) = det_spin(dets_p(:,idet))
            end do
        end if

        dets_p => NULL()

        if (write_determinants .and. parent) then
            fmt1 = int_fmt(ndets, padding=1)
            do i = 1, ndets
                write (det_unit,'('//fmt1//',4X)',advance='no') i
                call write_det(dets_list(:,i), det_unit, new_line=.true.)
            end do
        end if

    end subroutine enumerate_determinants

    function point_to_det(i) result(d)

        ! Return a variable of type dets which has components
        ! that point to the bit string, spin and wavevector (if
        ! doing a momentum space calculation) of a determinant.
        ! In:
        !    i: index of determinant.

        use system, only: system_type, hub_real

        type(det) :: d
        integer, intent(in) :: i

        d%f => dets_list(:,i)
        d%Ms => dets_Ms(i)
        if (system_type /= hub_real) d%k = dets_ksum(:,i)
    
    end function point_to_det

    pure function encode_det(occ_list) result(bit_list)

        ! In:
        !    occ_list(nel): integer list of occupied orbitals in the Slater determinant.
        ! Returns:
        !    bit_list(basis_length): a bit string representation of the occupied
        !        orbitals.   The first element contains the first i0_length basis
        !        functions, the second element the next i0_length and so on.  A basis
        !        function is occupied if the relevant bit is set.

        integer(i0) :: bit_list(basis_length)
        integer, intent(in) :: occ_list(nel)
        integer :: i, orb, bit_pos, bit_element

        bit_list = 0
        do i = 1, nel
            orb = occ_list(i)
            bit_pos = bit_lookup(1,orb)
            bit_element = bit_lookup(2,orb)
            bit_list(bit_element) = ibset(bit_list(bit_element), bit_pos)
        end do
        
    end function encode_det

    pure function decode_det(f) result(occ_list)

        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Returns:
        !    occ_list(nel): integer list of occupied orbitals in the Slater determinant.

        integer :: occ_list(nel)
        integer(i0), intent(in) :: f(basis_length)
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

    end function decode_det

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

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in), optional :: iunit
        logical, intent(in), optional :: new_line
        integer :: occ_list(nel), io, i

        if (present(iunit)) then
            io = iunit
        else
            io = 6
        end if

        occ_list = decode_det(f)

        write (io,'("(")', advance='no')
        do i = 1, nel-1
            write (io,'(i4,",")', advance='no') occ_list(i)
        end do
        write (io,'(i4,")")', advance='no') occ_list(i)
        if (present(new_line)) then
            if (new_line) write (io,'()')
        end if

    end subroutine write_det

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

    pure function get_excitation(f1,f2) result(excitation)

        ! In: 
        !    f1(basis_length): bit string representation of the Slater
        !        determinant.
        !    f2(basis_length): bit string representation of the Slater
        !        determinant.
        ! Returns:
        !    excitation: excit type containing the following information---
        !        excitation%nexcit: excitation level.
        !
        !    If the excitation is a single or double excitation then it also
        !    includes:
        ! 
        !        excitation%from_orbs(2): orbitals excited from in f1.
        !        excitation%to_orbs(2): orbitals excited to in f2.
        !        excitation%perm: true if an odd number of permutations are
        !            reqiured to align the determinants.
        !        The second element of from_orbs and to_orbs is zero for single
        !        excitations.

        use bit_utils

        type(excit) :: excitation
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer :: i, j, iexcit1, iexcit2, perm, iel1, iel2, shift
        logical :: test_f1, test_f2

        excitation = excit(0, 0, 0, .false.)

        if (any(f1/=f2)) then

            iexcit1 = 0
            iexcit2 = 0
            iel1 = 0
            iel2 = 0
            perm = 0

            ! Excitation level...
#ifdef _PGI
            ! Work round an *insane* bug in PGI where intrinsic bit operations
            ! return an integer(4) if the arguments are of a kind smaller than
            ! 4.  PGI gets it right if the kind is larger than 4, but that
            ! doesn't help us in this case...
            excitation%nexcit = sum(count_set_bits(int(ieor(f1,f2),i0)))/2
#else
            excitation%nexcit = sum(count_set_bits(ieor(f1,f2)))/2
#endif

            ! Finding the permutation to align the determinants is non-trivial.
            ! It turns out to be quite easy with bit operations.
            ! The idea is to do a "dumb" permutation where the determinants are 
            ! expressed in two sections: orbitals not involved in the excitation
            ! and those that are.  Each section is stored in ascending index
            ! order.
            ! To obtain such ordering requires (for each orbital that is
            ! involved in the excitation) a total of
            ! nel - iel - nexcit + iexcit
            ! where nel is the number of electrons, iel is the position of the 
            ! orbital within the list of occupied states in the determinant,
            ! nexcit is the total number of excitations and iexcit is the number
            ! of the "current" orbital involved in excitations.
            ! e.g. Consider (1, 2, 3, 4, 5) -> (1, 3, 5, 6, 7).
            ! (1, 2, 3, 4) goes to (1, 3, 2, 4).
            ! 2 is the first (iexcit=1) orbital found involved in the excitation
            ! and so requires 5 - 2 - 2 + 1 = 2 permutation to shift it to the
            ! first "slot" in the excitation "block" in the list of states.
            ! 4 is the second orbital found and requires 5 - 4 - 2 + 2 = 1
            ! permutation to shift it the end (last "slot" in the excitation
            ! block).
            ! Whilst the resultant number of permutations isn't necessarily the
            ! minimal number for the determinants to align, this is irrelevant
            ! as the Slater--Condon rules only care about whether the number of
            ! permutations are odd or even.
            shift = nel - excitation%nexcit 

            if (excitation%nexcit <= 2) then

                do i = 1, basis_length
                    if (f1(i) == f2(i)) cycle
                    do j = 0, i0_end

                        test_f1 = btest(f1(i),j)
                        test_f2 = btest(f2(i),j)

                        if (test_f2) iel2 = iel2 + 1

                        if (test_f1) then
                            iel1 = iel1 + 1
                            if (.not.test_f2) then
                                ! occupied in f1 but not in f2
                                iexcit1 = iexcit1 + 1
                                excitation%from_orb(iexcit1) = basis_lookup(j,i)
                                perm = perm + (shift - iel1 + iexcit1)
                            end if
                        else
                            if (test_f2) then
                                ! occupied in f1 but not in f2
                                iexcit2 = iexcit2 + 1
                                excitation%to_orb(iexcit2) = basis_lookup(j,i)
                                perm = perm + (shift - iel2 + iexcit2)
                            end if
                        end if

                    end do
                end do

                ! It seems that this test is faster than btest(perm,0)!
                excitation%perm = mod(perm,2) == 1

            end if
        end if

    end function get_excitation

end module determinants
