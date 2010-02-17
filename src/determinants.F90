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
    ! Refers to the wavevector of the i-th alpha spin-orbitals.
    integer, pointer :: ksum => NULL() 
end type det

! Store of determinant information.
! This will quickly become a memory issue, but for dealing with the FCI of small
! systems it is ok.

! Rather than creating an array of type(det), which leads to a serious
! memory overhead due to the need for pointers/allocatable arrays in the
! derived type, we instead create 3 separate variables.  A variable of type det
! can be used to point to the appropriate information, if needed.

! Bit list of the Slater determinant.  See note for f in det type.
! We only store determinants of the same Ms and (for momentum space
! calculations) same ksum at a time.
integer(i0), allocatable, target :: dets_list(:,:) ! (basis_length,ndets)

! Total spin of each Slater determinant stored in dets_list in units of electron spin (1/2).
integer, target :: dets_Ms 

! Sum of wavevectors of the occupied orbitals in each Slater determinant stored
! in det_list.  Sum is to within a reciprocal lattice vector.
! Refers to the wavevector of the i-th alpha spin-orbitals.
integer(i0), target :: dets_ksum

! Number of determinants stored in dets.
! This is the number of determinants enumerated in enumerate_determinants with
! the desired spin and momentum symmetry.
integer :: ndets

! Total size of determinant space.
integer :: tot_ndets

! Number of determinants of each symmetry.
integer, allocatable :: sym_space_size(:) ! (nsym)

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
            write (6,'(1X,a32,'//fmt1(3)//',/)') 'Total size of determinant space:', tot_ndets
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

        integer :: ierr

        deallocate(dets_list, stat=ierr)
        deallocate(bit_lookup, stat=ierr)
        deallocate(basis_lookup, stat=ierr)
        deallocate(sym_space_size, stat=ierr)

        close(det_unit, status='keep')

    end subroutine end_determinants

    subroutine set_spin_polarisation(Ms)

        ! Set the spin polarisation information stored in module-level
        ! variables:
        !    dets_Ms: spin of determinants that are being considered.
        !    nalpha, nbeta: number of alpha, beta electrons.
        !    nvirt_alpha, nvirt_beta: number of alpha, beta virtual spin-orbitals.
        ! In:
        !    Ms: spin of determinants that are being considered.

        use errors, only: stop_all

        integer, intent(in) :: Ms

        ! Find the number of determinants with the required spin.
        if (mod(Ms,2) /= mod(nel,2)) call stop_all('set_spin_polarisation','Required Ms not possible.')
         
        dets_Ms = Ms

        nbeta = (nel - Ms)/2
        nalpha = (nel + Ms)/2

        nvirt_alpha = nsites - nalpha
        nvirt_beta = nsites - nbeta

    end subroutine set_spin_polarisation

    subroutine find_sym_space_size()

        ! Finds the number of Slater determinants that can be formed from the
        ! basis functions belonging to each symmetry.
        ! Currently only crystal momentum symmetry is implemented.
        ! Note that this will overflow for large spaces (but we should be not
        ! doing FCI on such spaces anyway...).
        ! The spin polarisation of the system must be set first (by calling
        ! set_spin_polarisation).
        ! find_sym_space_size must be called first for each value of Ms before
        ! enumerating the determinant list.
        ! In:
        !   Ms: spin of determinants to be found. 

        use comb_m, only: binom
        use utils, only: int_fmt
        use bit_utils, only: first_perm, bit_permutation
        use symmetry, only: nsym, sym_table

        integer :: i, j, ierr, ibit
        integer :: nalpha_combinations, nbeta_combinations
        integer :: k_beta, k
        integer(i0) :: f_alpha, f_beta

        allocate(sym_space_size(nsym), stat=ierr)

        nbeta_combinations = binom(nbasis/2, nbeta)
        nalpha_combinations = binom(nbasis/2, nalpha)

        if (system_type == hub_real) then

            sym_space_size = nalpha_combinations*nbeta_combinations

        else

            ! Determinants are assigned a given symmetry by the sum of the
            ! wavevectors of the occupied basis functions.  This is because only
            ! doubly excitations are connected and D and D_{ij}^{ab} are only
            ! connected if k_i + k_j - k_a - k_b is a reciprocal lattice vector.
            ! Thus we can regard the sum of the wavevectors of the occupied
            ! spin-orbitals as a symmetry label.

            sym_space_size = 0

            do i = 1, nbeta_combinations

                ! Get beta orbitals.
                if (i == 1) then
                    f_beta = first_perm(nbeta)
                else
                    f_beta = bit_permutation(f_beta)
                end if

                k_beta = 1
                do ibit = 0, i0_end
                    if (btest(f_beta,ibit)) k_beta = sym_table(ibit+1, k_beta)
                end do

                do j = 1, nalpha_combinations

                    ! Get alpha orbitals.
                    if (j == 1) then
                        f_alpha = first_perm(nalpha)
                    else
                        f_alpha = bit_permutation(f_alpha)
                    end if

                    ! Symmetry of all orbitals.
                    k = k_beta
                    do ibit = 0, i0_end
                        if (btest(f_alpha,ibit)) k = sym_table(ibit+1, k)
                    end do

                    sym_space_size(k) = sym_space_size(k) + 1

                end do
            end do

        end if

        if (parent) then
            write (6,'(1X,a13,/,1X,13("-"),/)') 'Size of space'
            write (6,'(1X,a5,6X,a6)') 'Index','# dets'
            do i = 1, nsym
                write (6,'(1X,i4,i13)') i, sym_space_size(i)
            end do
            write (6,'()')
        end if

    end subroutine find_sym_space_size

    subroutine enumerate_determinants(ksum)
    
        ! Find the Slater determinants that can be formed from the
        ! basis functions.  The list of determinants is stored in the
        ! module level dets_list array.
        ! find_sym_space_size must be called first for each value of Ms.
        ! In:
        !   ksum: index of a wavevector.  Only determinants with the same
        !         wavevector (up to a reciprocal lattice vector) are stored.
        !         Ignored for the real space formulation of the Hubbard model.

        use comb_m, only: binom
        use errors, only: stop_all
        use utils, only: get_free_unit, int_fmt
        use bit_utils, only: first_perm, bit_permutation
        use symmetry, only: nsym, sym_table

        integer, intent(in) :: ksum

        integer :: i, j, idet, ierr, ibit
        integer :: nalpha_combinations, nbeta_combinations
        integer :: k_beta, k
        character(2) :: fmt1
        integer(i0) :: f_alpha, f_beta

        if (allocated(dets_list)) deallocate(dets_list, stat=ierr)

        nbeta_combinations = binom(nbasis/2, nbeta)
        nalpha_combinations = binom(nbasis/2, nalpha)

        ndets = sym_space_size(ksum)

        allocate(dets_list(basis_length, ndets), stat=ierr)

        ! Assume that we're not attempting to do FCI for more than
        ! a 2*i0_length spin orbitals, which is quite large... ;-)
        if (nbasis/2 > i0_length) then
            call stop_all('enumerate_determinants','Number of alpha spin functions longer than the an i0 integer.')
        end if

        idet = 0
        do i = 1, nbeta_combinations

            ! Get beta orbitals.
            if (i == 1) then
                f_beta = first_perm(nbeta)
            else
                f_beta = bit_permutation(f_beta)
            end if

            ! Symmetry of the beta orbitals.
            if (system_type /= hub_real) then
                k_beta = 1
                do ibit = 0, i0_end
                    if (btest(f_beta,ibit)) k_beta = sym_table(ibit+1, k_beta)
                end do
            end if

            do j = 1, nalpha_combinations

                ! Get alpha orbitals.
                if (j == 1) then
                    f_alpha = first_perm(nalpha)
                else
                    f_alpha = bit_permutation(f_alpha)
                end if

                ! Symmetry of all orbitals.
                if (system_type /= hub_real) then
                    k = k_beta
                    do ibit = 0, i0_end
                        if (btest(f_alpha,ibit)) k = sym_table(ibit+1, k)
                    end do
                end if

                if (system_type == hub_real .or. k == ksum) then

                    idet = idet + 1

                    ! Merge alpha and beta sets into determinant list.
                    ! Alpha orbitals are stored in the even bits, beta orbitals in
                    ! the odd bits (hence the conversion).
                    dets_list(:,idet) = 0
                    do ibit = 0, min(nbasis/2, i0_length/2-1)
                        if (btest(f_alpha,ibit)) dets_list(1,idet) = ibset(dets_list(1,idet), 2*ibit)
                        if (btest(f_beta,ibit)) dets_list(1,idet) = ibset(dets_list(1,idet), 2*ibit+1)
                    end do
                    do ibit = i0_length/2, max(nbasis/2, i0_end)
                        if (btest(f_alpha,ibit)) dets_list(2,idet) = ibset(dets_list(2,idet), 2*ibit-i0_length)
                        if (btest(f_beta,ibit)) dets_list(2,idet) = ibset(dets_list(2,idet), 2*ibit+1-i0_length)
                    end do

                end if

            end do
        end do

        dets_ksum = ksum

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
        d%Ms => dets_Ms
        if (system_type /= hub_real) d%ksum = dets_ksum
    
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

    pure subroutine decode_det_occ_unocc(f, occ_list, unocc_list_alpha, unocc_list_beta)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.  The unoccupied alpha and beta
        ! orbitals are given separately, as this is convenient for FCIQMC.
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    occ_list(nel): integer list of occupied spin-orbitals in the Slater determinant.
        !    unocc_list_alpha(nvirt_alpha): integer list of unoccupied alpha
        !        spin-orbitals in the Slater determinant.
        !    unocc_list_beta(nvirt_beta): integer list of unoccupied beta
        !        spin-orbitals in the Slater determinant.

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(out) :: occ_list(nel), unocc_list_alpha(nvirt_alpha), unocc_list_beta(nvirt_beta)
        integer :: i, j, iocc, iunocc_a, iunocc_b

        iocc = 0
        iunocc_a = 0
        iunocc_b = 0

        do i = 1, basis_length
            do j = 0, i0_end
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    occ_list(iocc) = basis_lookup(j, i)
                else
                    if (mod(j,2)==0) then
                        ! alpha state
                        iunocc_a = iunocc_a + 1
                        unocc_list_alpha(iunocc_a) = basis_lookup(j, i)
                    else
                        ! beta state 
                        iunocc_b = iunocc_b + 1
                        unocc_list_beta(iunocc_b) = basis_lookup(j, i)
                    end if
                end if
                ! Have we covered all basis functions?
                ! This avoids examining any "padding" at the end of f.
                if (iocc+iunocc_a+iunocc_b==nbasis) exit
            end do
        end do

    end subroutine decode_det_occ_unocc

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
