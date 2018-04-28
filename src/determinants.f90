module determinants

! Generation, inspection and manipulation of Slater determinants.

use const
use determinant_data
use parallel, only: parent

implicit none

contains

!--- Initialisation and finalisation of module-level variables ---

    subroutine init_determinants(sys, ex_level)

        ! Initialise determinant information: number of determinants, how
        ! they are stored as bit strings, lookup arrays for converting from
        ! integer list of orbitals to bit strings and vice versa.

        ! In:
        !    sys: system being studied.
        !    ex_level: truncation level being considered.

        use checking, only: check_allocate
        use system, only: sys_t

        use calc, only: ras, ras1, ras3, ras1_min

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: ex_level

        integer :: i, bit_pos, bit_element, ierr

        if (all(ras > 0)) then
            allocate(ras1(sys%basis%tot_string_len), stat=ierr)
            call check_allocate('ras1', sys%basis%tot_string_len, ierr)
            allocate(ras3(sys%basis%tot_string_len), stat=ierr)
            call check_allocate('ras3', sys%basis%tot_string_len, ierr)
            ras1 = 0
            ras3 = 0
            ras1_min = sys%nel - ex_level
            do i = 1, 2*ras(1) ! RAS is in *spatial* orbitals.
                bit_pos = sys%basis%bit_lookup(1,i)
                bit_element = sys%basis%bit_lookup(2,i)
                ras1(bit_element) = ibset(ras1(bit_element), bit_pos)
            end do
            do i = sys%basis%nbasis-2*ras(2)+1, sys%basis%nbasis
                bit_pos = sys%basis%bit_lookup(1,i)
                bit_element = sys%basis%bit_lookup(2,i)
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

        if (allocated(ras1)) then
            deallocate(ras1, stat=ierr)
            call check_deallocate('ras1',ierr)
            deallocate(ras3, stat=ierr)
            call check_deallocate('ras3',ierr)
        end if

    end subroutine end_determinants

!--- Initialisation and finalisation of det_info_t objects ---

    subroutine alloc_det_info_t(sys, det_info, allocate_bit_strings)

        ! Allocate the components of a det_info_t variable.

        ! In:
        !    sys: system to be studied (which defines the length of the
        !        det_info_t components).
        !    allocate_bit_strings (optional): if true (default), allocate the
        !        bit string attributes.  If false, then the bit string attributes
        !        can be used to point to already allocated bit strings.
        !        If set to false, the programmer *must* set allocated_bit_strings to
        !        false when calling dealloc_det_info_t.
        ! Out:
        !    det_info: det_info variable with components allocated to the
        !        appropriate sizes.

        use checking, only: check_allocate
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        logical, intent(in), optional :: allocate_bit_strings
        type(det_info_t), intent(inout) :: det_info
        logical :: alloc_f
        integer :: ierr

        ! Bit strings...
        if (present(allocate_bit_strings)) then
            alloc_f = allocate_bit_strings
        else
            alloc_f = .true.
        end if
        if (alloc_f) then
            allocate(det_info%f(sys%basis%tot_string_len), stat=ierr)
            call check_allocate('det_info%f',sys%basis%tot_string_len,ierr)
            allocate(det_info%f2(sys%basis%tot_string_len), stat=ierr)
            call check_allocate('det_info%f2',sys%basis%tot_string_len,ierr)
        end if

        ! Components for occupied basis functions...
        allocate(det_info%occ_list(sys%nel), stat=ierr)
        call check_allocate('det_info%occ_list',sys%nel,ierr)
        allocate(det_info%occ_list_alpha(sys%nel), stat=ierr)
        call check_allocate('det_info%occ_list_alpha',sys%nalpha,ierr)
        allocate(det_info%occ_list_beta(sys%nel), stat=ierr)
        call check_allocate('det_info%occ_list_beta',sys%nbeta,ierr)

        ! Components for unoccupied basis functions...
        allocate(det_info%unocc_list(sys%nvirt), stat=ierr)
        call check_allocate('det_info%unocc_list',sys%nvirt,ierr)
        allocate(det_info%unocc_list_alpha(sys%nvirt), stat=ierr)
        call check_allocate('det_info%unocc_list_alpha',sys%nvirt_alpha,ierr)
        allocate(det_info%unocc_list_beta(sys%nvirt), stat=ierr)
        call check_allocate('det_info%unocc_list_beta',sys%nvirt_beta,ierr)

        ! Components for symmetry summary of unoccupied basis functions...
        allocate(det_info%symunocc(2,sys%sym0_tot:sys%sym_max_tot), stat=ierr)
        call check_allocate('det_info%symunocc', 2*sys%nsym_tot, ierr)

        ! Weights for excitation generators for this det. 
        allocate(det_info%i_d_occ%weights(sys%nel), stat=ierr)
        call check_allocate('det_info%i_d_occ%weights', sys%nel, ierr)
        allocate(det_info%i_s_occ%weights(sys%nel), stat=ierr)
        call check_allocate('det_info%i_s_occ%weights', sys%nel, ierr)
        allocate(det_info%ia_s_weights_occ(sys%nvirt, sys%nel), stat=ierr)
        call check_allocate('det_info%ia_s_weights_occ', sys%nvirt*sys%nel, ierr)

        ! Orbitals that are different between this det and the ref det.
        allocate(det_info%ref_cdet_occ_list(sys%nel), stat=ierr)
        call check_allocate('det_info%ref_cdet_occ_list', sys%nel, ierr)

    end subroutine alloc_det_info_t

    subroutine dealloc_det_info_t(det_info, allocated_bit_strings)

        ! Deallocate the components of a det_info_t variable.

        ! In:
        !    allocated_bit_strings (optional): if true (default), the
        !        bit string attributes are allocated and must be deallocated.
        !        If false, then the bit string attributes were just used to
        !        point to already allocated bit strings and don't need to be
        !        deallocated.  This *must* correspond to the alloc_bit_strings
        !        argument given to alloc_det_info_t.
        ! Out:
        !    det_info: det_info variable with all components deallocated.

        use checking, only: check_deallocate

        logical, intent(in), optional :: allocated_bit_strings
        type(det_info_t), intent(inout) :: det_info
        integer :: ierr

        logical :: alloc_f

        if (present(allocated_bit_strings)) then
            alloc_f = allocated_bit_strings
        else
            alloc_f = .true.
        end if

        if (alloc_f) then
            deallocate(det_info%f, stat=ierr)
            call check_deallocate('det_info%f',ierr)
            deallocate(det_info%f2, stat=ierr)
            call check_deallocate('det_info%f2',ierr)
        end if
        deallocate(det_info%occ_list, stat=ierr)
        call check_deallocate('det_info%occ_list',ierr)
        deallocate(det_info%unocc_list, stat=ierr)
        call check_deallocate('det_info%unocc_list',ierr)
        deallocate(det_info%occ_list_alpha, stat=ierr)
        call check_deallocate('det_info%occ_list_alpha',ierr)
        deallocate(det_info%occ_list_beta, stat=ierr)
        call check_deallocate('det_info%occ_list_beta',ierr)
        deallocate(det_info%unocc_list_alpha, stat=ierr)
        call check_deallocate('det_info%unocc_list_alpha',ierr)
        deallocate(det_info%unocc_list_beta, stat=ierr)
        call check_deallocate('det_info%unocc_list_beta',ierr)
        deallocate(det_info%symunocc, stat=ierr)
        call check_deallocate('det_info%symunocc',ierr)
        deallocate(det_info%i_d_occ%weights, stat=ierr)
        call check_deallocate('det_info%i_d_occ%weights', ierr)
        deallocate(det_info%i_s_occ%weights, stat=ierr)
        call check_deallocate('det_info%i_s_occ%weights', ierr)
        deallocate(det_info%ia_s_weights_occ, stat=ierr)
        call check_deallocate('det_info%ia_s_weights_occ', ierr)
        deallocate(det_info%ref_cdet_occ_list, stat=ierr)
        call check_deallocate('det_info%ref_cdet_occ_list', ierr)

    end subroutine dealloc_det_info_t

!--- Encode determinant bit strings ---

    pure subroutine encode_det(basis_set, occ_list, bit_list)

        ! In:
        !    basis_set: information about the single-particle basis.
        !    occ_list(nel): integer list of occupied orbitals in the Slater determinant.
        ! Out:
        !    bit_list(tot_string_len): a bit string representation of the occupied
        !        orbitals.   The first element contains the first i0_length basis
        !        functions, the second element the next i0_length and so on.  A basis
        !        function is ocupied if the relevant bit is set.

        use basis_types, only: basis_t

        type(basis_t), intent(in) :: basis_set
        integer, intent(in) :: occ_list(:)
        integer(i0), intent(out) :: bit_list(basis_set%tot_string_len)
        integer :: i, orb, bit_pos, bit_element

        bit_list = 0
        do i = 1, size(occ_list)
            orb = occ_list(i)
            bit_pos = basis_set%bit_lookup(1,orb)
            bit_element = basis_set%bit_lookup(2,orb)
            bit_list(bit_element) = ibset(bit_list(bit_element), bit_pos)
        end do

    end subroutine encode_det

!--- Decode determinant bit strings ---

    pure subroutine decode_det(basis_set, f, occ_list, excit_gen_data)

        ! In:
        !    basis_set: information about the single-particle basis.
        !    f(:): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    occ_list(:): integer list of occupied orbitals in the Slater
        !        determinant. (min size: number of electrons.)

        ! This algorithm has a look over Nbits/256 rather than Nbits, and so
        ! is most likely dominated by O(N_el) scaling.

        use basis_types, only: basis_t
        use bit_table_256_m, only: bit_table_256
        use excit_gens, only: excit_gen_data_t

        type(basis_t), intent(in) :: basis_set
        integer(i0), intent(in) :: f(basis_set%tot_string_len)
        integer, intent(out) :: occ_list(:)
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data

        ! The lookup table contains the list of bits set for all possible integers contained in a given number of bits.
        ! Number of bits in integers in the lookup table (assume a power of 2!).
        integer, parameter :: field_size = ubound(bit_table_256, dim=1)
        ! Number of such bit chunks in integers of kind i0.
        integer, parameter :: nfields = i0_length/field_size
        ! Bit mask to extract a chunk containing field_size bits.
        integer(i0), parameter :: mask = 2**field_size - 1

        integer :: iel, ifield, nfound, nbits_seen
        integer(i0) :: offset, field

        ! WARNING: we assume that the basis functions 1,2,..., correspond to bits 0,1,...
        ! in the first integer of f and so on (i.e. basis_set%separate_strings is false).

        nfound = 0
        nbits_seen = 0
        outer: do iel = 1, basis_set%bit_string_len
            offset = 0
            do ifield = 1, nfields
                ! Inspect one byte at a time.
                field = iand(mask, ishft(f(iel), -offset))
                associate(in_field=>bit_table_256(0,field))
                    ! 1-based index in lookup table, which matches the orbitals indexing scheme.
                    occ_list(nfound+1:nfound+in_field) = bit_table_256(1:in_field, field) + nbits_seen
                    nfound = nfound + in_field
                end associate
                offset = offset + field_size
                nbits_seen = nbits_seen + field_size
            end do
            if (nfound == size(occ_list)) exit outer
        end do outer

    end subroutine decode_det

!--- Extract information from bit strings ---

    pure function spin_orb_list(basis_fns, orb_list) result(ms)

        ! In:
        !    basis_fns: list of single-particle basis functions.
        !    orb_list: list of orbitals (e.g. determinant).
        ! Returns:
        !    Ms: total spin of the determinant in units of electron spin (1/2).

        use basis_types, only: basis_fn_t

        integer :: ms
        type(basis_fn_t), intent(in) :: basis_fns(:)
        integer, intent(in) :: orb_list(:)

        integer :: i

        ms = 0
        do i = lbound(orb_list, dim=1), ubound(orb_list, dim=1)
            ms = ms + basis_fns(orb_list(i))%Ms
        end do

    end function spin_orb_list

!--- Output ---

    subroutine write_det(basis_set, nel, f, iunit, new_line)

        ! Write out a determinant as a list of occupied orbitals in the
        ! Slater determinant.
        ! In:
        !    basis_set: information about the single-particle basis.
        !    nel: number of electrons in system.
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        !    iunit (optional): io unit to which the output is written.
        !        Default: 6 (stdout).
        !    new_line (optional): if true, then a new line is written at
        !        the end of the list of occupied orbitals.  Default: no
        !        new line.

        use utils, only: int_fmt
        use basis_types, only: basis_t

        type(basis_t), intent(in) :: basis_set
        integer, intent(in) :: nel
        integer(i0), intent(in) :: f(basis_set%tot_string_len)
        integer, intent(in), optional :: iunit
        logical, intent(in), optional :: new_line
        integer :: occ_list(nel), io, i
        character(4) :: fmt1

        if (present(iunit)) then
            io = iunit
        else
            io = 6
        end if

        call decode_det(basis_set, f, occ_list)
        fmt1 = int_fmt(basis_set%nbasis,1)

        write (io,'("|")', advance='no')
        do i = 1, nel
            write (io,'('//fmt1//')', advance='no') occ_list(i)
        end do
        write (io,'(1X,">")', advance='no')
        if (present(new_line)) then
            if (new_line) write (io,'()')
        end if

    end subroutine write_det

    subroutine update_sys_spin_info(cdet, sys)

        ! Determine the spin polarisation from a given determinant and set
        ! system spin polarisation accordingly.

        ! In:
        !    cdet: det_info_t object with occ_list set, from which the total ms
        !       derives.
        ! In/Out:
        !    sys: sys_t object. On output spin polarisation (nalpha, nvirt ..)
        !       will be correctly set.

        use bit_utils, only: count_set_bits, count_even_set_bits
        use system, only: sys_t, heisenberg

        type(det_info_t), intent(in) :: cdet
        type(sys_t), intent(inout) :: sys

        select case (sys%system)
        case (heisenberg)
            sys%nel = sum(count_set_bits(cdet%f))
            sys%nvirt = sys%lattice%nsites - sys%nel
        case default
            sys%nalpha = sum(count_even_set_bits(cdet%f))
            sys%nbeta = sys%nel - sys%nalpha
            sys%nvirt_alpha = sys%basis%nbasis/2 - sys%nalpha
            sys%nvirt_beta = sys%basis%nbasis/2 - sys%nbeta
        end select

    end subroutine update_sys_spin_info

    pure function sum_sp_eigenvalues_occ_list(sys, occ_list) result(spe_sum)

        ! In:
        !    sys: system being studied.
        !    occ_list: list of occupied orbitals.
        ! Returns:
        !    Sum of the single particle energies.

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(sys%nel)

        integer :: iorb
        real(p) :: spe_sum

        spe_sum = 0.0_p

        do iorb = 1, sys%nel
            spe_sum = spe_sum + sys%basis%basis_fns(occ_list(iorb))%sp_eigv
        end do

    end function sum_sp_eigenvalues_occ_list

    pure function sum_sp_eigenvalues_bit_string(sys, f) result(spe_sum)

        ! In:
        !    sys: system being studied.
        !    f: bit-string representation of a determinant.
        ! Returns:
        !    Sum of the single particle energies.

        use system, only: sys_t

        real(p) :: spe_sum
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(:)
        integer :: occ(sys%nel)

        call decode_det(sys%basis, f, occ)
        spe_sum = sum_sp_eigenvalues_occ_list(sys, occ)

    end function sum_sp_eigenvalues_bit_string

end module determinants
