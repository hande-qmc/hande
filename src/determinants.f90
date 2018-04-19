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
        allocate(det_info%i_d_weights_occ(sys%nel), stat=ierr)
        call check_allocate('det_info%i_d_weights_occ', sys%nel, ierr)
        allocate(det_info%i_s_weights_occ(sys%nel), stat=ierr)
        call check_allocate('det_info%i_s_weights_occ', sys%nel, ierr)
        allocate(det_info%ia_s_weights_occ(sys%nvirt, sys%nel), stat=ierr)
        call check_allocate('det_info%ia_s_weights_occ', sys%nvirt*sys%nel, ierr)

        ! Orbitals that are different between this det and the ref det.
        allocate(det_info%diff_det_to_ref_orbs(sys%nel, sys%nel), stat=ierr)
        call check_allocate('det_info%diff_det_to_ref_orbs', sys%nel**2, ierr)

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
        deallocate(det_info%i_d_weights_occ, stat=ierr)
        call check_deallocate('det_info%i_d_weights_occ', ierr)
        deallocate(det_info%i_s_weights_occ, stat=ierr)
        call check_deallocate('det_info%i_s_weights_occ', ierr)
        deallocate(det_info%ia_s_weights_occ, stat=ierr)
        call check_deallocate('det_info%ia_s_weights_occ', ierr)
        deallocate(det_info%diff_det_to_ref_orbs, stat=ierr)
        call check_deallocate('det_info%diff_det_to_ref_orbs', ierr)

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

    pure subroutine decode_det_occ(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        !
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data

        call decode_det(sys%basis, f, d%occ_list)

    end subroutine decode_det_occ
    
    pure subroutine decode_det_spinocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals and integer lists containing occupied alpha and
        ! beta orbitals.
        ! In:
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha/_beta: integer list of occupied alpha/beta
        !            spin-orbitals in the Slater determinant.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        integer :: i, ialpha, ibeta
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data

        call decode_det(sys%basis, f, d%occ_list)

        ialpha = 1
        ibeta = 1

        do i = 1, sys%nel
            associate(orb=>d%occ_list(i))
                if (sys%basis%basis_fns(orb)%ms > 0) then
                    d%occ_list_alpha(ialpha) = orb
                    ialpha = ialpha + 1
                else
                    d%occ_list_beta(ibeta) = orb
                    ibeta = ibeta + 1
                end if
            end associate
        end do

    end subroutine decode_det_spinocc
    
    pure subroutine decode_det_occ_unocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.
        !
        ! We return the lists for alpha and beta electrons separately.
        !
        ! In:
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        unocc_list: integer list of unoccupied
        !            spin-orbitals in the Slater determinant.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, j, iocc, iunocc, orb, last_basis_ind

        ! A bit too much to do the chunk-based decoding of the occupied list and then fill
        ! in the remaining information.  We only use this in Hubbard model calculations in
        ! k-space, so for now just do a (slow) bit-wise inspection.

        iocc = 0
        iunocc = 0
        orb = 0

        do i = 1, sys%basis%bit_string_len - 1
            ! Manual unrolling allows us to avoid 2 mod statements
            ! and some branching.
            ! [todo] - note above still relevant?
            do j = 0, i0_end
                orb = orb + 1
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    d%occ_list(iocc) = orb
                else
                    iunocc = iunocc + 1
                    d%unocc_list(iunocc) = orb
                end if
            end do
        end do

        ! Deal with the last element in the determinant bit array separately.
        ! Note that decoding a bit string is surprisingly slow (or, more
        ! importantly, adds up when doing billions of times).
        ! Treating the last element as a special case rather than having an if
        ! statement in the above loop results a speedup of the Hubbard k-space
        ! FCIQMC calculations of 1.5%.
        last_basis_ind = sys%basis%nbasis - i0_length*(sys%basis%bit_string_len-1) - 1
        do j = 0, last_basis_ind, 2
            orb = orb + 1
            if (btest(f(i), j)) then
                iocc = iocc + 1
                d%occ_list(iocc) = orb
            else
                iunocc = iunocc + 1
                d%unocc_list(iunocc) = orb
            end if
        end do

    end subroutine decode_det_occ_unocc

    pure subroutine decode_det_occ_symunocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        ! In:
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        symunocc(2, sym0_tot:symmax_tot): number of unoccupied orbitals of each
        !            spin/symmetry.  The same indexing scheme is used for
        !            nbasis_sym_spin.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t


        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, ims, isym

        call decode_det(sys%basis, f, d%occ_list)

        d%symunocc = sys%read_in%pg_sym%nbasis_sym_spin

        do i = 1, sys%nel
            associate(orb=>d%occ_list(i))
                ims = (sys%basis%basis_fns(orb)%ms+3)/2
                isym = sys%basis%basis_fns(orb)%sym
            end associate
            d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
        end do

    end subroutine decode_det_occ_symunocc

    pure subroutine decode_det_spinocc_symunocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        ! In:
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha/_beta: integer list of occupied alpha/beta
        !            spin-orbitals in the Slater determinant.
        !        symunocc(2, sym0_tot:symmax_tot): number of unoccupied orbitals of each
        !            spin/symmetry.  The same indexing scheme is used for
        !            nbasis_sym_spin.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, ims, isym, ialpha, ibeta

        call decode_det(sys%basis, f, d%occ_list)

        d%symunocc = sys%read_in%pg_sym%nbasis_sym_spin
        ialpha = 1
        ibeta = 1

        do i = 1, sys%nel
            associate(orb=>d%occ_list(i))
                ims = (sys%basis%basis_fns(orb)%ms+3)/2
                isym = sys%basis%basis_fns(orb)%sym
                if (ims == 2) then
                    d%occ_list_alpha(ialpha) = orb
                    ialpha = ialpha + 1
                else
                    d%occ_list_beta(ibeta) = orb
                    ibeta = ibeta + 1
                end if
            end associate
            d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
        end do

    end subroutine decode_det_spinocc_symunocc
    
    pure subroutine decode_det_spinocc_spinunocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.
        !
        ! We return the lists for alpha and beta electrons separately.
        !
        ! In:
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha: integer list of occupied alpha
        !            spin-orbitals in the Slater determinant.
        !        occ_list_beta: integer list of occupied beta
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_alpha: integer list of unoccupied alpha
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_beta: integer list of unoccupied beta
        !            spin-orbitals in the Slater determinant.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, j, iocc, iocc_a, iocc_b, iunocc_a, iunocc_b, orb, last_basis_ind

        ! A bit too much to do the chunk-based decoding of the occupied list and then fill
        ! in the remaining information.  We only use this in Hubbard model calculations in
        ! k-space, so for now just do a (slow) bit-wise inspection.

        iocc = 0
        iocc_a = 0
        iocc_b = 0
        iunocc_a = 0
        iunocc_b = 0
        orb = 0

        do i = 1, sys%basis%bit_string_len - 1
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
        last_basis_ind = sys%basis%nbasis - i0_length*(sys%basis%bit_string_len-1) - 1
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
    
    pure subroutine decode_det_spinocc_spinsymunocc(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.
        !
        ! We return the lists for alpha and beta electrons separately.
        !
        ! In:
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha: integer list of occupied alpha
        !            spin-orbitals in the Slater determinant.
        !        occ_list_beta: integer list of occupied beta
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_alpha: integer list of unoccupied alpha
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_beta: integer list of unoccupied beta
        !            spin-orbitals in the Slater determinant.
        !        symunocc(2, sym0_tot:symmax_tot): number of unoccupied orbitals of each
        !            spin/symmetry.  The same indexing scheme is used for
        !            nbasis_sym_spin.

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, j, iocc, iocc_a, iocc_b, iunocc_a, iunocc_b, orb, last_basis_ind, isym, ims

        d%symunocc = sys%read_in%pg_sym%nbasis_sym_spin
        
        ! A bit too much to do the chunk-based decoding of the occupied list and then fill
        ! in the remaining information.

        iocc = 0
        iocc_a = 0
        iocc_b = 0
        iunocc_a = 0
        iunocc_b = 0
        orb = 0

        do i = 1, sys%basis%bit_string_len - 1
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
                    ims = 2
                    isym = sys%basis%basis_fns(orb)%sym
                    d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
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
                    ims = 1
                    isym = sys%basis%basis_fns(orb)%sym
                    d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
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
        last_basis_ind = sys%basis%nbasis - i0_length*(sys%basis%bit_string_len-1) - 1
        do j = 0, last_basis_ind, 2
            ! Test alpha orbital.
            orb = orb + 1
            if (btest(f(i), j)) then
                iocc = iocc + 1
                iocc_a = iocc_a + 1
                d%occ_list(iocc) = orb
                d%occ_list_alpha(iocc_a) = orb
                ims = 2
                isym = sys%basis%basis_fns(orb)%sym
                d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
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
                ims = 1
                isym = sys%basis%basis_fns(orb)%sym
                d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
            else
                iunocc_b = iunocc_b + 1
                d%unocc_list_beta(iunocc_b) = orb
            end if
        end do

    end subroutine decode_det_spinocc_spinsymunocc
    
    pure subroutine decode_det_ppMij(sys, f, d, excit_gen_data)

        ! WARNING: this decoder needs excit_gen_data!!!!

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.
        !
        ! We return the lists for alpha and beta electrons separately.
        ! Also return the weights to select i in a double excitation.
        !
        ! In:
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        occ_list_alpha: integer list of occupied alpha
        !            spin-orbitals in the Slater determinant.
        !        occ_list_beta: integer list of occupied beta
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_alpha: integer list of unoccupied alpha
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_beta: integer list of unoccupied beta
        !            spin-orbitals in the Slater determinant.
        !        i_d_weights_occ: weights to select i in a double excitation
        !            during excit gen call.
        !        i_d_weights_occ_tot: sum of i_d_weights_occ

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data

        call decode_det_spinocc_spinsymunocc(sys, f, d)
        call find_i_d_weights(sys%nel, excit_gen_data%excit_gen_pp%ppm_i_d_weights, d)

    end subroutine decode_det_ppMij
    
    pure subroutine decode_det_hb(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals and some weights for heat bath excit gen.
        !
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        i_d_weights_occ: weights to select i in a double excitation
        !            during excit gen call.
        !        i_d_weights_occ_tot: sum of i_d_weights_occ

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data

        call decode_det(sys%basis, f, d%occ_list)
        call find_i_d_weights(sys%nel, excit_gen_data%excit_gen_hb%i_weights, d)

    end subroutine decode_det_hb
    
    pure subroutine decode_det_hbu(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals and some weights for hbu excit gen.
        ! In:
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        symunocc(2, sym0_tot:symmax_tot): number of unoccupied orbitals of each
        !            spin/symmetry.  The same indexing scheme is used for
        !            nbasis_sym_spin.
        !        i_d_weights_occ: weights to select i in a double excitation
        !            during excit gen call.
        !        i_d_weights_occ_tot: sum of i_d_weights_occ

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t


        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data
        integer :: i, ims, isym

        call decode_det_occ_symunocc(sys, f, d)
        call find_i_d_weights(sys%nel, excit_gen_data%excit_gen_hb%i_weights, d)

    end subroutine decode_det_hbu
    
    pure subroutine decode_det_hbs(sys, f, d, excit_gen_data)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals and some weights for hbs excit gen.
        !
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(tot_string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        i_d_weights_occ: weights to select i in a double excitation
        !            during excit gen call.
        !        i_d_weights_occ_tot: sum of i_d_weights_occ

        use system, only: sys_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: d
        type(excit_gen_data_t), optional, intent(in) :: excit_gen_data

        call decode_det_occ_unocc(sys, f, d)
        call find_ia_single_weights(sys, d)

    end subroutine decode_det_hbs

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

!--- Helper Functions ---

    pure subroutine find_i_d_weights(nel, i_weights_precalc, d)

        ! Find weights to select i in a double excitation using pre-calculed heat bath weights.

        ! In:
        !   nel: number of electrons
        !   i_weights_precalc: precalculated weights that i could have. Need to reduce them to the ones
        !        that are actually occupied.
        ! In/Out:
        !   d: det_code_t object. Information about the determinant to decode.
        
        integer, intent(in) :: nel
        real(p), intent(in) :: i_weights_precalc(:)
        type(det_info_t), intent(inout) :: d

        integer :: pos_occ

        d%i_d_weights_occ_tot = 0.0_p
        do pos_occ = 1, nel
            d%i_d_weights_occ(pos_occ) = i_weights_precalc(d%occ_list(pos_occ))
            d%i_d_weights_occ_tot = d%i_d_weights_occ_tot + d%i_d_weights_occ(pos_occ)
        end do

    end subroutine find_i_d_weights
    
    pure subroutine find_ia_single_weights(sys, d)

        ! Create a random single excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! This calculates the weights for i and a as exact as possibly.
        ! For i, the weight is \sum_a H_ia, for a given i (a|i) it is H_ia.

        ! In:
        !    sys: system object being studied.
        ! In/Out:
        !    d: information about current determinant to decode

        use proc_pointers, only: abs_hmatel_ptr
        use system, only: sys_t
        use hamiltonian_data, only: hmatel_t
        use hamiltonian_molecular, only: slater_condon1_mol
        use hamiltonian_periodic_complex, only: slater_condon1_periodic_complex

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(inout) :: d
        
        type(hmatel_t) :: hmatel
        integer :: pos, i_ind, a_ind
        
        d%i_s_weights_occ_tot = 0.0_p
        do i_ind = 1, sys%nel
            d%i_s_weights_occ(i_ind) = 0.0_p
            do a_ind = 1, sys%nvirt
                ! [todo] - this reduces the computational time (by calculation ia_weights here)
                ! but more memory costs as after we have selected i, we don't need all elements in ia_weights. Need to balance
                ! these costs.
                ! [todo] - could consider rewritting with proc_pointer here but that would imply having to rewrite slater
                ! condon 1 mol / periodic complex. slater_condon1_excit_mol is not safe enough (no checks).
                if (sys%read_in%comp) then
                    hmatel%r = 0.0_p
                    hmatel%c = slater_condon1_periodic_complex(sys, d%occ_list, d%occ_list(i_ind), d%unocc_list(a_ind), &
                        .false.)
                else
                    hmatel%r = slater_condon1_mol(sys, d%occ_list, d%occ_list(i_ind), d%unocc_list(a_ind), .false.)
                    hmatel%c = cmplx(0.0_p, 0.0_p, p)
                end if
                d%ia_s_weights_occ(a_ind, i_ind) = abs_hmatel_ptr(hmatel)
                d%i_s_weights_occ(i_ind) = d%i_s_weights_occ(i_ind) + d%ia_s_weights_occ(a_ind, i_ind)
            end do
            d%i_s_weights_occ_tot = d%i_s_weights_occ_tot + d%i_s_weights_occ(i_ind)
        end do

    end subroutine find_ia_single_weights

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
