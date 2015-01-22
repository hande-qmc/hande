module determinants

! Generation, inspection and manipulation of Slater determinants.

use ccmc_data, only: cluster_t
use const
use parallel, only: parent

implicit none

! --- FCIQMC info ---

! A handy type for containing a lot of information about a determinant.
! This is convenient for passing around different amounts of info when
! we need consistent interfaces.
! Not all compenents are necessarily allocated: only those needed at the time.
type det_info_t
    ! bit representation of determinant.
    integer(i0), pointer :: f(:)  => NULL()  ! (string_len)
    integer(i0), pointer :: f2(:)  => NULL()  ! (string_len); for DMQMC
    ! List of occupied spin-orbitals.
    integer, pointer :: occ_list(:)  => NULL()  ! (nel)
    ! List of occupied alpha/beta spin-orbitals
    integer, pointer :: occ_list_alpha(:), occ_list_beta(:)
    ! List of unoccupied alpha/beta spin-orbitals
    integer, pointer :: unocc_list_alpha(:), unocc_list_beta(:)
    ! Number of unoccupied orbitals with each spin and symmetry.
    ! The first index maps to spin using (Ms+3)/2, where Ms=-1 is spin-down and
    ! Ms=1 is spin-up.
    integer, pointer :: symunocc(:,:) ! (2,sym0_tot:sym_max_tot)
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
end type det_info_t

contains

!--- Initialisation and finalisation of module-level variables ---

    subroutine init_determinants(sys)

        ! Initialise determinant information: number of determinants, how
        ! they are stored as bit strings, lookup arrays for converting from
        ! integer list of orbitals to bit strings and vice versa.

        ! In/Out:
        !    sys: system being studied.  On input system-level properties (as in
        !         init_system,  e.g. sys%nel, sys%lattice, etc) and basis-level
        !         properties (as in init_basis_t) have been set.

        use checking, only: check_allocate
        use utils, only: binom_r
        use utils, only: get_free_unit, int_fmt
        use calc, only: doing_calc, exact_diag, lanczos_diag, doing_calc, &
                        dmqmc_calc, ras, ras1, ras3, ras1_min, truncation_level
        use system, only: sys_t, heisenberg

        type(sys_t), intent(in) :: sys

        integer :: i, j, k, bit_pos, bit_element, ierr, site_index
        character(4) :: fmt1
        integer(int_64) :: tot_ndets

        tot_ndets = nint(binom_r(sys%basis%nbasis, sys%nel), int_64)
        fmt1 = int_fmt(tot_ndets, padding=1)
        if (parent .and. doing_calc(exact_diag+lanczos_diag)) &
                write (6,'(1X,a32,'//fmt1//')') 'Total size of determinant space:', tot_ndets

        if (all(ras > 0)) then
            allocate(ras1(sys%basis%string_len), stat=ierr)
            call check_allocate('ras1', sys%basis%string_len, ierr)
            allocate(ras3(sys%basis%string_len), stat=ierr)
            call check_allocate('ras3', sys%basis%string_len, ierr)
            ras1 = 0
            ras3 = 0
            ras1_min = sys%nel - truncation_level
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
            allocate(det_info%f(sys%basis%string_len), stat=ierr)
            call check_allocate('det_info%f',sys%basis%string_len,ierr)
            allocate(det_info%f2(sys%basis%string_len), stat=ierr)
            call check_allocate('det_info%f2',sys%basis%string_len,ierr)
        end if

        ! Components for occupied basis functions...
        allocate(det_info%occ_list(sys%nel), stat=ierr)
        call check_allocate('det_info%occ_list',sys%nel,ierr)
        allocate(det_info%occ_list_alpha(sys%nalpha), stat=ierr)
        call check_allocate('det_info%occ_list_alpha',sys%nalpha,ierr)
        allocate(det_info%occ_list_beta(sys%nbeta), stat=ierr)

        ! Components for unoccupied basis functions...
        call check_allocate('det_info%occ_list_beta',sys%nbeta,ierr)
        allocate(det_info%unocc_list_alpha(sys%nvirt_alpha), stat=ierr)
        call check_allocate('det_info%unocc_list_alpha',sys%nvirt_alpha,ierr)
        allocate(det_info%unocc_list_beta(sys%nvirt_beta), stat=ierr)
        call check_allocate('det_info%unocc_list_beta',sys%nvirt_beta,ierr)

        ! Components for symmetry summary of unoccupied basis functions...
        allocate(det_info%symunocc(2,sys%sym0_tot:sys%sym_max_tot), stat=ierr)
        call check_allocate('det_info%symunocc',size(det_info%symunocc),ierr)

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

    end subroutine dealloc_det_info_t

!--- Encode determinant bit strings ---

    pure subroutine encode_det(basis_set, occ_list, bit_list)

        ! In:
        !    basis_set: information about the single-particle basis.
        !    occ_list(nel): integer list of occupied orbitals in the Slater determinant.
        ! Out:
        !    bit_list(string_len): a bit string representation of the occupied
        !        orbitals.   The first element contains the first i0_length basis
        !        functions, the second element the next i0_length and so on.  A basis
        !        function is ocupied if the relevant bit is set.

        use basis_types, only: basis_t

        type(basis_t), intent(in) :: basis_set
        integer, intent(in) :: occ_list(:)
        integer(i0), intent(out) :: bit_list(basis_set%string_len)
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

    pure subroutine decode_det(basis_set, f, occ_list)

        ! In:
        !    basis_set: information about the single-particle basis.
        !    f(:): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    occ_list(:): integer list of occupied orbitals in the Slater
        !        determinant. (min size: number of electrons.)

        use basis_types, only: basis_t

        type(basis_t), intent(in) :: basis_set
        integer(i0), intent(in) :: f(basis_set%string_len)
        integer, intent(out) :: occ_list(:)
        integer :: i, j, iorb

        iorb = 1
        outer: do i = 1, basis_set%string_len
            do j = 0, i0_end
                if (btest(f(i), j)) then
                    occ_list(iorb) = basis_set%basis_lookup(j, i)
                    if (iorb == size(occ_list)) exit outer
                    iorb = iorb + 1
                end if
            end do
        end do outer

    end subroutine decode_det

    pure subroutine decode_det_occ(sys, f, d)

        ! Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        ! In:
        !    sys: system being studied (contains required basis information).
        !    f(string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        type(det_info_t), intent(inout) :: d
        integer :: i, j, iocc, iunocc_a, iunocc_b

        iocc = 0
        iunocc_a = 0
        iunocc_b = 0

        do i = 1, sys%basis%string_len
            do j = 0, i0_end
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    d%occ_list(iocc) = sys%basis%basis_lookup(j, i)
                end if
                if (iocc == sys%nel) exit
            end do
        end do

    end subroutine decode_det_occ

    pure subroutine decode_det_occ_spinunocc(sys, f, d)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.  The unoccupied alpha and beta
        ! orbitals are given separately, as this is convenient for FCIQMC.
        ! In:
        !    f(string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        unocc_list_alpha: integer list of unoccupied alpha
        !            spin-orbitals in the Slater determinant.
        !        unocc_list_beta: integer list of unoccupied beta
        !            spin-orbitals in the Slater determinant.

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        type(det_info_t), intent(inout) :: d
        integer :: i, j, iocc, iunocc_a, iunocc_b

        iocc = 0
        iunocc_a = 0
        iunocc_b = 0

        do i = 1, sys%basis%string_len
            do j = 0, i0_end
                if (btest(f(i), j)) then
                    iocc = iocc + 1
                    d%occ_list(iocc) = sys%basis%basis_lookup(j, i)
                else
                    if (mod(j,2)==0) then
                        ! alpha state (even bit index, odd basis function index)
                        iunocc_a = iunocc_a + 1
                        d%unocc_list_alpha(iunocc_a) = sys%basis%basis_lookup(j, i)
                    else
                        ! beta state (odd bit index, even basis function index)
                        iunocc_b = iunocc_b + 1
                        d%unocc_list_beta(iunocc_b) = sys%basis%basis_lookup(j, i)
                    end if
                end if
                ! Have we covered all basis functions?
                ! This avoids examining any "padding" at the end of f.
                if (iocc+iunocc_a+iunocc_b==sys%basis%nbasis) exit
            end do
        end do

    end subroutine decode_det_occ_spinunocc

    pure subroutine decode_det_spinocc_spinunocc(sys, f, d)

        ! Decode determinant bit string into integer lists containing the
        ! occupied and unoccupied orbitals.
        !
        ! We return the lists for alpha and beta electrons separately.
        !
        ! In:
        !    f(string_len): bit string representation of the Slater
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

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        type(det_info_t), intent(inout) :: d
        integer :: i, j, iocc, iocc_a, iocc_b, iunocc_a, iunocc_b, orb, last_basis_ind

        iocc = 0
        iocc_a = 0
        iocc_b = 0
        iunocc_a = 0
        iunocc_b = 0
        orb = 0

        do i = 1, sys%basis%string_len - 1
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
        last_basis_ind = sys%basis%nbasis - i0_length*(sys%basis%string_len-1) - 1
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

    pure subroutine decode_det_occ_symunocc(sys, f, d)

        !0 Decode determinant bit string into integer list containing the
        ! occupied orbitals.
        ! In:
        !    f(string_len): bit string representation of the Slater
        !        determinant.
        ! Out:
        !    d: det_info_t variable.  The following components are set:
        !        occ_list: integer list of occupied spin-orbitals in the
        !            Slater determinant.
        !        symunocc(2, sym0_tot:symmax_tot): number of unoccupied orbitals of each
        !            spin/symmetry.  The same indexing scheme is used for
        !            nbasis_sym_spin.

        use point_group_symmetry, only: nbasis_sym_spin
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        type(det_info_t), intent(inout) :: d
        integer :: i, j, iocc, iunocc_a, iunocc_b, orb, ims, isym

        iocc = 0
        iunocc_a = 0
        iunocc_b = 0

        d%symunocc = nbasis_sym_spin

        do i = 1, sys%basis%string_len
            do j = 0, i0_end
                if (btest(f(i), j)) then
                    orb = sys%basis%basis_lookup(j, i)
                    ims = (sys%basis%basis_fns(orb)%ms+3)/2
                    isym = sys%basis%basis_fns(orb)%sym
                    iocc = iocc + 1
                    d%occ_list(iocc) = orb
                    d%symunocc(ims, isym) = d%symunocc(ims, isym) - 1
                end if
                if (iocc == sys%nel) exit
            end do
        end do

    end subroutine decode_det_occ_symunocc

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
        !    f(string_len): bit string representation of the Slater
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
        integer(i0), intent(in) :: f(basis_set%string_len)
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

end module determinants
