module determinant_enumeration

! Module for complete enumeration of determinant Hilbert space.
! Used only in FCI and simple_fciqmc calculations.
! Rather a memory hog! ;-)

use const
use determinants

implicit none

!--- Info for FCI calculations ---

! Store of determinant information.
! This will quickly become a memory issue, but for dealing with the FCI of small
! systems it is ok.

! Rather than creating an array of derived types containing information about
! each determinant, which leads to a serious memory overhead due to the need for
! pointers/allocatable arrays in the derived type, we instead create 3 separate
! variables.

! Bit list of the Slater determinant.  See note for f in determinants module.
! We only store determinants of the same Ms and (for momentum space
! calculations) same ksum at a time.
integer(i0), allocatable :: dets_list(:,:) ! (basis_length,ndets)

! Number of determinants stored in dets.
! This is the number of determinants enumerated in enumerate_determinants with
! the desired spin and momentum symmetry.
integer :: ndets

! Total (exact) size of determinant space.
! Number of determinants of each symmetry.
integer, allocatable :: sym_space_size(:) ! (nsym)

! If true then the determinant list is written to determinant_file.
logical :: write_determinants = .false.
character(255) :: determinant_file = 'DETS'
integer :: det_unit

contains

!--- Initialisation and finalisation ---

    subroutine init_determinant_enumeration()

        ! Initialise determinant enumeration info,
        ! Note that memory allocation is done in the determinant enumeration
        ! routines.

        use utils, only: get_free_unit

        if (write_determinants) then
            det_unit = get_free_unit()
            open(det_unit, file=determinant_file, status='unknown')
        end if

    end subroutine init_determinant_enumeration

    subroutine end_determinant_enumeration()

        ! Clean up.

        use checking, only: check_deallocate

        integer :: ierr

        if (allocated(dets_list)) then
            deallocate(dets_list, stat=ierr)
            call check_deallocate('dets_list',ierr)
        end if
        if (allocated(sym_space_size)) then
            deallocate(sym_space_size, stat=ierr)
            call check_deallocate('sym_space_size',ierr)
        end if

        if (write_determinants) close(det_unit, status='keep')

    end subroutine end_determinant_enumeration

!--- Hilbert space enumeration ---

    subroutine enumerate_determinants(init, ref_sym, occ_list0)

        ! Find the Slater determinants that can be formed from the
        ! basis functions.

        ! On initialisation the number of determinants in each symmetry block
        ! for the specified spin polarisation is counted and store in
        ! the module-level variable sym_space_size.  Subseuent calls store
        ! determinants of the desired symmetry in the module-level variable
        ! dets_list array.

        ! In:
        !    init: if true, count the number of determinants rather than store
        !        them.  This must be done for each spin polarisation before
        !        subsequent calls to store the determinants.
        !    ref_sym(optional): index of an irreducible representation.  Only determinants
        !        with the same symmetry  are stored.  *MUST* be specified if
        !        init is false.
        !    occ_list0(optional): list of occupied orbitals in a determinant.
        !        Used as the reference determinant if a truncated CI space is
        !        being used.  *MUST* be specified if truncate_space is true.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use utils, only: binom_i, get_free_unit, int_fmt
        use bit_utils, only: first_perm, bit_permutation, decode_bit_string

        use basis, only: basis_length
        use calc, only: truncate_space, truncation_level
        use excitations, only: get_excitation_level
        use system, only: sym0, sym_max, nel, system_type, hub_real, heisenberg, ueg
        use symmetry, only: cross_product, symmetry_orb_list
        use ueg_system, only: ueg_basis_index

        logical, intent(in) :: init
        integer, intent(in), optional :: ref_sym
        integer, intent(in), optional :: occ_list0(nel)

        integer :: i, j, iel, idet, ierr
        integer :: nbeta_combinations
        integer, allocatable :: nalpha_combinations(:)
        integer :: sym_beta, sym, Ms, excit_level_alpha, excit_level_beta
        character(4) :: fmt1
        integer(i0) :: f(basis_length)
        integer, allocatable :: occ(:), comb(:,:), unocc(:)
        integer :: k(ndim), k_beta(ndim)
        type(det_info) :: d0

        if (init) then
            if (allocated(sym_space_size)) then
                deallocate(sym_space_size, stat=ierr)
                call check_deallocate('sym_space_size',ierr)
            end if
            allocate(sym_space_size(sym0:sym_max), stat=ierr)
            call check_allocate('sym_space_size',nsym,ierr)
            sym_space_size = 0
        else
            if (.not.present(ref_sym)) then
                call stop_all('enumerate_determinants', 'ref_sym must be specified when storing determinants.')
            end if
            if (allocated(dets_list)) then
                deallocate(dets_list, stat=ierr)
                call check_deallocate('dets_list',ierr)
            end if

            ndets = sym_space_size(ref_sym)

            allocate(dets_list(basis_length, ndets), stat=ierr)
            call check_allocate('dets_list',basis_length*ndets,ierr)
        end if

        call alloc_det_info(d0)

        ! For fermionic systems, we generate the alpha string and beta string
        ! separately and then combine them to form a determinant.
        ! If the Hilbert space is truncated, then only strings within the
        ! required excitation level are generated.
        ! For each beta string we generate, we then generate each
        ! (excitation-allowed) alpha string.  Hence if the Hilbert space is
        ! truncated, then the number of possilbe alpha strings is a function of
        ! the excitation level of the beta string.
        ! If the Hilbert space is not truncated, then all alpha strings can be
        ! combined with all beta strings.
        ! This scheme makes the generation of determinants in a truncated
        ! Hilbert space fast (as we only generate determinants within the
        ! truncated space and not those outside the truncated space but in the
        ! full Hilbert space).
        if (truncate_space) then
            if (.not.present(occ_list0)) then
                call stop_all('enumerate_determinants', 'Require reference for truncated CI spaces.')
            end if
            allocate(nalpha_combinations(0:nbeta), stat=ierr)
            nalpha_combinations = 0
            nbeta_combinations = 0
            ! Identical to the settings in the else statements if truncation_level == nel.
            do i = 0, min(truncation_level, nbeta)
                nbeta_combinations = nbeta_combinations + binom_i(nbeta, i)*binom_i(nvirt_beta, i)
                do j = 0, min(truncation_level-i, nalpha)
                    nalpha_combinations(i) = nalpha_combinations(i) + binom_i(nalpha, j)*binom_i(nvirt_alpha, j)
                end do
            end do
            call encode_det(occ_list0, d0%f)
            call decode_det_spinocc_spinunocc(d0%f, d0)
        else
            ! No matter what the excitation level of the beta string, all alpha
            ! combinations are allowed.
            allocate(nalpha_combinations(0:0), stat=ierr)
            nbeta_combinations = binom_i(nbasis/2, nbeta)
            nalpha_combinations = binom_i(nbasis/2, nalpha)
            truncation_level = nel
        end if
        call check_allocate('nalpha_combinations',size(nalpha_combinations),ierr)
        allocate(comb(2*truncation_level, 2), stat=ierr)
        call check_allocate('comb', size(comb), ierr)

        allocate(occ(nel), stat=ierr)
        call check_allocate('occ', nel, ierr)

        select case(system_type)

        case(heisenberg, chung_landau)

            ! See notes in system about how the Heisenberg model uses nel and
            ! nvirt.

            ! Just have spin up and spin down sites (no concept of unoccupied
            ! sites) so just need to arrange nel spin ups.  No need to
            ! interweave alpha and beta strings.

            if (init) then
                ! Easy to work out the number of spin kets ('determinants') in
                ! the Hilbert space.  No need to enumerate them all first.
                if (truncate_space) then
                    sym_space_size = 0
                    do i = 0, truncation_level
                        sym_space_size = sym_space_size + binom_i(nel, truncation_level)*binom_i(nsites-nel, truncation_level)
                    end do
                else
                    sym_space_size = binom_i(nsites, nel)
                end if
            else
                if (truncate_space) then
                    ! Need list of spin-down sites (ie 'unoccupied' bits).
                    allocate(unocc(nbasis-nel), stat=ierr)
                    call check_allocate('unocc', size(unocc), ierr)
                    call decode_bit_string(d0%f(1), unocc)
                    do i = 1, ndets
                        call next_string(i==1, .false., nbasis, nel, d0%occ_list, &
                                         unocc, truncation_level, comb(:,1),  &
                                         occ, excit_level_alpha)
                        call encode_det(occ, dets_list(:,i))
                    end do
                    deallocate(unocc, stat=ierr)
                    call check_deallocate('unocc', ierr)
                else
                    ! Easiest just to iterate through all possible bit permutations.
                    ! No symmetry to check!
                    ! Assume that we're not attempting to do FCI for more than
                    ! a i0_length sites , which is quite large... ;-)
                    if (nbasis > i0_length) then
                        call stop_all('enumerate_determinants','Number of spin functions longer than the an i0 integer.')
                    end if
                    dets_list(1,1) = first_perm(nel)
                    do i = 2, ndets
                        dets_list(1,i) = bit_permutation(dets_list(1,i-1))
                    end do
                end if
            end if

        case default

            if (init .and. system_type == hub_real .and. .not.truncate_space) then
                ! Closed form for size of Hilbert space.  No symmetry implemented!
                sym_space_size = nalpha_combinations*nbeta_combinations
            else
                ! Simply count the number of allowed determinants.
                ! CBA to find the closed form for the Hubbard model with local
                ! orbitals (ie no symmetry constraints).
                ! If not initialising (ie counting the number of allowed
                ! determinants) then we store the allowed determinants as they
                ! are generated.

                ! Determinants are assigned a given symmetry by the direct product
                ! of the representations spanned by the occupied orbitals.  For
                ! momentum space systems this amounts to the sum of the wavevectors
                ! of the occupied basis functions.  Thus we can regard the sum of
                ! the wavevectors of the occupied spin-orbitals as a symmetry label.

                idet = 0

                do i = 1, nbeta_combinations

                    ! Get beta orbitals.
                    call next_string(i==1, .false., nbasis/2, nbeta, d0%occ_list_beta, &
                                     d0%unocc_list_beta, truncation_level, comb(:,1),  &
                                     occ(1:nbeta), excit_level_beta)

                    ! Symmetry.
                    ! Have to treat the UEG as a special case as we don't consider
                    ! all possible symmetries for the UEG but rather only momenta
                    ! which exist in the basis set.
                    if (system_type == ueg) then
                        k_beta = 0
                        do iel = 1, nbeta
                            k_beta = k_beta + basis_fns(occ(iel))%l
                        end do
                    else
                        sym_beta = symmetry_orb_list(occ(1:nbeta))
                    end if

                    do j = 1, nalpha_combinations(excit_level_beta)

                        ! Get alpha orbitals.
                        call next_string(j==1, .true., nbasis/2, nalpha, d0%occ_list_alpha, &
                                         d0%unocc_list_alpha, truncation_level-excit_level_beta, &
                                         comb(:,2), occ(nbeta+1:), excit_level_alpha)

                        ! Symmetry of all orbitals.
                        if (system_type == ueg) then
                            k = k_beta
                            do iel = nbeta+1, nel
                                k = k + basis_fns(occ(iel))%l
                            end do
                            ! Symmetry label (convert from basis index).
                            sym = (ueg_basis_index(k,1)+1)/2
                        else
                            sym = cross_product(sym_beta, symmetry_orb_list(occ(nbeta+1:nel)))
                        end if

                        if (init) then
                            ! Count determinant.
                            if (sym >= lbound(sym_space_size,dim=1) .and. sym <= ubound(sym_space_size,dim=1)) then
                                ! Ignore symmetries outside the basis set (only affects
                                ! UEG currently).
                                sym_space_size(sym) = sym_space_size(sym) + 1
                            end if
                        else if (sym == ref_sym) then
                            ! Store determinant.
                            idet = idet + 1
                            call encode_det(occ, dets_list(:,idet))
                        end if

                    end do

                end do

            end if

        end select

        if (init .and. parent) then
            ! Output information about the size of the space.
            Ms = nalpha - nbeta
            write (6,'(1X,a25,/,1X,25("-"),/)') 'Size of determinant space'
            write (6,'(1X,a75,'//int_fmt(Ms,0)//',a1)') &
                     'The table below gives the number of determinants for each symmetry with Ms=', &
                     Ms,"."
            if (truncate_space) then
                write (6,'(1X,a24,'//int_fmt(truncation_level,1)//',1X,a54)') &
                    'Only determinants within', truncation_level, &
                    'excitations of the reference determinant are included.'
                write (6,'(1X,a29,1X)',advance='no') 'Reference determinant, |D0> ='
                call write_det(d0%f, new_line=.true.)
            end if
            write (6,'(/,1X,a14,6X,a6)') 'Symmetry index','# dets'
            do i = lbound(sym_space_size,dim=1), ubound(sym_space_size, dim=1)
                write (6,'(6X,i4,4X,i13)') i, sym_space_size(i)
            end do
            write (6,'()')
        else if (.not. init .and. write_determinants .and. parent) then
            ! Output the determinants.
            fmt1 = int_fmt(ndets, padding=1)
            do i = 1, ndets
                write (det_unit,'('//fmt1//',4X)',advance='no') i
                call write_det(dets_list(:,i), det_unit, new_line=.true.)
            end do
        end if

        deallocate(occ, stat=ierr)
        call check_deallocate('occ', ierr)
        deallocate(nalpha_combinations, stat=ierr)
        call check_deallocate('nalpha_combinations', ierr)
        deallocate(comb, stat=ierr)
        call check_deallocate('comb', ierr)
        call dealloc_det_info(d0)

    end subroutine enumerate_determinants

!--- Obtain the next basis function string ---

    subroutine next_string(init, alpha, norb_spin, nel_spin, occ_spin, unocc_spin, max_level, comb, string, excit_level)

        ! Generate the next string/list of orbitals occupied by electrons of
        ! a given spin.

        ! If there is no truncation of the Hilbert space, then all combinations
        ! of orbitals are generated.  If there is a truncation of the Hilbert
        ! space, then only combinations within the required number of
        ! excitations from a reference determinant are generated.

        ! In:
        !    init: set to true if the first call to next_string, false
        !        otherwise.  If true then comb is initialised appropriately.
        !    alpha: set true to produce a string of alpha orbitals; false for
        !        a string of beta orbitals.  Used only if max_level == nel_spin
        !        (ie no truncation of the space), otherwise orbitals are
        !        selected from occ_spin and unocc_spin.
        !    norb_spin: number of orbitals of the same spin.
        !    nel_spin: number of electrons of the same spin.
        !    occ_spin: list of occupied orbitals of the given spin in the
        !        reference determinant.  If max_level < nel_spin, this is used
        !        to select orbitals a spin string has in common with the
        !        reference.
        !    unocc_spin: list of unoccupied orbitals of the given spin in the
        !        reference determinant.  If max_level < nel_spin, this is used
        !        to select the orbitals a spin string does not have in common
        !        with the reference.
        !    max_level: maximum number of excitations from the reference to
        !        consider.  max_level == nel_spin implies no truncation of the
        !        space.
        ! In/Out:
        !    comb: internal state used to generate the next combination.  Do not
        !        change or set outside this routine.
        !    excit_level: number of excitations by which string and the
        !        reference differ.  Set to 0 if max_level == nel_spin and on
        !        initialisation otherwise.  Updated with the excitation level
        !        when generating a truncated space.  (Note that next_string
        !        generates all strings of a given excitation level before going
        !        to the next excitation level.)  Do not change or set outside
        !        this routine.
        ! Out:
        !    string: unique list of occupied spin-orbitals of the given spin.

        use errors, only: stop_all
        use utils, only: next_comb

        logical, intent(in) :: init, alpha
        integer, intent(in) :: norb_spin, nel_spin, max_level, occ_spin(:), unocc_spin(:)
        integer, intent(inout) :: comb(max(nel_spin, 2*max_level)), excit_level
        integer, intent(out) :: string(nel_spin)

        integer :: i, ierr

        if (max_level >= nel_spin) then
            ! Very simple: want all possible combinations.  Find a combination
            ! and hence create a string from it.
            excit_level = 0
            if (init) then
                forall(i=1:nel_spin) comb(i) = i-1
            else
                call next_comb(norb_spin, nel_spin, comb(1:nel_spin), ierr)
                if (ierr /= 0) &
                    call stop_all('next_string', 'Too many calls.  No combinations remaining.')
            end if
            ! Convert to spin orbitals.
            ! alpha (up) spin functions are odd.
            ! beta (down) spin functions are even.
            ! Elements in comb are in range [0,norb_spin-1]
            if (alpha) then
                string = 2*comb(:nel_spin) + 1
            else
                string = 2*(comb(:nel_spin)+ 1)
            end if
        else
            if (init) then
                ! Return reference as the first string.
                string = occ_spin
                excit_level = 0
                ! Prepare combination for first single excitation.
                ! CARE: this uses an undocumented feature of next_comb, where if
                ! the input comb is (-1), then next_comb returns (0).
                comb(1:2) = (/ 0, -1 /)
            else
                ! Must be doing (at least) single excitations now.
                if (excit_level == 0) excit_level = 1
                ! Select combination of valence orbitals to excite into.
                call next_comb(norb_spin-nel_spin, excit_level, comb(excit_level+1:2*excit_level), ierr)
                if (ierr == 1) then
                    ! Select next combination of core orbitals to excite from.
                    call next_comb(nel_spin, excit_level, comb(1:excit_level), ierr)
                    if (ierr == 0) then
                        ! Everything ok.
                    else if (excit_level == max_level) then
                        call stop_all('next_string', 'Too many calls.  No combinations remaining.')
                    else
                        ! No more combinations at this excitation level.  Go to the
                        ! next excitation level.
                        excit_level = excit_level + 1
                        forall (i=1:excit_level) comb(i) = i - 1
                    end if
                    ! Reset combination of the virtuals.
                    forall (i=1:excit_level) comb(i+excit_level) = i - 1
                end if
                ! Form string by exciting electrons from selected core orbitals
                ! into selected virtual orbitals.
                string = occ_spin
                forall (i=1:excit_level) string(comb(i)+1) = unocc_spin(comb(i+excit_level)+1)
            end if
        end if

    end subroutine next_string

end module determinant_enumeration
