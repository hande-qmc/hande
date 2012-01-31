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

        ! This finds the exact size of the space.  See estimate_hilbert_space
        ! for a Monte Carlo approach to estimating the size of the space (better
        ! for large systems where we can't do FCI).

        use checking, only: check_allocate, check_deallocate
        use utils, only: binom_i, int_fmt

        use basis, only: basis_length
        use calc, only: truncate_space, truncation_level
        use excitations, only: get_excitation_level
        use reference_determinant
        use system, only: sym0, sym_max, nel, system_type, hub_real, heisenberg, ueg
        use symmetry, only: cross_product, symmetry_orb_list
        use ueg_system, only: ueg_basis_index

        use bit_utils, only: first_perm, bit_permutation, decode_bit_string

        integer :: i, j, iel, ierr
        integer :: nalpha_combinations, nbeta_combinations
        integer :: sym_beta, sym
        integer(i0) :: f_alpha, f_beta, f0(basis_length,sym0:sym_max), f(basis_length)
        integer, allocatable :: occ(:)
        integer :: k(ndim), k_beta(ndim)

        if (allocated(sym_space_size)) then
            deallocate(sym_space_size, stat=ierr)
            call check_deallocate('sym_space_size',ierr)
        end if
        allocate(sym_space_size(sym0:sym_max), stat=ierr)
        call check_allocate('sym_space_size',nsym,ierr)

        nbeta_combinations = binom_i(nbasis/2, nbeta)
        nalpha_combinations = binom_i(nbasis/2, nalpha)

        allocate(occ(nel), stat=ierr)
        call check_allocate('occ', nel, ierr)

        select case(system_type)

        case(heisenberg)

            ! See notes in system about how the Heisenberg model uses nel and
            ! nvirt.
            if (truncate_space) then
                sym_space_size = binom_i(nsites-(nel-truncation_level),truncation_level)
            else
                sym_space_size = binom_i(nsites, nel)
            end if

        case default

            if (system_type == hub_real .and. .not.truncate_space) then
                sym_space_size = nalpha_combinations*nbeta_combinations
            else
                ! Simply count the number of allowed determinants.
                ! CBA to find the closed form for the Hubbard model with local
                ! orbitals (ie no symmetry constraints).

                ! Find reference det if required.
                if (truncate_space) then
                    do sym = sym0, sym_max
                        call set_reference_det(occ, .true., sym)
                        call encode_det(occ, f0(:,sym))
                    end do
                end if

                ! Determinants are assigned a given symmetry by the direct product
                ! of the representations spanned by the occupied orbitals.  For
                ! momentum space systems this amounts to the sum of the wavevectors
                ! of the occupied basis functions.  This is because only doubly
                ! excitations are connected and D and D_{ij}^{ab} are only connected
                ! if k_i + k_j - k_a - k_b is a reciprocal lattice vector.  Thus we
                ! can regard the sum of the wavevectors of the occupied
                ! spin-orbitals as a symmetry label.

                sym_space_size = 0

                do i = 1, nbeta_combinations

                    ! Get beta orbitals.
                    if (i == 1) then
                        f_beta = first_perm(nbeta)
                    else
                        f_beta = bit_permutation(f_beta)
                    end if

                    call decode_bit_string(f_beta, occ)
                    ! Convert to beta orbitals.
                    occ = 2*(occ+1)
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

                    do j = 1, nalpha_combinations

                        ! Get alpha orbitals.
                        if (j == 1) then
                            f_alpha = first_perm(nalpha)
                        else
                            f_alpha = bit_permutation(f_alpha)
                        end if

                        call decode_bit_string(f_alpha, occ(nbeta+1:nel))
                        ! Convert to alpha orbitals.
                        occ(nbeta+1:nel) = 2*occ(nbeta+1:nel)+1
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

                        if (sym >= lbound(sym_space_size,dim=1) .and. sym <= ubound(sym_space_size,dim=1)) then
                            ! Ignore symmetries outside the basis set (only affects
                            ! UEG currently).
                            if (truncate_space) then
                                ! Check it's within truncated Hilbert space.
                                call encode_det(occ, f)
                                if (get_excitation_level(f,f0(:,sym)) <= truncation_level) then
                                    sym_space_size(sym) = sym_space_size(sym) + 1
                                end if
                            else
                                sym_space_size(sym) = sym_space_size(sym) + 1
                            end if
                        end if

                    end do
                end do
            end if

        end select

        if (parent) then
            write (6,'(1X,a25,/,1X,25("-"),/)') 'Size of determinant space'
            write (6,'(1X,a75,'//int_fmt(dets_Ms,0)//',a1,/)') &
                     'The table below gives the number of determinants for each symmetry with Ms=', &
                     dets_Ms,"."
            write (6,'(1X,a14,6X,a6)') 'Symmetry index','# dets'
            do i = lbound(sym_space_size,dim=1), ubound(sym_space_size, dim=1)
                write (6,'(6X,i4,4X,i13)') i, sym_space_size(i)
            end do
            write (6,'()')
        end if

        deallocate(occ, stat=ierr)
        call check_deallocate('occ', ierr)

    end subroutine find_sym_space_size

    subroutine enumerate_determinants(ref_sym)

        ! Find the Slater determinants that can be formed from the
        ! basis functions.  The list of determinants is stored in the
        ! module level dets_list array.
        ! find_sym_space_size must be called first for each value of Ms.
        ! In:
        !   ref_sym: index of an irreducible representation.  Only determinants
        !   with the same symmetry  are stored.

        use checking, only: check_allocate, check_deallocate
        use utils, only: binom_i
        use errors, only: stop_all
        use utils, only: get_free_unit, int_fmt
        use bit_utils, only: first_perm, bit_permutation, decode_bit_string

        use calc, only: truncate_space, truncation_level
        use excitations, only: get_excitation_level
        use reference_determinant
        use symmetry, only: cross_product, symmetry_orb_list
        use ueg_system, only: ueg_basis_index

        integer, intent(in) :: ref_sym

        integer :: i, j, iel, idet, ierr
        integer :: nalpha_combinations, nbeta_combinations
        integer :: sym_beta, sym
        character(4) :: fmt1
        integer(i0) :: f_alpha, f_beta, f0(basis_length), f(basis_length)
        integer, allocatable :: occ(:)
        integer :: k(ndim), k_beta(ndim)

        allocate(occ(nel), stat=ierr)
        call check_allocate('occ', nel, ierr)

        if (allocated(dets_list)) then
            deallocate(dets_list, stat=ierr)
            call check_deallocate('dets_list',ierr)
        end if

        nbeta_combinations = binom_i(nbasis/2, nbeta)
        nalpha_combinations = binom_i(nbasis/2, nalpha)

        ndets = sym_space_size(ref_sym)

        allocate(dets_list(basis_length, ndets), stat=ierr)
        call check_allocate('dets_list',basis_length*ndets,ierr)

        ! Find reference det if required.
        if (truncate_space) then
            call set_reference_det(occ, .true., ref_sym)
            call encode_det(occ, f0)
        end if

        select case(system_type)

        case(heisenberg)

            ! Just have spin up and spin down sites (no concept of unoccupied
            ! sites) so just need to arrange nel spin ups.  No need to
            ! interweave alpha and beta strings.

            ! Assume that we're not attempting to do FCI for more than
            ! a i0_length sites , which is quite large... ;-)
            if (nbasis > i0_length) then
                call stop_all('enumerate_determinants','Number of spin functions longer than the an i0 integer.')
            end if

            if (truncate_space) then
                f = first_perm(nel)
                i = 0
                if (get_excitation_level(f,f0) <= truncation_level) then
                    i = i + 1
                    dets_list(:,i) = f
                end if
                do while(i<ndets)
                    f = bit_permutation(f(1))
                    if (get_excitation_level(f,f0) <= truncation_level) then
                        i = i + 1
                        dets_list(:,i) = f
                    end if
                end do
            else
                dets_list(1,1) = first_perm(nel)
                do i = 2, ndets
                    dets_list(1,i) = bit_permutation(dets_list(1,i-1))
                end do
            end if

        case default

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

                call decode_bit_string(f_beta, occ)
                ! Convert to beta orbitals.
                occ = 2*(occ+1)
                ! Symmetry of the beta orbitals.
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

                do j = 1, nalpha_combinations

                    ! Get alpha orbitals.
                    if (j == 1) then
                        f_alpha = first_perm(nalpha)
                    else
                        f_alpha = bit_permutation(f_alpha)
                    end if

                    call decode_bit_string(f_alpha, occ(nbeta+1:nel))
                    ! Convert to alpha orbitals.
                    occ(nbeta+1:nel) = 2*occ(nbeta+1:nel)+1
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

                    if (sym == ref_sym) then

                        if (truncate_space) then
                            call encode_det(occ, f)
                            if (get_excitation_level(f,f0) <= truncation_level) then
                                ! Within truncation level.
                                idet = idet + 1
                                if (separate_strings) then
                                    ! Merge alpha and beta sets into determinant list.
                                    dets_list(:,idet) = (/ f_alpha, f_beta /)
                                else
                                    dets_list(:,idet) = f
                                end if
                            end if
                        else
                            idet = idet + 1
                            if (separate_strings) then
                                ! Merge alpha and beta sets into determinant list.
                                dets_list(:,idet) = (/ f_alpha, f_beta /)
                            else
                                call encode_det(occ, f)
                                dets_list(:,idet) = f
                            end if
                        end if

                    end if

                end do
            end do

        end select

        dets_sym = ref_sym

        if (write_determinants .and. parent) then
            fmt1 = int_fmt(ndets, padding=1)
            do i = 1, ndets
                write (det_unit,'('//fmt1//',4X)',advance='no') i
                call write_det(dets_list(:,i), det_unit, new_line=.true.)
            end do
        end if

        deallocate(occ, stat=ierr)
        call check_deallocate('occ', ierr)

    end subroutine enumerate_determinants

end module determinant_enumeration
