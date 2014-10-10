module hilbert_space

! Estimating the size of the Hilbert space.

implicit none

integer :: nhilbert_cycles

private
public :: nhilbert_cycles, estimate_hilbert_space

contains

    subroutine estimate_hilbert_space(sys)

        ! Based on Appendix A in George Booth's thesis.

        ! Find the size of the Hilbert space which is of the same symmetry as
        ! the reference determinant.

        ! The idea is to randomly generate a determinant within a Hilbert space
        ! and then test whether the determinant is within the desired symmetry
        ! subspace.  The fraction of determinants within the subspace hence provides
        ! an estimate for the size of that subspace without requiring an explicit
        ! enumeration.

        ! See find_sym_space_size for a dumb but exact enumeration of the size
        ! of the space (which is needed for FCI calculations).

        ! In/Out:
        !    sys: system being studied.  Unaltered on output.

        use basis, only: write_basis_fn
        use calc, only: sym_in, ms_in, truncate_space, truncation_level, seed
        use const, only: dp
        use checking, only: check_allocate, check_deallocate
        use determinants, only: encode_det, det_info_t, alloc_det_info_t,  &
                                dealloc_det_info_t, decode_det_occ_spinunocc
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use fciqmc_data, only: occ_list0
        use reference_determinant, only: set_reference_det
        use symmetry, only: symmetry_orb_list
        use system
        use parallel
        use utils, only: binom_r, rng_init_info

        type(sys_t), intent(inout) :: sys

        integer :: icycle, i, ierr, b
        integer :: ref_sym, det_sym
        integer(i0) :: f(sys%basis%string_len), f0(sys%basis%string_len)
        integer :: occ_list(sys%nel)
        real(dp) :: space_size, naccept, nsuccess, weight
        real(dp), allocatable :: ptrunc_level(:)
#ifdef PARALLEL
        real(dp) :: proc_space_size(nprocs), sd_space_size
#endif

        type(dSFMT_t) :: rng
        type(sys_t) :: sys_bak
        type(det_info_t) :: det0

        if (parent) write (6,'(1X,a13,/,1X,13("-"),/)') 'Hilbert space'

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        select case(sys%system)

        case(heisenberg)

            ! Symmetry not currently implemented for the Heisenberg code.
            ! There is one spin per site, so it's just a case of how many ways
            ! there are to arrange the nalpha spins across the lattice%nsites (or
            ! equivalently the nbeta spins across the lattice%nsites).
            ! See comments in system for how nel and nvirt are used in the
            ! Heisenberg model.
            if (truncate_space) then
                space_size = binom_r(sys%lattice%nsites-(sys%nel-truncation_level),truncation_level)
            else
                space_size = binom_r(sys%lattice%nsites, sys%nel)
            end if
            if (parent) write (6,'(1X,a,g12.4,/)') 'Size of space is', space_size

        case default

            if ((sys%system == hub_real .or. sys%system == chung_landau) &
                                                .and. .not.truncate_space) then
                ! Symmetry not currently implemented for the real space lattice
                ! code.
                ! Just a case of how we arrange the alpha and beta electrons across
                ! the alpha orbitals and beta orbitals.  As we're dealing with the
                ! simplest possible lattice model, the number of orbitals of each
                ! spin is equal to the number of sites.
                associate(sg=>sys, sl=>sys%lattice)
                    if (parent) write (6,'(1X,a,g12.4,/)') 'Size of space is', &
                                    binom_r(sl%nsites, sg%nalpha)*binom_r(sl%nsites, sg%nbeta)
                end associate
            else

                ! Perform a Monte Carlo sampling of the space.

                if (sym_in < sys%sym_max) then
                    call set_reference_det(sys, occ_list0, .false., sym_in)
                else
                    call set_reference_det(sys, occ_list0, .false.)
                end if
                call encode_det(sys%basis, occ_list0, f0)

                ! Symmetry of the reference determinant.
                ref_sym = symmetry_orb_list(sys, occ_list0)

                if (truncate_space) then
                    ! Generate a determinant with a given truncation level up to
                    ! some maximum level (ie truncation_level) using a probability
                    ! proportional to the number of determinants at that truncation level
                    ! in the entire space relative to the size of the truncated Hilbert
                    ! space (excluding the reference).
                    allocate(ptrunc_level(truncation_level), stat=ierr)
                    call check_allocate('ptrunc_level', size(ptrunc_level), ierr)
                    space_size = 0
                    do i = 1, truncation_level
                        ! Number of possible determinants at excitation level i is (whilst
                        ! conserving the spin of the reference):
                        !   \sum_b C(N_b,b) C(M-N_b,b) C(N_a, i-b) C(M-N_a,i-b)
                        ! ie the number of ways of choosing the electrons and the number of ways
                        ! of choosing the holes.  N_a (N_b) is the number of alpha (beta)
                        ! electrons and M is the number of alpha/beta spin orbitals.
                        do b = max(0,i-sys%nalpha), min(sys%nbeta,i)
                            space_size = space_size &
                                             + (binom_r(sys%nbeta, b)*binom_r(sys%basis%nbasis/2-sys%nbeta,b) &
                                             *binom_r(sys%nalpha, i-b)*binom_r(sys%basis%nbasis/2-sys%nalpha,i-b))
                        end do
                        ptrunc_level(i) = space_size
                    end do
                    ptrunc_level = ptrunc_level / space_size
                    ! Need a completely decoded representation of the reference for efficient excitations.
                    call alloc_det_info_t(sys, det0)
                    call decode_det_occ_spinunocc(sys, f0, det0)
                    det0%f = f0
                else
                    ! Size of the complete Hilbert space in the desired symmetry block is given
                    ! by
                    !   C(nalpha_orbitals, nalpha_electrons)*C(nbeta_orbitals, nbeta_electrons).
                    space_size = binom_r(sys%basis%nbasis/2,sys%nalpha) * binom_r(sys%basis%nbasis/2,sys%nbeta)
                end if

                if (parent) then
                    write (6,'(1X,a34)',advance='no') 'Symmetry of reference determinant:'
                    if (sys%momentum_space) then
                        call write_basis_fn(sys, sys%basis%basis_fns(2*ref_sym), new_line=.true., print_full=.false.)
                    else
                        write (6,'(1X,i2)') ref_sym
                    end if
                end if

                naccept = 0
                nsuccess = 0
                do icycle = 1, nhilbert_cycles
                    ! Generate a random determinant.
                    if (truncate_space) then
                        ! More efficient sampling (especially if Hilbert space is large
                        ! and truncation level is low) if we just excite randomly from the
                        ! reference determinant.
                        call gen_random_det_truncate_space(rng, sys, truncation_level, det0, ptrunc_level, &
                                                           occ_list, weight)
                    else
                        call gen_random_det_full_space(rng, sys, f, occ_list)
                        ! gen_random_det_full_space always randomly generates
                        ! all determinants with uniform probability.
                        weight = 1
                    end if
                    ! The weight of the generated determinant is essentially how many
                    ! times it needs to be included such that it occurs with the same
                    ! probability as all other determinants.
                    nsuccess = nsuccess + weight
                    det_sym = symmetry_orb_list(sys, occ_list)
                    ! Is this the same symmetry as the reference determinant?
                    if (det_sym == ref_sym) naccept = naccept + weight
                end do

                ! space_size is the max number of excitations from the reference.  Account
                ! for the fraction of determinants with the correct symmetry.
                ! Note nsuccess == nhilbert_cycles if truncate_space is false.
                space_size = (space_size * naccept) / nsuccess
                if (truncate_space) then
                    ! We never allowed for the generation of the reference determinant
                    ! but, by definition, it is in the truncated Hilbert space so
                    ! explicitly include it.
                    space_size = space_size + 1
                end if

#ifdef PARALLEL
                ! If we did this on multiple processors then we can get an estimate
                ! of the error as well as a better mean...
                call mpi_gather(space_size, 1, mpi_real8, proc_space_size, 1, mpi_real8, root, mpi_comm_world, ierr)
                space_size = sum(proc_space_size)/nprocs
                sd_space_size = sqrt(sum((proc_space_size-space_size)**2))/(nprocs-1)
                if (parent) then
                    write (6,'(1X,a41,1X,es10.4)',advance='no') 'Monte-Carlo estimate of size of space is:', space_size
                    if (nprocs > 1) then
                        write (6,'(1X,a3,1X,es10.4)') '+/-', sd_space_size
                    else
                        write (6,'()')
                    end if
                end if
#else
                if (parent) write (6,'(1X,a41,1X,es10.4)') 'Monte-Carlo estimate of size of space is:', space_size
#endif

                if (parent) write (6,'()')

            end if

        end select

        if (truncate_space) then
            deallocate(ptrunc_level, stat=ierr)
            call check_deallocate('ptrunc_level', ierr)
        end if

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine estimate_hilbert_space

    subroutine gen_random_det_full_space(rng, sys, f, occ_list)

        ! Generate a random determinant with uniform probability in a (full) Hilbert space.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.  Unaltered on output.
        ! Out:
        !     f: bit string representation of random determinant.  Must have dimensions of
        !        (at least) string_len.
        !     occ_list: list of occupied orbital of random determinant.  Must have
        !         dimensions of (at least) sys%nel.

        use const, only: i0
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        integer(i0), intent(out) :: f(:)
        integer, intent(out) :: occ_list(:)

        integer :: iel, a, a_pos, a_el, b, b_pos, b_el

        ! Alpha electrons.
        f = 0_i0
        iel = 0
        if (sys%nalpha > 0) then
            do
                ! generate random number 1,3,5,...
                a = 2*int(get_rand_close_open(rng)*(sys%basis%nbasis/2))+1
                a_pos = sys%basis%bit_lookup(1,a)
                a_el = sys%basis%bit_lookup(2,a)
                if (.not.btest(f(a_el), a_pos)) then
                    ! found unoccupied alpha orbital.
                    f(a_el) = ibset(f(a_el), a_pos)
                    iel = iel + 1
                    occ_list(iel) = a
                    if (iel == sys%nalpha) exit
                end if
            end do
        end if
        ! Beta electrons.
        if (sys%nbeta > 0) then
            do
                ! generate random number 2,4,6,...
                b = 2*int(get_rand_close_open(rng)*(sys%basis%nbasis/2))+2
                b_pos = sys%basis%bit_lookup(1,b)
                b_el = sys%basis%bit_lookup(2,b)
                if (.not.btest(f(b_el), b_pos)) then
                    ! FOUND Unoccupied beta orbital.
                    f(b_el) = ibset(f(b_el), b_pos)
                    iel = iel + 1
                    occ_list(iel) = b
                    if (iel == sys%nel) exit
                end if
            end do
        end if

    end subroutine gen_random_det_full_space

    subroutine gen_random_det_truncate_space(rng, sys, truncation_level, det0, ptrunc_level, occ_list, weight)

        ! Generate a random excited determinant with in a truncated Hilbert space defined
        ! by a maximum number of excitations permitted from a single reference determinant.

        ! Note:
        ! * we never generate the reference determinant.
        ! * determinants are not generated uniformly within the same excitation level (see
        !   weight).

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.  Unaltered on output.
        !    truncation_level: max number of excitations permitted from the reference
        !        determinant.
        !    det0: det_info_t object for the reference determinant.  Must contain appropriately
        !        set occ_list, unocc_list_alpha and unocc_list_beta components.
        !    ptrunc_level: cumulative probability of generating a determinant at a given
        !        excitation level, i.e. ptrunc_level(i) is the probability of generating a
        !        determinant at or below excitation level i.  The probability of
        !        selecting a given excitation level must be proportional to the fraction
        !        of determinants in the space with that excitation level.
        ! Out:
        !    occ_list: list of occupied orbital of random determinant.  Must have
        !        dimensions of (at least) sys%nel.
        !    weight: the weight of the generated determinant relative to other
        !        determinants at the same excitation level.

        ! NOTE: the returned occ_list is *not* ordered.  Many procedures assume determinant
        ! lists are ordered by orbital index, so occ_list must be sorted before use in
        ! such routines.

        use errors, only: stop_all

        use const, only: i0, dp
        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t
        use utils, only: binom_r

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: truncation_level
        type(det_info_t), intent(in) :: det0
        real(dp), intent(in) :: ptrunc_level(:)
        integer, intent(out) :: occ_list(:)
        real(dp), intent(out) :: weight

        integer :: ilevel, nalpha_selected
        integer :: a, a_pos, tmp, nspin_sel
        real(dp) :: plevel
        integer :: unocc_alpha(sys%nvirt_alpha), unocc_beta(sys%nvirt_beta)

        occ_list = det0%occ_list
        unocc_alpha = det0%unocc_list_alpha
        unocc_beta = det0%unocc_list_beta
        plevel = get_rand_close_open(rng)
        nalpha_selected = 0
        do ilevel = 1, truncation_level
            ! Select an orbital and move it to the end (where it will be replaced with
            ! the excited orbital) so we don't attempt to pick it again.
            a_pos = int(get_rand_close_open(rng)*(sys%nel-(ilevel-1)))+1
            a = occ_list(a_pos)
            tmp = occ_list(sys%nel-(ilevel-1))
            occ_list(sys%nel-(ilevel-1)) = occ_list(a_pos)
            occ_list(a_pos) = tmp
            ! and exicted into, explicitly conserving spin.
            if (mod(a,2) == 1) then
                ! last entry in virtual array which has not yet been selected.
                nspin_sel = sys%nvirt_alpha - nalpha_selected
                if (nspin_sel < 1) &
                    call stop_all('gen_random_det_truncate_space', 'Cannot handle such a spin-constrained system')
                a_pos = int(get_rand_close_open(rng)*nspin_sel)+1
                a = unocc_alpha(a_pos)
                occ_list(sys%nel-(ilevel-1)) = a
                tmp = unocc_alpha(nspin_sel)
                unocc_alpha(nspin_sel) = a
                unocc_alpha(a_pos) = tmp
                nalpha_selected = nalpha_selected + 1
            else
                ! last entry in virtual array which has not yet been selected.
                nspin_sel = sys%nvirt_beta - (ilevel - nalpha_selected)
                if (nspin_sel < 1) &
                    call stop_all('gen_random_det_truncate_space', 'Cannot handle such a spin-constrained system')
                a_pos = int(get_rand_close_open(rng)*nspin_sel)+1
                a = unocc_beta(a_pos)
                occ_list(sys%nel-(ilevel-1)) = a
                tmp = unocc_beta(nspin_sel)
                unocc_beta(nspin_sel) = a
                unocc_beta(a_pos) = tmp
            end if
            if (plevel < ptrunc_level(ilevel)) exit
        end do

        ! We need to give equal value to each determinant.  We already take into account
        ! the different numbers of determinants at each excitation level by choosing to
        ! generate a determinant of a given level with probability proportional to the
        ! fraction of determinants at that excitation level in the Hilbert space.
        ! However, we need to take into account the ways a determinant could be generated
        ! within a given excitation level, e.g. if we excite 1,3 from the reference
        ! |1,2,3,4> to 5,7, we could have selected 1->5,3->7 or 1->7,3->5.  Hence
        ! |2,4,5,7> is twice as likely to come up as e.g. |1,2,5,6>, where the only the
        ! excitation 3->5,4->6 is possible as we conserve spin.
        ! For an arbitrary determinant at excitation level n, there are n_a! n_b! = n_a! (n-n_a)
        ! ways it could have been generated, where n_a (n_b) is the number of alpha (beta)
        ! electrons involved in the excitation.
        ! We need to divide the contribution through by this weight, such that all
        ! determinants within a given excitation level contribute equally.
        ! Weights within a given excitation level are relative and it is convenient if the
        ! weights are integer.  We thus choose an excitation involving all alpha (or all
        ! beta) electrons to have a weight of 1 (i.e. when n_a = n or n_b = n).  Hence by
        ! weighting all determinants by n!/(n_a! (n-n_a)!) we unbias for the
        ! non-uniformity in selection.
        weight = binom_r(ilevel, nalpha_selected)

    end subroutine gen_random_det_truncate_space

end module hilbert_space
