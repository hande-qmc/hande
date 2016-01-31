module hilbert_space

! Estimating the size of the Hilbert space.

implicit none

integer :: nhilbert_cycles

private
public :: nhilbert_cycles, estimate_hilbert_space, gen_random_det_full_space, gen_random_det_truncate_space

contains

    subroutine estimate_hilbert_space(sys, ex_level, nattempts, ncycles, occ_list0, rng_seed)

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
        ! In:
        !    ex_level: maximum excitation level relative to the reference
        !       determinant to include in the Hilbert space.  If negative or
        !       greater than the number of electrons, the entire Hilbert space
        !       is considered.
        !    nattempts: number of attempts, i.e. number of random determinants to
        !       generate, per Monte Carlo cycle.
        !    ncycles: number of blocks of nattempts to perform.  Statistics are
        !       estimated based upon the space size estimate from each block.
        !    occ_list0: reference determinant.  If not allocated, then a best
        !       guess is generated based upon the spin and symmetry quantum
        !       numbers.
        !    rng_seed: seed to initialise the random number generator.

        use basis, only: write_basis_fn
        use const, only: dp
        use checking, only: check_allocate, check_deallocate
        use determinants, only: encode_det, det_info_t, alloc_det_info_t,  &
                                dealloc_det_info_t, decode_det_spinocc_spinunocc
        use dSFMT_interface, only: dSFMT_t, dSFMT_init, dSFMT_end
        use json_out
        use reference_determinant, only: set_reference_det
        use symmetry, only: symmetry_orb_list
        use system
        use parallel
        use utils, only: binom_r

        type(sys_t), intent(inout) :: sys
        integer, intent(in) :: ex_level, nattempts
        integer, intent(inout), allocatable :: occ_list0(:)
        integer, intent(in) :: ncycles
        integer, intent(in), optional :: rng_seed

        integer :: truncation_level, icycle, i, ierr, a, n, ireport, ireport_ind, nattempts_local
        integer :: ref_sym, det_sym
        integer(i0) :: f(sys%basis%string_len), f0(sys%basis%string_len)
        integer :: occ_list(sys%nel)
        integer, parameter :: nreport_block = 100
        real(dp) :: full_space_size, space_size_mean, space_size_mean2, space_size_se, delta
        real(dp) :: naccept(nreport_block), naccept_sum(nreport_block)
        real(dp), allocatable :: ptrunc_level(:,:)

        type(dSFMT_t) :: rng
        type(sys_t) :: sys_bak
        type(det_info_t) :: det0
        logical :: truncate_space
        type(json_out_t) :: js

        if (parent) write (6,'(1X,a13,/,1X,13("-"),/)') 'Hilbert space'

        truncate_space = ex_level >= 0 .and. ex_level <= sys%nel
        if (truncate_space) then
            truncation_level = ex_level
        else
            truncation_level = sys%nel
        end if

        call dSFMT_init(rng_seed+iproc, 50000, rng)
        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, sys)

        if (parent) then
            call json_object_init(js, tag=.true.)
            call sys_t_json(js, sys)
            call json_write_key(js, 'ex_level', truncation_level)
            call json_write_key(js, 'nattempts', nattempts)
            call json_write_key(js, 'ncycles', ncycles)
            call json_write_key(js, 'occ_list', occ_list0)
            call json_write_key(js, 'rng_seed', rng_seed, .true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io,'()')
        end if

        select case(sys%system)

        case(heisenberg)

            ! Symmetry not currently implemented for the Heisenberg code.
            ! There is one spin per site, so it's just a case of how many ways
            ! there are to arrange the nalpha spins across the lattice%nsites (or
            ! equivalently the nbeta spins across the lattice%nsites).
            ! See comments in system for how nel and nvirt are used in the
            ! Heisenberg model.
            if (truncate_space) then
                full_space_size = binom_r(sys%lattice%nsites-(sys%nel-truncation_level),truncation_level)
            else
                full_space_size = binom_r(sys%lattice%nsites, sys%nel)
            end if
            if (parent) write (6,'(1X,a,g12.4,/)') 'Size of space is', full_space_size

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

                if (sys%symmetry < sys%sym_max) then
                    call set_reference_det(sys, occ_list0, .false., sys%symmetry)
                else
                    call set_reference_det(sys, occ_list0, .false.)
                end if
                call encode_det(sys%basis, occ_list0, f0)

                ! Symmetry of the reference determinant.
                ref_sym = symmetry_orb_list(sys, occ_list0)

                if (truncate_space) then
                    ! Generate a determinant with a given excitation level up to
                    ! some maximum level (ie truncation_level) and number of
                    ! alpha excitations using a probability proportional to the
                    ! number of determinants fulfilling those criteria relative
                    ! to the size of the truncated Hilbert space (excluding the
                    ! reference).
                    allocate(ptrunc_level(0:min(sys%nalpha, truncation_level), truncation_level), stat=ierr)
                    call check_allocate('ptrunc_level', truncation_level*(min(sys%nalpha, truncation_level)+1), ierr)
                    ptrunc_level = 0.0_dp
                    full_space_size = 0.0_dp
                    do i = 1, truncation_level
                        ! Number of possible determinants at excitation level i is (whilst
                        ! conserving the spin of the reference):
                        !   \sum_a C(N_a,a) C(M-N_a,a) C(N_b, i-a) C(M-N_b,i-a)
                        ! ie the number of ways of choosing the electrons and the number of ways
                        ! of choosing the holes.  N_a (N_b) is the number of alpha (beta)
                        ! electrons and M is the number of alpha/beta spin orbitals.
                        do a = max(0,i-sys%nbeta), min(sys%nalpha,i)
                            full_space_size = full_space_size + &
                                binom_r(sys%nalpha,a)*binom_r(sys%nvirt_alpha,a) &
                                *binom_r(sys%nbeta,i-a)*binom_r(sys%nvirt_beta,i-a)
                            ptrunc_level(a,i) = full_space_size
                        end do
                    end do
                    ptrunc_level = ptrunc_level / full_space_size
                    ! Need a completely decoded representation of the reference for efficient excitations.
                    call alloc_det_info_t(sys, det0)
                    call decode_det_spinocc_spinunocc(sys, f0, det0)
                    det0%f = f0
                else
                    ! Size of the complete Hilbert space in the desired symmetry block is given
                    ! by
                    !   C(nalpha_orbitals, nalpha_electrons)*C(nbeta_orbitals, nbeta_electrons).
                    full_space_size = binom_r(sys%basis%nbasis/2,sys%nalpha) * binom_r(sys%basis%nbasis/2,sys%nbeta)
                end if

                if (parent) then
                    write (6,'(1X,a34)',advance='no') 'Symmetry of reference determinant:'
                    if (sys%momentum_space) then
                        call write_basis_fn(sys, sys%basis%basis_fns(2*ref_sym), new_line=.true., print_full=.false.)
                    else
                        write (6,'(1X,i2)') ref_sym
                    end if
                    write (6,'(/,1X,"space size: estimate of the Hilbert space size from a single iteration.")')
                    write (6,'(1X,"mean: running estimate of the mean of the Hilbert space size.")')
                    write (6,'(1X,"std. err.: running estimate of the standard error in the estimate of the mean.")')
                    write (6,'(/,1X,"# iterations  space size    mean          std. err.")')
                end if

                nattempts_local = nattempts/nprocs
                if (iproc < mod(nattempts,nprocs)) nattempts_local = nattempts_local + 1

                naccept = 0.0_dp
                space_size_mean = 0.0_dp
                space_size_mean2 = 0.0_dp
                do ireport = 1, ncycles
                    ireport_ind = mod(ireport, nreport_block)
                    if (ireport_ind == 0) ireport_ind = nreport_block
                    do icycle = 1, nattempts_local
                        ! Generate a random determinant.
                        if (truncate_space) then
                            ! More efficient sampling (especially if Hilbert space is large
                            ! and truncation level is low) if we just excite randomly from the
                            ! reference determinant.
                            call gen_random_det_truncate_space(rng, sys, truncation_level, det0, ptrunc_level, &
                                                               occ_list)
                        else
                            call gen_random_det_full_space(rng, sys, f, occ_list)
                        end if
                        det_sym = symmetry_orb_list(sys, occ_list)
                        ! Is this the same symmetry as the reference determinant?
                        if (det_sym == ref_sym) naccept(ireport_ind) = naccept(ireport_ind) + 1
                    end do
                    if (ireport == ncycles .or. ireport_ind == nreport_block) then
#ifdef PARALLEL
                        call MPI_Reduce(naccept, naccept_sum, nreport_block, MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr)
#else
                        naccept_sum = naccept
#endif
                        ! space_size is the max number of excitations from the reference.  Account
                        ! for the fraction of determinants with the correct symmetry.
                        naccept_sum = (full_space_size * naccept_sum) / nattempts
                        if (truncate_space) then
                            ! We never allowed for the generation of the reference determinant but, by
                            ! definition,it is in the truncated Hilbert space so explicitly include it.
                            naccept_sum = naccept_sum + 1
                        end if

                        ! On-fly mean and standard error estimation, as proposed by Knuth/Welford.
                        ! See https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
                        do i = 1, ireport_ind
                            n = ireport - ireport_ind + i
                            delta = naccept_sum(i) - space_size_mean
                            space_size_mean = space_size_mean + delta/n
                            space_size_mean2 = space_size_mean2 + delta*(naccept_sum(i) - space_size_mean)
                            space_size_se = sqrt(space_size_mean2 / (n*(n-1)))
                            if (parent) write (6,'(1X,i12,3(2X,es12.6))') n, naccept_sum(i), space_size_mean, space_size_se
                        end do
                        naccept = 0.0_p
                    end if
                end do

                if (parent) write (6,'(/,1X,"Monte-Carlo estimate of size of space is: ",es12.6," +/- ",es12.6,/)') &
                                    space_size_mean, space_size_se

                if (truncate_space) call dealloc_det_info_t(det0)

            end if

        end select

        if (allocated(ptrunc_level)) then
            deallocate(ptrunc_level, stat=ierr)
            call check_deallocate('ptrunc_level', ierr)
        end if

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

        call dSFMT_end(rng)

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

    subroutine gen_random_det_truncate_space(rng, sys, truncation_level, det0, ptrunc_level, occ_list)

        ! Generate a random excited determinant with in a truncated Hilbert space defined
        ! by a maximum number of excitations permitted from a single reference determinant.

        ! Note:
        ! * we never generate the reference determinant.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.  Unaltered on output.
        !    truncation_level: max number of excitations permitted from the reference
        !        determinant.
        !    det0: det_info_t object for the reference determinant.  Must contain 
        !        appropriately set occ_list, occ_list_alpha, occ_list_beta,
        !        unocc_list_alpha and unocc_list_beta components.
        !    ptrunc_level: cumulative probability of generating a determinant at a given
        !        excitation level containing a given number of alpha
        !        excitations, i.e. ptrunc_level(na,i) is the probability of
        !        generating a determinant at excitation level i involving up to
        !        na alpha excitations or generating a determinant at a lower
        !        excitation level.
        ! Out:
        !    occ_list: list of occupied orbital of random determinant.  Must have
        !        dimensions of (at least) sys%nel.

        ! NOTE: the returned occ_list is *not* ordered.  Many procedures assume determinant
        ! lists are ordered by orbital index, so occ_list must be sorted before use in
        ! such routines.

        use const, only: dp
        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: truncation_level
        type(det_info_t), intent(in) :: det0
        real(dp), intent(in) :: ptrunc_level(0:,:)
        integer, intent(out) :: occ_list(:)

        integer :: ilevel, ialpha
        real(dp) :: plevel
        integer :: unocc_alpha(sys%nvirt_alpha), unocc_beta(sys%nvirt_beta)

        occ_list(:sys%nalpha) = det0%occ_list_alpha(:sys%nalpha)
        occ_list(sys%nalpha+1:) = det0%occ_list_beta(:sys%nbeta)
        unocc_alpha = det0%unocc_list_alpha(:sys%nvirt_alpha)
        unocc_beta = det0%unocc_list_beta(:sys%nvirt_beta)

        ! Select a given sector of the space (i.e. excitation level and number
        ! of alpha excitations) according to the size of the sector.
        plevel = get_rand_close_open(rng)
        outer: do ilevel = 1, truncation_level
            do ialpha = 0, min(ilevel, sys%nalpha)
                if (plevel < ptrunc_level(ialpha, ilevel)) exit outer
            end do
        end do outer

        ! Create a random excitation uniformly within the sector.
        call excit_spin_channel(rng, sys%nalpha, sys%nvirt_alpha, ialpha, occ_list(:sys%nalpha), unocc_alpha)
        call excit_spin_channel(rng, sys%nbeta, sys%nvirt_beta, ilevel-ialpha, occ_list(sys%nalpha+1:), unocc_beta)

        contains

            subroutine excit_spin_channel(rng, nspin, nvirt_spin, nexcit, occ_spin, unocc_spin)

                ! Excite a number of electrons of a given spin to virtual
                ! orbitals.

                ! In:
                !    nspin: number of electrons of the given spin.
                !    nvirt_spin: number of virtual orbitals of the given spin.
                !    nexcit: number of electrons to excit_t.
                ! In/Out:
                !    rng: as above.
                !    occ_spin: on input, list of occupied orbitals of the given
                !       spin from a reference determinant.  On output, nexcit
                !       orbitals have been switched with unocc_spin.
                !    unocc_spin: on input, list of occupied orbitals of the given
                !       spin (given occ_spin).  On output, nexcit orbitals have
                !       been switched with occ_spin.

                use errors, only: stop_all
                use dSFMT_interface, only: dSFMT_t, get_rand_close_open

                type(dSFMT_t), intent(inout) :: rng
                integer, intent(in) :: nspin, nvirt_spin, nexcit
                integer, intent(inout) :: occ_spin(:), unocc_spin(:)

                integer :: a, a_pos, tmp, nspin_sel, i

                do i = 1, nexcit
                    ! Select an orbital and move it to the end (where it will be
                    ! replaced with the excited orbital) so we don't attempt to
                    ! pick it again.
                    a_pos = int(get_rand_close_open(rng)*(nspin-(i-1)))+1
                    a = occ_spin(a_pos)
                    tmp = occ_spin(nspin-(i-1))
                    occ_spin(nspin-(i-1)) = a
                    occ_spin(a_pos) = tmp
                    ! and excited into, explicitly conserving spin.
                    ! last entry in virtual array which has not yet been selected.
                    nspin_sel = nvirt_spin - (i-1)
                    if (nspin_sel < 1) &
                        call stop_all('gen_random_det_truncate_space', 'Cannot handle such a spin-constrained system')
                    a_pos = int(get_rand_close_open(rng)*nspin_sel)+1
                    a = unocc_spin(a_pos)
                    occ_spin(nspin-(i-1)) = a
                    tmp = unocc_spin(nspin_sel)
                    unocc_spin(nspin_sel) = a
                    unocc_spin(a_pos) = tmp
                end do

            end subroutine excit_spin_channel

    end subroutine gen_random_det_truncate_space

end module hilbert_space
