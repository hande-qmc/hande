module excit_gen_ringium

! Module for excitation generator for ringium.

! As the basis functions used are angular momentum eigenfunctions, the allowed excitations must
! conserve angular momentum, so only double excitations are allowed.

use const

implicit none

contains

    subroutine gen_excit_ringium_no_renorm(rng, sys, pattempt_single, cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate the probability of
        ! selecting that excitation.

        ! In:
        !    sys: system object being studied.
        !    pattempt_single: Probability to attempt a single excitation.
        !        Unused but included for interface compatibility.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        ! Out:
        ! In/Out:
        !    rng: random number generator.
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant 

        use determinants, only: det_info_t
        use excitations, only: excit_t, find_excitation_permutation2
        use system, only: sys_t
        use hamiltonian_ringium, only: slater_condon2_ringium_excit
        use dSFMT_interface, only: dSFMT_t
        use excit_gen_ueg, only: choose_ij_k

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pattempt_single
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen, hmatel
        type(excit_t), intent(out) :: connection

        logical :: allowed
        integer :: ij_lz(1)
        integer :: ij_spin

        connection%nexcit = 2
        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_k(rng, sys, cdet%occ_list, connection%from_orb(1), connection%from_orb(2), ij_lz, ij_spin)

        ! 2. Select a random pair of spin orbitals to excite to.
        call choose_ab_ringium(rng, sys, cdet%f, ij_lz(1), connection%to_orb(1), connection%to_orb(2), allowed)

        if (allowed) then
            ! 3. Calculate the generation probability of the excitation.
            ! For one-band systems this depends only upon the orbitals excited from.
            pgen = calc_pgen_ringium(sys)

            call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

            ! 4. find the connecting matrix element.
            hmatel = slater_condon2_ringium_excit(sys, connection%from_orb(1), connection%from_orb(2), &
                        connection%to_orb(1), connection%to_orb(2), connection%perm)
        else
            pgen = 1.0_p
            hmatel = 0.0_p
        end if

    end subroutine gen_excit_ringium_no_renorm

    subroutine choose_ab_ringium(rng, sys, f, ij_lz, a, b, allowed)

        ! Choose a random pair of orbitals (a,b) into which electrons are excited.
        ! (a,b) are chosen such that the excitation (i,j)->(a,b) conserves angular momentum.

        ! In:
        !   sys: system being studied.
        !   f: bit string representation of Slater determinant.
        !   ij_lz: sum of the angular momenta of the i and j orbitals.
        ! In/Out:
        !   rng: random number generator.
        ! Out:
        !   a, b: Indices of the virtual orbitals selected for the excitation.
        !   allowed: if true, then the excitation is allowed.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer, intent(in) :: ij_lz
        integer, intent(out) :: a, b
        logical, intent(out) :: allowed

        integer :: lz_b, ind

        ! Select a random unoccupied orbital as a
        ! This does not take account of the fact that the required b orbital to conserve angular
        ! momentum might not be in the basis set - in that case we just reject the excitation.
        do
            a = 2 * int(get_rand_close_open(rng)*sys%basis%nbasis/2) + 1
            ! If a is unoccupied, then found orbital to excite to
            if (.not.btest(f(sys%basis%bit_lookup(2,a)), sys%basis%bit_lookup(1,a))) exit
        end do

        ! Once three orbitals of the excitation have been chosen, this determines the fourth.

        lz_b = ij_lz - sys%basis%basis_fns(a)%l(1)

        if (abs(lz_b) > sys%ringium%maxlz) then
            allowed = .false.
        else
            ! As all electrons are spin up, the basis function index for given lz > 0
            ! is 2lz+1 or 2lz-1 for -lz
            b = 2*abs(lz_b)
            if (lz_b >= 0) then
                b = b + 1
            else
                b = b -1
            end if
            ! Allowed if b is unoccupied (and not a)
            allowed = .not.btest(f(sys%basis%bit_lookup(2,b)), sys%basis%bit_lookup(1,b))
            if (a == b) allowed = .false.
        end if
        ! Ensure a < b (as assumed by other procedures)
        if (a > b) then
            ind = a
            a = b
            b = ind
        end if

    end subroutine choose_ab_ringium

    pure function calc_pgen_ringium(sys) result(pgen)

        ! Calculate the probability of an excitation.  Due to the (naive)
        ! algorithm used, all allowed excitations have the same probability.

        ! In:
        !   sys: system being studied
        ! Returns:
        !   probability of generating this excitation.

        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        real(p) :: pgen

        ! i and j are chosen uniformly from the nel*(nel-1)/2 possible pairs
        ! a is chosen randomly from the nspatial - nel unoccupied orbitals
        ! b is then fixed
        ! a and b could have been chosen in either order so
        ! p_gen = p(a|i,j) p(b|i,j,a) + p(b|i,j) p(a|i,j,b)
        !       = 4/(nel*(nel-1)*(nbasis/2-nel))

        pgen = 4.0/(sys%nel*(sys%nel-1)*(sys%basis%nbasis/2-sys%nel))

    end function calc_pgen_ringium

end module excit_gen_ringium
