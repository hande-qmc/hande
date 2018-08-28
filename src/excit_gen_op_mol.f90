module excit_gen_op_mol

! Module for random excitation generators and related routines for operators
! other than the Hamiltonian for generic systems, where the integrals are read
! in (i.e. molecular systems).

! See excit_gen_mol for more information about standard excitation generators.
! See hfs_fciqmc for more information about how to sample arbitrary operators
! within FCIQMC.

use const

implicit none

contains

!=== One-body operators ===

! If \hat{O_1} is a one-body operator, then we can largely re-use the single
! excitation generators.

! This is a hack until the excitation generators support generating excitations
! for operators of arbitrary symmetries.

    subroutine gen_excit_one_body_mol(rng, sys, excit_gen_data, cdet, pgen, connection, matel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the corresponding matrix element of
        ! a one body operator.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Data for excitation generator (not used) 
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    rng: random number generator.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    matel: < D | O_1 | D' >, the matrix element of a one-body operator
        !        between a determinant and a connected determinant in molecular
        !        systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinant_data, only: det_info_t
        use excitations, only: excit_t, find_excitation_permutation1
        use excit_gens, only: excit_gen_data_t
        use excit_gen_mol, only: choose_ia_mol, calc_pgen_single_mol
        use operators, only: one_body1_mol_excit
        use system, only: sys_t
        use hamiltonian_data

        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(inout) :: cdet
        real(p), intent(out) :: pgen
        type(excit_t), intent(out) :: connection
        type(hmatel_t), intent(out) :: matel
        logical, intent(out) :: allowed_excitation

        integer :: op_sym

        op_sym = sys%read_in%one_body_op_integrals%op_sym

        ! 1. Select orbital to excite from and orbital to excite into.
        call choose_ia_mol(rng, sys, op_sym, cdet%f, cdet%occ_list, cdet%symunocc, connection%from_orb(1), &
                           connection%to_orb(1), allowed_excitation)
        connection%nexcit = 1

        if (allowed_excitation) then
            ! 2. Probability of generating this excitation.
            pgen = calc_pgen_single_mol(sys, op_sym, cdet%occ_list, cdet%symunocc, connection%to_orb(1))

            ! 3. Parity of permutation required to line up determinants.
            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)

            ! 4. Find the connecting matrix element.
            matel%r = one_body1_mol_excit(sys, connection%from_orb(1), connection%to_orb(1), connection%perm)
        else
            ! We have a highly restrained system and this det has no single
            ! excitations at all.  To avoid reweighting pattempt_single and
            ! pattempt_double (an O(N^3) operation), we simply return a null
            ! excitation
            matel%r = 0.0_p
            pgen = 1.0_p
        end if

    end subroutine gen_excit_one_body_mol

    subroutine gen_excit_one_body_mol_no_renorm(rng, sys, excit_gen_data, cdet, pgen, connection, matel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the corresponding matrix element of
        ! a one body operator.

        ! This doesn't exclude the case where, having selected all orbitals
        ! involved in the excitation, the final orbital selected is already
        ! occupied and so cannot be excited into.  Whilst this is somewhat
        ! wasteful (generating excitations which can't be performed), there is
        ! a balance between the cost of generating forbidden excitations and the
        ! O(N) cost of renormalising the generation probabilities.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Data for excitation generator (not used) 
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    rng: random number generator.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    matel: < D | O_1 | D' >, the matrix element of a one-body operator
        !        between a determinant and a connected determinant in molecular
        !        systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinant_data, only: det_info_t
        use excitations, only: excit_t, find_excitation_permutation1
        use excit_gens, only: excit_gen_data_t
        use excit_gen_mol, only: find_ia_mol, calc_pgen_single_mol_no_renorm
        use operators, only: one_body1_mol_excit
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t
        use hamiltonian_data

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(inout) :: cdet
        real(p), intent(out) :: pgen
        type(excit_t), intent(out) :: connection
        type(hmatel_t), intent(out) :: matel
        logical, intent(out) :: allowed_excitation

        integer :: op_sym

        op_sym = sys%read_in%one_body_op_integrals%op_sym

        ! 1. Select orbital to excite from and orbital to excite into.
        call find_ia_mol(rng, sys, op_sym, cdet%f, cdet%occ_list, connection%from_orb(1), &
                         connection%to_orb(1), allowed_excitation)
        connection%nexcit = 1

        if (allowed_excitation) then
            ! 2. Probability of generating this excitation.
            pgen = calc_pgen_single_mol_no_renorm(sys, connection%to_orb(1))

            ! 3. Parity of permutation required to line up determinants.
            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)

            ! 4. Find the connecting matrix element.
            matel%r = one_body1_mol_excit(sys, connection%from_orb(1), connection%to_orb(1), connection%perm)
        else
            ! We have a highly restrained system and this det has no single
            ! excitations at all.  To avoid reweighting pattempt_single and
            ! pattempt_double (an O(N^3) operation), we simply return a null

            ! excitation
            matel%r = 0.0_p
            pgen = 1.0_p
        end if

    end subroutine gen_excit_one_body_mol_no_renorm

end module excit_gen_op_mol
