module excit_gen_mol_complex

! Module for random excitation generators and related routines for the molecular
! system (ie one read in from an FCIDUMP file) with complex integral values. Most
! related routines are reused from excit_gen_mol to avoid duplication.

! See top-level comments in spawning about the overall aim of the spawning step.

use const, only: i0, p

implicit none

contains

!--- Excitation generation ---

    subroutine gen_excit_mol_complex(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element. This
        ! function is specifically for use with systems with complex coefficients
        ! and hamiltonian elements.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Data for the excitation generator.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant in molecular systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation
        !        was not generated

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2
        use excit_gens, only: excit_gen_data_t
        use hamiltonian_molecular_complex, only: slater_condon1_mol_excit_complex, &
                                                slater_condon2_mol_excit_complex
        use system, only: sys_t

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use excit_gen_mol, only: choose_ia_mol, choose_ij_mol, choose_ab_mol, &
                                calc_pgen_single_mol, calc_pgen_double_mol

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        complex(p), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation

        integer :: ij_sym, ij_spin

        ! 1. Select single or double.

        if (get_rand_close_open(rng) < excit_gen_data%pattempt_single) then

            ! 2a. Select orbital to excite from and orbital to excite into.
            call choose_ia_mol(rng, sys, sys%read_in%pg_sym%gamma_sym, cdet%f, cdet%occ_list, cdet%symunocc, &
                               connection%from_orb(1), connection%to_orb(1), allowed_excitation)
            connection%nexcit = 1

            if (allowed_excitation) then
                ! 3a. Probability of generating this excitation.
                pgen = excit_gen_data%pattempt_single*calc_pgen_single_mol(sys, sys%read_in%pg_sym%gamma_sym, cdet%occ_list, &
                                                                   cdet%symunocc, connection%to_orb(1))

                ! 4a. Parity of permutation required to line up determinants.
                call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)

                ! 5a. Find the connecting matrix element.
                hmatel = slater_condon1_mol_excit_complex(sys, cdet%occ_list, connection%from_orb(1), connection%to_orb(1), &
                                                  connection%perm)
            else
                ! We have a highly restrained system and this det has no single
                ! excitations at all.  To avoid reweighting pattempt_single and
                ! pattempt_double (an O(N^3) operation), we simply return a null
                ! excitation
                hmatel = cmplx(0.0_p, 0.0_p, p)
                pgen = 1.0_p
            end if

        else

            ! 2b. Select orbitals to excite from and orbitals to excite into.
            call choose_ij_mol(rng, sys, cdet%occ_list, connection%from_orb(1), connection%from_orb(2), ij_sym, ij_spin)
            call choose_ab_mol(rng, sys, cdet%f, ij_sym, ij_spin, cdet%symunocc, connection%to_orb(1), &
                               connection%to_orb(2), allowed_excitation)
            connection%nexcit = 2

            if (allowed_excitation) then

                ! 3b. Probability of generating this excitation.
                pgen = excit_gen_data%pattempt_double*calc_pgen_double_mol(sys, ij_sym, connection%to_orb(1), &
                                                                                   connection%to_orb(2), ij_spin, cdet%symunocc)

                ! 4b. Parity of permutation required to line up determinants.
                ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

                ! 5b. Find the connecting matrix element.
                hmatel = slater_condon2_mol_excit_complex(sys, connection%from_orb(1), connection%from_orb(2), &
                                                  connection%to_orb(1), connection%to_orb(2), connection%perm)
            else
                ! Carelessly selected ij with no possible excitations.  Such
                ! events are not worth the cost of renormalising the generation
                ! probabilities.
                ! Return a null excitation.
                hmatel = cmplx(0.0_p, 0.0_p, p)
                pgen = 1.0_p
            end if

        end if

    end subroutine gen_excit_mol_complex

end module excit_gen_mol_complex
