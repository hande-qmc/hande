module importance_sampling

! Module for applying importance sampling transforms to the Hamiltonian, e.g. to
! stochastically sample \tilde{H}_{ij} = \psi^{(T})_i H_{ij} 1/\psi^{(T})_j
! rather than just H_{ij}.

! Instead of propogating using

! dc_i
! ---- = - \sum_j H_{ij} c_j
!  dt

! we instead propogate using

! df_i      \sum_j \psi^{(T)}_i H_{ij}      1        f_j
! ---- = -                             ------------
!  dt                                  \psi^{(T)}_j 

! where

! f_i = \psi^{(T)}_i c_i

! and \psi^{(T)}_i is the weight of the trial function on i.

! As we currently generate the excitations uniformly, this amounts to simply
! scaling H_{ij} appropriately after generating the excitation (as usual) but
! before testing whether or not to accept the spawning attempt (which depends
! upon the connecting Hamiltonian matrix element, which has changed).

! This module contains the appropriate functions required to transform the
! Hamiltonian matrix element according to the trial function being used.

! Note that the projected estimator must also be adjusted accordingly.

! In order to be accessible via a common function pointer, these subroutines
! should have the interface defined by i_trial_fn.

use const

implicit none

contains

    pure function importance_sampling_weight(trial, cdet, pop) result(weight)

        ! Calculate c_i from f_i.

        ! In:
        !    trial: guiding wavefunction in use, containing both wavefunction type and (if relevant)
        !       associated data.  See importance_sampling_data for supported wavefunctions.
        !    cdet: det_info_t object representing the i-th state (tensor product, determinant, etc.).
        !    pop: current population on i, giving a stochastic representation of f_i.
        ! Returns:
        !    stochastic representation of c_i.

        use determinants, only: det_info_t
        use importance_sampling_data

        real(p) :: weight
        type(trial_t), intent(in) :: trial
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop
        integer :: n

        select case(trial%guide)
        case(single_basis)
            weight = pop
        case(neel_singlet_guiding)
            n = nint(cdet%data(size(cdet%data)-1))
            weight = pop/trial%wfn_dat(n)
        end select

    end function importance_sampling_weight

    subroutine neel_trial_state(sys, cdet, connection, trial_func, hmatel)

        ! Apply the transformation to the Hamiltonian matrix element due to
        ! using the Neel singlet state as the trial function.

        ! In:
        !    sys: system being studied.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !    trial_func: importance sampling weights.
        ! In/Out:
        !    hmatel: on input, untransformed matrix element connecting two spin
        !        functions (kets).  On output, transformed matrix element,
        !        \Psi_i^T H_{ij} 1/\Psi_j^T.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        real(p), allocatable, intent(in) :: trial_func(:)
        real(p), intent(inout) :: hmatel

        integer :: up_spins_to, up_spins_from
        integer :: bit_position, bit_element

        ! Find the number of up spins on sublattice 1.
        ! WARNING: we assume the number of up spins is the penultimate element in cdet%data.
        up_spins_from = nint(cdet%data(size(cdet%data)-1))
        ! For the spin up which was flipped to create the connected
        ! basis function, find whether this spin was on sublattice 1 or 2.
        ! If it was on sublattice 1, the basis function we go to has 1 less
        ! up spin on sublattice 1, else it will have one more spin up here.
        bit_position = sys%basis%bit_lookup(1,connection%from_orb(1))
        bit_element = sys%basis%bit_lookup(2,connection%from_orb(1))
        if (btest(sys%heisenberg%lattice_mask(bit_element), bit_position)) then
            up_spins_to = up_spins_from-1
        else
            up_spins_to = up_spins_from+1
        end if

        ! For a given number of spins up on sublattice 1, n, the corresponding
        ! ampltidue of this basis function in the trial function is stored as
        ! trial_func(n), for this particular trial function. Hence we have:
        hmatel = (trial_func(up_spins_to)*hmatel)/trial_func(up_spins_from)

    end subroutine neel_trial_state

    subroutine dmqmc_weighting_fn(sys, cdet, connection, trial_func, hmatel)

        ! Apply a transformation to the Hamiltonian matrix element by
        ! reducing the probability of spawning to higher excitation levels
        ! between the two ends of the DMQMC bitstring. The exact trial function
        ! used is specified by the users upon input, and stored in the vector
        ! accumulated_probs.

        ! In:
        !    sys: system being studied.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !    trial_func: importance sampling weights.
        ! In/Out:
        !    hmatel: on input, untransformed matrix element connecting two spin
        !        functions (kets).  On output, transformed matrix element,
        !        \Psi^T(cdet%f,cdet%f2) H_{ij} 1/\Psi(f_new,cdet%f2).
        !        The factors which the Hamiltonian are multiplied by depend
        !        on the level which we come from, and go to, and so depend on
        !        the two ends of the bitstring we spawn from, and the new
        !        bitstring we spawn onto. Note this is different to more conventional
        !        importance sampling.

        use determinants, only: det_info_t
        use excitations, only: excit_t, get_excitation_level, create_excited_det
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        real(p), allocatable, intent(in) :: trial_func(:)
        real(p), intent(inout) :: hmatel
        integer(i0) :: f_new(sys%basis%string_len)
        integer :: excit_level_old, excit_level_new

        excit_level_old = get_excitation_level(cdet%f,cdet%f2)

        call create_excited_det(sys%basis, cdet%f, connection, f_new)

        excit_level_new = get_excitation_level(f_new,cdet%f2)

        hmatel = trial_func(excit_level_old)*hmatel*(1/trial_func(excit_level_new))

    end subroutine dmqmc_weighting_fn

    subroutine interaction_picture_reweighting_free(sys, cdet, connection, trial_func, hmatel)

        ! Apply transformation to the Hamiltonian matrix element so that
        ! excitation generations are performed with the weighting factors
        ! e^{-0.5*(beta-tau)(E_i-E_k)} H_{ik} = H_I(beta-tau) = e{-0.5(beta-tau) H^0} H e^{0.5(beta-tau) H^0}.
        ! The trial function is essentially meaningless, but we abuse its
        ! meaning to pass in beta-tau. Here E_i = \sum_{occ} \varepsilon_i,
        ! where \varepsilon_i are the single-particle eigenvalues.

        ! In:
        !    sys: system being studied.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !    trial_func: importance sampling weights, used to pass 0.5*(beta-tau)
        !        in the last element (sys%max_number_excitations+1).
        ! In/Out:
        !    hmatel: on input, untransformed matrix element connecting two spin
        !        functions (kets).  On output, transformed matrix element,
        !        is e^{-0.5*(beta-tau)(E_i-E_k)} H_{ik}.
        !        The factors which the Hamiltonian are multiplied by depend
        !        on the level which we come from, and go to, and so depend on
        !        the two ends of the bitstring we spawn from, and the new
        !        bitstring we spawn onto. Note this is not really importance sampling
        !        in the usual sense.

        use determinants, only: det_info_t
        use excitations, only: excit_t, get_excitation_level, create_excited_det
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        real(p), allocatable, intent(in) :: trial_func(:)
        real(p), intent(inout) :: hmatel

        integer :: iorb
        integer(i0) :: f_new(sys%basis%string_len)
        real(p) :: diff_ijab

        diff_ijab = 0.0_p

        call create_excited_det(sys%basis, cdet%f, connection, f_new)

        ! For the case when H^0 = \sum_i \varepsilon_i n_i, then the energy differences are simple to
        ! evaluate : E_i - E_k = (\varepsilon_a+\varepsilon_b) - (\varepsilon_i+\varepsilon_j).
        do iorb = 1, connection%nexcit
            diff_ijab = diff_ijab + sys%basis%basis_fns(connection%to_orb(iorb))%sp_eigv - &
                                    sys%basis%basis_fns(connection%from_orb(iorb))%sp_eigv
        end do
        hmatel = exp(-trial_func(sys%max_number_excitations+1)*diff_ijab) * hmatel

    end subroutine interaction_picture_reweighting_free

    subroutine interaction_picture_reweighting_hartree_fock(sys, cdet, connection, trial_func, hmatel)

        ! Apply transformation to the Hamiltonian matrix element so that
        ! excitation generations are performed with the weighting factors
        ! e^{-0.5*(beta-tau)(E_i-E_k)} H_{ik}. The trial function is essentially
        ! meaningless, but we abuse its meaning to pass in beta-tau. Here E_i = <D_i| H |D_i>.

        ! In:
        !    sys: system being studied.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !    trial_func: importance sampling weights, used to pass 0.5*(beta-tau)
        !        in the last element (sys%max_number_excitations+1).
        ! In/Out:
        !    hmatel: on input, untransformed matrix element connecting two spin
        !        functions (kets).  On output, transformed matrix element,
        !        is e^{-0.5*(beta-tau)(E_i-E_k)} H_{ik}.
        !        The factors which the Hamiltonian are multiplied by depend
        !        on the level which we come from, and go to, and so depend on
        !        the two ends of the bitstring we spawn from, and the new
        !        bitstring we spawn onto. Note this is not really importance sampling
        !        in the usual sense.

        use determinants, only: det_info_t
        use excitations, only: excit_t, get_excitation_level, create_excited_det
        use system, only: sys_t
        use proc_pointers, only: trial_dm_ptr
        use hamiltonian_ueg, only: exchange_energy_orb

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        real(p), allocatable, intent(in) :: trial_func(:)
        real(p), intent(inout) :: hmatel

        integer :: iorb, new, occ_list(sys%nel+2), a, b, i, j
        integer(i0) :: f_new(sys%basis%string_len)
        real(p) :: E_i, E_k, diff_ijab

        if (abs(hmatel) > depsilon) then
            call create_excited_det(sys%basis, cdet%f, connection, f_new)
            a = connection%to_orb(1)
            b = connection%to_orb(2)
            i = connection%from_orb(1)
            j = connection%from_orb(2)
            ! Order occ_list so that a and b orbitals appear first for convenience later.
            occ_list(1) = i
            occ_list(2) =  j
            new = 2
            do iorb = 1, sys%nel
                if (cdet%occ_list(iorb) /= i .and. cdet%occ_list(iorb) /= j) then
                    new = new + 1
                    occ_list(new) = cdet%occ_list(iorb)
                end if
            end do
            ! Set last two entries to be new orbitals.
            occ_list(sys%nel+1) = a
            occ_list(sys%nel+2) = b

            diff_ijab = 0.0_p

            ! Work out <D|H|D> - <D_{ij}^{ab}|H|D_{ij}^{ab}> as
            ! E_i - E_k = (\varepsilon_a+ex_a+\varepsilon_b+ex_b) -
            ! (\varepsilon_i+ex_i+\varepsilon_j+ex_j) - double counting.
            ! Where ex_i is the Hartree-Fock exchange potential of orbital i
            ! .i.e. ex_i ~ -\sum_{occ != i} 1/|k_occ-k_i|^2 (for the UEG).
            do iorb = 1, 2
                diff_ijab = diff_ijab - &
                    & (sys%basis%basis_fns(connection%from_orb(iorb))%sp_eigv + &
                    & exchange_energy_orb(sys, occ_list(:sys%nel), connection%from_orb(iorb))) + &
                    & (sys%basis%basis_fns(connection%to_orb(iorb))%sp_eigv + &
                    & exchange_energy_orb(sys, occ_list(3:), connection%to_orb(iorb)))
            end do
            ! Double counted the ab, and ij exchange terms so remove these explicitly.
            if (mod(a,2) == mod(b,2)) then
            diff_ijab = diff_ijab + &
                        & sys%ueg%exchange_int(sys%lattice%box_length(1), sys%basis, connection%to_orb(1), connection%to_orb(2))
            end if
            if (mod(i,2) == mod(j,2)) then
                diff_ijab = diff_ijab - &
                        & sys%ueg%exchange_int(sys%lattice%box_length(1), sys%basis, connection%from_orb(1), connection%from_orb(2))
            end if
            hmatel = exp(-trial_func(sys%max_number_excitations+1)*diff_ijab) * hmatel
        end if

    end subroutine interaction_picture_reweighting_hartree_fock

end module importance_sampling
