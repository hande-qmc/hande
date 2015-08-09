module canonical_energy_estimates

! Estimate the various energy estimates in the canonical ensemble.

use const

implicit none

enum, bind(c)
    ! Ensure energy estimates are always first.
    ! Total energy estimate from single-particle eigenvalues. Here H_0 = \sum_i epsilon_i n_i, where
    ! epsilon_i is a single particle energy eigenvalue (e.g. the kinetic energy
    ! for the ueg or the Hartree-Fock orbital energy for molecular systems) and
    ! n_i is the number operator for orbital i.
    enumerator :: ke_idx = 1
    ! Hartree-Fock energy in non-interacting ensemble. Here H_0 is again H_0 = \sum_i epsilon_i n_i,
    ! as above but we now evalate Tr (H \rho_0) / Tr (\rho_0).
    enumerator :: hf_idx
    ! Finite-T Hartree-Fock (sort of). Here we evaluate Tr (H \rho_HF) /
    ! Tr (\rho_HF), where \rho_HF = \sum_I e^{-\beta E_HF^I}, where E_HF^I are
    ! now many-particle Hartree-Fock energies of a given many-particle state I.
    enumerator :: hft_idx
    ! Finite-T reweighted estimate for partition function.
    enumerator :: hf_part_idx
    ! Index for checking for interaction with the calculation
    enumerator :: comms_found_idx
    ! last_idx-1 gives number of estimates.
    enumerator :: last_idx
end enum


contains

    subroutine estimate_canonical_energy(sys, fermi_temperature, beta, nsamples, ncycles, all_spin_sectors, rng_seed)

        ! From the Fermi factors calculated in the grand canonical ensemble we can
        ! estimate the total energy in the canonical ensemble by generating determinants
        ! with arbitrary particle number and only keeping those with <N> = nel.

        ! In:
        !    beta: target temperature.
        !    fermi_temperature: if true, rescale beta as the inverse reduced temperature:
        !        beta = 1/\Theta = T_F/T.  If false, then beta is in atomic units.
        !    nsamples: number of samples to use each cycle
        !    ncycles: number of Monte Carlo cycles to perform, over which the kinetic energy is
        !        estimated, along with an estimate of the standard error.
        !    rng_seed (optional): seed to initialise the random number generator.
        !       Default: seed based upon the hash of the time and calculation UUID.
        !    all_spin_sectors: if true then average over ms otherwise only
        !       provide estimates for specified spin sector.
        ! In/Out:
        !    sys: system being studied.

        use system
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use parallel
        use utils, only: rng_init_info

        use calc, only: GLOBAL_META, gen_seed
        use hamiltonian_ueg, only: exchange_energy_ueg
        use hamiltonian_molecular, only: double_counting_correction_mol
        use determinants, only: sum_sp_eigenvalues
        use interact, only: calc_interact, check_comms_file
        use proc_pointers, only: energy_diff_ptr
        use errors, only: stop_all
        use reference_determinant, only: set_reference_det

        type(sys_t), intent(inout) :: sys
        logical, intent(in) :: all_spin_sectors
        real(p), intent(in) :: beta
        logical, intent(in) :: fermi_temperature
        integer, intent(in) :: nsamples
        integer, intent(in), optional :: ncycles
        integer, intent(in), optional :: rng_seed

        real(dp) :: p_single(sys%basis%nbasis/2)
        real(dp) :: r
        integer :: occ_list(sys%nel), seed
        logical :: gen
        real(p) :: energy(hf_part_idx), beta_loc, hfx
        integer :: ierr, ireport, iorb
        integer(int_64) :: iaccept
        real(p) :: local_estimators(last_idx-1), estimators(last_idx-1)

        type(sys_t) :: sys_bak
        type (dSFMT_t) :: rng
        logical :: soft_exit, comms_found
        integer :: ngen, nalpha_allowed, nbeta_allowed
        real(p) :: energy_zero, ref_shift
        integer, allocatable :: occ_list0(:)

        if (parent) write (6,'(1X,a16,/,1X,16("-"),/)') 'Canonical energy'

        if (present(rng_seed)) then
            seed = rng_seed
        else
            seed = gen_seed(GLOBAL_META%uuid)
        end if

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)
        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, sys)

        beta_loc = beta
        if (fermi_temperature) then
            beta_loc = beta_loc / sys%ueg%ef
        end if

        if (sys%symmetry < sys%sym_max) then
            call set_reference_det(sys, occ_list0, .false., sys%symmetry)
        else
            call set_reference_det(sys, occ_list0, .false.)
        end if

        select case(sys%system)
        case (ueg)
            energy_diff_ptr => exchange_energy_ueg
            energy_zero = 0.0_p
        case (read_in)
            energy_diff_ptr => double_counting_correction_mol
            energy_zero = sys%read_in%Ecore
        case default
            call stop_all('estimate_canonical_energy', 'Not implemented for selected model Hamiltonian')
        end select
        ! Add a shift to the exponentials of the Boltzman factors to prevent
        ! numerical problems at lower temperatures.
        ref_shift = energy_diff_ptr(sys, occ_list0)

        if (parent) then
            write (6,'(1X,a90)') 'E_0: Estimate for total energy from single-particle Hamiltonian i.e. 1/Z_0 Tr(\rho_0 H_0).'
            write (6,'(1X,a66)') 'E_HF0: Estimate for Hartree-Fock-0 energy i.e. 1/Z_0 Tr(\rho_0 H).'
            write (6,'(1X,a91)') '\sum\rho_HF_{ii}H_{ii}: Estimate for numerator of "Hartree-Fock" energy i.e. Tr(\rho_HF H).'
            write (6,'(1X,a77)') '\sum\rho_HF_{ii}: Estimate for denominator of "Hatree-Fock" energy i.e. Z_HF.'
            write (6,'()')
        end if

        if (parent) write (6,'(1X,a12,19X,a3,17X,a5,4x,a22,6X,a16)') &
                    '# iterations', 'E_0', 'E_HF0', '\sum\rho_HF_{ii}H_{ii}', '\sum\rho_HF_{ii}'

        forall (iorb=1:sys%basis%nbasis:2) p_single(iorb/2+1) = 1.0_p / &
                                                          (1+exp(beta_loc*(sys%basis%basis_fns(iorb)%sp_eigv-sys%chem_pot)))

        if (all_spin_sectors) then
            ! If averaging over spin we need to allow for the number of
            ! alpha/beta orbitals to fluctuate between [0,sys%nel].
            nalpha_allowed = sys%nel
            nbeta_allowed = sys%nel
        else
            nalpha_allowed = sys%nalpha
            nbeta_allowed = sys%nbeta
        end if

        do ireport = 1, ncycles
            local_estimators = 0.0_p
            iaccept = 0 ! running number of samples this report cycle.
            do while (iaccept < nsamples)
                ngen = 0
                occ_list = 0
                if (nalpha_allowed > 0) call generate_allowed_orbital_list(sys, rng, p_single, nalpha_allowed, &
                                                                           1, occ_list, ngen)
                if (.not. all_spin_sectors .and. ngen /= sys%nalpha) then
                    cycle
                else if (ngen > nalpha_allowed) then
                    cycle
                end if
                if (nbeta_allowed > 0) call generate_allowed_orbital_list(sys, rng, p_single, sys%nel-ngen, &
                                                                          0, occ_list, ngen)
                if (ngen /= sys%nel) cycle
                iaccept = iaccept + 1
                ! Calculate Kinetic and Hartree-Fock exchange energies.
                energy(ke_idx) = sum_sp_eigenvalues(sys, occ_list)
                hfx = energy_diff_ptr(sys, occ_list)
                energy(hf_idx) = energy(ke_idx) + hfx + energy_zero
                ! We generate determinants with probability p(i1,..,iN) =
                ! 1/Z_0 \prod_{i} p(i1)X...Xp(iN), where Z_0 is the
                ! non-interacting canonical partition function, and p(i1) =
                ! e^{-\beta \varepsilon_i1}. We can instead
                ! calculate Z_HF = \sum_{i} e^{-\beta E_HF(i)} by reweighting,
                ! i.e.,
                ! p(i1,..,iN)_HF = 1/Z' e^{-beta(E_HF(i)-E_0(i))}p(i1,...,iN),
                ! where Z' = \sum_{i} e^{-\beta(E_HF(i)-E_0(i))}.
                energy(hf_part_idx) = exp(-beta_loc*(hfx-ref_shift))
                energy(hft_idx) = energy(hf_part_idx)*energy(hf_idx)
                local_estimators(ke_idx:hf_part_idx) = local_estimators(ke_idx:hf_part_idx) + energy
            end do

            if (check_comms_file()) local_estimators(comms_found_idx) = 1.0_p

#ifdef PARALLEL
            ! More efficient in parallel.
            call mpi_allreduce(local_estimators, estimators, last_idx-1, mpi_preal, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
            estimators = local_estimators
#endif
            ! Average over processors.
            estimators(ke_idx:hf_part_idx) = estimators(ke_idx:hf_part_idx) / (nprocs*nsamples)
            if (parent) write(6,'(3X,i10,5X,2(es17.10,5X),4X,2(es17.10,5X))') ireport, estimators(ke_idx), &
                                                             estimators(hf_idx), estimators(hft_idx), estimators(hf_part_idx)
            comms_found = abs(estimators(comms_found_idx)) > depsilon
            call calc_interact(comms_found, soft_exit)
            if (soft_exit) exit
        end do

        if (parent) write(6, '()')

        ! Return sys in unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine estimate_canonical_energy

    subroutine generate_allowed_orbital_list(sys, rng, porb, nselect, spin_factor, occ_list, ngen)

        ! Generate a list of orbitals according to their single
        ! particle GC orbital occupancy probabilities.

        ! In:
        !    sys: system being studied.
        !    porb: porb(i) gives the probabilty of selecting
        !        the orbital i.
        !    nselect: number of orbitals to select.
        !    spin_factor: integer to account for odd/even ordering of
        !        alpha/beta spin orbitals. Set to 1 for alpha spins, 0 for beta spins.
        ! In/Out:
        !    rng: random number generator.
        !    ngen: running total number of electrons calculated. Should be zerod
        !       upon first call to generate_allowed_orbital_list.
        ! Out:
        !    occ_list: array containing occupied orbitals.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        real(dp), intent(in) :: porb(:)
        integer, intent(in) :: nselect
        integer, intent(in) :: spin_factor
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: occ_list(:)
        integer, intent(inout) :: ngen

        integer :: iorb, iselect
        real(dp) :: r

        iselect = 0

        do iorb = 1, sys%basis%nbasis/2
            ! Select a random orbital.
            r = get_rand_close_open(rng)
            if (porb(iorb) > r) then
                iselect = iselect + 1
                ngen = ngen + 1
                if (iselect > nselect) then
                    ! Selected too many.
                    exit
                end if
                occ_list(ngen) = 2*iorb - spin_factor
            end if
        end do

    end subroutine generate_allowed_orbital_list

end module canonical_energy_estimates
