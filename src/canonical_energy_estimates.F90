module canonical_energy_estimates

! Estimate the various energy estimates in the canonical ensemble.
! Let H = T + V, then we can estimate thermodynamic values for <H>,
! where <X> = Tr(X\rho)/Tr(\rho), for the density matrix \rho using
! various level of approximations.
! Here we calculate <H> using two different density matrices, namely:
! 1) \rho_0 = \sum_I e^{-\beta \sum_{i}^{M} \varepsilon_i n_i} |D_I><D_I|,
! for many particle states |D_I>, and varepsilon_i is the single-particle orbital
! energy.
! 2) \rho_HF = \sum_I e^{-\beta E^HF_I} |D_I><D_I|, where E^HF_I = <D_I|H|D_I>.
! With these definitions we can the calculate the kinetic energy (T) and potential
! energies V) using a Monte Carlo sampling scheme.
! Note that for molecular systems we don't actually evaluate the kinetic energy
! directly but rather the sum of Hartree-Fock eigenvalues. The double counting
! correction is then added from the potential contribution.
! The folowing enumerators set the location in the estimates array.

use const

implicit none

enum, bind(c)
    ! Ensure energy estimates are always first.
    ! <T>_0.
    enumerator :: ke_idx = 1
    ! <V>_0.
    enumerator :: pe_idx
    ! <T>_HF.
    enumerator :: hf_ke_idx
    ! <V>_HF.
    enumerator :: hf_pe_idx
    ! Finite-T reweighted estimate for partition function = Tr(\rho_HF).
    enumerator :: hf_part_idx
    ! Index for checking for interaction with the calculation.
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
        !    rng_seed : seed to initialise the random number generator.
        !    all_spin_sectors: if true then average over ms otherwise only
        !       provide estimates for specified spin sector.
        ! In/Out:
        !    sys: system being studied.

        use system
        use dSFMT_interface, only: dSFMT_t, dSFMT_init, dSFMT_end
        use json_out
        use parallel

        use hamiltonian_ueg, only: exchange_energy_ueg
        use hamiltonian_molecular, only: double_counting_correction_mol
        use determinants, only: sum_sp_eigenvalues
        use interact, only: calc_interact, check_comms_file
        use proc_pointers, only: energy_diff_ptr
        use errors, only: stop_all
        use reference_determinant, only: set_reference_det
        use chem_pot, only: find_chem_pot

        type(sys_t), intent(inout) :: sys
        logical, intent(in) :: all_spin_sectors
        real(p), intent(in) :: beta
        logical, intent(in) :: fermi_temperature
        integer, intent(in) :: nsamples
        integer, intent(in) :: ncycles
        integer, intent(in) :: rng_seed

        real(dp) :: p_single(sys%basis%nbasis/2)
        integer :: occ_list(sys%nel)
        real(p) :: energy(hf_part_idx), beta_loc, mu
        integer :: ireport, iorb
        integer(int_64) :: iaccept
        real(p) :: local_estimators(last_idx-1), estimators(last_idx-1)

        type(sys_t) :: sys_bak
        type (dSFMT_t) :: rng
        logical :: soft_exit, comms_found
        integer :: ngen, nalpha_allowed, nbeta_allowed
        real(p) :: ref_shift
        integer, allocatable :: occ_list0(:)
        type(json_out_t) :: js
#ifdef PARALLEL
        integer :: ierr
#endif

        if (parent) write (6,'(1X,a16,/,1X,16("-"),/)') 'Canonical energy'

        call dSFMT_init(rng_seed+iproc, 50000, rng)
        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, sys)

        if (parent) then
            call json_object_init(js, tag=.true.)
            call sys_t_json(js, sys)
            call json_write_key(js, 'all_spin_sectors', all_spin_sectors)
            call json_write_key(js, 'beta', beta)
            call json_write_key(js, 'fermi_temperature', fermi_temperature)
            call json_write_key(js, 'nsamples', nsamples)
            call json_write_key(js, 'ncycles', ncycles)
            call json_write_key(js, 'rng_seed', rng_seed, terminal=.true.)
            call json_write_key(js, 'chem_pot', mu, terminal=.true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io,'()')
        end if

        beta_loc = beta
        if (fermi_temperature) then
            beta_loc = beta_loc / sys%ueg%ef
        end if
        mu = find_chem_pot(sys, beta_loc)

        if (sys%symmetry < sys%sym_max) then
            call set_reference_det(sys, occ_list0, .false., sys%symmetry)
        else
            call set_reference_det(sys, occ_list0, .false.)
        end if

        select case(sys%system)
        case (ueg)
            energy_diff_ptr => exchange_energy_ueg
        case (read_in)
            energy_diff_ptr => double_counting_correction_mol
        case default
            call stop_all('estimate_canonical_energy', 'Not implemented for selected model Hamiltonian')
        end select
        ! Add a shift to the exponentials of the Boltzman factors to prevent
        ! numerical problems at lower temperatures.
        ref_shift = energy_diff_ptr(sys, occ_list0)

        if (parent) then
            write (6,'(1X,a87)') '<T>_0: Estimate for kinetic energy in non-interacting ensemble i.e. 1/Z_0 Tr(\rho_0 T).'
            write (6,'(1X,a96)') '<V>_0: Estimate for potential energy in non-interacting ensemble energy i.e. 1/Z_0 Tr(\rho_0 V).'
            write (6,'(1X,a70)') 'Tr(T\rho_HF): Estimate for numerator of "Hartree-Fock" kinetic energy.'
            write (6,'(1X,a72)') 'Tr(V\rho_HF): Estimate for numerator of "Hartree-Fock" potential energy.'
            write (6,'(1X,a72)') 'Tr(\rho_HF): Estimate for denominator of "Hatree-Fock" energy i.e. Z_HF.'
            write (6,'()')
        end if

        if (parent) write (6,'(1X,a12,17X,a5,17X,a5,14x,a12,10X,a12,11X,a11)') &
                    '# iterations', '<T>_0', '<V>_0', 'Tr(T\rho_HF)', 'Tr(V\rho_HF)', 'Tr(\rho_HF)'

        forall (iorb=1:sys%basis%nbasis:2) p_single(iorb/2+1) = 1.0_p / &
                                                          (1+exp(beta_loc*(sys%basis%basis_fns(iorb)%sp_eigv-mu)))

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
                ! Add in the core contribution here for molecular systems.
                energy(pe_idx) = energy_diff_ptr(sys, occ_list)
                ! We generate determinants with probability p(i1,..,iN) =
                ! 1/Z_0 \prod_{i} p(i1)X...Xp(iN), where Z_0 is the
                ! non-interacting canonical partition function, and p(i1) =
                ! e^{-\beta \varepsilon_i1}. We can instead
                ! calculate Z_HF = \sum_{i} e^{-\beta E_HF(i)} by reweighting,
                ! i.e.,
                ! p(i1,..,iN)_HF = 1/Z' e^{-beta(E_HF(i)-E_0(i))}p(i1,...,iN),
                ! where Z' = \sum_{i} e^{-\beta(E_HF(i)-E_0(i))}.
                energy(hf_part_idx) = exp(-beta_loc*(energy(pe_idx)-ref_shift))
                energy(hf_ke_idx:hf_pe_idx) = energy(hf_part_idx)*energy(ke_idx:pe_idx)
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
            if (parent) write(6,'(3X,i10,5X,2(es17.10,5X),4X,3(es17.10,5X))') ireport, estimators(ke_idx), &
                                                             estimators(pe_idx), estimators(hf_ke_idx), &
                                                             estimators(hf_pe_idx), estimators(hf_part_idx)
            comms_found = abs(estimators(comms_found_idx)) > depsilon
            call calc_interact(comms_found, soft_exit)
            if (soft_exit) exit
        end do

        if (parent) write(6, '()')

        ! Return sys in unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

        call dSFMT_end(rng)

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
