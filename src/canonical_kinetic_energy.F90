module canonical_kinetic_energy

! Estimate the canonical free-electron total energy.

use const

implicit none

! Number of report loops use to estimate the energy.
integer :: nkinetic_cycles

enum, bind(c)
    ! Ensure energy estimates are always first.
    ! Index for kinetic energy.
    enumerator :: ke_idx = 1
    ! Index for Hartree-Fock energy.
    enumerator :: hf_idx
    ! Index for Variance in kinetic energy.
    enumerator :: ke_var_idx
    ! Index for Variance in Hartree-Fock energy.
    enumerator :: hf_var_idx
    ! Index for squared means needed for variance estimation in parallel.
    enumerator :: ke_sq_idx
    enumerator :: hf_sq_idx
    ! Index for checking for interaction with the calculation
    enumerator :: comms_found_idx
    ! last_idx-1 gives number of estimates.
    enumerator :: last_idx
end enum


contains

    subroutine estimate_kinetic_energy(sys)

        ! From the Fermi factors calculated in the grand canonical ensemble we can
        ! estimate the total energy in the canonical ensemble by generating determinants
        ! with arbitrary particle number and only keeping those with <N> = nel.

        ! In:
        !    sys: system being studied.

        use system
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use parallel
        use calc, only: ms_in, seed, fermi_temperature
        use fciqmc_data, only: init_beta, D0_population
        use utils, only: rng_init_info
        use hamiltonian_ueg, only: sum_sp_eigenvalues, potential_energy_ueg
        use interact, only: calc_interact, check_comms_file

        type(sys_t), intent(inout) :: sys

        real(dp) :: p_single(sys%basis%nbasis/2)
        real(dp) :: r
        integer :: occ_list(sys%nel)
        logical :: gen
        real(p) :: energy(2)
        real(p) :: delta(2), mean(2), std(2)
        integer :: ierr, ireport, iorb
        integer(int_64) :: iaccept
        real(p) :: local_estimators(last_idx-1), estimators(last_idx-1)

        type(sys_t) :: sys_bak
        type (dSFMT_t) :: rng
        logical :: soft_exit, comms_found

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)
        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        if (fermi_temperature) then
            init_beta = init_beta / sys%ueg%ef
        end if

        if (parent) then
            write (6,'(1X,a72)') 'E_0: Current estimate for thermal kinetic energy i.e. 1/Z Tr(\rho_0 H_0).'
            write (6,'(1X,a79)') 'E_HF: Current estimate for thermal "Hartree-Fock" energy i.e. 1/Z Tr(\rho_0 H).'
            write (6,'()')
        end if

        if (parent) write (6,'(1X,a12,6X,a3,19X,a7,15X,a4,18X,a10)') '# iterations', 'E_0', 'E_Error', 'E_HF', 'E_HF-Error'

        forall (iorb=1:sys%basis%nbasis:2) p_single(iorb/2+1) = 1.0_p / &
                                                          (1+exp(init_beta*(sys%basis%basis_fns(iorb)%sp_eigv-sys%ueg%chem_pot)))

        iaccept = 0 ! running total number of samples.
        delta = 0.0_p
        mean = 0.0_p
        local_estimators = 0.0_p
        do ireport = 1, nkinetic_cycles
            do while (iaccept < ireport*D0_population)
                if (sys%nalpha > 0) call generate_allowed_orbital_list(sys, rng, p_single, sys%nalpha, &
                                                                       1, occ_list(:sys%nalpha), gen)
                if (.not. gen) cycle
                if (sys%nbeta > 0) call generate_allowed_orbital_list(sys, rng, p_single, sys%nbeta, &
                                                                      0, occ_list(sys%nalpha+1:), gen)
                if (.not. gen) cycle
                iaccept = iaccept + 1
                ! Calculate Kinetic and Hartree-Fock energies.
                energy(ke_idx) = sum_sp_eigenvalues(sys, occ_list)
                energy(hf_idx) = energy(ke_idx) + potential_energy_ueg(sys, occ_list)
                ! Estimate mean and variance using Knuth's online algorithm.
                ! See: http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm
                delta = energy - local_estimators(ke_idx:hf_idx)
                ! Update means.
                local_estimators(ke_idx:hf_idx) = local_estimators(ke_idx:hf_idx) + delta/iaccept
                ! Mean squared.
                local_estimators(ke_sq_idx:hf_sq_idx) = local_estimators(ke_idx:hf_idx)**2.0_p
                ! Update variances (actually the variance is
                ! estimators(ke_var_idx:) / (iaccept-1), see link for details)
                local_estimators(ke_var_idx:hf_var_idx) =  local_estimators(ke_var_idx:hf_var_idx) + &
                                                        delta*(energy-local_estimators(ke_idx:hf_idx))
            end do

            if (check_comms_file()) local_estimators(comms_found_idx) = 1.0_p

#ifdef PARALLEL
            ! More efficient in parallel.
            call mpi_reduce(local_estimators, estimators, last_idx-1, mpi_preal, MPI_SUM, root, MPI_COMM_WORLD, ierr)
#else
            estimators = local_estimators
#endif
            ! Average means over processors.
            mean = estimators(ke_idx:hf_idx) / nprocs
            ! Find variance in average over processors (denoted var(\sum_i^nprocs)).
            ! Assuming the samples are uncorrelated and the number of samples
            ! contributing to each estimate is the same, then:
            ! var(\sum_i^nprocs) = 1/nprocs (\sum_i^nprocs var(i) + \sum_i^nprocs mean^2(i))
            !                      - mean(\sum_i^nprocs).
            ! This follows from the fact that (var + mean) = 1/m \sum_i^m x_i^2
            ! and mean = 1/m \sum_i^m x_i.
            ! The standard deviation of the mean is then sqrt(var/(nprocs*iaccept)).
            std = sqrt(((estimators(ke_var_idx:hf_var_idx)/(real(iaccept-1,p)) + &
                  estimators(ke_sq_idx:hf_sq_idx))/nprocs - mean**2.0_p)/(nprocs*real(iaccept, p)))
            if (parent) write(6,'(3X,i10,5X,4(es17.10,5X))') ireport, mean(ke_idx), std(ke_idx), mean(hf_idx), std(hf_idx)
            comms_found = abs(estimators(comms_found_idx)) > depsilon
            call calc_interact(comms_found, soft_exit)
            if (soft_exit) exit
        end do

        if (parent) write(6, '()')

    end subroutine estimate_kinetic_energy

    subroutine generate_allowed_orbital_list(sys, rng, porb, nselect, spin_factor, occ_list, gen)

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
        ! Out:
        !    occ_list: array containing occupied orbitals.
        !    gen: true if generation attempt was successful (i.e. nselect orbitals were actually selected).

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        real(dp), intent(in) :: porb(:)
        integer, intent(in) :: nselect
        integer, intent(in) :: spin_factor
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: occ_list(:)
        logical, intent(out) :: gen

        integer :: iorb, iselect
        real(dp) :: r

        iselect = 0
        occ_list = 0

        do iorb = 1, sys%basis%nbasis/2
            ! Select a random orbital.
            r = get_rand_close_open(rng)
            if (porb(iorb) > r) then
                iselect = iselect + 1
                if (iselect > nselect) then
                    ! Selected too many.
                    gen = .false.
                    exit
                end if
                occ_list(iselect) = 2*iorb - spin_factor
            end if
        end do
        if (iselect == nselect) then
            gen = .true.
        else
            gen = .false.
        end if

    end subroutine generate_allowed_orbital_list

end module canonical_kinetic_energy
