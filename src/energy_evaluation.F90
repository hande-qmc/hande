module energy_evaluation

! This module contains procedure for evaluating and estimating the energy of
! a system based upon the population dynamics of an FCIQMC calculation.

use const

implicit none

contains

    subroutine update_energy_estimators(ntot_particles_old)

        ! Update the shift and average the shift and projected energy
        ! estimators.

        ! Should be called every report loop in an FCIQMC/iFCIQMC calculation.

        ! Inout:
        !    ntot_particles_old: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop.

        use fciqmc_data, only: nparticles, sampling_size, target_particles, ncycles, rspawn,   &
                               proj_energy, shift, vary_shift, vary_shift_from,                &
                               vary_shift_from_proje, D0_population, fold_line
        use hfs_data, only: proj_hf_expectation
        use calc, only: doing_calc, hfs_fciqmc_calc, folded_spectrum

        use parallel

        integer(lint), intent(inout) :: ntot_particles_old(sampling_size)

        real(dp) :: ir(sampling_size+4), ir_sum(sampling_size+4)
        integer(lint) :: ntot_particles(sampling_size)
        integer :: ierr

        ! Need to sum the number of particles and the projected energy over
        ! all processors.
        ir(1:sampling_size) = nparticles
        ir(sampling_size+1) = proj_energy
        ir(sampling_size+2) = proj_hf_expectation
        ir(sampling_size+3) = D0_population
        ir(sampling_size+4) = rspawn

        ! Don't bother to optimise for running in serial.  This is a fast
        ! routine and is run only once per report loop anyway!
#ifdef PARALLEL
        call mpi_allreduce(ir, ir_sum, size(ir), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        ir_sum = ir
        ierr = 0 ! Prevent warning about unused variable in serial so -Werror can be used.
#endif

        ntot_particles = nint(ir_sum(1:sampling_size), lint)
        proj_energy = ir_sum(sampling_size+1)
        proj_hf_expectation = ir_sum(sampling_size+2)
        D0_population = ir_sum(sampling_size+3)
        rspawn = ir_sum(sampling_size+4)

        if (vary_shift) then
            call update_shift(ntot_particles_old(1), ntot_particles(1), ncycles)
            if (doing_calc(hfs_fciqmc_calc)) then
                call update_hf_shift(ntot_particles_old(1), ntot_particles(1), ntot_particles_old(2), &
                                     ntot_particles(2), ncycles)
            end if
        end if
        ntot_particles_old = ntot_particles
        if (ntot_particles(1) > target_particles .and. .not.vary_shift) then
            vary_shift = .true.
            if (vary_shift_from_proje) then
                if(doing_calc(folded_spectrum)) then
                  !if running a folded spectrum calculation, set the shift to
                  !instantaneously be the projected energy of the folded hamiltonian
                  shift = (proj_energy/D0_population - fold_line)**2
                else
                  ! Set shift to be instantaneous projected energy.
                  shift = proj_energy/D0_population
                endif
            else
                shift = vary_shift_from
            end if
        end if

        ! average energy quantities over report loop.
        proj_energy = proj_energy/ncycles
        D0_population = D0_population/ncycles
        ! Similarly for the HFS estimator
        proj_hf_expectation = proj_hf_expectation/ncycles
        ! average spawning rate over report loop and processor.
        rspawn = rspawn/(ncycles*nprocs)

    end subroutine update_energy_estimators

    subroutine update_shift(nparticles_old, nparticles, nupdate_steps)

        ! Update the shift according to:
        !  shift(beta) = shift(beta-A*tau) - xi*log(N_w(tau)/N_w(beta-A*tau))/(A*tau)
        ! where
        !  * shift(beta) is the shift at imaginary time beta;
        !  * A*tau is the amount of imaginary time between shift-updates (=# of
        !    Monte Carlo cycles between updating the shift);
        !  * xi is a damping factor (0.05-0.10 is appropriate) to prevent large fluctations;
        !  * N_w(beta) is the total number of particles at imaginary time beta.
        ! The running average of the shift is also updated.
        ! In:
        !    nparticles_old: N_w(beta-A*tau).
        !    nparticles: N_w(beta).

        use calc, only: doing_calc, folded_spectrum
        use fciqmc_data, only: shift, tau, shift_damping, dmqmc_factor

        integer(lint), intent(in) :: nparticles_old, nparticles
        integer, intent(in) :: nupdate_steps

        ! dmqmc_factor is included to account for a factor of 1/2 introduced into tau in
        ! DMQMC calculations. In all other calculation types, it is set to 1, and so can be ignored.
        shift = shift - log(real(nparticles,p)/nparticles_old)*shift_damping/(dmqmc_factor*tau*nupdate_steps)

    end subroutine update_shift

    subroutine update_hf_shift(nparticles_old, nparticles, nhf_particles_old, nhf_particles, nupdate_steps)

        use fciqmc_data, only: tau, shift_damping
        use hfs_data, only: hf_shift

        integer(lint), intent(in) :: nparticles_old, nparticles, nhf_particles_old, nhf_particles
        integer, intent(in) :: nupdate_steps

        hf_shift = hf_shift - &
                 (shift_damping/(tau*nupdate_steps)) &
                 *(real(nhf_particles,p)/nparticles - real(nhf_particles_old,p)/nparticles_old)

    end subroutine update_hf_shift

    pure subroutine update_proj_energy_hub_k(f0, cdet, pop, D0_pop_sum, proj_energy_sum)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for the Hubbard model in momentum space only.
        ! In:
        !    f0: reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info
        use excitations, only: excit, get_excitation
        use hamiltonian_hub_k, only: slater_condon2_hub_k

        integer(i0), intent(in) :: f0(:)
        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum

        type(excit) :: excitation
        real(p) :: hmatel

        excitation = get_excitation(cdet%f, f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon2_hub_k(excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            proj_energy_sum = proj_energy_sum + hmatel*pop
        end if

    end subroutine update_proj_energy_hub_k

    pure subroutine update_proj_energy_hub_real(f0, cdet, pop, D0_pop_sum, proj_energy_sum)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for the Hubbard model in real space only.
        ! In:
        !    f0: reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info
        use excitations, only: excit, get_excitation
        use hamiltonian_hub_real, only: slater_condon1_hub_real

        integer(i0), intent(in) :: f0(:)
        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum

        type(excit) :: excitation
        real(p) :: hmatel

        excitation = get_excitation(cdet%f, f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 1) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon1_hub_real(excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
            proj_energy_sum = proj_energy_sum + hmatel*pop
        end if

    end subroutine update_proj_energy_hub_real

    pure subroutine update_proj_energy_mol(f0, cdet, pop, D0_pop_sum, proj_energy_sum)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for molecular systems (i.e. those defined by an
        ! FCIDUMP file).
        !
        ! In:
        !    f0: reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use basis, only: basis_fns
        use determinants, only: decode_det, det_info
        use excitations, only: excit, get_excitation
        use hamiltonian_molecular, only: slater_condon1_mol_excit, slater_condon2_mol_excit
        use point_group_symmetry, only: cross_product_pg_basis
        use system, only: nel

        integer(i0), intent(in) :: f0(:)
        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum

        type(excit) :: excitation
        real(p) :: hmatel
        integer :: occ_list(nel), ij_sym, ab_sym

        excitation = get_excitation(cdet%f, f0)

        select case(excitation%nexcit)
        case (0)
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        case(1)
            ! Have a determinant connected to the reference determinant by
            ! a single excitation: add to projected energy.
            ! decode
            ! Is excitation symmetry allowed?
            if (basis_fns(excitation%from_orb(1))%Ms == basis_fns(excitation%to_orb(1))%Ms .and. &
                    basis_fns(excitation%from_orb(1))%sym == basis_fns(excitation%to_orb(1))%sym) then
                call decode_det(cdet%f, occ_list)
                hmatel = slater_condon1_mol_excit(occ_list, excitation%from_orb(1), excitation%to_orb(1), &
                                                  excitation%perm)
                proj_energy_sum = proj_energy_sum + hmatel*pop
            end if
        case(2)
            ! Have a determinant connected to the reference determinant by
            ! a double excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (basis_fns(excitation%from_orb(1))%Ms+basis_fns(excitation%from_orb(2))%Ms == &
                    basis_fns(excitation%to_orb(1))%Ms+basis_fns(excitation%to_orb(2))%Ms) then
                ij_sym = cross_product_pg_basis(excitation%from_orb(1), excitation%from_orb(2))
                ab_sym = cross_product_pg_basis(excitation%to_orb(1), excitation%to_orb(2))
                if (ij_sym == ab_sym) then
                    hmatel = slater_condon2_mol_excit(excitation%from_orb(1), excitation%from_orb(2), &
                                                      excitation%to_orb(1), excitation%to_orb(2),     &
                                                      excitation%perm)
                    proj_energy_sum = proj_energy_sum + hmatel*pop
                end if
            end if
        end select

    end subroutine update_proj_energy_mol

    subroutine update_proj_hfs_hub_k(idet, inst_proj_energy, inst_proj_hf_t1)

        ! Add the contribution of the current determinant to the projected
        ! energy in an identical way to update_proj_energy_hub_k.

        ! Also add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimator.

        ! This procedure is for the Hubbard model in momentum space only.

        ! In:
        !    idet: index of current determinant in the main walker list.
        ! In/Out:
        !    inst_proj_energy: running total of the \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !    This is updated if D_i is connected to D_0 (and isn't D_0).

        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population_cycle, proj_energy
        use excitations, only: excit, get_excitation
        use hamiltonian_hub_k, only: slater_condon2_hub_k
        use hfs_data, only: D0_hf_population

        integer, intent(in) :: idet
        real(p), intent(inout) :: inst_proj_energy, inst_proj_hf_t1
        type(excit) :: excitation
        real(p) :: hmatel

        excitation = get_excitation(walker_dets(:,idet), f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_population_cycle = D0_population_cycle + walker_population(1,idet)
            D0_hf_population = D0_hf_population + walker_population(2,idet)
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon2_hub_k(excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            inst_proj_energy = inst_proj_energy + hmatel*walker_population(1,idet)
            inst_proj_hf_t1 = inst_proj_hf_t1 + hmatel*walker_population(2,idet)
        end if

    end subroutine update_proj_hfs_hub_k

end module energy_evaluation
