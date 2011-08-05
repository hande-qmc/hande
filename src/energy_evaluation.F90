module energy_evaluation

! This module contains procedure for evaluating and estimating the energy of
! a system based upon the population dynamics of an FCIQMC calculation.

use const

implicit none

contains

    subroutine update_energy_estimators(ireport, ntot_particles_old)

        ! Update the shift and average the shift and projected energy
        ! estimators.

        ! Should be called every report loop in an FCIQMC/iFCIQMC calculation.

        ! In:
        !    ireport: index of the current report loop.
        ! Inout:
        !    ntot_particles_old: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop.

        use fciqmc_data, only: nparticles, sampling_size, target_particles, ncycles, rspawn,   &
                               proj_energy, av_proj_energy, av_D0_population, shift, av_shift, &
                               vary_shift, start_vary_shift, vary_shift_from,                  &
                               vary_shift_from_proje, fsfciqmc_vary_shift_from_proje,          &
                               D0_population, fold_line, fs_offset
        use hfs_data, only: proj_hf_expectation, av_proj_hf_expectation
        use calc, only: doing_calc, hfs_fciqmc_calc

        use parallel
        
        integer, intent(in) :: ireport
        integer, intent(inout) :: ntot_particles_old(sampling_size)

#ifdef PARALLEL
        real(dp) :: ir(sampling_size+4), ir_sum(sampling_size+4)
        integer :: ntot_particles(sampling_size), ierr

            ! Need to sum the number of particles and the projected energy over
            ! all processors.
            ir(1:sampling_size) = nparticles
            ir(sampling_size+1) = proj_energy
            ir(sampling_size+2) = proj_hf_expectation
            ir(sampling_size+3) = D0_population
            ir(sampling_size+4) = rspawn
            call mpi_allreduce(ir, ir_sum, size(ir), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
            ntot_particles = nint(ir_sum(1:sampling_size))
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
                start_vary_shift = ireport
                if(fsfciqmc_vary_shift_from_proje) then
                    !if running a folded spectrum calculation, set the shift to
                    !instantaneously be the projected energy of the folded hamiltonian
                    shift = (proj_energy/D0_population - fold_line)**2 + fs_offset
                else if (vary_shift_from_proje) then
                    ! Set shift to be instantaneous projected energy.
                    shift = proj_energy/D0_population
                else
                    shift = vary_shift_from
                end if
            end if
#else
            if (vary_shift) then
                call update_shift(ntot_particles_old(1), nparticles(1), ncycles)
                if (doing_calc(hfs_fciqmc_calc)) then
                    call update_hf_shift(ntot_particles_old(1), nparticles(1), ntot_particles_old(2), nparticles(2), ncycles)
                end if
            end if
            ntot_particles_old = nparticles
            if (nparticles(1) > target_particles .and. .not.vary_shift) then
                vary_shift = .true.
                start_vary_shift = ireport
                if(fsfciqmc_vary_shift_from_proje) then
                    !if running a folded spectrum calculation, set the shift to
                    !instantaneously be the projected energy of the folded hamiltonian
                    shift = (proj_energy/D0_population - fold_line)**2 + fs_offset
                else if (vary_shift_from_proje) then
                    ! Set shift to be instantaneous projected energy.
                    shift = proj_energy/D0_population
                else
                    shift = vary_shift_from
                end if
            end if
#endif

            ! Running average projected energy 
            ! Note that as proj_energy and D0_population are accumulated over
            ! the report loop, proj_energy/D0_population is 
            !   \sum_j <D_j|H|D_0> <N_j>/<N_0>,
            ! where <N_j> is the mean of the population on determinant j over
            ! the report loop.
            ! As a result, the running accumulation of the projected energy
            ! need only be divided by the the number of report loops in order to
            ! get an estimate of the average projected energy.
            av_proj_energy = av_proj_energy + proj_energy
            av_D0_population = av_D0_population + D0_population
            ! average energy quantities over report loop.
            proj_energy = proj_energy/ncycles
            D0_population = D0_population/ncycles
            ! Similarly for the HFS estimator
            av_proj_hf_expectation = av_proj_hf_expectation + proj_hf_expectation
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

        use fciqmc_data, only: shift, tau, shift_damping, av_shift

        integer, intent(in) :: nparticles_old, nparticles, nupdate_steps

        shift = shift - log(real(nparticles,p)/nparticles_old)*shift_damping/(tau*nupdate_steps)
        av_shift = av_shift + shift

    end subroutine update_shift

    subroutine update_hf_shift(nparticles_old, nparticles, nhf_particles_old, nhf_particles, nupdate_steps)

        use fciqmc_data, only: tau, shift_damping
        use hfs_data, only: hf_shift, av_hf_shift

        integer, intent(in) :: nparticles_old, nparticles, nhf_particles_old, nhf_particles, nupdate_steps

        hf_shift = hf_shift - &
                 (shift_damping/(tau*nupdate_steps)) &
                 *(real(nhf_particles,p)/nparticles - real(nhf_particles_old,p)/nparticles_old)
        av_hf_shift = av_hf_shift + hf_shift

    end subroutine update_hf_shift

    subroutine update_proj_energy_hub_k(idet)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i
        ! If the current determinant is the reference determinant, then
        ! N_0 is stored as D0_population.  This makes normalisation very
        ! efficient.
        ! This procedure is for the Hubbard model in momentum space only.
        ! In:
        !    idet: index of current determinant in the main walker list.

        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
        use excitations, only: excit, get_excitation
        use hamiltonian, only: slater_condon2_hub_k

        integer, intent(in) :: idet
        type(excit) :: excitation
        real(p) :: hmatel

        excitation = get_excitation(walker_dets(:,idet), f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_population = D0_population + walker_population(1,idet)
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to 
            ! projected energy.
            hmatel = slater_condon2_hub_k(excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            proj_energy = proj_energy + hmatel*walker_population(1,idet)
        end if

    end subroutine update_proj_energy_hub_k

    subroutine update_proj_energy_hub_real(idet)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i
        ! If the current determinant is the reference determinant, then
        ! N_0 is stored as D0_population.  This makes normalisation very
        ! efficient.
        ! This procedure is for the Hubbard model in real space only.
        ! In:
        !    idet: index of current determinant in the main walker list.

        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
        use excitations, only: excit, get_excitation
        use hamiltonian, only: slater_condon1_hub_real

        integer, intent(in) :: idet
        type(excit) :: excitation
        real(p) :: hmatel

        excitation = get_excitation(walker_dets(:,idet), f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_population = D0_population + walker_population(1,idet)
        else if (excitation%nexcit == 1) then
            ! Have a determinant connected to the reference determinant: add to 
            ! projected energy.
            hmatel = slater_condon1_hub_real(excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
            proj_energy = proj_energy + hmatel*walker_population(1,idet)
        end if

    end subroutine update_proj_energy_hub_real

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

        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
        use excitations, only: excit, get_excitation
        use hamiltonian, only: slater_condon2_hub_k
        use hfs_data, only: D0_hf_population

        integer, intent(in) :: idet
        real(p), intent(inout) :: inst_proj_energy, inst_proj_hf_t1
        type(excit) :: excitation
        real(p) :: hmatel

        excitation = get_excitation(walker_dets(:,idet), f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_population = D0_population + walker_population(1,idet)
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
