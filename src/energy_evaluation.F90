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
                               proj_energy, shift, vary_shift, vary_shift_from,                  &
                               vary_shift_from_proje, D0_population,                           &
                               fold_line
        use hfs_data, only: proj_hf_O_hpsip, proj_hf_H_hfpsip, hf_signed_pop, D0_hf_population, hf_shift
        use calc, only: doing_calc, hfs_fciqmc_calc, folded_spectrum

        use parallel

        integer(lint), intent(inout) :: ntot_particles_old(sampling_size)

        real(dp) :: ir(sampling_size+7), ir_sum(sampling_size+7)
        integer(lint) :: ntot_particles(sampling_size), new_hf_signed_pop
        integer :: ierr

        ! Need to sum the number of particles and the projected energy over
        ! all processors.
        ir(1:sampling_size) = nparticles
        ir(sampling_size+1) = proj_energy
        ir(sampling_size+2) = D0_population
        ir(sampling_size+3) = rspawn

        if (doing_calc(hfs_fciqmc_calc)) then
            ! HFS calculations also need to know \tilde{N} = \sum_i sign(N_j^(H)) N_j^(HF),
            ! where N_j^(H) is the population of Hamiltonian walkers on j and
            ! N_j^(HF) the population of Hellmann-Feynman walkers on j.
            ir(sampling_size+4) = calculate_hf_signed_pop()
            ir(sampling_size+5) = proj_hf_O_hpsip
            ir(sampling_size+6) = proj_hf_H_hfpsip
            ir(sampling_size+7) = D0_hf_population
        end if

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
        D0_population = ir_sum(sampling_size+2)
        rspawn = ir_sum(sampling_size+3)
        new_hf_signed_pop = nint(ir_sum(sampling_size+4), lint)
        proj_hf_O_hpsip = ir_sum(sampling_size+5)
        proj_hf_H_hfpsip = ir_sum(sampling_size+6)
        D0_hf_population = ir_sum(sampling_size+7)

        if (vary_shift) then
            call update_shift(ntot_particles_old(1), ntot_particles(1), ncycles)
            if (doing_calc(hfs_fciqmc_calc)) then
                call update_hf_shift(ntot_particles_old(1), ntot_particles(1), hf_signed_pop, &
                                     new_hf_signed_pop, ncycles)
            end if
        end if
        ntot_particles_old = ntot_particles
        hf_signed_pop = new_hf_signed_pop
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
                  hf_shift = proj_hf_O_hpsip/D0_population + proj_hf_H_hfpsip/D0_population &
                                                           - (proj_energy*D0_hf_population)/D0_population**2
                endif
            else
                shift = vary_shift_from
            end if
        end if

        ! average energy quantities over report loop.
        proj_energy = proj_energy/ncycles
        D0_population = D0_population/ncycles
        ! Similarly for the HFS estimator
        D0_hf_population = D0_hf_population/ncycles
        proj_hf_O_hpsip = proj_hf_O_hpsip/ncycles
        proj_hf_H_hfpsip = proj_hf_H_hfpsip/ncycles
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

        ! Update the Hellmann-Feynman shift, \tilde{S}.
        ! In:
        !    nparticles_old: N_w(beta-A*tau); total Hamiltonian population at beta-Atau.
        !    nparticles: N_w(beta); total Hamiltonian population at beta.
        !    nhf_particles_old: N_w(beta-A*tau); total Hellmann-Feynman (signed) population at beta-Atau.
        !    nhf_particles: N_w(beta); total Hellmann-Feynman (signed) population at beta.
        !
        ! WARNING:
        ! The Hellmann-Feynman signed population is not simply the sum over
        ! Hellmann-Feynman walkers but also involves the Hamiltonian walkers and
        ! *must* be calculated using calculate_hf_signed_pop.

        use fciqmc_data, only: tau, shift_damping
        use hfs_data, only: hf_shift

        integer(lint), intent(in) :: nparticles_old, nparticles, nhf_particles_old, nhf_particles
        integer, intent(in) :: nupdate_steps

        ! Given the definition of the shift, S, \tilde{S} \equiv \frac{dS}{d\alpha}|_{\alpha=0}.
        ! Hence \tilde{S}(\beta) =
        !           \tilde{S}(\beta-A\tau)
        !           - \frac{\xi}{A\tau} [ \frac{\tilde{N}_w(\beta)}{N_w(\beta)}
        !                                 - \frac{\tilde{N}_w(\beta-A\tau)}{N_w(\beta-A\tau)} ]
        ! where N_w(\beta) is the total population of (Hamiltonian) walkers at
        ! imaginary time \beta and \tilde{N}_w = \frac{dN_w}{d\alpha}|_{\alpha=0}.
        ! The latter quantity is calculated in calculate_hf_signed pop.

        hf_shift = hf_shift - &
                 (shift_damping/(tau*nupdate_steps)) &
                 *(real(nhf_particles,p)/nparticles - real(nhf_particles_old,p)/nparticles_old)

    end subroutine update_hf_shift

    function calculate_hf_signed_pop() result(hf_signed_pop)

        ! Find
        !    \sum_j sign(N_j(\beta)) \tilde{N}_j(\beta)
        ! where N_j(\beta) is the Hamiltonian population on j at imaginary time
        ! \beta and \tilde{N}_j(\beta) is the Hellmann-Feynman population on
        ! j at imaginary time \beta.

        use fciqmc_data, only: walker_population, tot_walkers
        use hfs_data, only: alpha0

        integer(lint) :: hf_signed_pop

        integer :: i

        hf_signed_pop = 0_lint
        do i = 1, tot_walkers
            if (walker_population(1,i) == 0) then
                if (alpha0 < 0) then
                    ! letting alpha->0_-
                    hf_signed_pop = hf_signed_pop - abs(walker_population(2,i))
                else
                    ! letting alpha->0_+
                    hf_signed_pop = hf_signed_pop + abs(walker_population(2,i))
                end if
            else
                hf_signed_pop = hf_signed_pop + sign(1, walker_population(1,i))*walker_population(2,i)
            end if
        end do

    end function calculate_hf_signed_pop

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
        use hamiltonian_hub_k, only: slater_condon2_hub_k

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
        use hamiltonian_hub_real, only: slater_condon1_hub_real

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

    subroutine update_proj_energy_mol(idet)

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
        ! This procedure is for molecular systems (i.e. those defined by an
        ! FCIDUMP file).
        !
        ! In:
        !    idet: index of current determinant in the main walker list.

        use basis, only: basis_fns
        use determinants, only: decode_det
        use excitations, only: excit, get_excitation
        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
        use hamiltonian_molecular, only: slater_condon1_mol_excit, slater_condon2_mol_excit
        use point_group_symmetry, only: cross_product_pg_basis
        use system, only: nel

        integer, intent(in) :: idet

        type(excit) :: excitation
        real(p) :: hmatel
        integer :: occ_list(nel), ij_sym, ab_sym

        excitation = get_excitation(walker_dets(:,idet), f0)

        select case(excitation%nexcit)
        case (0)
            ! Have reference determinant.
            D0_population = D0_population + walker_population(1,idet)
        case(1)
            ! Have a determinant connected to the reference determinant by
            ! a single excitation: add to projected energy.
            ! decode
            ! Is excitation symmetry allowed?
            if (basis_fns(excitation%from_orb(1))%Ms == basis_fns(excitation%to_orb(1))%Ms .and. &
                    basis_fns(excitation%from_orb(1))%sym == basis_fns(excitation%to_orb(1))%sym) then
                call decode_det(walker_dets(:,idet), occ_list)
                hmatel = slater_condon1_mol_excit(occ_list, excitation%from_orb(1), excitation%to_orb(1), &
                                                  excitation%perm)
                proj_energy = proj_energy + hmatel*walker_population(1,idet)
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
                    proj_energy = proj_energy + hmatel*walker_population(1,idet)
                end if
            end if
        end select

    end subroutine update_proj_energy_mol

    subroutine update_proj_hfs_hamiltonian_hub_k(idet)

        ! Add the contribution of the current determinant to the projected
        ! energy in an identical way to update_proj_energy_hub_k.

        ! Also add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimator.

        ! For debugging purposes, this procedure is for when we are sampling
        ! O=H.

        ! This procedure is for the Hubbard model in momentum space only.

        ! In:
        !    idet: index of current determinant in the main walker list.

        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
        use excitations, only: excit, get_excitation
        use hamiltonian_hub_k, only: slater_condon2_hub_k
        use hfs_data, only: D0_hf_population, proj_hf_O_hpsip, proj_hf_H_hfpsip

        integer, intent(in) :: idet
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
            proj_energy = proj_energy + hmatel*walker_population(1,idet)
            ! DEBUG/TESTING: For now, just using O=H
            ! In this case, \sum_j O_0j c_j = proj_energy
            proj_hf_O_hpsip = proj_energy
            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*walker_population(2,idet)
        end if

    end subroutine update_proj_hfs_hamiltonian_hub_k

    subroutine update_proj_hfs_diagonal_hub_k(idet)

        ! Add the contribution of the current determinant to the projected
        ! energy in an identical way to update_proj_energy_hub_k.

        ! Also add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimator.

        ! This procedure is for when we are sampling an operator, O, which is
        ! diagonal in the Slater determinant space.

        ! This procedure is for the Hubbard model in momentum space only.

        ! In:
        !    idet: index of current determinant in the main walker list.

        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
        use excitations, only: excit, get_excitation
        use hamiltonian_hub_k, only: slater_condon2_hub_k
        use hfs_data, only: D0_hf_population, proj_hf_O_hpsip, proj_hf_H_hfpsip

        integer, intent(in) :: idet
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
            proj_energy = proj_energy + hmatel*walker_population(1,idet)

            ! O is diagonal in the determinant basis.  As we are actually
            ! sampling O - <D0|O|D0>, this means that \sum_j O_j0 c_j = 0.

            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*walker_population(2,idet)
        end if

    end subroutine update_proj_hfs_diagonal_hub_k

end module energy_evaluation
