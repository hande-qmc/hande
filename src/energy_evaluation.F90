Module energy_evaluation

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
        use hfs_data, only: proj_hf_O_hpsip, proj_hf_H_hfpsip, hf_signed_pop, D0_hf_population, hf_shift
        use calc, only: doing_calc, hfs_fciqmc_calc, folded_spectrum

        use parallel

        real(dp), intent(inout) :: ntot_particles_old(sampling_size)

        real(dp) :: ir(sampling_size+7), ir_sum(sampling_size+7)
        real(dp) :: new_hf_signed_pop
        real(dp) :: ntot_particles(sampling_size)
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

        ntot_particles = ir_sum(1:sampling_size)
        proj_energy = ir_sum(sampling_size+1)
        D0_population = ir_sum(sampling_size+2)
        rspawn = ir_sum(sampling_size+3)
        if (doing_calc(hfs_fciqmc_calc)) then
            new_hf_signed_pop = ir_sum(sampling_size+4)
            proj_hf_O_hpsip = ir_sum(sampling_size+5)
            proj_hf_H_hfpsip = ir_sum(sampling_size+6)
            D0_hf_population = ir_sum(sampling_size+7)
        end if

        if (vary_shift) then
            call update_shift(shift(1), ntot_particles_old(1), ntot_particles(1), ncycles)
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

    subroutine update_shift(loc_shift, nparticles_old, nparticles, nupdate_steps)

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

        real(p), intent(inout) :: loc_shift
        real(dp), intent(in) :: nparticles_old, nparticles
        integer, intent(in) :: nupdate_steps

        ! dmqmc_factor is included to account for a factor of 1/2 introduced into tau in
        ! DMQMC calculations. In all other calculation types, it is set to 1, and so can be ignored.
        loc_shift = loc_shift - log(nparticles/nparticles_old)*shift_damping/(dmqmc_factor*tau*nupdate_steps)

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

        real(dp), intent(in) :: nparticles_old, nparticles, nhf_particles_old, nhf_particles
        integer, intent(in) :: nupdate_steps

        ! Given the definition of the shift, S, \tilde{S} \equiv \frac{dS}{d\alpha}|_{\alpha=0}.
        ! Hence \tilde{S}(\beta) =
        !           \tilde{S}(\beta-A\tau)
        !           - \frac{\xi}{A\tau} [ \frac{\tilde{N}_w(\beta)}{N_w(\beta)}
        !                                 - \frac{\tilde{N}_w(\beta-A\tau)}{N_w(\beta-A\tau)} ]
        ! where N_w(\beta) is the total population of (Hamiltonian) walkers at
        ! imaginary time \beta and \tilde{N}_w = \frac{dN_w}{d\alpha}|_{\alpha=0}.
        ! The latter quantity is calculated in calculate_hf_signed fpop.

        hf_shift = hf_shift - &
                 (shift_damping/(tau*nupdate_steps)) &
                 *(nhf_particles/nparticles - nhf_particles_old/nparticles_old)

    end subroutine update_hf_shift

    function calculate_hf_signed_pop() result(hf_signed_pop)

        ! Find
        !    \sum_j sign(N_j(\beta)) \tilde{N}_j(\beta)
        ! where N_j(\beta) is the Hamiltonian population on j at imaginary time
        ! \beta and \tilde{N}_j(\beta) is the Hellmann-Feynman population on
        ! j at imaginary time \beta.

        use fciqmc_data, only: walker_population, tot_walkers, real_factor, sampling_size
        use hfs_data, only: alpha0

        real(dp) :: hf_signed_pop
        real(dp) :: real_population(sampling_size)

        integer :: i

        hf_signed_pop = 0.0_dp
        do i = 1, tot_walkers
            real_population = real(abs(walker_population(:,i)),dp)/real_factor
            if (walker_population(1,i) == 0_int_p) then
                if (alpha0 < 0) then
                    ! letting alpha->0_-
                    hf_signed_pop = hf_signed_pop - real_population(2)
                else
                    ! letting alpha->0_+
                    hf_signed_pop = hf_signed_pop + real_population(2)
                end if
            else
                hf_signed_pop = hf_signed_pop + sign(1.0_dp, real_population(1))*&
                                                 real_population(2)
            end if
        end do

    end function calculate_hf_signed_pop

    pure subroutine update_proj_energy_hub_k(sys, f0, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for the Hubbard model in momentum space only.

        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        ! Out:
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info
        use excitations, only: excit, get_excitation
        use hamiltonian_hub_k, only: slater_condon2_hub_k
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit), intent(out) :: excitation
        real(p), intent(out) :: hmatel

        excitation = get_excitation(sys%nel, cdet%f, f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon2_hub_k(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            proj_energy_sum = proj_energy_sum + hmatel*pop
        end if

    end subroutine update_proj_energy_hub_k

    pure subroutine update_proj_energy_hub_real(sys, f0, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for the Hubbard model in real space only.

        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        ! Out:
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info
        use excitations, only: excit, get_excitation
        use hamiltonian_hub_real, only: slater_condon1_hub_real
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit), intent(out) :: excitation
        real(p), intent(out) :: hmatel

        excitation = get_excitation(sys%nel, cdet%f, f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 1) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon1_hub_real(sys, excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
            proj_energy_sum = proj_energy_sum + hmatel*pop
        end if

    end subroutine update_proj_energy_hub_real

    pure subroutine update_proj_energy_mol(sys, f0, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for molecular systems (i.e. those defined by an
        ! FCIDUMP file).

        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        ! Out:
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use basis, only: basis_fns
        use determinants, only: det_info
        use excitations, only: excit, get_excitation
        use hamiltonian_molecular, only: slater_condon1_mol_excit, slater_condon2_mol_excit
        use point_group_symmetry, only: cross_product_pg_basis
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit), intent(out) :: excitation
        real(p), intent(out) :: hmatel

        integer :: ij_sym, ab_sym

        excitation = get_excitation(sys%nel, cdet%f, f0)
        hmatel = 0.0_p

        select case(excitation%nexcit)
        case (0)
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        case(1)
            ! Have a determinant connected to the reference determinant by
            ! a single excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (basis_fns(excitation%from_orb(1))%Ms == basis_fns(excitation%to_orb(1))%Ms .and. &
                    basis_fns(excitation%from_orb(1))%sym == basis_fns(excitation%to_orb(1))%sym) then
                hmatel = slater_condon1_mol_excit(sys, cdet%occ_list, excitation%from_orb(1), excitation%to_orb(1), &
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
                    hmatel = slater_condon2_mol_excit(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                                      excitation%to_orb(1), excitation%to_orb(2),     &
                                                      excitation%perm)
                    proj_energy_sum = proj_energy_sum + hmatel*pop
                end if
            end if
        end select

    end subroutine update_proj_energy_mol

    pure subroutine update_proj_energy_ueg(sys, f0, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for the electron gas only.
        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        ! Out:
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info
        use excitations, only: excit, get_excitation
        use hamiltonian_ueg, only: slater_condon2_ueg
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit), intent(out) :: excitation
        real(p), intent(out) :: hmatel

        excitation = get_excitation(sys%nel, cdet%f, f0)
        hmatel = 0.0_p

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon2_ueg(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            proj_energy_sum = proj_energy_sum + hmatel*pop
        end if

    end subroutine update_proj_energy_ueg


    subroutine update_proj_hfs_hamiltonian(sys, f, fpop, f_hfpop, fdata, excitation, hmatel, &
                                           D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)

        ! Add the contribution of the current determinant to the projected
        ! energy in an identical way to update_proj_energy_hub_k.

        ! Also add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimators.

        ! For debugging purposes, this procedure is for when we are sampling
        ! O=H.

        ! This procedure is for the Hubbard model in momentum space only.

        ! In:
        !    sys: system being studied.  Unused.
        !    f(basis_length): bit string representation of the Slater determinant, D_i.
        !    fpop: Hamiltonian population on the determinant.
        !    f_hfpop: Hellmann-Feynman population on the determinant.
        !    fdata(:): additional information about the determinant (unused, for
        !       interface compatibility only).
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.
        ! In/Out:
        !    D0_hf_population: running total of the Hellmann-Feynman population
        !       on the reference.  Only updated if D_i *is* the reference determinant.
        !    proj_hf_O_hpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|O|D_0> N_i, where N_i is the Hamiltonian population on D_i.
        !    proj_hf_H_fhpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|H|D_0> \tilde{N}_i, where \tilde{N}_i is the
        !       Hellmann-Feynman population on D_i.

        use basis, only: basis_length
        use excitations, only: excit
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_hf_pop = D0_hf_pop + f_hfpop
        else if (excitation%nexcit <= 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected estimators.
            ! DEBUG/TESTING: For now, just using O=H
            ! In this case, \sum_j O_0j c_j = proj_energy
            proj_hf_O_hpsip = proj_hf_O_hpsip + hmatel*fpop
            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        end if

    end subroutine update_proj_hfs_hamiltonian

    subroutine update_proj_hfs_diagonal(sys, f, fpop, f_hfpop, fdata, excitation, hmatel, &
                                              D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)

        ! Add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimator.

        ! This procedure is for when we are sampling an operator, O, which is
        ! diagonal in the Slater determinant space.

        ! In:
        !    sys: system being studied.  Unused.
        !    f(basis_length): bit string representation of the Slater determinant, D_i
        !       (unused, for interface compatibility only).
        !    fpop: Hamiltonian population on the determinant (unused, for interface
        !       compatibility only).
        !    f_hfpop: Hellmann-Feynman population on the determinant.
        !    fdata(:): additional information about the determinant (unused, for
        !       interface compatibility only).
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.
        ! In/Out:
        !    D0_hf_population: running total of the Hellmann-Feynman population
        !       on the reference.  Only updated if D_i *is* the reference determinant.
        !    proj_hf_O_hpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|O|D_0> N_i, where N_i is the Hamiltonian population on D_i.
        !    proj_hf_H_fhpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|H|D_0> \tilde{N}_i, where \tilde{N}_i is the
        !       Hellmann-Feynman population on D_i.

        use basis, only: basis_length
        use excitations, only: excit
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_hf_pop = D0_hf_pop + f_hfpop
        else if (excitation%nexcit <= 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.

            ! O is diagonal in the determinant basis.  As we are actually
            ! sampling O - <D0|O|D0>, this means that \sum_j O_j0 c_j = 0.

            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        end if

    end subroutine update_proj_hfs_diagonal

    subroutine update_proj_hfs_double_occ_hub_k(sys, f, fpop, f_hfpop, fdata, excitation, hmatel, &
                                              D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)

        ! Add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimator.

        ! This procedure is for when we are sampling D, the double occupancy
        ! operator, in the Bloch (momentum) basis set for the Hubbard model.

        ! In:
        !    sys: system being studied.  Requires hubbard%u and lattice%nsites.
        !    f(basis_length): bit string representation of the Slater determinant, D_i
        !       (unused, for interface compatibility only).
        !    fpop: Hamiltonian population on the determinant.
        !    f_hfpop: Hellmann-Feynman population on the determinant.
        !    fdata(:): additional information about the determinant (unused, for
        !       interface compatibility only).
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.
        ! In/Out:
        !    D0_hf_population: running total of the Hellmann-Feynman population
        !       on the reference.  Only updated if D_i *is* the reference determinant.
        !    proj_hf_O_hpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|O|D_0> N_i, where N_i is the Hamiltonian population on D_i.
        !    proj_hf_H_fhpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|H|D_0> \tilde{N}_i, where \tilde{N}_i is the
        !       Hellmann-Feynman population on D_i.

        use basis, only: basis_length
        use excitations, only: excit
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip

        ! Note: two-electron operator.

        select case(excitation%nexcit)
        case(0)
            ! Have reference determinant.
            D0_hf_pop = D0_hf_pop + f_hfpop
        case(2)
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.

            !\hat{O}_0j = H_0j / (U L), where L is the number of sites.
            ! sampling \hat{O} - <D0|O|D0>, this means that \sum_j O_j0 c_j = 0.
            proj_hf_O_hpsip = proj_hf_O_hpsip + (hmatel/(sys%hubbard%u*sys%lattice%nsites))*fpop

            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        end select

    end subroutine update_proj_hfs_double_occ_hub_k

    subroutine update_proj_hfs_one_body_mol(sys, f, fpop, f_hfpop, fdata, excitation, hmatel, &
                                              D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)

        ! Add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimator.

        ! This procedure is for when we are sampling O_1, a one-body operator in
        ! a molecular system (i.e. where the integrals have been read in).

        ! In:
        !    sys: system being studied.  Unused.
        !    f(basis_length): bit string representation of the Slater determinant, D_i
        !       (unused, for interface compatibility only).
        !    fpop: Hamiltonian population on the determinant.
        !    f_hfpop: Hellmann-Feynman population on the determinant.
        !    fdata(:): additional information about the determinant (unused, for
        !       interface compatibility only).
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.
        ! In/Out:
        !    D0_hf_population: running total of the Hellmann-Feynman population
        !       on the reference.  Only updated if D_i *is* the reference determinant.
        !    proj_hf_O_hpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|O_1|D_0> N_i, where N_i is the Hamiltonian population on D_i.
        !    proj_hf_H_fhpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|H|D_0> \tilde{N}_i, where \tilde{N}_i is the
        !       Hellmann-Feynman population on D_i.

        use basis, only: basis_length
        use excitations, only: excit
        use operators, only: one_body1_mol
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip

        real(p) :: matel

        ! Note: one-electron operator.

        select case(excitation%nexcit)
        case(0)
            ! Have reference determinant.
            D0_hf_pop = D0_hf_pop + f_hfpop
        case(1)
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.

            ! \sum_j O_0j c_j
            matel = one_body1_mol(sys, excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
            proj_hf_O_hpsip = proj_hf_O_hpsip + matel*fpop

            ! \sum_j H_0j \tilde{c}_j
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        case(2)
            ! O is a one-body operator => no contributions from double
            ! excitations to \sum_j O_0j c_j.

            ! \sum_j H_0j \tilde{c}_j
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        end select

    end subroutine update_proj_hfs_one_body_mol

end module energy_evaluation
