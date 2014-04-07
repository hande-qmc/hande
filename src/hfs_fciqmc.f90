module hellmann_feynman_sampling

! Module for performing Hellmann--Feynman sampling in FCIQMC in order to obtain
! expectation values of arbitrary operators which do not commute with the
! Hamiltonian.

! This involves sampling a 'pumped' diffusion equation in an adjoint space to
! that used to sample the Hamiltonian.  Fortunately we can use essentially the
! same dynamics to sample the adjoint space.  The connection between the
! Hamiltonian space and the adjoint 'Hellmann--Feynman' space is defined by the
! operator being sampled.

! See documentation/theory/hellmann_feynman/hf.tex for details.

use const

implicit none

contains

    subroutine do_hfs_fciqmc(sys)

        ! Run the FCIQMC algorithm starting from the initial walker
        ! distribution and perform Hellmann--Feynman sampling in conjunction on
        ! the instantaneous wavefunction.

        ! This is implemented using F2003 function pointers.  This allows us to
        ! avoid many system dependent if blocks, which are constant for a given
        ! calculation.  Avoiding such branching is worth the extra verbosity
        ! (especially if procedures are written to be sufficiently modular that
        ! implementing a new system can reuse many existing routines) as it
        ! leads to much faster code.

        ! In:
        !    sys: system being studied.

        use parallel

        use annihilation, only: direct_annihilation
        use basis, only: basis_length
        use calc, only: seed, initiator_approximation
        use death, only: stochastic_hf_cloning
        use determinants, only:det_info, alloc_det_info, dealloc_det_info
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit
        use fciqmc_data, only: shift
        use hfs_data
        use interact, only: fciqmc_interact
        use qmc_common
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use utils, only: rng_init_info
        use proc_pointers
        use system, only: sys_t
        use restart_hdf5, only: restart_info_global, dump_restart_hdf5

        type(sys_t), intent(in) :: sys

        integer :: idet, ireport, icycle, iparticle, hf_initiator_flag, h_initiator_flag
        integer(lint) :: nattempts, nparticles_old(sampling_size)
        type(det_info) :: cdet

        integer :: nspawned, ndeath
        type(excit) :: connection
        type(dSFMT_t) :: rng
        real(p) :: hmatel
        type(excit), parameter :: null_excit = excit( 0, [0,0,0,0], [0,0,0,0], .false.)

        logical :: soft_exit

        real :: t1, t2

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

        ! Allocate det_info components.
        call alloc_det_info(sys, cdet, .false.)

        ! from restart
        nparticles_old = tot_nparticles

        ! Main fciqmc loop.

        if (parent) call write_fciqmc_report_header()
        call initial_fciqmc_status(sys)

        ! Initialise timer.
        call cpu_time(t1)

        do ireport = 1, nreport

            ! Zero report cycle quantities.
            proj_energy = 0.0_p
            proj_hf_O_hpsip = 0.0_p
            proj_hf_H_hfpsip = 0.0_p
            D0_population = 0.0_p
            D0_hf_population = 0.0_p
            rspawn = 0.0_p

            do icycle = 1, ncycles

                ! Reset the current position in the spawning array to be the
                ! slot preceding the first slot.
                qmc_spawn%head = qmc_spawn%head_start

                ! Number of spawning attempts that will be made.
                ! Each Hamiltonian particle gets a chance to spawn a Hamiltonian
                ! particle, clone/die, spawn a Hellmann-Feynman particle and clone
                ! itself into a Hellmann-Feynman particle.  Each H-F particle
                ! gets a chance to spawn and a chance to clone/die.
                ! This is used for accounting later, not for controlling the spawning.
                nattempts = 4*nparticles(1) + 2*nparticles(2)

                ! Reset death counter.
                ndeath = 0

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)
                    cdet%data => walker_data(:,idet)

                    call decoder_ptr(sys, cdet%f, cdet)

                    ! It is much easier to evaluate projected values at the
                    ! start of the FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    call update_proj_energy_ptr(sys, f0, cdet, real(walker_population(1,idet),p),  &
                                                D0_population_cycle, proj_energy, connection, hmatel)
                    call update_proj_hfs_ptr(sys, cdet%f, walker_population(1,idet),&
                                             walker_population(2,idet), cdet%data,  &
                                             connection, hmatel, D0_hf_population,  &
                                             proj_hf_O_hpsip, proj_hf_H_hfpsip)

                    ! Is this determinant an initiator?
                    ! A determinant can be an initiator in the Hamiltonian space
                    ! or the Hellmann-Feynman space or both.
                    ! The initiator_flag attribute of det_info is checked and passed to the
                    ! annihilation routine in the appropriate create_spawned_particle_*
                    ! routine, so we must set cdet%initiator_flag
                    ! appropriately...
                    call set_parent_flag_ptr(walker_population(1,idet), cdet%f, h_initiator_flag)
                    call set_parent_flag_ptr(walker_population(2,idet), cdet%f, hf_initiator_flag)
                    cdet%initiator_flag = h_initiator_flag

                    do iparticle = 1, abs(walker_population(1,idet))

                        ! Attempt to spawn Hamiltonian walkers..
                        call spawner_ptr(rng, sys, qmc_spawn, cdet, walker_population(1,idet), gen_excit_ptr, nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) call create_spawned_particle_ptr(cdet, connection, nspawned, 1, qmc_spawn)

                        ! Attempt to spawn Hellmann--Feynman walkers from
                        ! Hamiltonian walkers.
                        call spawner_hfs_ptr(rng, sys, cdet, walker_population(1,idet), &
                                             gen_excit_hfs_ptr, nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) call create_spawned_particle_ptr(cdet, connection, nspawned, 2, qmc_spawn)

                    end do

                    cdet%initiator_flag = hf_initiator_flag

                    do iparticle = 1, abs(walker_population(2,idet))

                        ! Attempt to spawn Hellmann--Feynman walkers from
                        ! Hellmann--Feynman walkers.
                        call spawner_ptr(rng, sys, qmc_spawn, cdet, walker_population(2,idet), gen_excit_ptr, nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) call create_spawned_particle_ptr(cdet, connection, nspawned, 2, qmc_spawn)

                    end do

                    ! Now deal with the diagonal events.

                    ! *IMPORTANT*
                    ! My mother told me that it is important to play nice, be
                    ! fair and not beat up younger brothers.  Monte Carlo is
                    ! a particularly fickle playmate and one must follow her
                    ! advice in order to avoid biased (a polite way of saying
                    ! wrong) results.

                    ! In this instance, we must give all particles the same
                    ! opportunities.  In particular:
                    ! * Hamiltonian walkers that existed at the beginning of the
                    !   loop *must* have an opportunity to create Hellmann--Feynman 
                    !   versions of themselves.
                    ! * Hamiltonian walkers that existed at the beginning of the
                    !   loop *must* have an opportunity to die/clone themselves.
                    ! * Hellmann--Feynman walkers that existed at the beginning
                    !   of the loop must have an opportunity to die/clone
                    !   themselves.

                    ! As the death/clone routines act in place, we must ensure
                    ! the above requirements are met and (e.g.) walkers that are
                    ! created don't get an additional death/cloning opportunity.

                    ! Clone or die: Hellmann--Feynman walkers.
                    call death_ptr(rng, walker_data(1,idet), shift(1), walker_population(2,idet), nparticles(2), ndeath)

                    ! Clone Hellmann--Feynman walkers from Hamiltonian walkers.
                    ! Not in place, must set initiator flag.
                    cdet%initiator_flag = h_initiator_flag
                    call stochastic_hf_cloning(rng, walker_data(2,idet), walker_population(1,idet), nspawned)
                    if (nspawned /= 0) call create_spawned_particle_ptr(cdet, null_excit, nspawned, 2, qmc_spawn)

                    ! Clone or die: Hamiltonian walkers.
                    call death_ptr(rng, walker_data(1,idet), shift(1), walker_population(1,idet), nparticles(1), ndeath)

                end do

                ! Add the spawning rate (for the processor) to the running
                ! total.
                rspawn = rspawn + spawning_rate(ndeath, nattempts)

                call direct_annihilation(sys, initiator_approximation)

            end do

            ! Update the energy estimators (shift & projected energy).
            call update_energy_estimators(nparticles_old)

            call cpu_time(t2)

            ! t1 was the time at the previous iteration, t2 the current time.
            ! t2-t1 is thus the time taken by this report loop.
            if (parent) call write_fciqmc_report(ireport, nparticles_old, t2-t1, .false.)

            ! Write restart file if required.
!            if (mod(ireport,write_restart_file_every_nreports) == 0) &
!                call dump_restart(mc_cycles_done+ncycles*ireport, nparticles_old)

            ! cpu_time outputs an elapsed time, so update the reference timer.
            t1 = t2

            call fciqmc_interact(soft_exit)
            if (soft_exit) exit

        end do

        if (parent) then
            call write_fciqmc_final(ireport)
            write (6,'()')
        end if

        call load_balancing_report()

        if (soft_exit) then
            mc_cycles_done = mc_cycles_done + ncycles*ireport
        else
            mc_cycles_done = mc_cycles_done + ncycles*nreport
        end if

        if (dump_restart_file) then
            call dump_restart_hdf5(restart_info_global, mc_cycles_done, nparticles_old)
            if (parent) write (6,'()')
        end if

        call dealloc_det_info(cdet, .false.)

    end subroutine do_hfs_fciqmc

end module hellmann_feynman_sampling
