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

    subroutine do_hfs_fciqmc(sys, qmc_in)

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
        ! In/Out:
        !    qmc_in: input options relating to QMC methods.

        use parallel

        use annihilation, only: direct_annihilation
        use calc, only: seed
        use death, only: stochastic_death, stochastic_hf_cloning
        use determinants, only:det_info_t, alloc_det_info_t, dealloc_det_info_t
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit_t, get_excitation
        use fciqmc_data, only: real_factor
        use hfs_data
        use interact, only: calc_interact, check_comms_file
        use qmc_common
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use utils, only: rng_init_info
        use proc_pointers
        use system, only: sys_t
        use restart_hdf5, only: restart_info_global, dump_restart_hdf5
        use qmc_data, only: qmc_in_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in

        integer :: idet, ireport, icycle, iparticle, hf_initiator_flag, h_initiator_flag
        integer(int_64) :: nattempts
        real(p) :: nparticles_old(sampling_size)
        real(p) :: real_population(sampling_size)
        type(det_info_t) :: cdet

        integer(int_p) :: nspawned, ndeath
        integer :: nspawn_events
        type(excit_t) :: connection
        type(dSFMT_t) :: rng
        real(p) :: hmatel
        type(excit_t), parameter :: null_excit = excit_t( 0, [0,0], [0,0], .false.)

        logical :: soft_exit, comms_found

        real :: t1, t2

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

        ! Allocate det_info_t components.
        call alloc_det_info_t(sys, cdet, .false.)

        ! from restart
        nparticles_old = tot_nparticles

        ! Main fciqmc loop.

        if (parent) call write_fciqmc_report_header()
        call initial_fciqmc_status(sys, qmc_in)

        ! Initialise timer.
        call cpu_time(t1)

        do ireport = 1, qmc_in%nreport

            ! Zero report cycle quantities.
            proj_energy = 0.0_p
            proj_hf_O_hpsip = 0.0_p
            proj_hf_H_hfpsip = 0.0_p
            D0_population = 0.0_p
            D0_hf_population = 0.0_p
            rspawn = 0.0_p

            do icycle = 1, qmc_in%ncycles

                ! Reset the current position in the spawning array to be the
                ! slot preceding the first slot.
                qmc_spawn%head = qmc_spawn%head_start

                ! Number of spawning attempts that will be made.
                ! Each Hamiltonian particle gets a chance to spawn a Hamiltonian
                ! particle, clone/die, spawn a Hellmann-Feynman particle and clone
                ! itself into a Hellmann-Feynman particle.  Each H-F particle
                ! gets a chance to spawn and a chance to clone/die.
                ! This is used for accounting later, not for controlling the spawning.
                nattempts = nint(4*nparticles(1) + 2*nparticles(2))

                ! Reset death counter.
                ndeath = 0_int_p

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)
                    cdet%data => walker_data(:,idet)

                    call decoder_ptr(sys, cdet%f, cdet)

                    ! Extract the real sign from the encoded sign.
                    real_population = real(walker_population(1,idet),p)/real_factor

                    ! It is much easier to evaluate projected values at the
                    ! start of the FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    connection = get_excitation(sys%nel, sys%basis, cdet%f, f0)
                    call update_proj_energy_ptr(sys, f0, cdet, real_population(1),  &
                                                D0_population, proj_energy, connection, hmatel)
                    ! [todo] - JSS: pass real populations through to HFS projected energy update
                    call update_proj_hfs_ptr(sys, cdet%f, int(walker_population(1,idet)),&
                                             int(walker_population(2,idet)), cdet%data,  &
                                             connection, hmatel, D0_hf_population,  &
                                             proj_hf_O_hpsip, proj_hf_H_hfpsip)

                    ! Is this determinant an initiator?
                    ! A determinant can be an initiator in the Hamiltonian space
                    ! or the Hellmann-Feynman space or both.
                    ! The initiator_flag attribute of det_info_t is checked and passed to the
                    ! annihilation routine in the appropriate create_spawned_particle_*
                    ! routine, so we must set cdet%initiator_flag
                    ! appropriately...
                    call set_parent_flag_ptr(real_population(1), qmc_in%initiator_pop, cdet%f, 1, h_initiator_flag)
                    call set_parent_flag_ptr(real_population(2), qmc_in%initiator_pop, cdet%f, 1, hf_initiator_flag)
                    cdet%initiator_flag = h_initiator_flag

                    do iparticle = 1, abs(walker_population(1,idet))

                        ! Attempt to spawn Hamiltonian walkers..
                        call spawner_ptr(rng, sys, qmc_in%tau, qmc_spawn%cutoff, real_factor, cdet, walker_population(1,idet), &
                                         gen_excit_ptr, nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0_int_p) &
                            call create_spawned_particle_ptr(sys%basis, cdet, connection, nspawned, 1, qmc_spawn)

                        ! Attempt to spawn Hellmann--Feynman walkers from
                        ! Hamiltonian walkers.
                        ! [todo] - JSS: real populations for HFS spawner.
                        call spawner_hfs_ptr(rng, sys, qmc_in%tau, qmc_spawn%cutoff, real_factor, cdet, walker_population(1,idet), &
                                             gen_excit_hfs_ptr, nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0_int_p) &
                            call create_spawned_particle_ptr(sys%basis, cdet, connection, nspawned, 2, qmc_spawn)

                    end do

                    cdet%initiator_flag = hf_initiator_flag

                    do iparticle = 1, abs(walker_population(2,idet))

                        ! Attempt to spawn Hellmann--Feynman walkers from
                        ! Hellmann--Feynman walkers.
                        call spawner_ptr(rng, sys, qmc_in%tau, qmc_spawn%cutoff, real_factor, cdet, walker_population(2,idet), &
                                         gen_excit_ptr, nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0_int_p) &
                            call create_spawned_particle_ptr(sys%basis, cdet, connection, nspawned, 2, qmc_spawn)

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
                    call stochastic_death(rng, qmc_in%tau, walker_data(1,idet), shift(1), walker_population(2,idet), &
                                           nparticles(2), ndeath)

                    ! Clone Hellmann--Feynman walkers from Hamiltonian walkers.
                    ! Not in place, must set initiator flag.
                    cdet%initiator_flag = h_initiator_flag
                    ! [todo] - JSS: real populations for HFS spawner.
                    call stochastic_hf_cloning(rng, qmc_in%tau, walker_data(2,idet), walker_population(1,idet), nspawned)
                    if (nspawned /= 0) call create_spawned_particle_ptr(sys%basis, cdet, null_excit, nspawned, 2, qmc_spawn)

                    ! Clone or die: Hamiltonian walkers.
                    call stochastic_death(rng, qmc_in%tau, walker_data(1,idet), shift(1), walker_population(1,idet), &
                                           nparticles(1), ndeath)

                end do

                call direct_annihilation(sys, rng, qmc_in, nspawn_events)

            end do

            ! Test for a comms file so MPI communication can be combined with
            ! energy_estimators communication
            comms_found = check_comms_file()
            ! Update the energy estimators (shift & projected energy).
            call update_energy_estimators(qmc_in, nspawn_events, nparticles_old, comms_found)

            call cpu_time(t2)

            ! t1 was the time at the previous iteration, t2 the current time.
            ! t2-t1 is thus the time taken by this report loop.
            if (parent) call write_fciqmc_report(qmc_in, ireport, nparticles_old, t2-t1, .false.)

            ! Write restart file if required.
!            if (mod(ireport,write_restart_file_every_nreports) == 0) &
!                call dump_restart(mc_cycles_done+qmc_in%ncycles*ireport, nparticles_old)

            ! cpu_time outputs an elapsed time, so update the reference timer.
            t1 = t2

            call calc_interact(comms_found, soft_exit, qmc_in)
            if (soft_exit) exit

        end do

        if (parent) write (6,'()')
        call load_balancing_report(qmc_spawn%mpi_time)

        if (soft_exit) then
            mc_cycles_done = mc_cycles_done + qmc_in%ncycles*ireport
        else
            mc_cycles_done = mc_cycles_done + qmc_in%ncycles*qmc_in%nreport
        end if

        if (dump_restart_file) then
            call dump_restart_hdf5(restart_info_global, mc_cycles_done, nparticles_old)
            if (parent) write (6,'()')
        end if

        call dealloc_det_info_t(cdet, .false.)

    end subroutine do_hfs_fciqmc

end module hellmann_feynman_sampling
