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

    subroutine do_hfs_fciqmc(sys, qmc_in, restart_in, load_bal_in, reference_in)

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
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.
        !    reference_in: current reference determinant.  If not set (ie
        !       components allocated) then a best guess is made based upon the
        !       desired spin/symmetry.
        ! In/Out:
        !    qmc_in: input options relating to QMC methods.

        use parallel
        use checking, only: check_allocate

        use annihilation, only: direct_annihilation
        use death, only: stochastic_death, stochastic_hf_cloning
        use determinants, only: alloc_det_info_t, dealloc_det_info_t, sum_fock_values_occ_list
        use determinant_data, only: det_info_t
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit_t, get_excitation
        use qmc_io, only: write_qmc_report_header, write_qmc_report
        use hfs_data
        use ifciqmc
        use interact, only: calc_interact, check_comms_file
        use qmc_common
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use proc_pointers
        use qmc, only: init_qmc
        use system, only: sys_t
        use restart_hdf5, only: dump_restart_hdf5, restart_info_t, init_restart_info_t
        use qmc_data, only: qmc_in_t, restart_in_t, load_bal_in_t, qmc_state_t, annihilation_flags_t
        use logging, only: logging_t
        use reference_determinant, only: reference_t
        use hamiltonian_data
        use energy_evaluation, only: get_sanitized_projected_energy

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(in) :: reference_in

        integer :: idet, ireport, icycle, hf_initiator_flag, h_initiator_flag, ierr
        integer(int_p) :: iparticle
        integer(int_64) :: nattempts
        real(dp), allocatable :: nparticles_old(:)
        real(p), allocatable :: real_population(:)
        type(det_info_t) :: cdet

        integer(int_p) :: nspawned, ndeath, dummy
        integer :: nspawn_events
        type(excit_t) :: connection
        type(dSFMT_t) :: rng
        type(hmatel_t) :: hmatel
        type(excit_t), parameter :: null_excit = excit_t( 0, [0,0], [0,0], .false.)
        type(qmc_state_t), target :: qs
        type(annihilation_flags_t) :: annihilation_flags
        type(restart_info_t) :: ri
        character(36) :: uuid_restart
        type(logging_t) :: logging_info

        logical :: soft_exit, comms_found, error, restart_proj_est

        real :: t1, t2
        integer :: iunit, restart_version_restart

        iunit = 6

        if (parent) then
            write (iunit,'(1X,"FCIQMC (with Hellmann-Feynman sampling")')
            write (iunit,'(1X,"--------------------------------------",/)')
        end if

        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, 6, annihilation_flags, qs, uuid_restart, &
                        restart_version_restart)

        allocate(nparticles_old(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('nparticles_old', qs%psip_list%nspaces, ierr)
        allocate(real_population(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('real_population', qs%psip_list%nspaces, ierr)

        call dSFMT_init(qmc_in%seed+iproc, 50000, rng)

        ! Allocate det_info_t components.
        call alloc_det_info_t(sys, cdet, .false.)

        ! from restart
        nparticles_old = qs%psip_list%tot_nparticles

        ! Main fciqmc loop.

        if (parent) call write_qmc_report_header(qs%psip_list%nspaces)
        restart_proj_est = (restart_in%read_restart .and. restart_version_restart >= 2)
        if (.not.restart_proj_est) call initial_ci_projected_energy(sys, qs, .false., nparticles_old)
        call initial_qmc_status(sys, qmc_in, qs, nparticles_old, .false.)

        ! Initialise timer.
        call cpu_time(t1)

        call init_restart_info_t(ri, write_id=restart_in%write_id)

        do ireport = 1, qmc_in%nreport

            qs%estimators%proj_energy_old = get_sanitized_projected_energy(qs)
            ! Zero report cycle quantities.
            qs%estimators%proj_energy = 0.0_p
            qs%estimators%proj_hf_O_hpsip = 0.0_p
            qs%estimators%proj_hf_H_hfpsip = 0.0_p
            qs%estimators%D0_population = 0.0_p
            qs%estimators%D0_hf_population = 0.0_p
            qs%spawn_store%rspawn = 0.0_p

            do icycle = 1, qmc_in%ncycles

                ! Reset the current position in the spawning array to be the
                ! slot preceding the first slot.
                qs%spawn_store%spawn%head = qs%spawn_store%spawn%head_start

                ! Number of spawning attempts that will be made.
                ! Each Hamiltonian particle gets a chance to spawn a Hamiltonian
                ! particle, clone/die, spawn a Hellmann-Feynman particle and clone
                ! itself into a Hellmann-Feynman particle.  Each H-F particle
                ! gets a chance to spawn and a chance to clone/die.
                ! This is used for accounting later, not for controlling the spawning.
                nattempts = nint(4*qs%psip_list%nparticles(1) + 2*qs%psip_list%nparticles(2))

                ! Reset death counter.
                ndeath = 0_int_p

                do idet = 1, qs%psip_list%nstates ! loop over walkers/dets

                    cdet%f = qs%psip_list%states(:,idet)
                    cdet%data => qs%psip_list%dat(:,idet)

                    call decoder_ptr(sys, cdet%f, cdet)
                    if (qs%propagator%quasi_newton) &
                        cdet%fock_sum = sum_fock_values_occ_list(sys, qs%propagator%sp_fock, cdet%occ_list) - qs%ref%fock_sum

                    ! Extract the real sign from the encoded sign.
                    real_population = real(qs%psip_list%pops(1,idet),p)/qs%psip_list%pop_real_factor

                    ! It is much easier to evaluate projected values at the
                    ! start of the FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    connection = get_excitation(sys%nel, sys%basis, cdet%f, qs%ref%f0)
                    call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, real_population,  &
                                                qs%estimators(1), connection, hmatel)
                    ! [todo] - JSS: pass real populations through to HFS projected energy update
                    call update_proj_hfs_ptr(sys, cdet%f, int(qs%psip_list%pops(1,idet)),&
                                             int(qs%psip_list%pops(2,idet)), cdet%data,  &
                                             connection, hmatel, qs%estimators(1)%D0_hf_population,  &
                                             qs%estimators(1)%proj_hf_O_hpsip, qs%estimators(1)%proj_hf_H_hfpsip)

                    ! Is this determinant an initiator?
                    ! A determinant can be an initiator in the Hamiltonian space
                    ! or the Hellmann-Feynman space or both.
                    ! The initiator_flag attribute of det_info_t is checked and passed to the
                    ! annihilation routine in the appropriate create_spawned_particle_*
                    ! routine, so we must set cdet%initiator_flag
                    ! appropriately...
                    call set_parent_flag([real_population(1)], qmc_in%initiator_pop, 1, .true., &
                                            h_initiator_flag)
                    call set_parent_flag([real_population(2)], qmc_in%initiator_pop, 1, .true., &
                                            hf_initiator_flag)
                    cdet%initiator_flag = h_initiator_flag

                    do iparticle = 1, abs(qs%psip_list%pops(1,idet))

                        ! Attempt to spawn Hamiltonian walkers..
                        call spawner_ptr(rng, sys, qs, qs%spawn_store%spawn%cutoff, qs%psip_list%pop_real_factor, &
                                         cdet, qs%psip_list%pops(1,idet), gen_excit_ptr, qs%trial%wfn_dat, &
                                         logging_info, nspawned, dummy, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0_int_p) then
                            associate(spawn=>qs%spawn_store%spawn)
                                call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, connection, nspawned, 1, spawn)
                            end associate
                        end if

                        ! Attempt to spawn Hellmann--Feynman walkers from
                        ! Hamiltonian walkers.
                        ! [todo] - JSS: real populations for HFS spawner.
                        call spawner_hfs_ptr(rng, sys, qs, qs%spawn_store%spawn%cutoff, qs%psip_list%pop_real_factor, &
                                             cdet, qs%psip_list%pops(1,idet), gen_excit_hfs_ptr, qs%trial%wfn_dat, logging_info, &
                                             nspawned, dummy, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0_int_p) then
                            associate(spawn=>qs%spawn_store%spawn)
                                call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, connection, nspawned, 2, spawn)
                            end associate
                        end if
                    end do

                    cdet%initiator_flag = hf_initiator_flag

                    do iparticle = 1, abs(qs%psip_list%pops(2,idet))

                        ! Attempt to spawn Hellmann--Feynman walkers from
                        ! Hellmann--Feynman walkers.
                        call spawner_ptr(rng, sys, qs, qs%spawn_store%spawn%cutoff, qs%psip_list%pop_real_factor, &
                                         cdet, qs%psip_list%pops(2,idet), gen_excit_ptr, qs%trial%wfn_dat, logging_info, &
                                         nspawned, dummy, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0_int_p) then
                            associate(spawn=>qs%spawn_store%spawn)
                                call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, connection, nspawned, 2, spawn)
                             end associate
                         end if

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
                    call stochastic_death(rng, sys, qs, cdet%fock_sum, qs%psip_list%dat(2,idet), qs%shift(1), &
                                          qs%estimators(1)%proj_energy_old, logging_info, qs%psip_list%pops(2,idet), &
                                          qs%psip_list%nparticles(2), ndeath)

                    ! Clone Hellmann--Feynman walkers from Hamiltonian walkers.
                    ! Not in place, must set initiator flag.
                    cdet%initiator_flag = h_initiator_flag
                    ! [todo] - JSS: real populations for HFS spawner.
                    call stochastic_hf_cloning(rng, qs%tau, qs%shift(2), qs%psip_list%dat(2,idet), &
                                               qs%psip_list%pops(1,idet), nspawned)
                    if (nspawned /= 0) then
                        associate(spawn=>qs%spawn_store%spawn)
                            call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, null_excit, nspawned, 2, spawn)
                        end associate
                    end if

                    ! Clone or die: Hamiltonian walkers.
                    call stochastic_death(rng, sys, qs, cdet%fock_sum, qs%psip_list%dat(1,idet), qs%shift(1), &
                                          qs%estimators(1)%proj_energy_old, logging_info, qs%psip_list%pops(1,idet), &
                                          qs%psip_list%nparticles(1), ndeath)

                end do

                call direct_annihilation(sys, rng, qs%ref, annihilation_flags, qs%psip_list, qs%spawn_store%spawn, nspawn_events)

            end do

            ! Test for a comms file so MPI communication can be combined with
            ! energy_estimators communication
            comms_found = check_comms_file()
            ! Update the energy estimators (shift & projected energy).
            error = qs%spawn_store%spawn%error .or. qs%psip_list%error
            call update_energy_estimators(qmc_in, qs, nspawn_events, nparticles_old, load_bal_in=load_bal_in, &
                                          comms_found=comms_found, error=error)
            if (error) exit

            call cpu_time(t2)

            ! t1 was the time at the previous iteration, t2 the current time.
            ! t2-t1 is thus the time taken by this report loop.
            if (parent) call write_qmc_report(qmc_in, qs, ireport, nparticles_old, t2-t1, .false., .false.)

            ! Write restart file if required.
!            if (mod(ireport,write_restart_file_every_nreports) == 0) &
!                call dump_restart(mc_cycles_done+qmc_in%ncycles*ireport, nparticles_old)

            ! cpu_time outputs an elapsed time, so update the reference timer.
            t1 = t2

            call calc_interact(comms_found, iunit, soft_exit, qs)
            if (soft_exit) exit

        end do

        if (parent) write (iunit,'()')
        call load_balancing_report(qs%psip_list%nparticles, qs%psip_list%nstates, qmc_in%use_mpi_barriers, &
                                   qs%spawn_store%spawn%mpi_time)

        if (soft_exit .or. error) then
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*ireport
        else
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*qmc_in%nreport
        end if

        if (restart_in%write_restart) then
            call dump_restart_hdf5(qs, qs%mc_cycles_done, nparticles_old, sys%basis%nbasis, .false., sys%basis%info_string_len)
            if (parent) write (iunit,'()')
        end if

        call dealloc_det_info_t(cdet, .false.)

    end subroutine do_hfs_fciqmc

end module hellmann_feynman_sampling
