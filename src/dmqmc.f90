module dmqmc

! Main loop for performing DMQMC calculations, similar to the
! file fciqmc.f90 which carries out the main loop in FCIQMC
! calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine do_dmqmc()

        ! Run DMQMC calculation. We run from a beta=0 to a value of beta
        ! specified by the user and then repeat this main loop beta_loops
        ! times, to accumulate statistics for each value for beta.

        use parallel
        use annihilation, only: direct_annihilation
        use basis, only: basis_length, bit_lookup, nbasis
        use death, only: stochastic_death
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        use dmqmc_procedures, only: random_distribution_heisenberg
        use dmqmc_procedures, only: update_sampling_weights, output_and_alter_weights
        use dmqmc_estimators, only: update_dmqmc_estimators, call_dmqmc_estimators
        use dmqmc_estimators, only: call_rdm_procedures
        use excitations, only: excit, get_excitation_level
        use fciqmc_restart, only: dump_restart, write_restart_file_every_nreports
        use qmc_common
        use interact, only: fciqmc_interact
        use system, only: nel
        use calc, only: seed, doing_dmqmc_calc, dmqmc_energy, initiator_approximation
        use calc, only: dmqmc_staggered_magnetisation, dmqmc_energy_squared
        use system, only: system_type, heisenberg
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use utils, only: int_fmt
        use errors, only: stop_all

        integer :: idet, ireport, icycle, iparticle, iteration, ireplica
        integer :: beta_cycle
        integer(lint) :: nparticles_old(sampling_size)
        integer(lint) :: nattempts
        type(det_info) :: cdet1, cdet2
        integer :: nspawned, ndeath
        type(excit) :: connection
        integer :: spawning_end
        logical :: soft_exit
        real :: t1, t2
        type(dSFMT_t) :: rng

        ! Allocate det_info components. We need two cdet objects
        ! for each 'end' which may be spawned from in the DMQMC algorithm.
        call alloc_det_info(cdet1, .false.)
        call alloc_det_info(cdet2, .false.)

        ! Main DMQMC loop.
        if (parent) call write_fciqmc_report_header()
        ! Initialise timer.
        call cpu_time(t1)

        initial_shift = shift
        ! When we accumulate data throughout a run, we are actually accumulating
        ! results from the psips distribution from the previous iteration.
        ! For example, in the first iteration, the trace calculated will be that
        ! of the initial distribution, which corresponds to beta=0. Hence, in the
        ! output we subtract one from the iteration number, and run for one more
        ! report loop, asimplemented in the line of code below.
        nreport = nreport+1

        do beta_cycle = 1, beta_loops
            ! Reset the current position in the spawning array to be the
            ! slot preceding the first slot.
            qmc_spawn%head = qmc_spawn%head_start
            tot_walkers = 0
            shift = initial_shift
            nparticles = 0
            if (allocated(reduced_density_matrix)) reduced_density_matrix = 0.0_p
            if (dmqmc_vary_weights) dmqmc_accumulated_probs = 1.0_p
            if (dmqmc_find_weights) excit_distribution = 0
            vary_shift = .false.

            if (beta_cycle .ne. 1 .and. parent) then
                write (6,'(a32,i7)') &
                       " # Resetting beta... Beta loop =", beta_cycle
                write (6,'(a52,'//int_fmt(seed,1)//',a1)') &
                    " # Resetting random number generator with a seed of:", seed+iproc+beta_cycle-1, "."
            end if
            ! Reset the random number generator with seed = seed + 1 (each
            ! iteration)
            call dSFMT_init(seed+iproc+beta_cycle-1, 50000, rng)

            ! Need to place psips randomly along the diagonal at the
            ! start of every iteration. Pick orbitals randomly, each
            ! with equal probability, so that when electrons are placed
            ! on these orbitals they will have the correct spin and symmetry.
            ! Initial particle distribution.
            select case(system_type)
            case(heisenberg)
                call random_distribution_heisenberg(rng)
            case default
                call stop_all('init_proc_pointers','DMQMC not implemented for this system.')
            end select

            call direct_annihilation(initiator_approximation)

            nparticles_old = nint(D0_population)

            do ireport = 1, nreport
                ! Zero report cycle quantities.
                rspawn = 0.0_p
                trace = 0.0_p
                estimator_numerators = 0.0_p
                if (calculate_excit_distribution) excit_distribution = 0

                do icycle = 1, ncycles
                    qmc_spawn%head = qmc_spawn%head_start
                    iteration = (ireport-1)*ncycles + icycle

                    ! Number of spawning attempts that will be made.
                    ! Each particle and each end gets to attempt to
                    ! spawn onto a connected determinant and a chance
                    ! to die/clone.
                    nattempts = 4*nparticles(1)*sampling_size

                    ! Reset death counter.
                    ndeath = 0

                    do idet = 1, tot_walkers ! loop over walkers/dets

                        ! f points to the bitstring that is spawning, f2 to the
                        ! other bit string.
                        cdet1%f => walker_dets(:basis_length,idet)
                        cdet1%f2 => walker_dets((basis_length+1):(2*basis_length),idet)
                        cdet2%f => walker_dets((basis_length+1):(2*basis_length),idet)
                        cdet2%f2 => walker_dets(:basis_length,idet)

                        ! Decode and store the the relevant information for
                        ! both bitstrings. Both of these bitstrings are required
                        ! to refer to the correct element in the density matrix.
                        call decoder_ptr(cdet1%f, cdet1)
                        call decoder_ptr(cdet2%f, cdet2)

                        ! Call wrapper function which calls all requested estimators
                        ! to be updated, and also always updates the trace separately.
                        if (icycle == 1) call call_dmqmc_estimators(idet, iteration)

                        do ireplica = 1, sampling_size
                            do iparticle = 1, abs(walker_population(ireplica,idet))
                                ! Spawn from the first end.
                                spawning_end = 1
                                ! Attempt to spawn.
                                call spawner_ptr(rng, cdet1, walker_population(ireplica,idet), gen_excit_ptr, nspawned, connection)
                                ! Spawn if attempt was successful.
                                if (nspawned /= 0) then
                                    call create_spawned_particle_dm_ptr(cdet1%f, cdet2%f, connection, nspawned, spawning_end, &
                                                                        ireplica, qmc_spawn)
                                end if

                                ! Now attempt to spawn from the second end.
                                spawning_end = 2
                                call spawner_ptr(rng, cdet2, walker_population(ireplica,idet), gen_excit_ptr, nspawned, connection)
                                if (nspawned /= 0) then
                                    call create_spawned_particle_dm_ptr(cdet2%f, cdet1%f, connection, nspawned, spawning_end, &
                                                                        ireplica, qmc_spawn)
                                end if
                            end do

                            ! Clone or die.
                            ! We have contributions to the clone/death step from both ends of the
                            ! current walker. We do both of these at once by using walker_data(:,idet)
                            ! which, when running a DMQMC algorithm, stores the average of the two diagonal
                            ! elements corresponding to the two indicies of the density matrix (the two ends).
                            call stochastic_death(rng, walker_data(ireplica,idet), walker_population(ireplica,idet), &
                                                  nparticles(ireplica), ndeath)
                        end do
                    end do

                    ! Add the spawning rate (for the processor) to the running total.
                    rspawn = rspawn + spawning_rate(ndeath, nattempts)

                    ! Perform the annihilation step where the spawned walker list is merged with the
                    ! main walker list, and walkers of opposite sign on the same sites are annihilated.
                    call direct_annihilation(initiator_approximation)

                    ! If doing importance sampling *and* varying the weights of the trial function, call a routine
                    ! to update these weights and alter the number of psips on each excitation level accordingly.
                    if (dmqmc_vary_weights .and. iteration <= finish_varying_weights) call update_sampling_weights(rng)

                end do

                ! If averaging the shift to use in future beta loops, add contirubtion from this report.
                if (average_shift_until > 0) shift_profile(ireport) = shift_profile(ireport) + shift

                ! Update the shift and desired thermal quantites.
                call update_dmqmc_estimators(nparticles_old, ireport)

                call cpu_time(t2)

                ! t1 was the time at the previous iteration, t2 the current time.
                ! t2-t1 is thus the time taken by this report loop.
                if (parent) call write_fciqmc_report(ireport, nparticles_old, t2-t1, .false.)
                ! Write restart file if required.
                if (mod(ireport,write_restart_file_every_nreports) == 0) &
                    call dump_restart(mc_cycles_done+ncycles*ireport, nparticles_old)

                ! cpu_time outputs an elapsed time, so update the reference timer.
                t1 = t2

                call fciqmc_interact(soft_exit)
                if (soft_exit) exit

            end do

            if (soft_exit) exit

            ! If have just finished last beta loop of accumulating the shift, then perform
            ! the averaging and set average_shift_until to -1. This tells the shift update
            ! algorithm to use the values for shift stored in shift_profile.
            if (beta_cycle == average_shift_until) then
                shift_profile = shift_profile/average_shift_until
                average_shift_until = -1
            end if

            ! Calculate and output all requested estimators based on the reduced dnesity matrix.
            if (doing_reduced_dm) call call_rdm_procedures(beta_cycle)
            ! Calculate and output new weights based on the psip distirubtion in the previous loop.
            if (dmqmc_find_weights) call output_and_alter_weights()

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

        if (dump_restart_file) call dump_restart(mc_cycles_done, nparticles_old, vspace=.true.)

        call dealloc_det_info(cdet1, .false.)
        call dealloc_det_info(cdet2, .false.)

    end subroutine do_dmqmc

end module dmqmc
