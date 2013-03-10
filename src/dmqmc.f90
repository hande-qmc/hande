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
        use dmqmc_estimators, only: update_dmqmc_estimators, call_dmqmc_estimators
        use excitations, only: excit, get_excitation_level
        use fciqmc_restart, only: dump_restart, write_restart_file_every_nreports
        use qmc_common
        use interact, only: fciqmc_interact
        use system, only: nel
        use calc, only: seed, doing_dmqmc_calc, dmqmc_energy
        use calc, only: dmqmc_staggered_magnetisation, dmqmc_energy_squared
        use dSFMT_interface, only: dSFMT_init
        use utils, only: int_fmt

        integer :: idet, ireport, icycle, iparticle
        integer :: beta_cycle
        integer(lint) :: nparticles_old(sampling_size)
        integer(lint) :: nattempts, nparticles_start_report
        type(det_info) :: cdet1, cdet2
        integer :: nspawned, ndeath
        type(excit) :: connection
        integer :: spawning_end
        logical :: soft_exit
        real :: t1, t2

        ! Allocate det_info components. We need two cdet objects
        ! for each 'end' which may be spawned from in the DMQMC algorithm.
        call alloc_det_info(cdet1)
        call alloc_det_info(cdet2)

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
            spawning_head = spawning_block_start
            tot_walkers = 0
            shift = initial_shift
            nparticles = 0
            vary_shift = .false.

            ! Need to place psips randomly along the diagonal at the
            ! start of every iteration. Pick orbitals randomly, each
            ! with equal probability, so that when electrons are placed
            ! on these orbitals they will have the correct spin and symmetry.
            call dmqmc_initial_distribution_ptr()

            call direct_annihilation()

            if (beta_cycle .ne. 1 .and. parent) then
                write (6,'(a32,i7)') &
                       " # Resetting beta... Beta loop =", beta_cycle
                ! Reset the random number generator with seed = seed + 1
                seed = seed + 1
                call dSFMT_init(seed + iproc)
                write (6,'(a52,'//int_fmt(seed,1)//',a1)') &
                    " # Resetting random number generator with a seed of:", seed, "."
            end if

            nparticles_old = nint(D0_population)

            do ireport = 1, nreport

                ! Zero report cycle quantities.
                rspawn = 0.0_p
                trace = 0
                estimator_numerators = 0
                nparticles_start_report = nparticles_old(1)

                do icycle = 1, ncycles
                    spawning_head = spawning_block_start

                    ! Number of spawning attempts that will be made.
                    ! Each particle and each end gets to attempt to
                    ! spawn onto a connected determinant and a chance
                    ! to die/clone.
                    nattempts = 4*nparticles(1)

                    ! Reset death counter
                    ndeath = 0

                    do idet = 1, tot_walkers ! loop over walkers/dets
                        cdet1%f = walker_dets(:basis_length,idet)
                        cdet2%f = walker_dets((basis_length+1):(2*basis_length),idet)

                        ! Decode and store the the relevant information for
                        ! both bitstrings. Both of these bitstrings are required
                        ! to refer to the correct element in the density matrix.
                        call decoder_ptr(cdet1%f, cdet1)
                        call decoder_ptr(cdet2%f, cdet2)

                        ! Call wrapper function which calls all requested estimators
                        ! to be updated, and also always updates the trace separately.
                        if (icycle == 1) call call_dmqmc_estimators(idet)

                        do iparticle = 1, abs(walker_population(1,idet))
                            ! Spawn from the first end.
                            spawning_end = 1
                            ! Attempt to spawn.
                            call spawner_ptr(cdet1, walker_population(1,idet), nspawned, connection)
                            ! Spawn if attempt was successful.
                            if (nspawned /= 0) then
                                call create_spawned_particle_dm_ptr(cdet1%f, cdet2%f, connection, nspawned, spawning_end)
                            end if

                            ! Now attempt to spawn from the second end.
                            spawning_end = 2
                            call spawner_ptr(cdet2, walker_population(1,idet), nspawned, connection)
                            if (nspawned /= 0) then
                                call create_spawned_particle_dm_ptr(cdet2%f, cdet1%f, connection, nspawned, spawning_end)
                            end if
                        end do

                        ! Clone or die.
                        ! We have contirbutions to the clone/death step from both ends of the
                        ! current walker. We do both of these at once by using walker_data(1,idet)
                        ! which, when running a DMQMC algorithm, stores the average of the two diagonal
                        ! elements corresponding to the two indicies of the density matrix (the two ends).
                        call stochastic_death(walker_data(1,idet), walker_population(1,idet), nparticles(1), ndeath)
                    end do

                    ! Add the spawning rate (for the processor) to the running
                    ! total.
                    rspawn = rspawn + spawning_rate(ndeath, nattempts)

                    ! Perform the annihilation step where the spawned walker list is merged with
                    ! the main walker list, and walkers of opposite sign on the same sites are
                    ! annihilated.
                    call direct_annihilation()

                end do

                old_shift=shift
                ! Update the shift and desired thermal quantites.
                call update_dmqmc_estimators(nparticles_old)

                call cpu_time(t2)

                ! t1 was the time at the previous iteration, t2 the current time.
                ! t2-t1 is thus the time taken by this report loop.
                if (parent) call write_fciqmc_report(ireport, nparticles_start_report, t2-t1)
                ! Write restart file if required.
                if (mod(ireport,write_restart_file_every_nreports) == 0) &
                    call dump_restart(mc_cycles_done+ncycles*ireport, nparticles_old(1))

                ! cpu_time outputs an elapsed time, so update the reference timer.
                t1 = t2

                call fciqmc_interact(soft_exit)
                if (soft_exit) exit

            end do

        end do

        if (parent) then
            call write_fciqmc_final(ireport)
            write (6,'()')
        end if

        call load_balancing_report()

        if (dump_restart_file) call dump_restart(mc_cycles_done+ncycles*nreport, nparticles_old(1))

        call dealloc_det_info(cdet1)
        call dealloc_det_info(cdet2)

    end subroutine do_dmqmc

end module dmqmc
