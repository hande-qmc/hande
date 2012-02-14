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
        use dmqmc_estimators, only: call_rdm_procedures
        use excitations, only: excit, get_excitation_level
        use fciqmc_common
        use fciqmc_restart, only: dump_restart
        use interact, only: fciqmc_interact
        use system, only: nel
        use calc, only: seed, doing_dmqmc_calc, dmqmc_energy
        use calc, only: dmqmc_staggered_magnetisation, dmqmc_energy_squared
        use dSFMT_interface, only: dSFMT_init
        use utils, only: int_fmt

        integer :: idet, ireport, icycle, iparticle, iteration
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
            spawning_head = spawning_block_start
            tot_walkers = 0
            shift = initial_shift
            nparticles = 0
            if (allocated(reduced_density_matrix)) reduced_density_matrix = 0
            vary_shift = .false.

            ! Need to place psips randomly along the diagonal at the
            ! start of every iteration. Pick orbitals randomly, each
            ! with equal probability, so that when electrons are placed
            ! on these orbitals they will have the correct spin and symmetry.
            call dmqmc_initial_distribution_ptr()

            call direct_annihilation()

            if (beta_cycle .ne. 1) then
                ! Reset the random number generator with seed = seed + nprocs
                seed = seed + nprocs
                call dSFMT_init(seed + iproc)
                if (parent) then
                    write (6,'(a32,i7)') &
                        " # Resetting beta... Beta loop =", beta_cycle
                    write (6,'(a52,'//int_fmt(seed,1)//',a1)') &
                        " # Resetting random number generator with a seed of:", seed, "."
                end if
            end if

            nparticles_old = nint(D0_population)

            do ireport = 1, nreport
                ! Zero report cycle quantities.
                rspawn = 0.0_p
                trace = 0
                estimator_numerators = 0
                if (allocated(excit_distribution)) excit_distribution = 0
                nparticles_start_report = nparticles_old(1)

                do icycle = 1, ncycles
                    spawning_head = spawning_block_start
                    iteration = (ireport-1)*ncycles + icycle

                    ! Number of spawning attempts that will be made.
                    ! Each particle and each end gets to attempt to
                    ! spawn onto a connected determinant and a chance
                    ! to die/clone.
                    nattempts = 4*nparticles(1)

                    ! Reset death counter
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

                ! Update the shift and desired thermal quantites.
                call update_dmqmc_estimators(nparticles_old)

                call cpu_time(t2)

                ! t1 was the time at the previous iteration, t2 the current time.
                ! t2-t1 is thus the time taken by this report loop.
                if (parent) call write_fciqmc_report(ireport, nparticles_start_report, t2-t1)

                ! cpu_time outputs an elapsed time, so update the reference timer.
                t1 = t2

                call fciqmc_interact(soft_exit)
                if (soft_exit) exit

            end do
            
            if (soft_exit) then
                exit
            else if (doing_reduced_dm) then
                call call_rdm_procedures()
            end if
        end do

        if (parent) then
            call write_fciqmc_final(ireport)
            write (6,'()')
        end if

        call load_balancing_report()

        if (dump_restart_file) call dump_restart(mc_cycles_done+ncycles*nreport, nparticles_old(1))

        call dealloc_det_info(cdet1, .false.)
        call dealloc_det_info(cdet2, .false.)

    end subroutine do_dmqmc

end module dmqmc
