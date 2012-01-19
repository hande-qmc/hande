module dmqmc

! Main loop for performing DMQMC calculations, similar to the
! file fciqmc.f90 which carries out the main loop in FCIQMC
! calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine dmqmc_main()

        ! Wrapper around dmqmc calculation procedures to set the appropriate procedures
        ! that are to be called for the current dmqmc calculation.

        use system, only: system_type, hub_k, hub_real, heisenberg, hub_k_coulomb, hubt
        use hamiltonian, only: slater_condon0_hub_k, slater_condon0_hub_real
        use hamiltonian, only: diagonal_element_heisenberg
        use heisenberg_estimators, only: update_proj_energy_heisenberg_basic
        use dmqmc_estimators, only: dmqmc_energy_heisenberg, dmqmc_stag_mag_heisenberg
        use dmqmc_estimators, only: dmqmc_energy_hub_real, dmqmc_energy_squared_heisenberg
        use dmqmc_estimators, only: dmqmc_correlation_function_heisenberg
        use dmqmc_procedures, only: random_distribution_heisenberg, random_distribution_entire_space
        use determinants, only: decode_det_spinocc_spinunocc, decode_det_occ
        use energy_evaluation, only: update_proj_energy_hub_k, update_proj_hfs_hub_k
        use spawning, only: spawn_hub_k, spawn_hub_real
        use spawning, only: create_spawned_particle_density_matrix, create_spawned_particle_truncated_density_matrix
        use spawning, only: spawn_heisenberg
        use calc, only: dmqmc_calc, doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_energy_squared, dmqmc_correlation

        use excitations, only: enumerate_all_excitations_hub_k, enumerate_all_excitations_hub_real

        real(dp) :: hub_matel

        ! set function pointers
        select case(system_type)
        case (hub_k)
            decoder_ptr => decode_det_spinocc_spinunocc
            update_proj_energy_ptr => update_proj_energy_hub_k
            spawner_ptr => spawn_hub_k
            sc0_ptr => slater_condon0_hub_k
            hub_matel = hub_k_coulomb
        case (hub_real)
            decoder_ptr => decode_det_occ
            if (doing_dmqmc_calc(dmqmc_energy)) update_dmqmc_energy_ptr => dmqmc_energy_hub_real
            spawner_ptr => spawn_hub_real
            sc0_ptr => slater_condon0_hub_real
            hub_matel = hubt
            dmqmc_initial_distribution_ptr => random_distribution_entire_space
        case (heisenberg)
            ! Only need occupied orbitals list, as for the real Hubbard case.
            decoder_ptr => decode_det_occ
            if (doing_dmqmc_calc(dmqmc_energy)) update_dmqmc_energy_ptr => dmqmc_energy_heisenberg
            if (doing_dmqmc_calc(dmqmc_energy_squared)) update_dmqmc_energy_squared_ptr => &
                                                           dmqmc_energy_squared_heisenberg
            if (doing_dmqmc_calc(dmqmc_correlation)) update_dmqmc_correlation_ptr => &
                                                           dmqmc_correlation_function_heisenberg
            if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) &
                         update_dmqmc_stag_mag_ptr => dmqmc_stag_mag_heisenberg
            spawner_ptr => spawn_heisenberg
            sc0_ptr => diagonal_element_heisenberg
            dmqmc_initial_distribution_ptr => random_distribution_heisenberg
        end select

        if (truncate_space) then
            create_spawned_particle_dm_ptr => create_spawned_particle_truncated_density_matrix
        else
            create_spawned_particle_dm_ptr => create_spawned_particle_density_matrix
        end if

        call do_dmqmc()

    end subroutine dmqmc_main

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
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit, get_excitation_level
        use fciqmc_common
        use fciqmc_restart, only: dump_restart
        use interact, only: fciqmc_interact
        use system, only: nel
        use calc, only: seed, doing_dmqmc_calc, dmqmc_energy
        use calc, only: dmqmc_staggered_magnetisation, dmqmc_energy_squared
        use dSFMT_interface, only: dSFMT_init
        use utils, only: int_fmt

        integer :: i, idet, ireport, icycle, iparticle
        integer :: beta_cycle
        integer(lint) :: nparticles_old(sampling_size)
        integer(lint) :: nparticles_start_report
        type(det_info) :: cdet1, cdet2
        integer :: nspawned, nattempts, ndeath
        type(excit) :: connection
        integer :: bit_pos, bit_element
        integer :: walker_pop_1, walker_pop_2, spawning_end
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

            nparticles_old = D0_population

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
                call update_dmqmc_estimators(ireport, nparticles_old)

                call cpu_time(t2)

                ! t1 was the time at the previous iteration, t2 the current time.
                ! t2-t1 is thus the time taken by this report loop.
                if (parent) call write_fciqmc_report(ireport, nparticles_start_report, t2-t1)

                ! cpu_time outputs an elapsed time, so update the reference timer.
                t1 = t2

                call fciqmc_interact(ireport, soft_exit)
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
