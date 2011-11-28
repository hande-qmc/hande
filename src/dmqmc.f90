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
        use dmqmc_estimators, only: call_dmqmc_estimators
        use determinants, only: decode_det_spinocc_spinunocc, decode_det_occ
        use energy_evaluation, only: update_proj_energy_hub_k, update_proj_hfs_hub_k
        use energy_evaluation, only: update_proj_energy_hub_real
        use spawning, only: spawn_hub_k, spawn_hub_real, create_spawned_particle_density_matrix
        use spawning, only: spawn_heisenberg
        use calc, only: dmqmc_calc

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
            update_proj_energy_ptr => update_proj_energy_hub_real
            spawner_ptr => spawn_hub_real
            sc0_ptr => slater_condon0_hub_real
            hub_matel = hubt
        case (heisenberg)
            ! Only need occupied orbitals list, as for the real Hubbard case.
            decoder_ptr => decode_det_occ
            update_dmqmc_energy_ptr => call_dmqmc_estimators
            spawner_ptr => spawn_heisenberg
            sc0_ptr => diagonal_element_heisenberg
        end select

        create_spawned_particle_dm_ptr => create_spawned_particle_density_matrix

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
        use dmqmc_procedures, only: random_distribution_heisenberg, initial_dmqmc_status
        use dmqmc_estimators, only: update_dmqmc_estimators
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit
        use fciqmc_common
        use fciqmc_restart, only: dump_restart
        use interact, only: fciqmc_interact
        use system, only: nel
        use calc, only: seed
        use dSFMT_interface, only: dSFMT_init
        use utils, only: int_fmt
        
        integer :: i, idet, ireport, icycle, iparticle
        integer :: beta_cycle, nparticles_old(sampling_size)
        type(det_info) :: cdet1, cdet2
        integer :: nspawned, nattempts, ndeath
        type(excit) :: connection
        integer :: bit_pos, bit_element
        integer :: walker_pop_1, walker_pop_2, spawning_end
        integer :: beta_index
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
 
        do beta_cycle = 1, beta_loops
        ! Reset the current position in the spawning array to be the
        ! slot preceding the first slot.
        spawning_head = spawning_block_start
        tot_walkers = 0
        nparticles = 0
        shift = 0
        av_shift =0
        start_vary_shift = 0
        vary_shift = .false.
 
        ! Need to place psips randomly along the diagonal at the
        ! start of every iteration. Pick orbitals randomly, each
        ! with equal probability, so that when electrons are placed
        ! one these orbitals they will have the correct spin and symmetry.
        call random_distribution_heisenberg()

        call direct_annihilation()

        if (beta_cycle .ne. 1 .and. parent) then
           write (6,'(a32,i7)') &
                   " # Resetting beta... Beta loop =", beta_cycle
        ! Reset the random number generator with seed = seed + 1
           seed = seed + 1
           call dSFMT_init(seed + iproc)   
           write (6,'(a52,'//int_fmt(seed,1)//',a1)') ' # Resetting random number generator with a seed of:', seed, '.'
        end if
        call initial_dmqmc_status()

        nparticles_old = dmqmc_npsips

        do ireport = 1, nreport

            ! Zero report cycle quantities.
            rspawn = 0.0_p
            beta_index = 0
            trace = 0
            total_trace = 0
            thermal_energy = 0
            total_thermal_energy = 0

            do icycle = 1, ncycles
                beta_index = beta_index + 1
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
                    cdet1%idet = idet
                    cdet2%idet = idet

                    ! Decode and store the the relevant information for
                    ! both bitstrings. Both of these bitstrings are required
                    ! to refer to the correct element in the density matrix.
                    call decoder_ptr(cdet1%f, cdet1)
                    call decoder_ptr(cdet2%f, cdet2)

                    ! Add the contributions from this particular density matrix
                    ! element to the quantites required to calculate the desired 
                    ! properties, such as the trace of the density matrix.
                    call update_dmqmc_energy_ptr(idet, beta_index)

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
                    ! current walker. We do both of these at once by using
                    ! (walker_energies(1,idet)+walker_energies(2,idet))/2 as the correct
                    ! energy. 
                    call stochastic_death((walker_energies(1,idet)+walker_energies(2,idet))/2, &
                             walker_population(1,idet), nparticles(1), ndeath)
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
            call update_dmqmc_estimators(ireport, nparticles_old)

            call cpu_time(t2)

            ! t1 was the time at the previous iteration, t2 the current time.
            ! t2-t1 is thus the time taken by this report loop.
            if (parent) call write_fciqmc_report(ireport, nparticles_old(1), t2-t1)

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
