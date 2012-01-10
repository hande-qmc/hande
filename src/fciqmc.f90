module fciqmc

! Module for performing optimised (hopefully!) full configuration interaction
! quantum monte carlo (FCIQMC) calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine fciqmc_main()

        ! Wrapper around fciqmc calculation procedures to set the appropriate procedures
        ! that are to be called for the current fciqmc calculation.

        use system, only: system_type, hub_k, hub_real, heisenberg, hub_k_coulomb, hubt, staggered_magnetic_field
        use system, only: trial_function, single_basis, neel_singlet
        use system, only: guiding_function, no_guiding, neel_singlet_guiding
        use hellmann_feynman_sampling
        use hamiltonian, only: slater_condon0_hub_k, slater_condon0_hub_real
        use hamiltonian, only: diagonal_element_heisenberg, diagonal_element_heisenberg_staggered
        use heisenberg_estimators, only: update_proj_energy_heisenberg_basic, update_proj_energy_&
                                                                          &heisenberg_neel_singlet
        use heisenberg_estimators, only: update_proj_energy_heisenberg_positive
        use determinants, only: decode_det_spinocc_spinunocc, decode_det_occ
        use energy_evaluation, only: update_proj_energy_hub_k, update_proj_hfs_hub_k, update_proj_energy_hub_real
        use spawning, only: spawn_hub_k, spawn_hub_real, create_spawned_particle, create_spawned_particle_initiator, &
                            spawn_hub_k_no_renorm, spawn_hub_real_no_renorm
        use spawning, only: spawn_heisenberg, spawn_heisenberg_importance_sampling
        use death, only: stochastic_death

        use calc, only: initiator_fciqmc, hfs_fciqmc_calc, ct_fciqmc_calc, fciqmc_calc, folded_spectrum, doing_calc

        use ct_fciqmc, only: do_ct_fciqmc
        use excitations, only: enumerate_all_excitations_hub_k, enumerate_all_excitations_hub_real
        use ifciqmc, only: init_ifciqmc, set_parent_flag, set_parent_flag_dummy

        use folded_spectrum_utils
        use spawning, only: gen_excit_hub_k, gen_excit_hub_real

        real(dp) :: hub_matel

        ! set function pointers
        select case(system_type)
        case (hub_k)
            decoder_ptr => decode_det_spinocc_spinunocc
            update_proj_energy_ptr => update_proj_energy_hub_k
            if (no_renorm) then
                spawner_ptr => spawn_hub_k_no_renorm
            else
                spawner_ptr => spawn_hub_k
            end if
            sc0_ptr => slater_condon0_hub_k
            hub_matel = hub_k_coulomb
            spawner_ptr => spawn_hub_k
            death_ptr => stochastic_death
            if(doing_calc(folded_spectrum)) gen_excit_ptr => gen_excit_hub_k
        case (hub_real)
            decoder_ptr => decode_det_occ
            update_proj_energy_ptr => update_proj_energy_hub_real
            if (no_renorm) then
                spawner_ptr => spawn_hub_real_no_renorm
            else
                spawner_ptr => spawn_hub_real
            end if
            sc0_ptr => slater_condon0_hub_real
            hub_matel = hubt
            spawner_ptr => spawn_hub_real
            death_ptr => stochastic_death
            if(doing_calc(folded_spectrum)) gen_excit_ptr => gen_excit_hub_real
        case (heisenberg)
            ! Only need occupied orbitals list, as for the real Hubbard case.
            decoder_ptr => decode_det_occ
            death_ptr => stochastic_death
            ! Set which trial wavefunction to use for the energy estimator.
            select case(trial_function)
            case (single_basis)
                update_proj_energy_ptr => update_proj_energy_heisenberg_basic
            case (neel_singlet)
                update_proj_energy_ptr => update_proj_energy_heisenberg_neel_singlet
            end select
            ! Set which guiding wavefunction to use, if requested.
            select case(guiding_function)
            case (no_guiding)
                spawner_ptr => spawn_heisenberg
            case (neel_singlet_guiding)
                spawner_ptr => spawn_heisenberg_importance_sampling
            end select             
            ! Set whether the staggered magnetisation is to be calculated.
            if (abs(staggered_magnetic_field) > 0.0_p) then
                sc0_ptr => diagonal_element_heisenberg_staggered
            else
                sc0_ptr => diagonal_element_heisenberg
            end if
        end select

        if(doing_calc(folded_spectrum)) call init_folded_spectrum()

        if (doing_calc(initiator_fciqmc)) then
            call init_ifciqmc()
            set_parent_flag_ptr => set_parent_flag
            create_spawned_particle_ptr => create_spawned_particle_initiator
        else
            set_parent_flag_ptr => set_parent_flag_dummy
            create_spawned_particle_ptr => create_spawned_particle
        end if

        if (doing_calc(ct_fciqmc_calc)) then
            call do_ct_fciqmc(hub_matel)
        else
            if (doing_calc(hfs_fciqmc_calc)) then
                call init_hellmann_feynman_sampling()
                call do_hfs_fciqmc(update_proj_hfs_hub_k)
            else
                call do_fciqmc()
            end if
        end if

    end subroutine fciqmc_main

    subroutine do_fciqmc()

        ! Run the FCIQMC or initiator-FCIQMC algorithm starting from the initial walker
        ! distribution using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        use parallel
  
        use annihilation, only: direct_annihilation
        use basis, only: basis_length, bit_lookup, nbasis
        use calc, only: folded_spectrum, doing_calc
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit
        use interact, only: fciqmc_interact
        use fciqmc_restart, only: dump_restart
        use system, only: nel
        use spawning, only: create_spawned_particle_initiator
        use fciqmc_common
        use ifciqmc, only: set_parent_flag
        use folded_spectrum_utils, only: cdet_excit

        integer :: i, idet, ireport, icycle, iparticle
        integer(lint) :: nparticles_old(sampling_size)
        type(det_info) :: cdet

        integer :: nspawned, nattempts, ndeath
        type(excit) :: connection

        integer :: bit_pos, bit_element

        logical :: soft_exit

        real :: t1, t2

        ! Allocate det_info components.
        call alloc_det_info(cdet)
        if (doing_calc(folded_spectrum)) call alloc_det_info(cdet_excit)

        ! from restart
        nparticles_old = nparticles_old_restart

        ! Main fciqmc loop.
        if (parent) call write_fciqmc_report_header()
        call initial_fciqmc_status()
        ! Initialise timer.
        call cpu_time(t1)
        
        do ireport = 1, nreport

            ! Zero report cycle quantities.
            proj_energy = 0.0_p
            rspawn = 0.0_p
            D0_population = 0.0_p

            do icycle = 1, ncycles

                ! Reset the current position in the spawning array to be the
                ! slot preceding the first slot.
                spawning_head = spawning_block_start

                ! Number of spawning attempts that will be made.
                ! Each particle gets to attempt to spawn onto a connected
                ! determinant and a chance to die/clone.
                nattempts = 2*nparticles(1)

                ! Reset death counter
                ndeath = 0

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)
                    cdet%idet = idet

                    call decoder_ptr(cdet%f, cdet)

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the i-FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    call update_proj_energy_ptr(idet)

                    ! Is this determinant an initiator?
                    call set_parent_flag_ptr(walker_population(1,idet), cdet%f, cdet%initiator_flag)

                    do iparticle = 1, abs(walker_population(1,idet))
                        
                        ! Attempt to spawn.
                        call spawner_ptr(cdet, walker_population(1,idet), nspawned, connection)

                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) then
                            call create_spawned_particle_ptr(cdet, connection, nspawned, spawned_pop)
                        end if

                    end do

                    ! Clone or die.
                    call death_ptr(walker_data(1,idet), walker_population(1,idet), nparticles(1), ndeath)

                end do

                ! Add the spawning rate (for the processor) to the running
                ! total.
                rspawn = rspawn + spawning_rate(ndeath, nattempts)

                ! D0_population is communicated in the direct_annihilation
                ! algorithm for efficiency.
                call direct_annihilation()

            end do
            
            ! Update the energy estimators (shift & projected energy).
            call update_energy_estimators(ireport, nparticles_old)

            call cpu_time(t2)

            ! t1 was the time at the previous iteration, t2 the current time.
            ! t2-t1 is thus the time taken by this report loop.
            if (parent) call write_fciqmc_report(ireport, nparticles_old(1), t2-t1)

            ! cpu_time outputs an elapsed time, so update the reference timer.
            t1 = t2

            call fciqmc_interact(ireport, soft_exit)
            if (soft_exit) exit
            if (mod(ireport, select_ref_det_every_nreports) == 0) call select_ref_det()

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

        if (dump_restart_file) call dump_restart(mc_cycles_done, nparticles_old(1))

        call dealloc_det_info(cdet)
        if (doing_calc(folded_spectrum)) call dealloc_det_info(cdet_excit)

    end subroutine do_fciqmc

end module fciqmc
