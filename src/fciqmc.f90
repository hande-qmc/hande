module fciqmc

! Module for performing optimised (hopefully!) full configuration interaction
! quantum monte carlo (FCIQMC) calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine do_fciqmc()

        ! Run the FCIQMC or initiator-FCIQMC algorithm starting from the initial walker
        ! distribution using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        use parallel

        use annihilation, only: direct_annihilation
        use basis, only: basis_length, nbasis
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

        integer :: idet, ireport, icycle, iparticle
        integer(lint) :: nparticles_old(sampling_size)
        type(det_info) :: cdet

        integer :: nspawned, nattempts, ndeath
        type(excit) :: connection

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
                ! This is used for accounting later, not for controlling the spawning.
                nattempts = 2*nparticles(1)

                ! Reset death counter
                ndeath = 0

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)
                    cdet%data => walker_data(:,idet)

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
