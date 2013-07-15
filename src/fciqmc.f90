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
        use calc, only: folded_spectrum, doing_calc, seed
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        use excitations, only: excit
        use fciqmc_restart, only: dump_restart, write_restart_file_every_nreports
        use system, only: nel
        use spawning, only: create_spawned_particle_initiator
        use qmc_common
        use ifciqmc, only: set_parent_flag
        use folded_spectrum_utils, only: cdet_excit
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use utils, only: rng_init_info

        integer :: idet, ireport, icycle, iparticle
        integer(lint) :: nattempts, nparticles_old(sampling_size)
        type(det_info) :: cdet
        type(dSFMT_t) :: rng

        integer :: nspawned, ndeath
        type(excit) :: connection

        logical :: soft_exit

        real :: t1

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

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
            call init_report_loop()

            do icycle = 1, ncycles

                call init_mc_cycle(nattempts, ndeath)

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)
                    cdet%data => walker_data(:,idet)

                    call decoder_ptr(cdet%f, cdet)

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the i-FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    call update_proj_energy_ptr(cdet, real(walker_population(1,idet),p))

                    ! Is this determinant an initiator?
                    call set_parent_flag_ptr(walker_population(1,idet), cdet%f, cdet%initiator_flag)

                    do iparticle = 1, abs(walker_population(1,idet))

                        ! Attempt to spawn.
                        call spawner_ptr(rng, cdet, walker_population(1,idet), nspawned, connection)

                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) then
                            call create_spawned_particle_ptr(cdet, connection, nspawned, spawned_pop)
                        end if

                    end do

                    ! Clone or die.
                    call death_ptr(rng, walker_data(1,idet), walker_population(1,idet), nparticles(1), ndeath)

                end do

                call direct_annihilation()

                call end_mc_cycle(ndeath, nattempts)

            end do

            call end_report_loop(ireport, nparticles_old, t1, soft_exit)

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

        if (dump_restart_file) call dump_restart(mc_cycles_done, nparticles_old(1))

        call dealloc_det_info(cdet)
        if (doing_calc(folded_spectrum)) call dealloc_det_info(cdet_excit)

    end subroutine do_fciqmc

end module fciqmc
