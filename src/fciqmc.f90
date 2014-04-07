module fciqmc

! Module for performing optimised (hopefully!) full configuration interaction
! quantum monte carlo (FCIQMC) calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine do_fciqmc(sys)

        ! Run the FCIQMC or initiator-FCIQMC algorithm starting from the initial walker
        ! distribution using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        ! In:
        !    sys: system being studied.

        use parallel

        use annihilation, only: direct_annihilation
        use basis, only: basis_length, nbasis
        use calc, only: folded_spectrum, doing_calc, seed, initiator_approximation
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        use excitations, only: excit
        use spawning, only: create_spawned_particle_initiator
        use qmc_common
        use ifciqmc, only: set_parent_flag
        use folded_spectrum_utils, only: cdet_excit
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use utils, only: rng_init_info
        use system, only: sys_t
        use restart_hdf5, only: restart_info_global, dump_restart_hdf5

        type(sys_t), intent(in) :: sys

        integer :: idet, ireport, icycle, iparticle
        integer(lint) :: nattempts, nparticles_old(sampling_size)
        type(det_info) :: cdet
        type(dSFMT_t) :: rng
        type(fciqmc_bloom_stats_t) :: bloom_stats

        integer :: nspawned, ndeath
        type(excit) :: connection
        real(p) :: hmatel

        logical :: soft_exit

        real :: t1

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

        ! Allocate det_info components.
        call alloc_det_info(sys, cdet, .false.)
        ! Folded spectrum *needs* the bit strings to be allocated as it needs
        ! be able to manipulate the bit string to create excited states.
        if (doing_calc(folded_spectrum)) call alloc_det_info(sys, cdet_excit)

        ! from restart
        nparticles_old = tot_nparticles

        ! Main fciqmc loop.
        if (parent) call write_fciqmc_report_header()
        call initial_fciqmc_status(sys)
        ! Initialise timer.
        call cpu_time(t1)

        do ireport = 1, nreport

            ! Zero report cycle quantities.
            call init_report_loop()

            do icycle = 1, ncycles

                call init_mc_cycle(nattempts, ndeath)

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f => walker_dets(:,idet)
                    cdet%data => walker_data(:,idet)

                    call decoder_ptr(sys, cdet%f, cdet)

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the i-FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    call update_proj_energy_ptr(sys, f0, cdet, real(walker_population(1,idet),p), D0_population_cycle, &
                                                proj_energy, connection, hmatel)

                    ! Is this determinant an initiator?
                    call set_parent_flag_ptr(walker_population(1,idet), cdet%f, cdet%initiator_flag)

                    do iparticle = 1, abs(walker_population(1,idet))

                        ! Attempt to spawn.
                        call spawner_ptr(rng, sys, cdet, walker_population(1,idet), gen_excit_ptr, nspawned, connection)

                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) then
                            if(abs(nspawned) >= bloom_stats%n_bloom) then
                                if(bloom_stats%nwarnings < bloom_stats%nverbose_warnings ) then
                                    call print_bloom_warning(nspawned) 
                                end if
                                bloom_stats%nwarnings = bloom_stats%nwarnings + 1
                            end if
                            call create_spawned_particle_ptr(cdet, connection, nspawned, 1, qmc_spawn)
                        end if

                    end do

                    ! Clone or die.
                    call death_ptr(rng, walker_data(1,idet), shift(1), walker_population(1,idet), nparticles(1), ndeath)

                end do

                call direct_annihilation(sys, initiator_approximation)

                call end_mc_cycle(ndeath, nattempts)

            end do

            ! Calculate the number of warnings which occurred only on the last iteration
            bloom_stats%nwarnings_last = bloom_stats%nwarnings - bloom_stats%nwarnings_last 
            if(bloom_stats%nwarnings_last > 0 .and. tau_search .and. .not. vary_shift) then
                call decrease_tau()
            end if
            bloom_stats%nwarnings_last = bloom_stats%nwarnings

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

        if (dump_restart_file) then
            call dump_restart_hdf5(restart_info_global, mc_cycles_done, nparticles_old)
            if (parent) write (6,'()')
        end if

        call dealloc_det_info(cdet, .false.)
        if (doing_calc(folded_spectrum)) call dealloc_det_info(cdet_excit)

    end subroutine do_fciqmc

    subroutine print_bloom_warning(nspawned)

       ! Print a nice warning so the user knows a bloom has occurred
        use utils, only: int_fmt
        integer, intent(in) :: nspawned

        write(6, '(1X, "# Warning a Blooming event occured a single psips spawned"'//int_fmt(nspawned,1)//&
            '" children.")'), nspawned 

    end subroutine print_bloom_warning

end module fciqmc
