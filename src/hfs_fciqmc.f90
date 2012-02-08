module hellmann_feynman_sampling

! Module for performing Hellmann--Feynman sampling in FCIQMC.

! See documentation/theory/hellmann_feynman/hf.tex for details.

use const

use hfs_data

use proc_pointers

implicit none

contains

    subroutine do_hfs_fciqmc()

        ! Run the FCIQMC algorithm starting from the initial walker
        ! distribution and perform Hellmann--Feynman sampling in conjunction on
        ! the instantaneous wavefunction.

        ! This is implemented using F2003 function pointers.  This allows us to
        ! avoid many system dependent if blocks, which are constant for a given
        ! calculation.  Avoiding such branching is worth the extra verbosity
        ! (especially if procedures are written to be sufficiently modular that
        ! implementing a new system can reuse many existing routines) as it
        ! leads to much faster code.

        use parallel

        use annihilation, only: direct_annihilation
        use basis, only: basis_length
        use death, only: stochastic_hf_cloning
        use determinants, only:det_info, alloc_det_info
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit
        use interact, only: fciqmc_interact
        use fciqmc_restart, only: dump_restart, write_restart_file_every_nreports
        use spawning, only: create_spawned_particle
        use fciqmc_common

        integer :: idet, ireport, icycle, iparticle
        integer(lint) :: nattempts, nparticles_old(sampling_size)
        type(det_info) :: cdet

        integer :: nspawned, ndeath
        type(excit) :: connection

        logical :: soft_exit

        real :: t1, t2

        ! Allocate det_info components.
        call alloc_det_info(cdet)

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
            proj_hf_O_hpsip = 0.0_p
            proj_hf_H_hfpsip = 0.0_p
            D0_population = 0.0_p
            D0_hf_population = 0.0_p
            rspawn = 0.0_p

            do icycle = 1, ncycles

                ! Reset the current position in the spawning array to be the
                ! slot preceding the first slot.
                spawning_head = spawning_block_start

                ! Number of spawning attempts that will be made.
                ! Each Hamiltonian particle gets a chance to spawn a Hamiltonian
                ! particle, clone/die, spawn a Hellmann-Feynman particle and clone
                ! itself into a Hellmann-Feynman particle.  Each H-F particle
                ! gets a chance to spawn and a chance to clone/die.
                ! This is used for accounting later, not for controlling the spawning.
                nattempts = 4*nparticles(1) + 2*nparticles(2)

                ! Reset death counter.
                ndeath = 0

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)

                    call decoder_ptr(cdet%f, cdet)

                    ! It is much easier to evaluate projected values at the
                    ! start of the FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    call update_proj_energy_ptr(idet)

                    do iparticle = 1, abs(walker_population(1,idet))

                        ! Attempt to spawn Hamiltonian walkers..
                        call spawner_ptr(cdet, walker_population(1,idet), nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) call create_spawned_particle(cdet, connection, nspawned, spawned_pop)

                        ! Attempt to spawn Hellmann--Feynman walkers from
                        ! Hamiltonian walkers.
                        call spawner_hfs_ptr(cdet, walker_population(1,idet), nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) call create_spawned_particle(cdet, connection, nspawned, spawned_hf_pop)

                    end do

                    do iparticle = 1, abs(walker_population(2,idet))

                        ! Attempt to spawn Hellmann--Feynman walkers from
                        ! Hellmann--Feynman walkers.
                        call spawner_ptr(cdet, walker_population(2,idet), nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) call create_spawned_particle(cdet, connection, nspawned, spawned_hf_pop)

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
                    call death_ptr(walker_data(1,idet), walker_population(2,idet), nparticles(2), ndeath)

                    ! Clone Hellmann--Feynman walkers from Hamiltonian walkers.
                    call stochastic_hf_cloning(walker_data(2,idet), walker_population(1,idet), &
                                               walker_population(2,idet), nparticles(2))

                    ! Clone or die: Hamiltonian walkers.
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
            call update_energy_estimators(nparticles_old)

            call cpu_time(t2)

            ! t1 was the time at the previous iteration, t2 the current time.
            ! t2-t1 is thus the time taken by this report loop.
            if (parent) call write_fciqmc_report(ireport, nparticles_old(1), t2-t1)

            if (parent) then
                ! DEBUG output only.
                ! TODO: include with std. output.
                write (17,'(i6,6f18.8,2i8)') ireport, hf_shift, proj_energy, proj_hf_O_hpsip, proj_hf_H_hfpsip, &
                             D0_population, D0_hf_population, nparticles_old
                call flush(17)
            end if

            ! Write restart file if required.
            if (mod(ireport,write_restart_file_every_nreports) == 0) &
                call dump_restart(mc_cycles_done+ncycles*ireport, nparticles_old(1))

            ! cpu_time outputs an elapsed time, so update the reference timer.
            t1 = t2

            call fciqmc_interact(soft_exit)
            if (soft_exit) exit

        end do

        if (parent) then
            call write_fciqmc_final(ireport)
            write (6,'()')
        end if

        call load_balancing_report()

        if (dump_restart_file) call dump_restart(mc_cycles_done+ncycles*nreport, nparticles_old(1))

    end subroutine do_hfs_fciqmc

end module hellmann_feynman_sampling
