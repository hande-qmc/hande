module ccmc

! Module  for performing coupled cluster Monte Carlo (CCMC) calculations).

! Due to the similarities with FCIQMC, we can re-use lots of the same routines
! (especially the spawning, death and annihilation).  As a result, the structure
! of do_ccmc is remarkably similar to the other do_*mc routines.

implicit none

contains

    subroutine do_ccmc()

        ! Run the CCMC algorithm starting from the initial walker distribution
        ! using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        use parallel

        use annihilation, only: direct_annihilation
        use basis, only: basis_length, nbasis
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        use excitations, only: excit
        use energy_evaluation, only: update_energy_estimators
        use fciqmc_data
        use fciqmc_common
        use fciqmc_restart, only: dump_restart
        use interact, only: fciqmc_interact
        use proc_pointers

        integer :: ireport, icycle
!        integer :: iparticle ! not used yet
        integer(lint) :: iattempt, nattempts, nparticles_old(sampling_size)
        type(det_info) :: cdet

        integer :: ndeath
!        integer :: nspawned ! not used yet
!        type(excit) :: connection  ! not used yet

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

                ! Allow one spawning & death attempt for each walker on the
                ! processor.
                do iattempt = 1, nparticles(1)
                    ! TODO: select cluster size
                    ! TODO: find cluster
                    ! TODO: projected estimator.
                    ! TODO: evolve.
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

            ! cpu_time outputs an elapsed time, so update the reference timer.
            t1 = t2

            call fciqmc_interact(soft_exit)
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

    end subroutine do_ccmc

end module ccmc
