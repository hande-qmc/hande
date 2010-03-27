module fciqmc

! Module for performing optimised (hopefully!) full configuration interaction
! quantum monte carlo (FCIQMC) calculations.

use fciqmc_data
implicit none

contains

    subroutine init_fciqmc()

        ! Initialisation for fciqmc calculations.
        ! Setup the spin polarisation for the system, initialise the RNG,
        ! allocate the required memory for the list of walkers and set the
        ! initial walker.

        use errors, only: stop_all
        use parallel, only: nprocs, parent
        use utils, only: int_fmt

        use basis, only: basis_length
        use calc, only: sym_in, ms_in
        use determinants, only: encode_det, set_spin_polarisation, write_det
        use hamiltonian, only: get_hmatel_real, slater_condon0_hub_real, slater_condon0_hub_k
        use fciqmc_restart, only: read_restart
        use system, only: nel, nalpha, nbeta, system_type, hub_real, hub_k

        integer :: ierr
        integer :: i

        if (nprocs > 1) call stop_all('init_fciqmc','Not (yet!) a parallel algorithm.')

        if (parent) write (6,'(1X,a6,/,1X,6("-"),/)') 'FCIQMC'

        ! Allocate main walker lists.
        allocate(walker_dets(basis_length,walker_length), stat=ierr)
        allocate(walker_population(walker_length), stat=ierr)
        allocate(walker_energies(walker_length), stat=ierr)

        ! Allocate spawned walker lists.
        allocate(spawned_walker_dets(basis_length,spawned_walker_length), stat=ierr)
        allocate(spawned_walker_population(spawned_walker_length), stat=ierr)
        allocate(spawning_head(0:nprocs-1), stat=ierr)

        ! Set spin variables.
        call set_spin_polarisation(ms_in)

        ! Set initial walker population.
        ! occ_list could be set and allocated in the input.
        allocate(f0(basis_length), stat=ierr)
        if (restart) then
            if (.not.allocated(occ_list0)) allocate(occ_list0(nel), stat=ierr)
            call read_restart()
        else
            tot_walkers = 1
            walker_population(tot_walkers) = D0_population

            ! Set the reference determinant to be the spin-orbitals with the lowest
            ! kinetic energy which satisfy the spin polarisation.
            ! Note: this is for testing only!  The symmetry input is currently
            ! ignored.
            if (.not.allocated(occ_list0)) then
                allocate(occ_list0(nel), stat=ierr)
                forall (i=1:nalpha) occ_list0(i) = 2*i-1
                forall (i=1:nbeta) occ_list0(i+nalpha) = 2*i
            end if

            walker_dets(:,tot_walkers) = encode_det(occ_list0)

            walker_energies(tot_walkers) = 0.0_p

            ! Reference det
            f0 = walker_dets(:,tot_walkers)
            ! Energy of reference determinant.
            select case(system_type)
            case(hub_k)
                H00 = slater_condon0_hub_k(f0)
            case(hub_real)
                H00 = slater_condon0_hub_real(f0)
            end select
        end if

        if (parent) then
            write (6,'(1X,a29,1X)',advance='no') 'Reference determinant, |D0> ='
            call write_det(walker_dets(:,tot_walkers), new_line=.true.)
            write (6,'(1X,a16,f20.12)') 'E0 = <D0|H|D0> =',H00
            write (6,'(1X,a44,'//int_fmt(walker_population(tot_walkers),1)//')') &
                              'Initial population on reference determinant:',walker_population(tot_walkers)
            write (6,'(/,1X,a68,/)') 'Note that FCIQMC calculates the correlation energy relative to |D0>.'
        end if
        
    end subroutine init_fciqmc

    subroutine fciqmc_main()

        ! Wrapper around do_fciqmc to set the appropriate procedures that are to
        ! be called for the current fciqmc calculation.
        ! This is a bit hacky, but avoids lots of branching due to if blocks
        ! within the fciqmc algorithm.

        use system, only: system_type, hub_k, hub_real
        use hamiltonian, only: slater_condon0_hub_k, slater_condon0_hub_real
        use determinants, only: decode_det_spinocc_spinunocc
        use energy_evaluation, only: update_proj_energy_hub_k, update_proj_energy_hub_real
        use spawning, only: spawn_hub_k, spawn_hub_real

        select case(system_type)
        case(hub_k)
            call do_fciqmc(decode_det_spinocc_spinunocc, update_proj_energy_hub_k, spawn_hub_k, slater_condon0_hub_k)
        case(hub_real)
            call do_fciqmc(decode_det_spinocc_spinunocc, update_proj_energy_hub_real, spawn_hub_real, slater_condon0_hub_real)
        end select

    end subroutine fciqmc_main

    subroutine do_fciqmc(decoder, update_proj_energy, spawner, sc0)

        ! Run the FCIQMC algorithm starting from the initial walker
        ! distribution.

        use parallel, only: parent
  
        use annihilation, only: direct_annihilation
        use basis, only: basis_length
        use death, only: stochastic_death
        use determinants, only: det_info
        use fciqmc_restart, only: dump_restart
        use system, only: nel, nalpha, nbeta, nvirt_alpha, nvirt_beta
        use energy_evaluation, only: update_shift

        ! It seems this interface block cannot go in a module when we're passing
        ! subroutines around as arguments.  Bummer.
        ! If only procedure pointers were more commonly implemented...
        interface
            subroutine decoder(f,d)
                use basis, only: basis_length
                use const, only: i0
                use determinants, only: det_info
                implicit none
                integer(i0), intent(in) :: f(basis_length)
                type(det_info), intent(inout) :: d
            end subroutine decoder
            subroutine update_proj_energy(idet, inst_proj_energy)
                use const, only: p
                implicit none
                integer, intent(in) :: idet
                real(p), intent(inout) :: inst_proj_energy
            end subroutine update_proj_energy
            subroutine spawner(d, parent_sign)
                use determinants, only: det_info
                implicit none
                type(det_info), intent(in) :: d
                integer :: parent_sign
            end subroutine spawner
            function sc0(f) result(hmatel)
                use basis, only: basis_length
                use const, only: i0, p
                implicit none
                real(p) :: hmatel
                integer(i0), intent(in) :: f(basis_length)
            end function sc0
        end interface

        integer :: ierr
        integer :: idet, ireport, icycle, iparticle, nparticles, nparticles_old
        type(det_info) :: cdet
        real(p) :: inst_proj_energy

! DEBUG CHECK ONLY.
!        integer :: sum1, sum2

        ! Allocate det_info components.
        allocate(cdet%f(basis_length), stat=ierr)
        allocate(cdet%occ_list(nel), stat=ierr)
        allocate(cdet%occ_list_alpha(nalpha), stat=ierr)
        allocate(cdet%occ_list_beta(nbeta), stat=ierr)
        allocate(cdet%unocc_list_alpha(nvirt_alpha), stat=ierr)
        allocate(cdet%unocc_list_beta(nvirt_beta), stat=ierr)

        ! from restart
        nparticles_old = nparticles_old_restart

        ! Main fciqmc loop.

        if (parent) call write_fciqmc_report_header()

        do ireport = 1, nreport

            ! Zero averaged projected energy.
            proj_energy = 0.0_p

            do icycle = 1, ncycles

                ! Zero instantaneous projected energy.
                inst_proj_energy = 0.0_p

                ! Reset the current position in the spawning array to be the
                ! slot preceding the first slot.
                spawning_head = 0

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)

                    call decoder(cdet%f, cdet)

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    call update_proj_energy(idet, inst_proj_energy)

                    do iparticle = 1, abs(walker_population(idet))
                        
                        ! spawn
                        call spawner(cdet, walker_population(idet))

                    end do

                    ! Clone or die.
                    call stochastic_death(idet)

                end do

! DEBUG CHECK ONLY.
!                sum1 = sum(walker_population(:tot_walkers)) + sum(spawned_walker_population(:spawning_head))
                call direct_annihilation(sc0)
! DEBUG CHECK ONLY.
!                sum2 = sum(walker_population(:tot_walkers))
!                if (sum1 /= sum2) then
!                    write (6,*) 'huh?!', sum1, sum2
!                    stop
!                end if

                ! normalise projected energy and add to running total.
                proj_energy = proj_energy + inst_proj_energy/D0_population

            end do

            ! Update the shift
            nparticles = sum(abs(walker_population(:tot_walkers))) ! This can be done more efficiently by counting as we go...
            if (vary_shift) then
                call update_shift(nparticles_old, nparticles, ncycles)
            end if
            nparticles_old = nparticles
            if (nparticles > target_particles .and. .not.vary_shift) then
                vary_shift = .true.
                start_vary_shift = ireport
            end if

            ! Running average projected energy 
            av_proj_energy = av_proj_energy + proj_energy
            ! average projected energy over report loop.
            proj_energy = proj_energy/ncycles

            if (parent) call write_fciqmc_report(ireport, nparticles)

        end do

        if (parent) then
            call write_fciqmc_final()
            write (6,'()')
        end if

        if (dump_restart_file) call dump_restart(mc_cycles_done+ncycles*nreport, nparticles_old)

    end subroutine do_fciqmc

end module fciqmc
