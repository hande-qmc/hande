module fciqmc

! Module for performing optimised (hopefully!) full configuration interaction
! quantum monte carlo (FCIQMC) calculations.

use fciqmc_data
implicit none

contains

    subroutine fciqmc_main()

        ! Wrapper around do_fciqmc and do_ifciqmc to set the appropriate procedures
        ! that are to be called for the current fciqmc calculation.
        ! This is a bit hacky, but avoids lots of branching due to if blocks
        ! within the fciqmc algorithm.

        use system, only: system_type, hub_k, hub_real, hub_k_coulomb, hubt
        use hellmann_feynman_sampling

        use hamiltonian, only: slater_condon0_hub_k, slater_condon0_hub_real
        use determinants, only: decode_det_spinocc_spinunocc, decode_det_occ
        use energy_evaluation, only: update_proj_energy_hub_k, update_proj_hfs_hub_k, update_proj_energy_hub_real
        use spawning, only: spawn_hub_k, spawn_hub_real

        use calc, only: initiator_fciqmc, hfs_fciqmc_calc, ct_fciqmc_calc, doing_calc

        use ct_fciqmc, only: do_ct_fciqmc
        use excitations, only: enumerate_all_excitations_hub_k, enumerate_all_excitations_hub_real

        if (doing_calc(initiator_fciqmc)) then
            select case(system_type)
            case(hub_k)
                call do_ifciqmc(decode_det_spinocc_spinunocc, update_proj_energy_hub_k, spawn_hub_k, slater_condon0_hub_k)
            case(hub_real)
                call do_ifciqmc(decode_det_occ, update_proj_energy_hub_real, spawn_hub_real, slater_condon0_hub_real)
            end select
        else if (doing_calc(ct_fciqmc_calc)) then
            select case(system_type)
            case(hub_k)
                call do_ct_fciqmc(decode_det_spinocc_spinunocc, update_proj_energy_hub_k, enumerate_all_excitations_hub_k,&
                                  slater_condon0_hub_k, hub_k_coulomb)
            case(hub_real)
                call do_ct_fciqmc(decode_det_occ, update_proj_energy_hub_real, enumerate_all_excitations_hub_real,&
                                  slater_condon0_hub_real, hubt)
            end select
        else
            select case(system_type)
            case(hub_k)
                if (doing_calc(hfs_fciqmc_calc)) then
                    call init_hellmann_feynman_sampling()
                    call do_hfs_fciqmc(decode_det_spinocc_spinunocc, update_proj_hfs_hub_k, spawn_hub_k, slater_condon0_hub_k)
                else
                    call do_fciqmc(decode_det_spinocc_spinunocc, update_proj_energy_hub_k, spawn_hub_k, slater_condon0_hub_k)
                end if
            case(hub_real)
                call do_fciqmc(decode_det_occ, update_proj_energy_hub_real, spawn_hub_real, slater_condon0_hub_real)
            end select
        end if

    end subroutine fciqmc_main

    subroutine do_fciqmc(decoder, update_proj_energy, spawner, sc0)

        ! Run the FCIQMC algorithm starting from the initial walker
        ! distribution.

        ! This is implemented by abusing fortran's ability to pass procedures as
        ! arguments (if only function pointers (F2003) were implemented in more
        ! compilers!).  This allows us to avoid many system dependent if blocks,
        ! which are constant for a given calculation.  Avoiding such branching
        ! is worth the extra verbosity (especially if procedures are written to
        ! be sufficiently modular that implementing a new system can reuse many
        ! existing routines) as it leads to much faster code.

        ! (This idea is now being "borrowed" for use in neci.  Bah...)

        ! In:
        !    decoder: relevant subroutine to decode/extract the necessary
        !        information from the determinant bit string.  See the
        !        determinants module.
        !    update_proj_energy: relevant subroutine to update the projected
        !        energy.  See the energy_evaluation module.
        !    spawner: relevant subroutine to attempt to spawn a walker from an
        !        existing walker.  See the spawning module.
        !    sc0: relevant function to evaluate the diagonal Hamiltonian matrix
        !    elements, <D|H|D>.  See the hamiltonian module.

        use parallel
  
        use annihilation, only: direct_annihilation
        use basis, only: basis_length
        use death, only: stochastic_death
        use determinants, only:det_info, alloc_det_info 
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit
        use interact, only: fciqmc_interact
        use fciqmc_restart, only: dump_restart
        use spawning, only: create_spawned_particle
        use fciqmc_common

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
            subroutine update_proj_energy(idet)
                use const, only: p
                implicit none
                integer, intent(in) :: idet
            end subroutine update_proj_energy
            subroutine spawner(d, parent_sign, nspawned, connection)
                use determinants, only: det_info
                use excitations, only: excit
                implicit none
                type(det_info), intent(in) :: d
                integer, intent(in) :: parent_sign
                integer, intent(out) :: nspawned
                type(excit), intent(out) :: connection
            end subroutine spawner
            function sc0(f) result(hmatel)
                use basis, only: basis_length
                use const, only: i0, p
                implicit none
                real(p) :: hmatel
                integer(i0), intent(in) :: f(basis_length)
            end function sc0
        end interface

        integer :: idet, ireport, icycle, iparticle, nparticles_old(sampling_size)
        type(det_info) :: cdet

        integer :: nspawned, nattempts
        type(excit) :: connection


        logical :: soft_exit

        real :: t1, t2

        ! Allocate det_info components.
        call alloc_det_info(cdet)

        ! from restart
        nparticles_old = nparticles_old_restart

        ! Main fciqmc loop.

        if (parent) call write_fciqmc_report_header()
        call initial_fciqmc_status(update_proj_energy)

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
                nattempts = nparticles(1)

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)

                    call decoder(cdet%f, cdet)

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    call update_proj_energy(idet)

                    do iparticle = 1, abs(walker_population(1,idet))
                        
                        ! Attempt to spawn.
                        call spawner(cdet, walker_population(1,idet), nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) call create_spawned_particle(cdet, connection, nspawned, spawned_pop)

                    end do

                    ! Clone or die.
                    call stochastic_death(walker_energies(1,idet), walker_population(1,idet), nparticles(1))

                end do

                ! Add the spawning rate (for the processor) to the running
                ! total.
                rspawn = rspawn + spawning_rate(nattempts)

                ! D0_population is communicated in the direct_annihilation
                ! algorithm for efficiency.
                call direct_annihilation(sc0)

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

!            call dump_restart(ireport*ncycles, nparticles_old(1))

        end do

        if (parent) then
            call write_fciqmc_final(ireport)
            write (6,'()')
        end if

        call load_balancing_report()

        if (dump_restart_file) call dump_restart(mc_cycles_done+ncycles*nreport, nparticles_old(1))

    end subroutine do_fciqmc

    subroutine do_ifciqmc(decoder, update_proj_energy, spawner, sc0)

        ! Run the initiator-FCIQMC algorithm starting from the initial walker
        ! distribution.

        ! See notes about the implementation of this using function pointers
        ! (F77 style rather than F2003, sadly!) in do_fciqmc.

        ! In:
        !    decoder: relevant subroutine to decode/extract the necessary
        !        information from the determinant bit string.  See the
        !        determinants module.
        !    update_proj_energy: relevant subroutine to update the projected
        !        energy.  See the energy_evaluation module.
        !    spawner: relevant subroutine to attempt to spawn a walker from an
        !        existing walker.  See the spawning module.
        !    sc0: relevant function to evaluate the diagonal Hamiltonian matrix
        !    elements, <D|H|D>.  See the hamiltonian module.

        use parallel
  
        use annihilation, only: direct_annihilation_initiator
        use basis, only: basis_length, bit_lookup, nbasis
        use death, only: stochastic_death
        use determinants, only: det_info, alloc_det_info
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit
        use interact, only: fciqmc_interact
        use fciqmc_restart, only: dump_restart
        use system, only: nel
        use spawning, only: create_spawned_particle_initiator
        use fciqmc_common

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
            subroutine update_proj_energy(idet)
                use const, only: p
                implicit none
                integer, intent(in) :: idet
            end subroutine update_proj_energy
            subroutine spawner(d, parent_sign, nspawned, connection)
                use determinants, only: det_info
                use excitations, only: excit
                implicit none
                type(det_info), intent(in) :: d
                integer, intent(in) :: parent_sign
                integer, intent(out) :: nspawned
                type(excit), intent(out) :: connection
            end subroutine spawner
            function sc0(f) result(hmatel)
                use basis, only: basis_length
                use const, only: i0, p
                implicit none
                real(p) :: hmatel
                integer(i0), intent(in) :: f(basis_length)
            end function sc0
        end interface

        integer :: i, idet, ireport, icycle, iparticle, nparticles_old(sampling_size)
        type(det_info) :: cdet

        integer :: nspawned, nattempts
        type(excit) :: connection

        integer :: parent_flag
        integer(i0) :: cas_mask(basis_length), cas_core(basis_length)
        integer :: bit_pos, bit_element

        logical :: soft_exit

        real :: t1, t2

        ! Allocate det_info components.
        call alloc_det_info(cdet)

        ! The complete active space (CAS) is given as (N_cas,N_active), where
        ! N_cas is the number of electrons in the N_active orbitals.
        ! The N-N_cas electrons occupy the lowest energy orbitals ("core"
        ! orbitals) for all determinants within the CAS.
        ! The 2M-N_core-N_active highest energy orbitals are inactive and are
        ! not occupied in any determinants within the CAS.
        ! Create a mask which has bits set for all core electrons and a mask
        ! which has bits set for all inactive orbitals.
        cas_mask = 0
        cas_core = 0
        ! Set core obitals.
        do i = 1, nel - CAS(1)
            bit_pos = bit_lookup(1,i)
            bit_element = bit_lookup(2,i)
            cas_mask = ibset(cas_mask(bit_element), bit_pos)
            cas_core = ibset(cas_core(bit_element), bit_pos)
        end do
        ! Set inactive obitals.
        do i = nel - CAS(1) + 2*CAS(2) + 1, nbasis
            bit_pos = bit_lookup(1,i)
            bit_element = bit_lookup(2,i)
            cas_mask = ibset(cas_mask(bit_element), bit_pos)
        end do
        ! Thus ANDing a determinant with cas_mask gives the electrons in the
        ! core or inactive orbitals.  The determinant is only in the CAS if the
        ! result is identical to the cas_core mask (i.e. all the core orbitals
        ! are filled and no electrons are in the inactive orbitals).

        ! from restart
        nparticles_old = nparticles_old_restart

        ! Main fciqmc loop.

        if (parent) call write_fciqmc_report_header()
        call initial_fciqmc_status(update_proj_energy)

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
                nattempts = nparticles(1)

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f = walker_dets(:,idet)

                    call decoder(cdet%f, cdet)

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the i-FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    call update_proj_energy(idet)

                    ! Is this determinant an initiator?
                    if (abs(walker_population(1,idet)) > initiator_population) then
                        ! Has a high enough population to be an initiator.
                        parent_flag = 0
                    else if (all(iand(cdet%f,cas_mask) == cas_core)) then
                        ! Is in the complete active space.
                        parent_flag = 0
                    else
                        ! Isn't an initiator.
                        parent_flag = 1
                    end if

                    do iparticle = 1, abs(walker_population(1,idet))
                        
                        ! Attempt to spawn.
                        call spawner(cdet, walker_population(1,idet), nspawned, connection)
                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) then
                            call create_spawned_particle_initiator(cdet, parent_flag, connection, nspawned, spawned_pop)
                        end if

                    end do

                    ! Clone or die.
                    call stochastic_death(walker_energies(1,idet), walker_population(1,idet), nparticles(1))

                end do

                ! Add the spawning rate (for the processor) to the running
                ! total.
                rspawn = rspawn + spawning_rate(nattempts)

                ! D0_population is communicated in the direct_annihilation
                ! algorithm for efficiency.
                call direct_annihilation_initiator(sc0)

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

!            call dump_restart(ireport*ncycles, nparticles_old(1))

        end do

        if (parent) then
            call write_fciqmc_final(ireport)
            write (6,'()')
        end if

        call load_balancing_report()

        if (dump_restart_file) call dump_restart(mc_cycles_done+ncycles*nreport, nparticles_old(1))

    end subroutine do_ifciqmc

end module fciqmc
