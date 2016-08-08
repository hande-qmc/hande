module fciqmc_data

! Data for fciqmc calculations and procedures which manipulate fciqmc and only
! fciqmc data.

use const
use spawn_data, only: spawn_t

implicit none

contains

    !--- Output procedures ---

    subroutine write_fciqmc_report_header(ntypes, comp)

        ! In:
        !    ntypes: number of particle types being sampled.
        ! In (optional):
        !    comp: Print out information of complex estimators.

        use calc, only: doing_calc, hfs_fciqmc_calc

        integer, intent(in) :: ntypes
        logical, optional, intent(in) :: comp
        logical :: comp_set

        comp_set = .false.
        if (present(comp)) comp_set = comp

        ! Data table info.
        write (6,'(1X,"Information printed out every QMC report loop:",/)')
        write (6,'(1X,"Shift: the energy offset calculated at the end of the report loop.")')
        write (6,'(1X,"H_0j: <D_0|H|D_j>, Hamiltonian matrix element.")')
        write (6,'(1X,"N_j: population of Hamiltonian particles on determinant D_j.")')
        if (doing_calc(hfs_fciqmc_calc)) then
            write (6,'(1X,"O_0j: <D_0|O|D_j>, operator matrix element.")')
            write (6,'(1X,a67)') "N'_j: population of Hellmann--Feynman particles on determinant D_j."
            write (6,'(1X,"# HF psips: current total population of Hellmann--Feynman particles.")')
        end if

        write (6,'(1X,"# H psips: current total population of Hamiltonian particles.")')
        write (6,'(1X,"# states: number of many-particle states occupied.")')
        write (6,'(1X,"# spawn_events: number of successful spawning events across all processors.")')
        write (6,'(1X,"R_spawn: average rate of spawning across all processors.")')
        write (6,'(1X,"time: average time per Monte Carlo cycle.",/)')
        write (6,'(1X,"Note that all particle populations are averaged over the report loop.",/)')

        if (comp_set) then
            write (6,'(1X,a13,(2X,a17),2(2X,a20),2(2X,a17))', advance='no') &
                     "# iterations ", "Shift            ", "Re{\sum H_0j N_j}  ", "Im{\sum H_0j N_j}  ", "Re{N_0}  ", "Im{N_0}  "
            write (6,'(4X,a9,8X)', advance='no') "# H psips"
        else
            write (6,'(1X,a13,3(2X,a17))', advance='no') &
                     "# iterations ", "Shift            ", "\sum H_0j N_j    ", "N_0              "
            if (doing_calc(hfs_fciqmc_calc)) then
                write (6,'(6(2X,a17))', advance='no') &
                    "H.F. Shift       ","\sum O_0j N_j    ","\sum H_0j N'_j   ","N'_0             ", &
                    "# H psips        ","# HF psips       "
            else
                write (6,'(4X,a9,8X)', advance='no') "# H psips"
            end if
        end if
        write (6,'(3X,"# states  # spawn_events  R_spawn    time")')

    end subroutine write_fciqmc_report_header

    subroutine write_fciqmc_report(qmc_in, qs, ireport, ntot_particles, elapsed_time, comment, non_blocking_comm, comp)

        ! Write the report line at the end of a report loop.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    qs: QMC state (containing shift and various estimators).
        !    ireport: index of the report loop.
        !    ntot_particles: total number of particles in main walker list.
        !    elapsed_time: time taken for the report loop.
        !    comment: if true, then prefix the line with a #.
        !    non_blocking_comm: true if using non-blocking communications
        ! In (optional):
        !    comp: if true, doing calculation with real and imaginary walkers
        !       so need to print extra parameters.

        use calc, only: doing_calc, hfs_fciqmc_calc
        use qmc_data, only: qmc_in_t, qmc_state_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        integer, intent(in) :: ireport
        real(dp), intent(in) :: ntot_particles(:)
        real, intent(in) :: elapsed_time
        logical, intent(in) :: comment, non_blocking_comm
        logical, intent(in), optional :: comp

        logical :: comp_set
        integer :: mc_cycles, i, j, ntypes

        ntypes = size(ntot_particles)

        ! For non-blocking communications we print out the nth report loop
        ! after the (n+1)st iteration. Adjust mc_cycles accordingly
        if (.not. non_blocking_comm) then
            mc_cycles = ireport*qmc_in%ncycles
        else
            mc_cycles = (ireport-1)*qmc_in%ncycles
        end if

        if (present(comp)) then
            comp_set = comp
        else
            comp_set = .false.
        end if

        if (comment) then
            write (6,'(1X,"#",1X)', advance='no')
        else
            write (6,'(3X)', advance='no')
        end if

        ! See also the format used in inital_fciqmc_status if this is changed.

        if (doing_calc(hfs_fciqmc_calc)) then
            write (6,'(i10,2X,6(es17.10,2X),es17.10,4X,es17.10,1X,es17.10)', advance = 'no') &
                                             qs%mc_cycles_done+mc_cycles, qs%shift(1),   &
                                             qs%estimators%proj_energy, qs%estimators%D0_population, &
                                             qs%shift(2), qs%estimators%proj_hf_O_hpsip, qs%estimators%proj_hf_H_hfpsip, &
                                             qs%estimators%D0_hf_population, &
                                             ntot_particles
        else if (comp_set) then
            write (6,'(i10,2X,es17.10,2X,2(es20.10,2X),2(es17.10,2X),2X,es17.10)', advance='no') &
                                         qs%mc_cycles_done+mc_cycles, qs%shift(1),   &
                                         real(qs%estimators%proj_energy_comp, p), aimag(qs%estimators%proj_energy_comp), &
                                         real(qs%estimators%D0_population_comp, p), aimag(qs%estimators%D0_population_comp), &
                                         ntot_particles(1) + ntot_particles(2)
        else
            write (6,'(i10,2X,2(es17.10,2X),es17.10,4X,es17.10)', advance='no') &
                                             qs%mc_cycles_done+mc_cycles, qs%shift(1),   &
                                             qs%estimators%proj_energy, qs%estimators%D0_population, &
                                             ntot_particles
        end if
        write (6,'(2X,i10,4X,i12,2X,f7.4,2X,f7.3)') qs%estimators%tot_nstates, qs%estimators%tot_nspawn_events, &
                                             qs%spawn_store%rspawn, elapsed_time/qmc_in%ncycles
    end subroutine write_fciqmc_report

    subroutine end_fciqmc(reference, psip_list, spawn)

        ! Deallocate fciqmc data arrays.

        ! In/Out (optional):
        !   reference: reference state. On exit, allocatable components are deallocated.
        !   psip_list: main particle_t object.  On exit, allocatable components are deallocated.
        !   spawn: spawn_t object.  On exit, allocatable components are deallocated.

        use checking, only: check_deallocate
        use spawn_data, only: dealloc_spawn_t
        use load_balancing, only: dealloc_parallel_t
        use qmc_data, only: particle_t
        use reference_determinant, only: reference_t, dealloc_reference_t

        type(reference_t), intent(inout), optional :: reference
        type(particle_t), intent(inout), optional :: psip_list
        type(spawn_t), intent(inout), optional :: spawn

        integer :: ierr

        if (present(reference)) then
            call dealloc_reference_t(reference)
        end if
        if (present(psip_list)) then
            if (allocated(psip_list%nparticles)) then
                deallocate(psip_list%nparticles, stat=ierr)
                call check_deallocate('psip_list%nparticles',ierr)
            end if
            if (allocated(psip_list%states)) then
                deallocate(psip_list%states, stat=ierr)
                call check_deallocate('psip_list%states',ierr)
            end if
            if (allocated(psip_list%pops)) then
                deallocate(psip_list%pops, stat=ierr)
                call check_deallocate('psip_list%pops',ierr)
            end if
            if (allocated(psip_list%dat)) then
                deallocate(psip_list%dat, stat=ierr)
                call check_deallocate('psip_list%dat',ierr)
            end if
            if (allocated(psip_list%nparticles_proc)) then
                deallocate(psip_list%nparticles_proc, stat=ierr)
                call check_deallocate('psip_list%nparticles_proc', ierr)
            end if
        end if
        if (present(spawn)) call dealloc_spawn_t(spawn)

    end subroutine end_fciqmc

end module fciqmc_data
