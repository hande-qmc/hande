module hande_top_level

! A very coarse interface to HANDE.

implicit none

contains

    subroutine init_hande(start_cpu_time, start_wall_time)

        ! Initialise HANDE (minimum global state).
        ! Print out information about the compiled executable.

        ! Out:
        !     start_cpu_time: cpu_time at the start of the calculation.
        !     start_wall_time: system_clock at the start of the calculation.

        use report, only: environment_report, comm_global_uuid, VCS_VERSION, GLOBAL_UUID
        use calc, only: init_calc_defaults
        use parallel, only: init_parallel, parallel_report, iproc, nprocs, nthreads, parent
        use qmc_data, only: qmc_in_global

        real, intent(out) :: start_cpu_time
        integer, intent(out) :: start_wall_time

        call init_parallel()

        call cpu_time(start_cpu_time)
        call system_clock(start_wall_time)

        if (parent) then
            write (6,'(/,a8,/)') 'HANDE'
            call environment_report()
        end if
        call comm_global_uuid()

        call init_calc_defaults(VCS_VERSION, GLOBAL_UUID, qmc_in_global%seed)

        if ((nprocs > 1 .or. nthreads > 1) .and. parent) call parallel_report()

    end subroutine init_hande

    subroutine init_calc(sys, reference)

        ! Initialise the calculation, read input options and initialse the system and
        ! basis functions to be used.

        ! In/Out:
        !     sys: system to be studied.  On input sys has default values.  On
        !          output, its values have been updated according to the input file
        !          and allocatable components have been appropriately allocated.
        !    reference: current reference determinant. If specified as input, it has been
        !          allocated and set on exit.

        use parse_input, only: read_input, check_input, distribute_input
        use system
        use basis, only: init_model_basis_fns
        use basis_types, only: copy_basis_t, dealloc_basis_t, init_basis_strings, print_basis_metadata
        use determinants, only: init_determinants
        use determinant_enumeration, only: init_determinant_enumeration
        use dmqmc_data, only: dmqmc_in_global, dmqmc_estimates_global
        use excitations, only: init_excitations
        use parallel, only: parent
        use real_lattice, only: init_real_space
        use momentum_symmetry, only: init_momentum_symmetry
        use point_group_symmetry, only: print_pg_symmetry_info
        use read_in_system, only: read_in_integrals
        use calc
        use ueg_system, only: init_ueg_proc_pointers
        use qmc_data, only: qmc_in_global, semi_stoch_in_global, fciqmc_in_global, &
                            ccmc_in_global, restart_in_global, load_bal_in_global
        use qmc_data, only: reference_t


        type(sys_t), intent(inout) :: sys
        type(reference_t), intent(inout) :: reference
        
        if (parent) call read_input(sys, qmc_in_global, fciqmc_in_global, ccmc_in_global, &
                                    semi_stoch_in_global, restart_in_global, reference, &
                                    load_bal_in_global, dmqmc_in_global, dmqmc_estimates_global%rdm_info)

        call distribute_input(sys, qmc_in_global, fciqmc_in_global, ccmc_in_global, &
                              semi_stoch_in_global, restart_in_global, load_bal_in_global, &
                              reference, dmqmc_in_global, dmqmc_estimates_global%rdm_info)

        call init_system(sys)
        ! Note: can't set ex_level to be full space if not set until *after* sys%nel is defined.
        ! It's not set until init_system for the Heisenberg model.
        if (reference%ex_level < 0) reference%ex_level = sys%nel

        call check_input(sys, qmc_in_global, fciqmc_in_global, ccmc_in_global, semi_stoch_in_global, &
                         restart_in_global, reference, load_bal_in_global, dmqmc_in_global)

        ! Initialise basis functions.
        if (sys%system == read_in) then
            call read_in_integrals(sys, cas_info=sys%cas)
        else
            call init_model_basis_fns(sys)
        end if

        call init_basis_strings(sys%basis)
        call print_basis_metadata(sys%basis, sys%nel, sys%system == heisenberg)
        call init_determinants(sys)
        call init_determinant_enumeration()

        call init_excitations(sys%basis)

        ! System specific.
        select case(sys%system)
        case(ueg)
            call init_momentum_symmetry(sys)
            call init_ueg_proc_pointers(sys%lattice%ndim, sys%ueg)
        case(hub_k)
            call init_momentum_symmetry(sys)
        case(hub_real, heisenberg, chung_landau)
            call init_real_space(sys)
        case(read_in)
            call print_pg_symmetry_info(sys)
        end select

    end subroutine init_calc

    subroutine run_calc(sys, reference)

        ! Run the calculation based upon the input options.

        ! In/Out:
        !    sys: system to be studied.  Note: sys may be altered during the
        !    calculation procedure but should be unaltered on exit of each
        !    calculation procedure.
        !    reference: reference determinant. If given as input, that is used, otherwise it is set
        !    during initialisation.

        use calc
        use diagonalisation, only: diagonalise
        use hilbert_space, only: estimate_hilbert_space
        use canonical_kinetic_energy, only: estimate_kinetic_energy
        use parallel, only: iproc, parent
        use qmc_data, only: qmc_in_global, semi_stoch_in_global, fciqmc_in_global, &
                            ccmc_in_global, restart_in_global, load_bal_in_global, &
                            annihilation_flags_global
        use qmc_data, only: reference_t
        use dmqmc_data, only: dmqmc_in_global, dmqmc_estimates_global
        use simple_fciqmc, only: do_simple_fciqmc
        use system, only: sys_t

        type(sys_t), intent(inout) :: sys
        type(reference_t), intent(inout) :: reference

        if (doing_calc(exact_diag+lanczos_diag)) call diagonalise(sys, reference)

        if (doing_calc(mc_hilbert_space)) then
            call estimate_hilbert_space(sys, reference, qmc_in_global%seed)
        end if

        if (doing_calc(mc_canonical_kinetic_energy)) then
            call estimate_kinetic_energy(sys, qmc_in_global, dmqmc_in_global)
        end if

        if (doing_calc(fciqmc_calc+hfs_fciqmc_calc+ct_fciqmc_calc+dmqmc_calc+ccmc_calc)) then
            if (doing_calc(simple_fciqmc_calc)) then
                call do_simple_fciqmc(sys, qmc_in_global, restart_in_global, reference)
            else
                call do_qmc(sys, qmc_in_global, fciqmc_in_global, ccmc_in_global, semi_stoch_in_global, &
                            restart_in_global, reference, load_bal_in_global, dmqmc_in_global, &
                            dmqmc_estimates_global, annihilation_flags_global)
            end if
        end if

    end subroutine run_calc

    subroutine end_calc(sys, reference)

        ! Clean up time!

        ! In/Out:
        !     sys: main system object.  All allocatable components are
        !          deallocated on exit.
        !     reference: defines the reference state. Allocatable components are
        !          deallocated on exit.

        use basis_types, only: dealloc_basis_t
        use calc
        use system, only: sys_t, end_lattice_system
        use determinants, only: end_determinants
        use excitations, only: end_excitations
        use diagonalisation, only: end_hamil
        use fciqmc_data, only: end_fciqmc
        use real_lattice, only: end_real_space
        use momentum_symmetry, only: end_momentum_symmetry
        use molecular_integrals, only: end_one_body_t, end_two_body_t
        use qmc_data, only: fciqmc_in_global
        use report, only: end_report
        use ueg_system, only: end_ueg_indexing
        use qmc_data, only: reference_t

        type(sys_t), intent(inout) :: sys
        type(reference_t), intent(inout) :: reference

        ! Deallocation routines.
        ! NOTE:
        !   end_ routines should surround every deallocate statement with a test
        !   that the array is allocated.
        call end_lattice_system(sys%lattice, sys%k_lattice, sys%real_lattice)

        call end_one_body_t(sys%read_in%one_e_h_integrals)
        call end_one_body_t(sys%read_in%one_body_op_integrals)
        call end_two_body_t(sys%read_in%coulomb_integrals)

        call dealloc_basis_t(sys%basis)
        call end_excitations(sys%basis%excit_mask)
        call end_momentum_symmetry()
        call end_ueg_indexing(sys%ueg)
        call end_determinants()
        call end_hamil()
        call end_real_space(sys%heisenberg)
        call end_fciqmc(fciqmc_in_global%non_blocking_comm, reference)

    end subroutine end_calc

    subroutine end_hande(start_cpu_time, start_wall_time)

        ! Terminate HANDE: print timing report and terminate MPI stack...

        ! In:
        !     start_cpu_time: cpu_time at the start of the calculation.
        !     start_wall_time: system_clock at the start of the calculation.

        use parallel, only: parent, end_parallel
        use report, only: end_report

        real, intent(in) :: start_cpu_time
        integer, intent(in) :: start_wall_time
        real :: end_cpu_time, wall_time
        integer :: end_wall_time, count_rate, count_max

        ! Calculation time.
        call cpu_time(end_cpu_time)
        call system_clock(end_wall_time, count_rate, count_max)
        if (end_wall_time < start_wall_time) then
            ! system_clock returns the time modulo count_max.
            ! Have ticked over to the next "block" (assume only one as this
            ! happens roughly once every 1 2/3 years with gfortran!)
            end_wall_time = end_wall_time + count_max
        end if
        wall_time = real(end_wall_time-start_wall_time)/count_rate
        if (parent) call end_report(wall_time, end_cpu_time-start_cpu_time)

        call end_parallel()

    end subroutine end_hande

! --- QMC wrapper ---

    subroutine do_qmc(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, reference, load_bal_in, &
                      dmqmc_in, dmqmc_estimates, annihilation_flags)

        ! Initialise and run stochastic quantum chemistry procedures.

        ! In/Out:
        !    sys: system being studied.  This should(!) be returned unaltered on
        !         output from each procedure, but might be varied during the
        !         run if needed.
        !    qmc_in: input options relating to QMC methods.
        !    fciqmc_in: input options relating to FCIQMC.
        !    reference: the reference determinant.
        !    dmqmc_in: input options relating to DMQMC.
        !    dmqmc_estimates: type containing all DMQMC estimates.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In:
        !    ccmc_in: input options relating to CCMC.
        !    semi_stoch_in: Input options for the semi-stochastic adaptation.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.

        use calc

        use ccmc, only: do_ccmc
        use ct_fciqmc, only: do_ct_fciqmc
        use dmqmc, only: do_dmqmc
        use fciqmc, only: do_fciqmc
        use hellmann_feynman_sampling, only: do_hfs_fciqmc

        use qmc_data, only: qmc_in_t, fciqmc_in_t, ccmc_in_t, semi_stoch_in_t
        use qmc_data, only: restart_in_t, reference_t, load_bal_in_t, annihilation_flags_t
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t
        use qmc, only: init_proc_pointers
        use system, only: sys_t, copy_sys_spin_info, set_spin_polarisation, heisenberg, ueg, read_in
        use parallel, only: nprocs

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(fciqmc_in_t), intent(inout) :: fciqmc_in
        type(ccmc_in_t), intent(inout) :: ccmc_in
        type(semi_stoch_in_t), intent(in) :: semi_stoch_in
        type(restart_in_t), intent(in) :: restart_in
        type(reference_t), intent(inout) :: reference
        type(load_bal_in_t), intent(inout) :: load_bal_in
        type(dmqmc_in_t), intent(inout) :: dmqmc_in
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        type(annihilation_flags_t), intent(inout) :: annihilation_flags

        real(p) :: hub_matel
        type(sys_t) :: sys_bak

        ! Initialise procedure pointers.
        call init_proc_pointers(sys, qmc_in, dmqmc_in, reference)

        ! Set spin variables.
        ! [todo] - handle all_spin_sectors more gracefully.  It should probably be handled by DMQMC-specific code but must be done before set_spin_polarisation.
        if (dmqmc_in%all_spin_sectors) then
            if (sys%system /= read_in .and. sys%system /= ueg) ms_in = sys%lattice%nsites
            if (sys%system == heisenberg) sys%max_number_excitations = sys%lattice%nsites/2
        end if
        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        ! Calculation-specifc initialisation and then run QMC calculation.

        if (doing_calc(dmqmc_calc)) then
            call do_dmqmc(sys, qmc_in, dmqmc_in, dmqmc_estimates, restart_in, load_bal_in, reference)
        else if (doing_calc(ct_fciqmc_calc)) then
            call do_ct_fciqmc(sys, qmc_in, restart_in, reference, load_bal_in, annihilation_flags, hub_matel)
        else if (doing_calc(ccmc_calc)) then
            call do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference)
        else
            ! Doing FCIQMC calculation (of some sort) using the original
            ! timestep algorithm.
            if (doing_calc(hfs_fciqmc_calc)) then
                call do_hfs_fciqmc(sys, qmc_in, restart_in, load_bal_in, reference)
            else
                call do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, load_bal_in, reference)
            end if
        end if

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine do_qmc

end module hande_top_level
