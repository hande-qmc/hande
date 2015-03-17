module hande_top_level

! A very coarse interface to HANDE.

implicit none

contains

    subroutine init_calc(sys, reference, start_cpu_time, start_wall_time)

        ! Initialise the calculation.
        ! Print out information about the compiled executable,
        ! read input options and initialse the system and basis functions
        ! to be used.

        ! In/Out:
        !     sys: system to be studied.  On input sys has default values.  On
        !          output, its values have been updated according to the input file
        !          and allocatable components have been appropriately allocated.
        !    reference: current reference determinant. If specified as input, it has been
        !          allocated and set on exit.
        ! Out:
        !     start_cpu_time: cpu_time at the start of the calculation.
        !     start_wall_time: system_clock at the start of the calculation.

        use report, only: environment_report, comm_global_uuid, VCS_VERSION, GLOBAL_UUID
        use parse_input, only: read_input, check_input, distribute_input
        use system
        use basis, only: init_model_basis_fns
        use basis_types, only: copy_basis_t, dealloc_basis_t, init_basis_strings, print_basis_metadata
        use calc
        use determinants, only: init_determinants
        use determinant_enumeration, only: init_determinant_enumeration
        use excitations, only: init_excitations
        use parallel, only: init_parallel, parallel_report, iproc, nprocs, nthreads, parent
        use real_lattice, only: init_real_space
        use momentum_symmetry, only: init_momentum_symmetry
        use point_group_symmetry, only: print_pg_symmetry_info
        use qmc_data, only: qmc_in, semi_stoch_in, fciqmc_in, ccmc_in, restart_in, load_bal_in
        use qmc_data, only: reference_t
        use dmqmc_data, only: dmqmc_in
        use read_in_system, only: read_in_integrals
        use ueg_system, only: init_ueg_proc_pointers

        type(sys_t), intent(inout) :: sys
        type(reference_t), intent(inout) :: reference
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

        call init_calc_defaults(VCS_VERSION, GLOBAL_UUID, qmc_in%seed)

        if ((nprocs > 1 .or. nthreads > 1) .and. parent) call parallel_report()

        if (parent) call read_input(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, reference, &
                                    load_bal_in, dmqmc_in)

        call distribute_input(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, &
                              reference, dmqmc_in)

        call init_system(sys, dmqmc_in)

        call check_input(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, reference, load_bal_in, &
                         dmqmc_in)

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
        use qmc, only: do_qmc
        use hilbert_space, only: estimate_hilbert_space
        use canonical_kinetic_energy, only: estimate_kinetic_energy
        use parallel, only: iproc, parent
        use qmc_data, only: qmc_in, semi_stoch_in, fciqmc_in, ccmc_in, restart_in, load_bal_in
        use qmc_data, only: reference_t
        use dmqmc_data, only: dmqmc_in
        use simple_fciqmc, only: do_simple_fciqmc
        use system, only: sys_t

        type(sys_t), intent(inout) :: sys
        type(reference_t), intent(inout) :: reference

        if (doing_calc(exact_diag+lanczos_diag)) call diagonalise(sys, reference)

        if (doing_calc(mc_hilbert_space)) then
            call estimate_hilbert_space(sys, reference, qmc_in%seed)
        end if

        if (doing_calc(mc_canonical_kinetic_energy)) then
            call estimate_kinetic_energy(sys, qmc_in)
        end if

        if (doing_calc(fciqmc_calc+hfs_fciqmc_calc+ct_fciqmc_calc+dmqmc_calc+ccmc_calc)) then
            if (doing_calc(simple_fciqmc_calc)) then
                call do_simple_fciqmc(sys, qmc_in, restart_in, reference)
            else
                call do_qmc(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, reference, load_bal_in, &
                            dmqmc_in)
            end if
        end if

    end subroutine run_calc

    subroutine end_calc(sys, reference, start_cpu_time, start_wall_time)

        ! Clean up time!

        ! In:
        !     start_cpu_time: cpu_time at the start of the calculation.
        !     start_wall_time: system_clock at the start of the calculation.
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
        use parallel, only: parent, end_parallel
        use real_lattice, only: end_real_space
        use momentum_symmetry, only: end_momentum_symmetry
        use molecular_integrals, only: end_one_body_t, end_two_body_t
        use qmc_data, only: fciqmc_in
        use report, only: end_report
        use ueg_system, only: end_ueg_indexing
        use qmc_data, only: reference_t

        type(sys_t), intent(inout) :: sys
        type(reference_t), intent(inout) :: reference
        real, intent(in) :: start_cpu_time
        integer, intent(in) :: start_wall_time
        real :: end_cpu_time, wall_time
        integer :: end_wall_time, count_rate, count_max

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
        call end_fciqmc(fciqmc_in%non_blocking_comm, reference)

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

    end subroutine end_calc

end module hande_top_level
