module hande_top_level

! A very coarse interface to HANDE.

implicit none

contains

    subroutine init_calc(start_cpu_time, start_wall_time)

        ! Initialise the calculation.
        ! Print out information about the compiled executable,
        ! read input options and initialse the system and basis functions
        ! to be used.

        ! Out:
        !     start_cpu_time: cpu_time at the start of the calculation.
        !     start_wall_time: system_clock at the start of the calculation.

        use report, only: environment_report
        use parse_input, only: read_input, check_input, distribute_input
        use system
        use basis, only: init_model_basis_fns
        use determinants, only: init_determinants
        use determinant_enumeration, only: init_determinant_enumeration
        use excitations, only: init_excitations
        use parallel, only: init_parallel, parallel_report, iproc, nprocs, nthreads, parent
        use hubbard_real, only: init_real_space
        use momentum_symmetry, only: init_momentum_symmetry
        use point_group_symmetry, only: print_pg_symmetry_info
        use read_in_system, only: read_in_integrals
        use calc
        use ueg_system, only: init_ueg_proc_pointers

        real, intent(out) :: start_cpu_time
        integer, intent(out) :: start_wall_time

        call init_parallel()

        call cpu_time(start_cpu_time)
        call system_clock(start_wall_time)

        if (parent) then
            write (6,'(/,a8,/)') 'Hubbard'
            call environment_report()
        end if

        if ((nprocs > 1 .or. nthreads > 1) .and. parent) call parallel_report()

        if (parent) call read_input()

        call distribute_input()

        call init_system(sys_global)

        call check_input()

        ! Initialise basis functions.
        if (sys_global%system == read_in) then
            call read_in_integrals(cas_info=sys_global%cas)
        else
            call init_model_basis_fns()
        end if

        call init_determinants()
        call init_determinant_enumeration()

        call init_excitations()

        ! System specific.
        select case(sys_global%system)
        case(ueg)
            call init_momentum_symmetry()
            call init_ueg_proc_pointers()
        case(hub_k)
            call init_momentum_symmetry()
        case(hub_real, heisenberg, chung_landau)
            call init_real_space()
        case(read_in)
            call print_pg_symmetry_info()
        end select

    end subroutine init_calc

    subroutine run_calc()

        ! Run the calculation based upon the input options.

        use calc
        use diagonalisation, only: diagonalise
        use qmc, only: do_qmc
        use hilbert_space, only: estimate_hilbert_space
        use parallel, only: iproc, parent
        use simple_fciqmc, only: do_simple_fciqmc, init_simple_fciqmc
        use utils, only: int_fmt

        if (doing_calc(exact_diag+lanczos_diag)) call diagonalise()

        if (doing_calc(mc_hilbert_space)) then
            call estimate_hilbert_space()
        end if

        if (doing_calc(fciqmc_calc+hfs_fciqmc_calc+ct_fciqmc_calc+dmqmc_calc+ccmc_calc)) then
            if (doing_calc(simple_fciqmc_calc)) then
                call init_simple_fciqmc()
                call do_simple_fciqmc()
            else 
                call do_qmc()
            end if
        end if

    end subroutine run_calc

    subroutine end_calc(start_cpu_time, start_wall_time)

        ! Clean up time!

        ! In:
        !     start_cpu_time: cpu_time at the start of the calculation.
        !     start_wall_time: system_clock at the start of the calculation.

        use calc
        use system
        use basis, only: end_basis_fns
        use determinants, only: end_determinants
        use excitations, only: end_excitations
        use diagonalisation, only: end_hamil
        use fciqmc_data, only: end_fciqmc
        use ifciqmc, only: end_ifciqmc
        use parallel, only: parent, end_parallel
        use hubbard_real, only: end_real_space
        use momentum_symmetry, only: end_momentum_symmetry
        use report, only: end_report

        real, intent(in) :: start_cpu_time
        integer, intent(in) :: start_wall_time
        real :: end_cpu_time, wall_time
        integer :: end_wall_time, count_rate, count_max

        ! Deallocation routines.
        ! NOTE:
        !   end_ routines should surround every deallocate statement with a test
        !   that the array is allocated.
        call end_lattice_system(sys_global%lattice)
        call end_basis_fns()
        call end_momentum_symmetry()
        call end_determinants()
        call end_excitations()
        call end_hamil()
        call end_real_space()
        call end_fciqmc()
        call end_ifciqmc()

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
