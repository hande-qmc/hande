module lua_hande_calc_utils

! Subroutines utilised within lua_hande_calc to initialise
! hande calculations without directly interacting with the
! lua API.

implicit none

contains

    subroutine init_output_unit(output_in, sys, io_unit)

        ! Initialises file if filename specified in input and
        ! returns the corresponding io unit.
        ! If not specified will just direct all output to stdout.
        ! If not writing to stdout will write message to stdout.

        ! In:
        !   output_in: derived type containing input parameters
        !       related to calculation output.
        !   sys: system object. Used to reprint all required
        !       information if requested by user.
        ! Out:
        !   io_unit: io unit allocated to given filename. If writing
        !       to stdout, set to 6.

        use qmc_data, only: output_in_t
        use report, only: environment_report
        use calc, only: calc_type, get_calculation_string
        use system, only: sys_t
        use parallel

#ifdef PARALLEL
        integer :: ierr
#endif
        type(output_in_t), intent(in) :: output_in
        type(sys_t), intent(in) :: sys
        integer, intent(out) :: io_unit

        if (output_in%out_filename == 'stdout') then
            ! Filename is still the default.
            io_unit = 6
        else
            ! Different filename; write note to stdout for clarity
            ! and initialise io_unit.

            if (parent) then
                write (6,'(1X,"Calculation")')
                write (6,'(1X,"-----------")')
                write (6,'(1X)')
                write (6,'(1X,"Writing ",a," calculation output to",1X,a,"...")') trim(get_calculation_string(calc_type)), &
                                    trim(output_in%out_filename)


                open(newunit=io_unit, file=output_in%out_filename, &
                        status='unknown')

                write (io_unit,'(/,a8,/)') 'HANDE'
                call environment_report(io=io_unit)
                if (nthreads > 1 .or. nprocs > 1) call parallel_report(io_unit)
                if (output_in%reprint_sys_info) then
                    ! Want to reprint all system information in output file
                    ! for ease of use.
                    write (io_unit,'(1X,"Reprinting all system information as requested.")')
                    write (io_unit,'(1X)')
                    call reprint_sys_info(sys, io_unit)
                else
                    write (io_unit,'(1X,"System information previously written to stdout.")')
                    write (io_unit,'(1X)')
                end if
            end if
#ifdef PARALLEL
            call mpi_bcast(io_unit, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
#endif
        end if

    end subroutine init_output_unit

    subroutine end_output_unit(filename, io_unit)

        ! Close calculation output file with the date and time.

        ! In:
        !   filename: filename calculation is writing to.
        !   io_unit: io unit to be closed.

        use report, only: write_date_time_close
        use calc, only: calc_type, get_calculation_string
        use parallel, only: parent

        character(255), intent(in) :: filename
        integer, intent(in) :: io_unit
        integer :: date_values(8)

        if (parent) then
            if (filename /= 'stdout') then
                call date_and_time(VALUES=date_values)

                call write_date_time_close(io_unit, date_values)

                write (6,'(1X,"Finished writing ",a," calculation output to",1X,a)') trim(get_calculation_string(calc_type)), &
                                    trim(filename)
                write (6,'(1x)')
            end if
        end if

    end subroutine end_output_unit

    subroutine reprint_sys_info(sys, io_unit)

        ! Reprints all information from system usually written during
        ! initialisation. For use when writing calculations to separate
        ! files.

        use system, only: sys_t, heisenberg, hub_k, read_in
        use basis, only: write_basis_fn_header, write_basis_fn
        use basis_types, only: print_basis_metadata
        use parallel, only: parent
        use momentum_symmetry, only: print_hubbard_k_symmetry_info
        use point_group_symmetry, only: print_pg_symmetry_info
        use momentum_sym_read_in, only: print_mom_sym_info

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: io_unit
        integer :: i

        if (parent) then
            call write_basis_fn_header(sys, iunit=io_unit)
            do i = 1, sys%basis%nbasis
                call write_basis_fn(sys, sys%basis%basis_fns(i), ind=i, iunit=io_unit, new_line=.true.)
            end do
            write (io_unit,'(/,1X,a8,f18.12)') 'E_core =', sys%read_in%Ecore
            call print_basis_metadata(sys%basis, sys%nel, sys%system == heisenberg, io_unit=io_unit)
            select case(sys%system)
            case(hub_k)
                call print_hubbard_k_symmetry_info(sys, io_unit)
            case(read_in)
                if (sys%momentum_space) then
                    call print_mom_sym_info(sys, io_unit)
                else
                    call print_pg_symmetry_info(sys, io_unit)
                end if
            end select
        end if

    end subroutine reprint_sys_info

end module lua_hande_calc_utils
