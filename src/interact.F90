module interact

! Module for interacting with running calculations.

implicit none

character(*), parameter :: comms_file = "HANDE.COMM"

contains

    subroutine calc_interact(comms_found, soft_exit, qmc_in, qs)

        ! Read HANDE.COMM if it exists in the working directory of any
        ! processor and set the variables according to the options defined in
        ! HANDE.COMM.

        ! In:
        !    comms_found: true if the file HANDE.COMM exists on nay processor.
        ! Out:
        !    softexit: true if SOFTEXIT is defined in HANDE.COMM, in which case
        !        any calculation should exit immediately and go to the
        !        post-processing steps.
        ! In/Out:
        !    qmc_in (optional): Input options relating to QMC methods.
        !    qs (optional): QMC calculation state. The shift and/or timestep may be updated.

        use aotus_module, only: open_config_chunk
        use aot_table_module, only: aot_get_val
        use aot_vector_module, only: aot_get_val
        use flu_binding, only: flu_State

        use utils, only: get_free_unit, read_file_to_buffer
        use parallel

        use qmc_data, only: qmc_in_t, qmc_state_t

        logical, intent(in) :: comms_found
        logical, intent(out) :: soft_exit
        type(qmc_in_t), optional, intent(inout) :: qmc_in
        type(qmc_state_t), optional, intent(inout) :: qs

        logical :: comms_exists, comms_read, eof
        integer :: proc, i, j, ierr, lua_err, iunit
        integer, allocatable :: ierr_arr(:)
#ifdef PARALLEL
        integer :: buf_len
#endif
#if ! defined(__GNUC__) || __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ > 7))
        character(:), allocatable :: buffer
#else
        character(1024**2) :: buffer
#endif
        type(flu_State) :: lua_state
        character(255) :: err_str

        ! Note that all output in this subroutine *must* be prepended with #.
        ! This enables the blocking script to ignore these lines whilst doing
        ! data analysis.

        soft_exit = .false.

        if (comms_found) then
            ! Check if file is on *this* process
            comms_exists = check_comms_file()

            ! Read in the HANDE.COMM file.
            ! This should be a very rare event, so we don't worry too much
            ! about optimised communications in this section.
            if (parent) then
                write (6,'(1X,"#",1X,62("-"))')
                write (6,'(1X,"#",1X,a21)') comms_file//' detected.'
                write (6,'(1X,"#",/,1X,"#",1X,a24,/,1X,"#")') 'Contents of '//comms_file//':'
                ! Flush output from parent processor so that processor which
                ! has the HANDE.COMM file can print out the contents without
                ! mixing the output.
                flush(6)
            end if
            ! Quick pause to ensure output is all done by this point.
#ifdef PARALLEL
            call mpi_barrier(mpi_comm_world, ierr)
#endif
            ! Slightly tricky bit: need to take into account multi-core
            ! machines where multiple processors can share the same disk and so
            ! be picking up the same HANDE.COMM file.  We want to ensure that
            ! only one processor reads it in (avoid race conditions!).
            ! Solution: loop over processors and place a blocking comms call at
            ! the end of each iteration.
            ! proc will end up holding the processor id that read in
            ! HANDE.COMM.
            comms_read = .false.
            do proc = 0, nprocs-1
                if (proc == iproc .and. comms_exists) then
                    ! Read in file.
                    iunit = get_free_unit()
                    open(iunit, file=comms_file, status='old')
                    call read_file_to_buffer(buffer, in_unit=iunit)
                    ! Don't want to keep HANDE.COMM around to be detected again on
                    ! the next Monte Carlo iteration.
                    close(iunit, status="delete")
                    comms_read = .true.
                end if
#ifdef PARALLEL
                call mpi_bcast(comms_read, 1, mpi_logical, proc, mpi_comm_world, ierr)
#endif
                if (comms_read) exit
            end do

#ifdef PARALLEL
            call mpi_bcast(buf_len, 1, MPI_INTEGER, proc, mpi_comm_world, ierr)
#if ! defined(__GNUC__) || __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ > 7))
            if (.not.parent) allocate(character(len=buf_len) :: buffer)
#endif
            call mpi_bcast(buffer, buf_len, MPI_CHARACTER, proc, mpi_comm_world, ierr)
#endif

            ! Only print out the HANDE.COMM file from the parent processor;
            ! prepend each line with # to ease data extraction.
            if (parent) then
                i = 1
                do
                    j = index(buffer(i:), new_line(buffer))
                    if (j == 0) exit
                    write (6,'(a1,a)', advance='no') '#', trim(buffer(i:j))
                    i = j+1
                end do
                write (6,'(1X, "#", a)') trim(buffer(i:))
            end if

            ! Now each processor has the HANDE.COMM script, attempt to execute it...
            call open_config_chunk(lua_state, buffer, lua_err, err_str)

            if (lua_err == 0) then
                ! ... and get variables from global state.
                call aot_get_val(soft_exit, ierr, lua_state, key='softexit')
                if (present(qs)) then
                    call aot_get_val(qs%tau, ierr, lua_state, key='tau')
                    call aot_get_val(qs%shift, ierr_arr, size(qs%shift), lua_state, key='shift')
                end if
                if (present(qmc_in)) then
                    call aot_get_val(qmc_in%target_particles, ierr, lua_state, key='target_population')
                    if (qmc_in%target_particles < 0 .and. present(qs)) qs%vary_shift = .true.
                end if
            end if

            if (parent) then
                if (lua_err == 0) then
                    write (6,'(1X,"#",/,1X,"#",1X,a)') "From now on we use the information provided in "//comms_file//"."
                else
                    write (6,'(1X,"# aotus/lua error code:", i3)') ierr
                    write (6,'(1X,"# error message:", a)') trim(err_str)
                    write (6,'(1X,"# Ignoring variables in ",a)') trim(comms_file)
                end if
                write (6,'(1X,"#",1X,62("-"))')
            end if

        end if

    end subroutine calc_interact

    subroutine check_interact(comms_found)

        ! Checks if there is a HANDE.COMM file present to interact with the calculation

        ! In/Out:
        !   comms_found: on entry, whether HANDE.COMM exists on this processor; on exit whether it
        !   exists on any

        use parallel

        logical, intent(inout) :: comms_found

        logical :: comms_found_any
        integer :: ierr

#ifdef PARALLEL
        call mpi_allreduce(comms_found, comms_found_any, 1, mpi_logical, mpi_lor, mpi_comm_world, ierr)
        comms_found = comms_found_any
#endif

        end subroutine check_interact

        function check_comms_file()

            ! Test whether HANDE.COMM is present on this processor

            logical :: check_comms_file

            inquire(file=comms_file, exist=check_comms_file)
            
        end function check_comms_file

end module interact
