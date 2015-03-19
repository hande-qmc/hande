module interact

! Module for interacting with running calculations.

implicit none

character(*), parameter :: comms_file = "HANDE.COMM"

contains

    subroutine calc_interact(comms_found, soft_exit, qmc_in)

        ! Read HANDE.COMM if it exists in the working directory of any
        ! processor and set the variables according to the options defined in
        ! HANDE.COMM.

        ! Out:
        !    softexit: true if SOFTEXIT is defined in HANDE.COMM, in which case
        !        any calculation should exit immediately and go to the
        !        post-processing steps.
        ! In/Out:
        !    qmc_in (optional): Input options relating to QMC methods.

        use input, line=>char
        use utils, only: get_free_unit
        use parallel

        use fciqmc_data, only: vary_shift, shift
        use qmc_data, only: qmc_in_t, walker_global

        logical, intent(in) :: comms_found
        logical, intent(out) :: soft_exit
        type(qmc_in_t), optional, intent(inout) :: qmc_in

        logical :: comms_exists, comms_read, eof
        integer :: proc, i
#ifdef PARALLEL
        integer :: ierr
#endif

        character(100) :: w

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
                    ir = get_free_unit()
                    open(ir, file=comms_file, status='old')
                    ! Will do our own echoing as want to prepend lines with '#'
                    call input_options(echo_lines=.false., skip_blank_lines=.true.)

                    do ! loop over lines in HANDE.COMM.
                        call read_line(eof)
                        if (eof) exit
                        write (6,'(1X,"#",1X,a)') trim(line)

                        call readu(w)
                        select case(w)
                        case('SOFTEXIT')
                            ! Exit calculation immediately.
                            soft_exit = .true.
                        case('TAU')
                            ! Change timestep.
                            if (present(qmc_in)) call readf(qmc_in%tau)
                        case('VARYSHIFT_TARGET')
                            call readf(qmc_in%target_particles)
                            if (qmc_in%target_particles < 0) then
                                ! start varying the shift now.
                                vary_shift = .true.
                            end if
                        case('SHIFT')
                            call readf(shift(1))
                            do i = 2, walker_global%sampling_size
                                shift(i) = shift(1)
                            end do
                        case default
                            write (6, '(1X,"#",1X,a24,1X,a)') 'Unknown keyword ignored:', trim(w)
                        end select

                    end do ! end reading of HANDE.COMM.

                    ! Don't want to keep HANDE.COMM around to be detected
                    ! again on the next iteration.
                    close(ir, status="delete")
                    comms_read = .true.
                end if
#ifdef PARALLEL
                call mpi_bcast(comms_read, 1, mpi_logical, proc, mpi_comm_world, ierr)
#endif
                if (comms_read) exit
            end do

#ifdef PARALLEL
            ! If in parallel need to broadcast data.
            call mpi_bcast(soft_exit, 1, mpi_logical, proc, mpi_comm_world, ierr)
            if (present(qmc_in)) call mpi_bcast(qmc_in%tau, 1, mpi_preal, proc, mpi_comm_world, ierr)
            call mpi_bcast(shift, walker_global%sampling_size, mpi_preal, proc, mpi_comm_world, ierr)
            call mpi_bcast(qmc_in%target_particles, 1, mpi_preal, proc, mpi_comm_world, ierr)
            call mpi_bcast(vary_shift, 1, mpi_logical, proc, mpi_comm_world, ierr)
#endif

            if (parent) write (6,'(1X,"#",/,1X,"#",1X,a59,/,1X,"#",1X,62("-"))')  &
                   "From now on we use the information provided in "//comms_file//"."

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
