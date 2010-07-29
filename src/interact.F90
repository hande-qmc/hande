module interact

! Module for interacting with running FCIQMC calculations.

implicit none

contains

    subroutine fciqmc_interact(ireport, soft_exit)

        ! Read FCIQMC.COMM if it exists in the working directory of any
        ! processor and set the variables according to the options defined in
        ! FCIQMC.COMM.

        ! In:
        !    ireport: index of the current report loop.
        ! Out:
        !    softexit: true if SOFTEXIT is defined in FCIQMC.COMM, in which case
        !        any fciqmc calculation should exit immediately and go to the
        !        post-processing steps.

        use input, line=>char
        use utils, only: get_free_unit
        use parallel

        use fciqmc_data, only: target_particles, tau, av_shift, av_proj_energy, &
                               av_D0_population, vary_shift, start_vary_shift,  &
                               start_averaging_from

        integer, intent(in) :: ireport
        logical, intent(out) :: soft_exit

        logical :: comms_exists, comms_found, comms_read, eof
        integer :: proc
        integer :: ierr

        character(100) :: w

        ! Note that all output in this subroutine *must* be prepended with #.
        ! This enables the blocking script to ignore these lines whilst doing
        ! data analysis.

        soft_exit = .false.

        inquire(file='FCIQMC.COMM', exist=comms_exists)

#ifdef PARALLEL
        call mpi_allreduce(comms_exists, comms_found, 1, mpi_logical, mpi_lor, mpi_comm_world, ierr)
#else
        comms_found = comms_exists
#endif

        if (comms_found) then
            ! Read in the FCIQMC.COMM file.
            ! This should be a very rare event, so we don't worry too much
            ! about optimised communications in this section.
            if (parent) then
                write (6,'(1X,"#",1X,62("-"))')
                write (6,'(1X,"#",1X,a21)') 'FCIQMC.COMM detected.'
                write (6,'(1X,"#",/,1X,"#",1X,a24,/,1X,"#")') 'Contents of FCIQMC.COMM:'
                ! Flush output from parent processor so that processor which
                ! has the FCIQMC.COMM file can print out the contents without 
                ! mixing the output.
                call flush(6)
            end if
            ! Quick pause to ensure output is all done by this point.
#ifdef PARALLEL
            call mpi_barrier(mpi_comm_world, ierr)
#endif
            ! Slightly tricky bit: need to take into account multi-core
            ! machines where multiple processors can share the same disk and so
            ! be picking up the same FCIQMC.COMM file.  We want to ensure that
            ! only one processor reads it in (avoid race conditions!).
            ! Solution: loop over processors and place a blocking comms call at
            ! the end of each iteration.
            ! proc will end up holding the processor id that read in
            ! FCIQMC.COMM.
            comms_read = .false.
            do proc = 0, nprocs-1
                if (proc == iproc .and. comms_exists) then
                    ! Read in file.
                    ir = get_free_unit()
                    open(ir, file='FCIQMC.COMM', status='old')
                    ! Will do our own echoing as want to prepend lines with '#'
                    call input_options(echo_lines=.false., skip_blank_lines=.true.)

                    do ! loop over lines in FCIQMC.COMM.
                        call read_line(eof)
                        if (eof) exit
                        write (6,'(1X,"#",1X,a)') trim(line)
                        
                        call readu(w)
                        select case(w)
                        case('SOFTEXIT')
                            ! Exit FCIQMC immediately.
                            soft_exit = .true.
                        case('TAU')
                            ! Change timestep.
                            call readf(tau)
                        case('VARYSHIFT_TARGET')
                            call readi(target_particles)
                            if (target_particles < 0) then
                                ! start varying the shift now.
                                vary_shift = .true.
                                start_vary_shift = ireport
                            end if
                        case('ZERO_MEANS')
                            av_proj_energy = 0.0_p
                            av_D0_population = 0.0_p
                            av_shift = 0.0_p
                            start_averaging_from = ireport
                        case default
                            write (6, '(1X,"#",1X,a24,1X,a)') 'Unknown keyword ignored:', trim(w)
                        end select

                    end do ! end reading of FCIQMC.COMM.

                    ! Don't want to keep FCIQMC.COMM around to be detected
                    ! again on the next FCIQMC iteration.
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
            call mpi_bcast(av_proj_energy, 1, mpi_preal, proc, mpi_comm_world, ierr)
            call mpi_bcast(av_shift, 1, mpi_preal, proc, mpi_comm_world, ierr)
            call mpi_bcast(tau, 1, mpi_preal, proc, mpi_comm_world, ierr)
            call mpi_bcast(target_particles, 1, mpi_integer, proc, mpi_comm_world, ierr)
            call mpi_bcast(vary_shift, 1, mpi_logical, proc, mpi_comm_world, ierr)
            call mpi_bcast(start_vary_shift, 1, mpi_integer, proc, mpi_comm_world, ierr)
            call mpi_bcast(start_averaging_from, 1, mpi_integer, proc, mpi_comm_world, ierr)
#endif

            if (parent) write (6,'(1X,"#",/,1X,"#",1X,a59,/,1X,"#",1X,62("-"))')  & 
                   "From now on we use the information provided in FCIQMC.COMM."
        end if

    end subroutine fciqmc_interact

end module interact
