module fciqmc_restart

! Module for dumping out and restarting FCIQMC files.

use fciqmc_data
implicit none

character(*), parameter :: restart_file_stem = 'restart'

! If negative, then this is set in the input and we need to read from
! a specific restart file.
integer :: read_restart_number = 0

! If negative, then this is set in the input and we need to write to
! a specific restart file.
integer :: write_restart_number = 0

contains

    subroutine dump_restart(nmc_cycles, nparticles_old)

        ! Write out the main walker list to file.

        use parallel, only: parent
        use utils, only: get_unique_filename, get_free_unit

        integer, intent(in) :: nmc_cycles, nparticles_old
        character(255) :: restart_file
        integer :: io, i
        integer, parameter :: restart_version = 1

        if (parent) then
            io = get_free_unit()
            if (write_restart_number < 0) then
                call get_unique_filename(restart_file_stem, .true., write_restart_number, restart_file)
            else
                call get_unique_filename(restart_file_stem, .true., 0, restart_file)
            end if

            write (6,'(1X,a23,1X,a,a1,/)') 'Writing restart file to',trim(restart_file),'.'

            open(io, file=restart_file)

            write (io,*) 'Restart version'
            write (io,*) restart_version
            write (io,*) '# number of cycles'
            write (io,*) nmc_cycles
            write (io,*) '# shift'
            write (io,*) nparticles_old, shift, vary_shift
            write (io,*) '# reference determinant'
            write (io,*) f0, occ_list0, D0_population, H00
            write (io,*) '# number of unique walkers'
            write (io,*) tot_walkers
            write (io,*) '# walkers'
            do i = 1, tot_walkers
                write (io,*) walker_dets(:,i), walker_population(i), walker_energies(i)
            end do
            close(io)
        end if

    end subroutine dump_restart

    subroutine read_restart()

        ! Read in the main walker list from file.

        use errors, only: stop_all
        use parallel, only: parent
        use utils, only: get_unique_filename, get_free_unit

        character(255) :: restart_file, junk
        integer :: io, i
        logical :: exists
        integer :: restart_version

        if (parent) then
            io = get_free_unit()
            if (read_restart_number < 0) then
                call get_unique_filename(restart_file_stem, .false., read_restart_number, restart_file)
            else
                call get_unique_filename(restart_file_stem, .false., 0, restart_file)
            end if

            write (6,'(1X,a25,1X,a,a1,/)') 'Reading from restart file',trim(restart_file),'.'

            inquire(file=restart_file, exist=exists)

            if (.not.exists) then
                call stop_all('read_restart','restart file '//trim(restart_file)//' does not exist.')
            end if

            open(io, file=restart_file)

            read (io,*) junk
            read (io,*) restart_version
            read (io,*) junk
            read (io,*) mc_cycles_done
            read (io,*) junk
            read (io,*) nparticles_old_restart, shift, vary_shift
            read (io,*) junk
            read (io,*) f0, occ_list0, D0_population, H00
            read (io,*) junk
            read (io,*) tot_walkers
            read (io,*) junk
            do i = 1, tot_walkers
                read (io,*) walker_dets(:,i), walker_population(i), walker_energies(i)
            end do
            close(io)
        end if

    end subroutine read_restart

end module fciqmc_restart
