module fciqmc_restart

! Module for dumping out and restarting FCIQMC files.

use parallel
use utils, only: get_unique_filename, get_free_unit

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

        integer, intent(in) :: nmc_cycles, nparticles_old
        character(255) :: restart_file
        integer :: io
        integer, parameter :: restart_version = 1
#ifdef PARALLEL
        integer :: nwalkers(0:nprocs-1), ierr, stat(MPI_STATUS_SIZE), i
        integer, parameter :: comm_tag = 123
        character(255) :: junk

        ! Total number of walkers on each processor.
        call mpi_gather(tot_walkers, 1, mpi_integer, nwalkers, 1, mpi_integer, root, mpi_comm_world, ierr)

#endif

        if (parent) then
            io = get_free_unit()
            if (write_restart_number < 0) then
                call get_unique_filename(restart_file_stem, .true., write_restart_number, restart_file)
            else
                call get_unique_filename(restart_file_stem, .true., 0, restart_file)
            end if

            write (6,'(1X,a23,1X,a,a1,/)') 'Writing restart file to',trim(restart_file),'.'

            open(io, file=restart_file)

            write (io,*) '# restart version'
            write (io,*) restart_version
            write (io,*) '# number of cycles'
            write (io,*) nmc_cycles
            write (io,*) '# shift'
            write (io,*) nparticles_old, shift, vary_shift
            write (io,*) '# reference determinant'
            write (io,*) f0, occ_list0, D0_population, H00
            write (io,*) '# number of unique walkers'
#ifdef PARALLEL
            write (io,*) sum(nwalkers)
            ! Write out walkers on parent processor to restart file.
            write (io,*) '# walker info'
            call write_walkers(tot_walkers, io)

            ! Communicate with all other processors.
            do i = 1, nprocs-1
                ! Receive walker infor from all other processors.
                call mpi_recv(walker_population, nwalkers(i), mpi_integer, i, comm_tag, mpi_comm_world, stat, ierr)
                call mpi_recv(walker_dets, nwalkers(i), mpi_det_integer, i, comm_tag, mpi_comm_world, stat, ierr)
                call mpi_recv(walker_energies, nwalkers(i), mpi_preal, i, comm_tag, mpi_comm_world, stat, ierr)
                ! Write out walkers from all other processors.
                call write_walkers(nwalkers(i), io)
            end do

            ! Read "self" info back in.
            call flush(io)
            rewind(io)
            do
                ! Read restart file until we've found the start of the
                ! walker information.
                read (io,'(a255)') junk
                call flush(6)
                if (index(junk,'walker info') /= 0) exit
            end do
            ! The next tot_walkers lines contain the walker info that came
            ! from the root processor.
            do i = 1, tot_walkers
                read (io, *) walker_dets(:,i), walker_population(i), walker_energies(i)
            end do
#else
            write (io,*) tot_walkers
            write (io,*) '# walker info'
            call write_walkers(tot_walkers, io)
#endif
            close(io)

        else
#ifdef PARALLEL
            ! Send walker info to root processor.
            call mpi_send(walker_population, tot_walkers, mpi_integer, root, comm_tag, mpi_comm_world, ierr)
            call mpi_send(walker_dets, tot_walkers, mpi_det_integer, root, comm_tag, mpi_comm_world, ierr)
            call mpi_send(walker_energies, tot_walkers, mpi_preal, root, comm_tag, mpi_comm_world, ierr)
#endif
        end if

        contains

            subroutine write_walkers(my_nwalkers, iunit)

                integer, intent(in) :: my_nwalkers, iunit
                integer :: iwalker

                do iwalker = 1, my_nwalkers
                    write (iunit,*) walker_dets(:,iwalker), walker_population(iwalker), walker_energies(iwalker)
                end do

            end subroutine write_walkers

    end subroutine dump_restart

    subroutine read_restart()

        ! Read in the main walker list from file.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use hashing, only: murmurhash_bit_string

        use basis, only: basis_length
        use system, only: nel

        character(255) :: restart_file, junk
        integer :: io, i
        logical :: exists
        integer :: restart_version
#ifdef PARALLEL
        integer :: global_tot_walkers, pop, iread, nread, ierr, dest
        real(p), allocatable :: scratch_energies(:)
        integer :: send_counts(0:nprocs-1), send_displacements(0:nprocs-1)
        real(p) :: energy
        integer(i0) :: det(basis_length)
        integer :: spawn_max(0:nprocs-1)
        logical :: done
#endif

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
        end if

        ! Just need to read in the walker information now.
#ifdef PARALLEL
        global_tot_walkers = tot_walkers
        tot_walkers = 0
        ! Read in walkers to spawning arrays.
        ! Restart file might have been produced with a different number of
        ! processors, thus need to hash walkers again to choose which
        ! processor to send each determinant to.
        ! Use the spawning arrays as scratch space.
        allocate(scratch_energies(spawned_walker_length), stat=ierr)
        call check_allocate('scratch_energies',spawned_walker_length,ierr)
        ! spawning_head_start gives the first slot in the spawning array for
        ! each processor.  Also want the last slot in the spawning array for
        ! each processor.
        forall (i=0:nprocs-2) spawn_max(i) = spawning_head(i+1) - 1
        spawn_max(nprocs-1) = spawned_walker_length
        iread = 1
        do
            ! read in a "block" of walkers.
            if (parent) then
                spawning_head = spawning_block_start
                do i = iread, global_tot_walkers
                    read (io,*) det, pop, energy
                    dest = modulo(murmurhash_bit_string(det, basis_length), nprocs)
                    spawning_head(dest) = spawning_head(dest) + 1
                    spawned_walkers(:basis_length, spawning_head(dest)) = det
                    spawned_walkers(basis_length+1, spawning_head(dest)) = pop
                    scratch_energies(spawning_head(dest)) = energy
                    ! Filled up spawning/scratch arrays?
                    if (any(spawning_head(:nprocs-1)-spawn_max == 0)) exit
                end do
                iread = i
                done = iread == global_tot_walkers + 1
            end if

            ! update the number of walkers on this processor from the number of
            ! walkers just read in.
            send_counts = spawning_head(:nprocs-1) - spawning_block_start(:nprocs-1)
            send_displacements = spawning_block_start(:nprocs-1)
            call mpi_scatter(send_counts, 1, mpi_integer, nread, 1, mpi_integer, root, mpi_comm_world, ierr)
            ! send walkers to their appropriate processor.
            call mpi_scatterv(scratch_energies, send_counts, send_displacements, mpi_preal, &
                              walker_energies(tot_walkers+1:), nread, mpi_preal, root,      &
                              mpi_comm_world, ierr)
            send_counts = send_counts*spawned_size
            send_displacements = send_displacements*spawned_size
            ! Easy to scatter into a different array.  Helpfully we already need
            ! spawned_walkers_recvd to be allocated for parallel calculations.
            ! :-)
            call mpi_scatterv(spawned_walkers, send_counts, send_displacements, mpi_det_integer, &
                              spawned_walkers_recvd, spawned_size*nread, mpi_det_integer, root, mpi_comm_world, ierr) 
            ! Transfer from spawned arrays to main walker arrays.
            do i = 1, nread
                walker_dets(:,i+tot_walkers) = spawned_walkers_recvd(:basis_length,i)
                walker_population(i+tot_walkers) = spawned_walkers_recvd(basis_length+1,i)
            end do
            tot_walkers = tot_walkers + nread

            call mpi_bcast(done, 1, mpi_logical, root, mpi_comm_world, ierr)
            if (done) exit
        end do
        deallocate(scratch_energies, stat=ierr)
        call check_deallocate('scratch_energies',ierr)
        ! Finally, need to broadcast the other information read in.
        call mpi_bcast(restart_version, 1, mpi_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(mc_cycles_done, 1, mpi_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(nparticles_old_restart, 1, mpi_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(shift, 1, mpi_preal, root, mpi_comm_world, ierr)
        call mpi_bcast(vary_shift, 1, mpi_logical, root, mpi_comm_world, ierr)
        call mpi_bcast(f0, basis_length, mpi_det_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(occ_list0, nel, mpi_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(D0_population, 1, mpi_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(H00, 1, mpi_preal, root, mpi_comm_world, ierr)
        ! Evaluate quantities based upon data read from restart file.
        if (nprocs > 1) then
            D0_proc = modulo(murmurhash_bit_string(f0, basis_length), nprocs)
        else
            D0_proc = iproc
        end if
#else
        do i = 1, tot_walkers
            read (io,*) walker_dets(:,i), walker_population(i), walker_energies(i)
        end do
#endif

        if (parent) close(io)

    end subroutine read_restart

end module fciqmc_restart
