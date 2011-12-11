module fciqmc_restart

! Module for dumping out and restarting FCIQMC files.

#include "cdefs.h"

use parallel
use utils, only: get_unique_filename, get_free_unit
use const, only: p, i0
use fciqmc_data
use basis, only: basis_length
use system, only: nel

implicit none

character(*), parameter :: restart_file_stem = 'restart'

! If negative, then this is set in the input and we need to read from
! a specific restart file.
integer :: read_restart_number = 0

! If negative, then this is set in the input and we need to write to
! a specific restart file.
integer :: write_restart_number = 0

! specifies if the restart file (to read --> in, to write --> out)
! is in binary (.true.) or ASCII (.false.) format;
! The latter requires substantially more space ( 1 byte per digit of output
! as opposed to 1 integer per integer of output!) but is human-readable
logical :: binary_fmt_in = .true., binary_fmt_out = .true. 

! An attempt to do generic programming in Fortran: These functions all print
! a variable of a specific type to a given unit, either in binary or ascii
! format
interface write_out
#if DET_SIZE != 32 
    ! for non-32 bit integers, the i0 kind is distinct from the default kind;
    ! thus we need distinct functions to deal with them
    module procedure write_out_int
    module procedure write_out_int_arr
#endif
    module procedure write_out_int_i0
    module procedure write_out_int_arr_i0
    module procedure write_out_float
    module procedure write_out_char
    module procedure write_out_logical
    module procedure write_out_i0arr_i_r
    module procedure write_out_i_r_l
    module procedure write_out_i0arr_iarr_2r
end interface write_out

interface read_in 
#if DET_SIZE != 32
    module procedure read_in_int
    module procedure read_in_int_arr
#endif
    module procedure read_in_int_i0
    module procedure read_in_int_arr_i0
    module procedure read_in_float
    module procedure read_in_char
    module procedure read_in_logical
    module procedure read_in_i0arr_i_r
    module procedure read_in_i_r_l
    module procedure read_in_i0arr_iarr_2r
end interface read_in

contains 

    subroutine dump_restart(nmc_cycles, nparticles_old)

        ! Write out the main walker list to file.

        integer, intent(in) :: nmc_cycles, nparticles_old
        character(255) :: restart_file
        integer :: io
        integer, parameter :: restart_version = 1
#ifdef PARALLEL
        integer :: nwalkers(0:nprocs-1), ierr, stat(MPI_STATUS_SIZE), i, scratch
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

            write (6,'(1X,"#",1X,a23,1X,a,a1)') 'Writing restart file to',trim(restart_file),'.'

            if (binary_fmt_out) then
                open(io,file=restart_file,form='unformatted')
            else
                open(io, file=restart_file)
            end if

            call write_out(io,'# restart version')
            call write_out(io,restart_version) 
            call write_out(io,'# number of cycles')
            call write_out(io,nmc_cycles)
            call write_out(io,'#shift')
            call write_out(io,nparticles_old,shift,vary_shift)
            call write_out(io,'#reference determinant')
            call write_out(io,f0,occ_list0,D0_population,H00)
            call write_out(io,'# number of unique walkers')
#ifdef PARALLEL
            call write_out(io,sum(nwalkers))
            ! Write out walkers on parent processor to restart file.
            call write_out(io,'# walker info')
            call write_walkers(io,tot_walkers)
            
            ! if writing in binary, there is no character marker telling us
            ! where the root processor's walkers are stored - thus use scratch
            ! so that we can get root's walkers back
            if (binary_fmt_out) then  
                scratch = get_free_unit()
                open(scratch,status="scratch",form='unformatted')
                call write_walkers(scratch,tot_walkers)
            end if

            ! Communicate with all other processors.
            do i = 1, nprocs-1
                ! Receive walker infor from all other processors.
                ! This overwrites the root processor's walkers
                call mpi_recv(walker_population, nwalkers(i), mpi_integer, i, comm_tag, mpi_comm_world, stat, ierr)
                call mpi_recv(walker_dets, nwalkers(i), mpi_det_integer, i, comm_tag, mpi_comm_world, stat, ierr)
                call mpi_recv(walker_data, nwalkers(i), mpi_preal, i, comm_tag, mpi_comm_world, stat, ierr)
                ! Write out walkers from all other processors.
                call write_walkers(io, nwalkers(i)) 
            end do

            ! we need to read back from scratch if in binary format
            if (binary_fmt_out) then
                call flush(scratch)
                rewind(scratch)
                ! best to preserve the binary_fmt state in case of 
                ! alteration of program mechaincs later
                if (binary_fmt_in) then
                    call read_walkers(scratch,tot_walkers)
                else
                    binary_fmt_in = .true.
                    call read_walkers(scratch,tot_walkers)
                    binary_fmt_in = .false.
                end if
                !no longer need the scratchfile
                close(scratch)
            else
                ! we can read from the input file and no need for scratch
                ! Read "self" info back in.
                call flush(io)
                rewind(io)
                do
                    ! Read restart file until we've found the start of the
                    ! walker information.
                    call read_in(io,junk,'(a255)')
                    if (index(junk,'walker info') /= 0) exit
                end do
                ! The next tot_walkers lines contain the walker info that came
                ! from the root processor.
                call read_walkers(io,tot_walkers)
            end if
#else
            call write_out(io,tot_walkers)
            call write_out(io,'# walker info')
            call write_walkers(io,tot_walkers)
#endif
            close(io)

        else
#ifdef PARALLEL
            ! Send walker info to root processor.
            call mpi_send(walker_population, tot_walkers, mpi_integer, root, comm_tag, mpi_comm_world, ierr)
            call mpi_send(walker_dets, tot_walkers, mpi_det_integer, root, comm_tag, mpi_comm_world, ierr)
            call mpi_send(walker_data, tot_walkers, mpi_preal, root, comm_tag, mpi_comm_world, ierr)
#endif
        end if

    end subroutine dump_restart

    subroutine read_restart()

        ! Read in the main walker list from file.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use hashing, only: murmurhash_bit_string

        character(255) :: restart_file, junk
        integer :: io
        logical :: exists
        integer :: restart_version
#ifdef PARALLEL
        integer :: global_tot_walkers, pop, iread, nread, ierr, dest, i
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
            
            if (binary_fmt_in) then
                open(io,file=restart_file,form='unformatted')
            else
                open(io, file=restart_file)
            end if
            call read_in(io,junk)
            call read_in(io,restart_version)
            call read_in(io,junk)
            call read_in(io,mc_cycles_done)
            call read_in(io,junk)
            call read_in(io,nparticles_old_restart,shift,vary_shift)
            call read_in(io,junk)
            call read_in(io,f0,occ_list0,D0_population,H00)
            call read_in(io,junk)
            call read_in(io,tot_walkers)
            call read_in(io,junk)
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
        iread = 0
        do
            ! read in a "block" of walkers.
            if (parent) then
                spawning_head = spawning_block_start
                do i = iread+1, global_tot_walkers
                    call read_in(io,det,pop,energy)
                    dest = modulo(murmurhash_bit_string(det, basis_length), nprocs)
                    spawning_head(dest) = spawning_head(dest) + 1
                    spawned_walkers(:basis_length, spawning_head(dest)) = det
                    spawned_walkers(spawned_pop, spawning_head(dest)) = pop
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
            call mpi_scatterv(scratch_energies, send_counts, send_displacements, mpi_preal,   &
                              walker_data(:,tot_walkers+1:), nread, mpi_preal, root,      &
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
                walker_population(1,i+tot_walkers) = spawned_walkers_recvd(basis_length+1,i)
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
        call mpi_bcast(D0_population, 1, mpi_preal, root, mpi_comm_world, ierr)
        call mpi_bcast(H00, 1, mpi_preal, root, mpi_comm_world, ierr)
        ! Evaluate quantities based upon data read from restart file.
        if (nprocs > 1) then
            D0_proc = modulo(murmurhash_bit_string(f0, basis_length), nprocs)
        else
            D0_proc = iproc
        end if
#else
        call read_walkers(io, tot_walkers)
#endif

        if (parent) close(io)

    end subroutine read_restart

    subroutine write_walkers(iunit, my_nwalkers)

        integer, intent(in) :: my_nwalkers, iunit
        integer :: iwalker

        do iwalker = 1, my_nwalkers
            call write_out(iunit,walker_dets(:,iwalker),&
                           walker_population(1,iwalker),&
                           walker_data(1,iwalker))
        end do

    end subroutine write_walkers

    subroutine read_walkers(iunit, my_nwalkers)

        integer, intent(in) :: my_nwalkers, iunit
        integer :: iwalker

        do iwalker = 1, my_nwalkers
            call read_in(iunit,walker_dets(:,iwalker),&
                         walker_population(1,iwalker),&
                         walker_data(1,iwalker))
        end do

    end subroutine read_walkers
    
#if DET_SIZE != 32
    
    subroutine write_out_int(wunit, a, fmt_string)
        
        ! write one integer to wunit with optional format string

        integer, intent(in) :: a, wunit
        character(*), intent(in), optional :: fmt_string
        
        if (binary_fmt_out) then
            write(wunit) a
        else if(present(fmt_string)) then
            write(wunit,fmt=fmt_string) a
        else 
            write(wunit,*) a
        end if

    end subroutine write_out_int

    subroutine write_out_int_arr(wunit, a, fmt_string)
    
        !print out an array of integers to wunit with optional format string

        integer, intent(in) :: wunit, a(:)
        character(*), intent(in), optional :: fmt_string
        
        if (binary_fmt_out) then
            write(wunit) a
        else if(present(fmt_string)) then
            write(wunit,fmt=fmt_string) a
        else 
            write(wunit,*) a
        end if
    end subroutine write_out_int_arr

#endif     
    
    subroutine write_out_int_i0(wunit, a, fmt_string)
        
        ! for non 32-bit integers, the i0 type is distinct and thus needs
        ! its own write_out procedure

        integer(i0), intent(in) :: a
        integer, intent(in) :: wunit
        character(*), intent(in), optional :: fmt_string
        
        if (binary_fmt_out) then
            write(wunit) a
        else if(present(fmt_string)) then
            write(wunit,fmt=fmt_string) a
        else 
            write(wunit,*) a
        end if

    end subroutine write_out_int_i0

    subroutine write_out_int_arr_i0(wunit, a, fmt_string)
        
        ! the i0 integer array needs its own procedure for non 32-bit 
        ! integers
     
        integer, intent(in) ::  wunit
        integer(i0), intent(in) :: a(:)
        character(*), intent(in), optional :: fmt_string
        
        if (binary_fmt_out) then
            write(wunit) a
        else if(present(fmt_string)) then
            write(wunit,fmt=fmt_string) a
        else 
            write(wunit,*) a
        end if

    end subroutine write_out_int_arr_i0

    subroutine write_out_float(wunit, a,fmt_string)

        ! write a p-precision floating point number to wunit

        real(p), intent(in) :: a
        integer, intent(in) :: wunit
        character(*), intent(in), optional :: fmt_string
        
        if (binary_fmt_out) then
            write(wunit) a
        else if(present(fmt_string)) then
            write(wunit,fmt=fmt_string) a
        else 
            write(wunit,*) a
        end if
    end subroutine write_out_float

    subroutine write_out_char(wunit, a, fmt_string)

        ! Write out a character to wunit in the optional format given
        ! by fmt_string. For binary output we do not want character output
        ! as they only serve as comments for the human reader, and thus if
        ! binary_fmt_out is set, this procedure does nothing

        character(*), intent(in) :: a
        integer, intent(in) :: wunit
        character(*), intent(in), optional :: fmt_string
        
        ! no character data for the binary format output file
        if (.not. binary_fmt_out) then
            if (present(fmt_string)) then
                write(wunit,fmt=fmt_string) a
            else
                write(wunit,*) a
            end if
        end if
    end subroutine write_out_char

    subroutine write_out_logical(wunit, a, fmt_string)
        
        ! write a logical variable to wunit with an optional format

        logical, intent(in):: a
        integer, intent(in) :: wunit
        character(*), intent(in),optional :: fmt_string

        if (binary_fmt_out) then
            write(wunit) a
        else if (present(fmt_string)) then
            write(wunit,fmt=fmt_string) a
        else
            write(wunit,*) a
        end if
    end subroutine write_out_logical

    ! we need specific routines for writing more than 1 entry per line
    ! due to lack of proper generic programming in Fortran 90.
    ! For binary format these procedures will do the same task as multiple
    ! calls to the above procedures, however for ASCII output they are 
    ! necessary to have output all on one line

    subroutine write_out_i0arr_i_r(wunit, i0arr, i, r, fmt_string)
        
        ! Write an i0 array, integer and real out on 1 line in optional format 

        integer, intent(in) :: i, wunit
        integer(i0), intent(in) :: i0arr(:)
        real(p), intent(in) :: r
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_out) then
            write(wunit) i0arr, i, r
        else if (present(fmt_string)) then
            write(wunit,fmt=fmt_string) i0arr, i, r
        else
            write(wunit,*) i0arr, i, r
        end if
    end subroutine write_out_i0arr_i_r

    subroutine write_out_i_r_l(wunit, i, r, l, fmt_string)
        
        ! write out an integer, real(p) and logical all on 1 line 
        ! in optional format
        
        integer, intent(in) :: i, wunit
        real(p), intent(in) :: r
        logical, intent(in) :: l
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_out) then
            write(wunit) i, r, l
        else if (present(fmt_string)) then
            write(wunit,fmt=fmt_string) i, r, l
        else
            write(wunit,*) i, r, l
        end if
    end subroutine write_out_i_r_l

    subroutine write_out_i0arr_iarr_2r(wunit, i0arr, iarr, r1, r2, fmt_string)
        
        ! write out i0 array, integer array and 2 real
        ! variables all on 1 line
        

        integer, intent(in) :: wunit, iarr(:)
        integer(i0), intent(in) :: i0arr(:)
        real(p), intent(in) :: r1, r2
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_out) then
            write(wunit) i0arr, iarr, r1, r2
        else if (present(fmt_string)) then
            write(wunit,fmt=fmt_string) i0arr, iarr, r1, r2
        else
            write(wunit,*) i0arr, iarr, r1, r2
        end if
    end subroutine write_out_i0arr_iarr_2r

#if DET_SIZE != 32

    subroutine read_in_int(runit, a, fmt_string)
        
        ! read in a single integer from runit

        integer, intent(out) :: a
        integer, intent(in) :: runit
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_in) then
            read(runit) a 
        else if (present(fmt_string)) then
            read(runit,fmt=fmt_string) a
        else
            read(runit,*) a
        end if
    end subroutine read_in_int

    subroutine read_in_int_arr(runit, a, fmt_string)

        ! read in an integer array from runit 
        
        integer, intent(in) :: runit
        integer, intent(out) :: a(:)
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_in) then
            read(runit) a
        else if (present(fmt_string)) then
            read(runit,fmt=fmt_string) a
        else
            read(runit,*) a
        end if
    end subroutine read_in_int_arr

#endif

    subroutine read_in_int_i0(runit, a, fmt_string)
        
        ! Read in a single integer of i0 kind from runit

        integer(i0), intent(out) :: a
        integer, intent(in) :: runit
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_in) then
            read(runit) a 
        else if (present(fmt_string)) then
            read(runit,fmt=fmt_string) a
        else
            read(runit,*) a
        end if
    end subroutine read_in_int_i0

    subroutine read_in_int_arr_i0(runit, a, fmt_string)

        ! Read in an array of i0-kind integers from runit
        
        integer, intent(in) :: runit
        integer(i0), intent(out) :: a(:)
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_in) then
            read(runit) a
        else if (present(fmt_string)) then
            read(runit,fmt=fmt_string) a
        else
            read(runit,*) a
        end if
    end subroutine read_in_int_arr_i0

    subroutine read_in_float(runit, a, fmt_string)

        ! Read in a p-kind real from runit

        real(p), intent(out) :: a
        integer, intent(in) :: runit
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_in) then
            read(runit) a 
        else if (present(fmt_string)) then
            read(runit,fmt=fmt_string) a
        else
            read(runit,*) a
        end if
    end subroutine read_in_float

    subroutine read_in_char(runit, a, fmt_string)

        ! Read in a string with optional format from runit. 
        ! As there is no character data in the binary restart files, 
        ! if binary_fmt_in is set, this procedure does nothing

        character(*), intent(out) :: a
        integer, intent(in) :: runit
        character(*), intent(in), optional :: fmt_string

        if (.not. binary_fmt_in) then
            if (present(fmt_string)) then
                read(runit,fmt=fmt_string) a
            else
                read(runit,*) a
            end if
        end if
    end subroutine read_in_char

    subroutine read_in_logical(runit, a, fmt_string)
        
        ! Read in a single logical variable from runit 

        integer, intent(in) :: runit
        logical, intent(out) :: a
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_in) then
            read(runit) a
        else if (present(fmt_string)) then
            read(runit,fmt=fmt_string) a
        else
            read(runit,*) a
        end if
    end subroutine read_in_logical

    subroutine read_in_i0arr_i_r(runit, i0arr, i, r, fmt_string)

        ! Read in an array of i0-kind integers, an integer and a real 
        ! from runit. If we are reading an ASCII file then this is all
        ! read from a single line

        integer, intent(in) :: runit
        integer, intent(out) :: i
        integer(i0), intent(out) :: i0arr(:)
        real(p), intent(out) :: r
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_in) then
            read(runit) i0arr, i, r
        else if (present(fmt_string)) then
            read(runit,fmt=fmt_string) i0arr, i, r
        else
            read(runit,*) i0arr, i, r
        end if
    end subroutine read_in_i0arr_i_r

    subroutine read_in_i_r_l(runit, i, r, l, fmt_string)
        
        ! Read in an integer, a real and a logial variable from runit

        integer, intent(in) :: runit
        integer, intent(out) :: i
        real(p), intent(out) :: r
        logical, intent(out) :: l
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_in) then
            read(runit) i, r, l
        else if (present(fmt_string)) then
            read(runit,fmt=fmt_string) i, r, l
        else
            read(runit,*) i, r, l
        end if
    end subroutine read_in_i_r_l

    subroutine read_in_i0arr_iarr_2r(runit, i0arr, iarr, r1, r2,fmt_string)
       
        ! Read in an array of i0-kind integers, an array of integers 
        ! and 2 reals/ If we are reading an ASCII format file then this is
        ! all read from a single line

        integer, intent(in) :: runit
        integer(i0), intent(out) :: i0arr(:)
        integer, intent(out) :: iarr(:)
        real(p), intent(out) :: r1, r2
        character(*), intent(in), optional :: fmt_string

        if (binary_fmt_in) then
            read(runit) i0arr, iarr, r1, r2
        else if (present(fmt_string)) then
            read(runit,fmt=fmt_string) i0arr, iarr, r1, r2
        else
            read(runit,*) i0arr, iarr, r1, r2
        end if
    end subroutine read_in_i0arr_iarr_2r

end module fciqmc_restart
