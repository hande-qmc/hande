module fciqmc_restart

! Module for dumping out and restarting FCIQMC files.

#include "cdefs.h"

use const, only: p, i0, lint

implicit none

character(*), parameter :: restart_file_stem = 'restart'

! If negative, then this is set in the input and we need to read from
! a specific restart file.
integer :: read_restart_number = 0

! If negative, then this is set in the input and we need to write to
! a specific restart file.
integer :: write_restart_number = 0

! Dump restart files every n report_loops
integer :: write_restart_file_every_nreports = huge(0)

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
    module procedure write_out_int_i0_arr
    module procedure write_out_int_lint_arr
    module procedure write_out_real
    module procedure write_out_logical
    module procedure write_out_char_str
    module procedure write_out_int_arr_i0_int_arr_real_arr
end interface write_out

interface read_in
#if DET_SIZE != 32
    ! for non-32 bit integers, the i0 kind is distinct from the default kind;
    ! thus we need distinct functions to deal with them
    module procedure read_in_int
    module procedure read_in_int_arr
#endif
    module procedure read_in_int_i0
    module procedure read_in_int_i0_arr
    module procedure read_in_int_lint_arr
    module procedure read_in_real
    module procedure read_in_char_str
    module procedure read_in_logical
    ! For determinant information
    module procedure read_in_int_arr_i0_int_arr_real_arr
end interface read_in

contains

!--- write/read restart file ---

    subroutine dump_restart(nmc_cycles, nparticles_old)

        ! Write out the main walker list to file.
        use hfs_data, only: hf_shift
        use fciqmc_data, only: sampling_size, info_size, shift, occ_list0, tot_walkers, &
                               nparticles, walker_population, walker_dets, walker_data

        use parallel
        use utils, only: get_unique_filename, get_free_unit

        integer, intent(in) :: nmc_cycles
        integer(lint) :: nparticles_old(:)
        character(255) :: restart_file
        integer :: io
        integer, parameter :: restart_version = 2
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

            call write_out(io,'# restart version', binary_fmt_out)
            call write_out(io,restart_version, binary_fmt_out)
            call write_out(io,'# number of cycles', binary_fmt_out)
            call write_out(io,nmc_cycles, binary_fmt_out)
            call write_out(io,'# sampling size', binary_fmt_out)
            call write_out(io, sampling_size, binary_fmt_out)
            call write_out(io,'# info size', binary_fmt_out)
            call write_out(io, info_size, binary_fmt_out)
            call write_out(io,'# nparticles', binary_fmt_out)
            call write_out(io,nparticles_old, binary_fmt_out)
            call write_out(io,'# shift', binary_fmt_out)
            call write_out(io,shift, binary_fmt_out)
            call write_out(io,'# Hellmann--Feynman shift', binary_fmt_out)
            call write_out(io,hf_shift, binary_fmt_out)
            call write_out(io,'# reference determinant: orbital list', binary_fmt_out)
            call write_out(io,occ_list0, binary_fmt_out)
            call write_out(io,'# number of unique walkers', binary_fmt_out)
#ifdef PARALLEL
            call write_out(io,sum(nwalkers), binary_fmt_out)
            ! Write out walkers on parent processor to restart file.
            call write_out(io,'# walker info', binary_fmt_out)
            call write_walkers(io,tot_walkers, binary_fmt_out)

            ! The walker arrays on root are overwritten by data from other
            ! processors (so only root does i/o).  Write the root walkers to
            ! a scratch file so we can easily read them back in later.
            ! Best to store this in binary, not matter what the user input
            ! options are.
            scratch = get_free_unit()
            open(scratch,status="scratch",form='unformatted')
            call write_walkers(scratch,tot_walkers, .true.)

            ! Communicate with all other processors.
            do i = 1, nprocs-1
                ! Receive walker infor from all other processors.
                ! This overwrites the root processor's walkers
                call mpi_recv(walker_population, sampling_size*nwalkers(i), mpi_integer, i, comm_tag, mpi_comm_world, stat, ierr)
                call mpi_recv(walker_dets, basis_length*nwalkers(i), mpi_det_integer, i, comm_tag, mpi_comm_world, stat, ierr)
                call mpi_recv(walker_data, (sampling_size+info_size)*nwalkers(i), mpi_preal, i, comm_tag, mpi_comm_world, stat, ierr)
                ! Write out walkers from all other processors.
                call write_walkers(io, nwalkers(i), binary_fmt_out)
            end do

            ! We need to read back from scratch as we overwrote walker data on
            ! root.
            call flush(scratch)
            rewind(scratch)
            call read_walkers(scratch,tot_walkers, sampling_size, info_size, .true.)
            ! No longer need the scratchfile
            close(scratch)
#else
            call write_out(io,tot_walkers, binary_fmt_out)
            call write_out(io,'# walker info', binary_fmt_out)
            call write_walkers(io,tot_walkers, binary_fmt_out)
#endif
            close(io)

        else
#ifdef PARALLEL
            ! Send walker info to root processor.
            call mpi_send(walker_population, sampling_size*tot_walkers, mpi_integer, root, comm_tag, mpi_comm_world, ierr)
            call mpi_send(walker_dets, basis_length*tot_walkers, mpi_det_integer, root, comm_tag, mpi_comm_world, ierr)
            call mpi_send(walker_data, (sampling_size+info_size)*tot_walkers, mpi_preal, root, comm_tag, mpi_comm_world, ierr)
#endif
        end if

    end subroutine dump_restart

    subroutine read_restart()

        ! Read in the main walker list from file.

        use basis, only: basis_length
        use determinants, only: encode_det
        use errors, only: stop_all
        use fciqmc_data, only: sampling_size, info_size, shift, occ_list0, tot_walkers,    &
                               mc_cycles_done, nparticles, walker_dets, walker_population, &
                               walker_data, spawned_walkers, spawned_walkers_recvd,        &
                               spawning_head, spawned_size, spawning_block_start,          &
                               spawned_walker_length
        use hashing, only: murmurhash_bit_string
        use hfs_data, only: O00, hf_shift
        use proc_pointers, only: op0_ptr
        use system, only: nel

        use checking, only: check_allocate, check_deallocate
        use parallel
        use utils, only: get_unique_filename, get_free_unit

        character(255) :: restart_file, junk
        integer :: io, i
        logical :: exists
        integer :: restart_version
        integer :: restart_sampling_size, restart_info_size
        integer(i0) :: det(basis_length)
#ifdef PARALLEL
        integer :: global_tot_walkers, pop(sampling_size), iread, nread, ierr, dest
        real(p), allocatable :: scratch_data(:,:)
        integer :: send_counts(0:nprocs-1), send_displacements(0:nprocs-1)
        real(p) :: tmp_data(sampling_size+info_size)
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

            ! Read in data.  Format is a description line followed by a line
            ! containing the corresponding data item.  Binary-formatted restart
            ! files contain only the data items.
            call read_in(io,junk,binary_fmt_in)
            call read_in(io,restart_version,binary_fmt_in)
            call read_in(io,junk,binary_fmt_in)
            call read_in(io,mc_cycles_done,binary_fmt_in)
            call read_in(io,junk, binary_fmt_in)
            call read_in(io, restart_sampling_size, binary_fmt_in)
            call read_in(io,junk, binary_fmt_in)
            call read_in(io, restart_info_size, binary_fmt_in)
            call read_in(io,junk,binary_fmt_in)
            ! In case sampling_size > restart_sampling_size
            nparticles = 0
            call read_in(io,nparticles(:restart_sampling_size),binary_fmt_in)
            call read_in(io,junk,binary_fmt_in)
            call read_in(io,shift,binary_fmt_in)
            call read_in(io,junk,binary_fmt_in)
            call read_in(io,hf_shift, binary_fmt_in)
            call read_in(io,junk,binary_fmt_in)
            call read_in(io,occ_list0,binary_fmt_in)
            call read_in(io,junk,binary_fmt_in)
            call read_in(io,tot_walkers,binary_fmt_in)
            call read_in(io,junk,binary_fmt_in)

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
        allocate(scratch_data(sampling_size+info_size,spawned_walker_length), stat=ierr)
        call check_allocate('scratch_data',size(scratch_data),ierr)
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
                    call read_in(io,det,pop,tmp_data,binary_fmt_in)
                    dest = modulo(murmurhash_bit_string(det, basis_length), nprocs)
                    spawning_head(dest) = spawning_head(dest) + 1
                    ! zero spawned array in case some elements were not set
                    ! (e.g. restart file from standard FCIQMC calculation but
                    ! now doing Hellmann--Feynman sampling)
                    spawned_walkers(:, spawning_head(dest)) = 0
                    scratch_data(:, spawning_head(dest)) = 0.0_p
                    spawned_walkers(:basis_length, spawning_head(dest)) = det
                    spawned_walkers(basis_length+1:basis_length+restart_sampling_size, spawning_head(dest)) = pop
                    scratch_data(:restart_sampling_size+restart_info_size,spawning_head(dest)) = tmp_data
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
            call mpi_scatterv(scratch_data, (sampling_size+info_size)*send_counts,               &
                                 (sampling_size+info_size)*send_displacements, mpi_preal,        &
                                 walker_data(:,tot_walkers+1:), (sampling_size+info_size)*nread, &
                                 mpi_preal, root, mpi_comm_world, ierr)
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
                walker_population(:,i+tot_walkers) = spawned_walkers_recvd(basis_length+1:basis_length+sampling_size,i)
            end do
            tot_walkers = tot_walkers + nread

            call mpi_bcast(done, 1, mpi_logical, root, mpi_comm_world, ierr)
            if (done) exit
        end do
        deallocate(scratch_data, stat=ierr)
        call check_deallocate('scratch_data',ierr)
        ! Finally, need to broadcast the other information read in.
        call mpi_bcast(restart_version, 1, mpi_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(mc_cycles_done, 1, mpi_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(restart_sampling_size, 1, mpi_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(restart_info_size, 1, mpi_integer, root, mpi_comm_world, ierr)
        call mpi_bcast(nparticles, size(nparticles), mpi_integer8, root, mpi_comm_world, ierr)
        call mpi_bcast(shift, 1, mpi_preal, root, mpi_comm_world, ierr)
        call mpi_bcast(hf_shift, 1, mpi_preal, root, mpi_comm_world, ierr)
        call mpi_bcast(occ_list0, nel, mpi_integer, root, mpi_comm_world, ierr)
#else
        call read_walkers(io, tot_walkers, restart_sampling_size, restart_info_size, binary_fmt_in)
#endif

        ! Fill in any gaps in walker_data
        if (info_size /= restart_info_size) then
            ! Seriously wrong...
            ! Doing (e.g.) a standard Heisenberg calculation based upon
            ! a restart file generated using a Neel trial state!
            call stop_all('read_restart','Info size in restart file does not match calculation info_size.')
        end if
        select case(sampling_size-restart_sampling_size)
        case(0)
            ! Restarting from same kind of calculation.  Have all the data
            ! required.
        case(1)
            ! HFS calculation started from a standard FCIQMC calculation.  Need
            ! to fill in <D|O|D> - <D0|O|D0>.
            call encode_det(occ_list0, det)
            O00 = op0_ptr(det)
            do i = 1, tot_walkers
                walker_data(2,i) = op0_ptr(walker_dets(:,i)) - O00
            end do
        case default
            ! No idea...
            call stop_all('read_restart','Sampling size in restart file not compatible with calculation.')
        end select

        if (parent) close(io)

    end subroutine read_restart

!--- Helper routines: read/write walker lists ---

    subroutine write_walkers(iunit, my_nwalkers, binary)

        ! Read in a list of walker data.

        ! In:
        !    iunit: file unit to read from
        !    my_nwalkers: number of walkers to read in
        !    binary: true if file is in binary format.

        use fciqmc_data, only: walker_population, walker_data, walker_dets

        integer, intent(in) :: my_nwalkers, iunit
        logical, intent(in) :: binary
        integer :: iwalker

        do iwalker = 1, my_nwalkers
            call write_out(iunit,walker_dets(:,iwalker),&
                           walker_population(:,iwalker),&
                           walker_data(:,iwalker), binary)
        end do

    end subroutine write_walkers

    subroutine read_walkers(iunit, my_nwalkers, restart_sampling_size, restart_info_size, binary)

        ! Read in a list of walker data.

        ! In:
        !    iunit: file unit to read from
        !    my_nwalkers: number of walkers to read in
        !    restart_sampling_size: sampling_size (see fciqmc_data) of QMC
        !        calculation which produced the restart file.
        !    restart_info_size: info_size (see fciqmc_data) of QMC
        !        calculation which produced the restart file.
        !    binary: true if file is in binary format.

        use fciqmc_data, only: walker_population, walker_data, walker_dets

        integer, intent(in) :: my_nwalkers, iunit, restart_sampling_size, restart_info_size
        logical, intent(in) :: binary
        integer :: iwalker

        do iwalker = 1, my_nwalkers
            walker_population(:,iwalker) = 0
            walker_data(:,iwalker) = 0
            call read_in(iunit,walker_dets(:,iwalker),&
                         walker_population(:restart_sampling_size,iwalker),&
                         walker_data(:restart_sampling_size+restart_info_size,iwalker), binary)
        end do

    end subroutine read_walkers

!--- Utility functions: writing out data ---

#if DET_SIZE != 32

    subroutine write_out_int(iunit, a, binary, fmt_string)

        ! Write an integer.

        ! In:
        !    iunit: file unit to write to.
        !    a: integer to write out.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the output of a using this string.

        integer, intent(in) :: a, iunit
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           write (iunit) a
        else if(present(fmt_string)) then
           write (iunit,fmt=fmt_string) a
        else
           write (iunit,*) a
        end if

    end subroutine write_out_int

    subroutine write_out_int_arr(iunit, a, binary, fmt_string)

        ! Write an integer array.

        ! In:
        !    iunit: file unit to write to.
        !    a: integer array to write out.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the output of a using this string.

        integer, intent(in) :: a(:), iunit
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           write (iunit) a
        else if(present(fmt_string)) then
           write (iunit,fmt=fmt_string) a
        else
           write (iunit,*) a
        end if

    end subroutine write_out_int_arr

#endif

    subroutine write_out_int_i0(iunit, a, binary, fmt_string)

        ! Write an integer of kind i0.

        ! In:
        !    iunit: file unit to write to.
        !    a: integer of kind i0 to write out.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the output of a using this string.

        integer, intent(in) :: iunit
        integer(i0), intent(in) :: a
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           write (iunit) a
        else if(present(fmt_string)) then
           write (iunit,fmt=fmt_string) a
        else
           write (iunit,*) a
        end if

    end subroutine write_out_int_i0

    subroutine write_out_int_i0_arr(iunit, a, binary, fmt_string)

        ! Write an integer array of kind i0.

        ! In:
        !    iunit: file unit to write to.
        !    a: integer array of kind i0 to write out.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the output of a using this string.

        integer, intent(in) :: iunit
        integer(i0), intent(in) :: a(:)
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           write (iunit) a
        else if(present(fmt_string)) then
           write (iunit,fmt=fmt_string) a
        else
           write (iunit,*) a
        end if

    end subroutine write_out_int_i0_arr

    subroutine write_out_int_lint_arr(iunit, a, binary, fmt_string)

        ! Write an integer array of kind lint.

        ! In:
        !    iunit: file unit to write to.
        !    a: integer array of kind lint to write out.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the output of a using this string.

        integer, intent(in) :: iunit
        integer(lint), intent(in) :: a(:)
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           write (iunit) a
        else if(present(fmt_string)) then
           write (iunit,fmt=fmt_string) a
        else
           write (iunit,*) a
        end if

    end subroutine write_out_int_lint_arr

    subroutine write_out_real(iunit, a, binary, fmt_string)

        ! Write a real of kind p.

        ! In:
        !    iunit: file unit to write to.
        !    a: real of kind p to write out.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the output of a using this string.

        integer, intent(in) :: iunit
        real(p), intent(in) :: a
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           write (iunit) a
        else if(present(fmt_string)) then
           write (iunit,fmt=fmt_string) a
        else
           write (iunit,*) a
        end if

    end subroutine write_out_real

    subroutine write_out_logical(iunit, a, binary, fmt_string)

        ! Write a logical.

        ! In:
        !    iunit: file unit to write to.
        !    a: logical to write out.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the output of a using this string.

        integer, intent(in) :: iunit
        logical, intent(in) :: a
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           write (iunit) a
        else if(present(fmt_string)) then
           write (iunit,fmt=fmt_string) a
        else
           write (iunit,*) a
        end if

    end subroutine write_out_logical

    subroutine write_out_char_str(iunit, a, suppress_output, fmt_string)

        ! Write a character string.

        ! In:
        !    iunit: file unit to write to.
        !    a: character string to write out.
        !    suppress_output: if true (e.g. when iunit is a binary file) then no
        !        output is printed out.
        !    fmt_string (optional): if present and binary is false, then format
        !        the output of a using this string.

        character(*), intent(in) :: a
        integer, intent(in) :: iunit
        logical, intent(in) :: suppress_output
        character(*), intent(in), optional :: fmt_string

        ! Tend to not want character data for the binary format output file.
        if (.not. suppress_output) then
            if (present(fmt_string)) then
               write (iunit,fmt=fmt_string) a
            else
               write (iunit,*) a
            end if
        end if

    end subroutine write_out_char_str

    subroutine write_out_int_arr_i0_int_arr_real_arr(iunit, i0arr, intarr, rarr, binary, fmt_string)

        ! Write an i0 array, integer array and real array out on 1 line.

        ! We need specific routines for writing more than 1 entry per line
        ! due to lack of proper generic programming in Fortran 90.
        ! For binary format these procedures will do the same task as multiple
        ! calls to the above procedures, however for ASCII output they are
        ! necessary to have output all on one line

        ! In:
        !    iunit: file unit to write to.
        !    i0arr: integer array of kind i0 to write out.
        !    intarr: integer array to write out.
        !    rarr: real array of kind of to write out.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the output of a using this string.

        integer, intent(in) :: iunit
        integer(i0), intent(in) :: i0arr(:)
        integer, intent(in) :: intarr(:)
        real(p), intent(in) :: rarr(:)
        logical, intent(in) :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           write (iunit) i0arr, intarr, rarr
        else if (present(fmt_string)) then
           write (iunit,fmt=fmt_string) i0arr, intarr, rarr
        else
           write (iunit,*) i0arr, intarr, rarr
        end if

    end subroutine write_out_int_arr_i0_int_arr_real_arr

!--- Utility functions: reading data ---

#if DET_SIZE != 32

    subroutine read_in_int(iunit, a, binary, fmt_string)

        ! Read an integer.

        ! In:
        !    iunit: file unit to read to.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the input of a using this string.
        ! Out:
        !    a: integer to read in.

        integer, intent(in) :: iunit
        integer, intent(out) :: a
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           read (iunit) a
        else if(present(fmt_string)) then
           read (iunit,fmt=fmt_string) a
        else
           read (iunit,*) a
        end if

    end subroutine read_in_int

    subroutine read_in_int_arr(iunit, a, binary, fmt_string)

        ! Read an integer array.

        ! In:
        !    iunit: file unit to read to.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the input of a using this string.
        ! Out:
        !    a: integer array to read in.

        integer, intent(in) :: iunit
        integer, intent(out) :: a(:)
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           read (iunit) a
        else if(present(fmt_string)) then
           read (iunit,fmt=fmt_string) a
        else
           read (iunit,*) a
        end if

    end subroutine read_in_int_arr

#endif

    subroutine read_in_int_i0(iunit, a, binary, fmt_string)

        ! Read an integer of kind i0.

        ! In:
        !    iunit: file unit to read to.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the input of a using this string.
        ! Out:
        !    a: integer of kind i0 to read in.

        integer, intent(in) :: iunit
        integer(i0), intent(out) :: a
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           read (iunit) a
        else if(present(fmt_string)) then
           read (iunit,fmt=fmt_string) a
        else
           read (iunit,*) a
        end if

    end subroutine read_in_int_i0

    subroutine read_in_int_i0_arr(iunit, a, binary, fmt_string)

        ! Read an integer array of kind i0.

        ! In:
        !    iunit: file unit to read to.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the input of a using this string.
        ! Out:
        !    a: integer array of kind i0 to read in.

        integer, intent(in) :: iunit
        integer(i0), intent(out) :: a(:)
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           read (iunit) a
        else if(present(fmt_string)) then
           read (iunit,fmt=fmt_string) a
        else
           read (iunit,*) a
        end if

    end subroutine read_in_int_i0_arr

    subroutine read_in_int_lint_arr(iunit, a, binary, fmt_string)

        ! Read an integer array of kind lint.

        ! In:
        !    iunit: file unit to read to.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the input of a using this string.
        ! Out:
        !    a: integer array of kind lint to read in.

        integer, intent(in) :: iunit
        integer(lint), intent(out) :: a(:)
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           read (iunit) a
        else if(present(fmt_string)) then
           read (iunit,fmt=fmt_string) a
        else
           read (iunit,*) a
        end if

    end subroutine read_in_int_lint_arr

    subroutine read_in_real(iunit, a, binary, fmt_string)

        ! Read a real of kind p.

        ! In:
        !    iunit: file unit to read to.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the input of a using this string.
        ! Out:
        !    a: real of kind p to read in.

        integer, intent(in) :: iunit
        real(p), intent(out) :: a
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           read (iunit) a
        else if(present(fmt_string)) then
           read (iunit,fmt=fmt_string) a
        else
           read (iunit,*) a
        end if

    end subroutine read_in_real

    subroutine read_in_logical(iunit, a, binary, fmt_string)

        ! Read a logical.

        ! In:
        !    iunit: file unit to read to.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the input of a using this string.
        ! Out:
        !    a: logical to read in.

        integer, intent(in) :: iunit
        logical, intent(out) :: a
        logical :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           read (iunit) a
        else if(present(fmt_string)) then
           read (iunit,fmt=fmt_string) a
        else
           read (iunit,*) a
        end if

    end subroutine read_in_logical

    subroutine read_in_char_str(iunit, a, suppress_input, fmt_string)

        ! Read a character string.

        ! In:
        !    iunit: file unit to read to.
        !    a: character string to read in.
        !    suppress_input: if true (e.g. when iunit is a binary file) then no
        !        input is read in.
        !    fmt_string (optional): if present and binary is false, then format
        !        the input of a using this string.

        integer, intent(in) :: iunit
        character(*), intent(out) :: a
        logical, intent(in) :: suppress_input
        character(*), intent(in), optional :: fmt_string

        ! Tend to not want character data for the binary format output file.
        if (.not. suppress_input) then
            if (present(fmt_string)) then
               read (iunit,fmt=fmt_string) a
            else
               read (iunit,*) a
            end if
        end if

    end subroutine read_in_char_str

    subroutine read_in_int_arr_i0_int_arr_real_arr(iunit, i0arr, intarr, rarr, binary, fmt_string)

        ! Read an i0 array, integer array and real array out on 1 line.
        ! We need specific routines for writing more than 1 entry per line
        ! due to lack of proper generic programming in Fortran 90.
        ! For binary format these procedures will do the same task as multiple
        ! calls to the above procedures, however for ASCII output they are
        ! necessary to have output all on one line
        ! In:
        !    iunit: file unit to write to.
        !    binary: if true, iunit is a binary file and no format statement is used.
        !    fmt_string (optional): if present and binary is false, then format
        !        the input of a using this string.
        ! Out:
        !    i0arr: integer array of kind i0 to read in.
        !    intarr: integer array to read in.
        !    rarr: real array of kind of to read in.

        integer, intent(in) :: iunit
        integer(i0), intent(out) :: i0arr(:)
        integer, intent(out) :: intarr(:)
        real(p), intent(out) :: rarr(:)
        logical, intent(in) :: binary
        character(*), intent(in), optional :: fmt_string

        if (binary) then
           read (iunit) i0arr, intarr, rarr
        else if (present(fmt_string)) then
           read (iunit,fmt=fmt_string) i0arr, intarr, rarr
        else
           read (iunit,*) i0arr, intarr, rarr
        end if

    end subroutine read_in_int_arr_i0_int_arr_real_arr

end module fciqmc_restart
