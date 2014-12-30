module restart_hdf5

    ! Restart functionality based on the HDF5 library.  Note: this is only
    ! for QMC (ie FCIQMC, DMQMC or CCMC) calculations).
    ! Due to the use of HDF5, the format is pretty flexible (i.e. we can easily add or
    ! remove data items, though changing the data structure of the existing output is more
    ! problematic from a backward-compatibility viewpoint.

    ! We save things we absolutely need to restart the calculation, some useful
    ! metadata (to make it possible to figure out where the restart file came
    ! from) and some small data items to make life easier and avoid recomputing
    ! them.

    ! WARNING: We use some of the Fortran 2003 interfaces so HDF5 must be
    ! compiled with them enabled (i.e. --enable-fortran --enable-fortran2003 in
    ! the configure line).

    ! See HDF5 documentation and tutorials (http://www.hdfgroup.org/HDF5/ and
    ! http://www.hdfgroup.org/HDF5/Tutor/).
    ! It's a bit hard to get going (not all examples are correct/helpful/self-explanatory!)
    ! but fortunately we restrict ourselves to just a simple usage and there is not enough
    ! space to regurgitate/add thorough explanation to the full HDF5 documentation.

    ! The HDF5 structure we use is:
    ! (Note: order does not matter as HDF5 requires explicitly statement of the group and
    ! dataspace names for each read and write operation.)

    ! /                                # ROOT/
    !
    !  metadata/
    !           restart version        # Version of restart module used to produce the restart file.
    !           hande version          # git sha1 hash.  For info only (not currently used on read-in).
    !           date                   # For info only (not currently used on read-in).
    !           uuid                   # UUID of calculation.  For info only (not currently used on read-in).
    !           calc_type              # Calculation type (as given by a parameter in calc).
    !           nprocs                 # Number of processors used in calculation.
    !           i0_length              # Number of bits in an i0 integer.
    !
    !  qmc/
    !      psips/
    !            determinants          # List of occupied determinant bit strings.
    !            populations           # population(s) on each determinant
    !            data                  # data associated with each determinant
    !            total population      # total population for each particle type.
    !            received_list         # list of walkers sent from last iteration of non-blocking calculation.
    !            processor map         # processor map used for load balancing
    !      state/
    !            shift                 # shift (energy offset/population control)
    !            ncycles               # number of Monte Carlo cycles performed
    !      reference/
    !                reference determinant # reference determinant
    !                reference population  # population on reference
    !                Hilbert space reference determinant # reference determinant
    !                                  # defining Hilbert space (see comments in
    !                                  # fciqmc_data for details).
    !
    !  rng/                            # Not used yet.

    ! where XXX/ indicates a group called XXX, YYY indicates a dataset called
    ! YYY and a nested structure indicates group membership and # is used to
    ! denote a comment.

    ! NOTE: We strive to maintain backwards compatibility.  Reading in HDF5
    ! files should always test that a data item exists before reading it; if it
    ! doesn't then handle the situation gracefully (e.g. by falling back to
    ! default values).

    implicit none

    private
    public :: dump_restart_hdf5, read_restart_hdf5, restart_info_global, restart_info_global_shift

    type restart_info_t
        ! If write_id is negative, then it was set by the user in the input file.  Set
        ! Y=-ID-1 to undo the transformation in parse_input, where Y is in the restart
        ! filename below.  If non-negative, generate Y such that the restart filename is unique.
        integer :: write_id ! ID number to write to.
        ! As for write_id but if positive find the highest possible value of Y  out of the
        ! existing restart files (assuming they have sequential values of Y).
        integer :: read_id  ! ID number to write to.
        ! Number of QMC iterations between writing out a restart file.
        integer :: write_restart_freq
        ! Stem to use for creating restart filenames (of the format restart_stem.Y.pX.H5,
        ! where X is the processor rank and Y is a common positive integer related to
        ! write_id/read_id.
        character(8), private :: restart_stem = 'HANDE.RS'
    end type restart_info_t

    ! Global restart info store until we have a calc type which is passed
    ! around...
    type(restart_info_t) :: restart_info_global = restart_info_t(0,0,huge(0))
    ! Global restart info to store the restart information about when the shift turns
    ! on.
    type(restart_info_t) :: restart_info_global_shift = restart_info_t(0,0,huge(0))

    ! Version id of the restart file *produced*.  Please increment if you add
    ! anything to dump_restart_hdf5!
    ! Note that the restart version is not currently used anywhere but might be helpful
    ! when writing post-processing utilities which act upon restart files.
    integer, parameter :: restart_version = 1

    ! Group names...
    character(*), parameter :: gmetadata = 'metadata',  &
                               gqmc = 'qmc',            &
                               gpsips = 'psips',        &
                               gstate = 'state',        &
                               gref = 'reference',      &
                               grng = 'rng'

    ! Dataspace names...
    character(*), parameter :: drestart = 'restart version',        &
                               dcalc = 'calc type',                 &
                               dhande = 'hande version',            &
                               ddate = 'date',                      &
                               duuid = 'uuid',                      &
                               dnprocs = 'nprocs',                  &
                               di0_length = 'i0_length',            &
                               ddets = 'determinants',              &
                               dpops = 'populations',               &
                               ddata = 'data',                      &
                               dshift = 'shift',                    &
                               dncycles = 'ncycles',                &
                               dtot_pop = 'total population',       &
                               dproc_map = 'processor map',         &
                               dspawn = 'received_list',            &
                               dnspawn = 'nspawn',                  &
                               dref = 'reference determinant',      &
                               dref_pop = 'reference population @ t-1', &
                               dhsref = 'Hilbert space reference determinant'

    contains

#ifndef DISABLE_HDF5
        subroutine init_restart_hdf5(ri, write_mode, filename, kinds)

            ! Initialise restart functionality:

            ! * print information line;
            ! * create restart filename;
            ! * create HDF5 types.

            ! NOTE: HDF5 library must be opened (h5open_f) before init_restart_hdf5 is
            ! called and not closed between calling init_restart_hdf5 and operating on
            ! the restart file to ensure the HDF5 types match those calculated here.

            ! In:
            !    ri: restart information. ri%restart_stem and ri%write_id/ri%read_id (as
            !        appropriate) are used.
            !    write_mode: true for writing out a restart file, false for reading one in.
            ! Out:
            !    filename: name of the restart file.
            !    kinds: derived type containing HDF5 types which correspond to the
            !        non-standard integer and real kinds used in HANDE.

            use hdf5_helper, only: hdf5_kinds_t, hdf5_kinds_init

            use parallel, only: nprocs, iproc, parent
            use utils, only: int_fmt, get_unique_filename

            type(restart_info_t), intent(in) :: ri
            logical, intent(in) :: write_mode
            character(*), intent(out) :: filename
            type(hdf5_kinds_t), intent(out) :: kinds

            character(10) :: proc_suf
            integer :: id, ierr
            logical :: exists

            if (write_mode) then
                id = ri%write_id
            else
                id = ri%read_id
            end if

            ! Figure out filename: restart_stem.Y.pX.H5, where Y is related to id and X is the processor rank.
            write (proc_suf,'(".p",'//int_fmt(iproc,0)//',".H5")') iproc
            call get_unique_filename(trim(ri%restart_stem), trim(proc_suf), write_mode, min(id,0), filename)

            ! New HDF5 files have a '.H5' suffix. However, older HANDE restart
            ! files do not have this. Therefore, if the above file does not
            ! exist, try without '.H5'.
            inquire(file=filename, exist=exists)
            if ((.not. write_mode) .and. (.not. exists)) then
                write (proc_suf,'(".p",'//int_fmt(iproc,0)//')') iproc
                call get_unique_filename(trim(ri%restart_stem), trim(proc_suf), write_mode, min(id,0), filename)
            end if

            if (parent) then
                if (write_mode) then
                    write (6,'(1X,"#",1X,"Writing restart file to",1X,a)', advance='no') trim(filename)
                else
                    write (6,'(1X,"Reading restart file from",1X,a)', advance='no') trim(filename)
                end if
                if (nprocs > 1) then
                    write (6,'(1X, "family.")')
                else
                    write (6,'(".")')
                end if
            end if

            call hdf5_kinds_init(kinds)

        end subroutine init_restart_hdf5
#endif

        subroutine dump_restart_hdf5(ri, ncycles, total_population)

            ! Write out a restart file.

            ! In:
            !    ri: restart information.  ri%restart_stem and ri%write_id are used.
            !    ncycles: number of Monte Carlo cycles performed.
            !    total_population: the total population of each particle type.

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_write
#endif
            use const
            use, intrinsic :: iso_c_binding
            use parallel, only: nprocs, iproc, parent, nthreads
            use utils, only: get_unique_filename, int_fmt

            use fciqmc_data, only: walker_dets, walker_population, walker_data, &
                                   shift, f0, hs_f0, tot_walkers,               &
                                   D0_population_cycle, par_info, received_list
            use calc, only: calc_type, non_blocking_comm, GLOBAL_META
            use errors, only: warning

            type(restart_info_t), intent(in) :: ri
            integer, intent(in) :: ncycles
            real(dp), intent(in) :: total_population(:)
#ifndef DISABLE_HDF5
            character(255) :: restart_file


            ! HDF5 kinds
            type(hdf5_kinds_t) :: kinds
            ! HDF5 handles
            integer(hid_t) :: file_id, group_id, subgroup_id

            integer :: date_time(8)
            character(19) :: date_str
            integer :: ierr
            type(c_ptr) :: ptr
            ! Shape of data (sub-)array to be written out.
            integer(HSIZE_T) :: dshape2(2)
            ! Temporary variables so for copying data to which we can also call c_ptr on.
            ! This allows us to use the same array functions for writing out (the small
            ! amount of) scalar data we have to write out.
            real(dp), allocatable, target :: tmp_pop(:)
            real(p), target :: tmp(1)


            ! Initialise HDF5 and open file.
            call h5open_f(ierr)
            call init_restart_hdf5(ri, .true., restart_file, kinds)
            ! NOTE: if file exists (ie user requested we re-use an existing file), then it is overwritten.
            call h5fcreate_f(restart_file, H5F_ACC_TRUNC_F, file_id, ierr)

            ! --- metadata group ---
            call h5gcreate_f(file_id, gmetadata, group_id, ierr)
            call h5gopen_f(file_id, gmetadata, group_id, ierr)

                call hdf5_write(group_id, dhande, GLOBAL_META%git_sha1)

                call hdf5_write(group_id, duuid, GLOBAL_META%uuid)

                call date_and_time(values=date_time)
                ! Print out current time and date as HH:MM:SS DD/MM/YYYY.
                write (date_str,'(2(i0.2,":"),i0.2,1X,2(i0.2,"/"),i4)') date_time(5:7), date_time(3:1:-1)
                call hdf5_write(group_id, ddate, date_str)

                call hdf5_write(group_id, dnprocs, nprocs)

                call hdf5_write(group_id, di0_length, i0_length)

                call hdf5_write(group_id, drestart, restart_version)

                call hdf5_write(group_id, dcalc, calc_type)

            call h5gclose_f(group_id, ierr)

            ! --- qmc group ---
            call h5gcreate_f(file_id, gqmc, group_id, ierr)
            call h5gopen_f(file_id, gqmc, group_id, ierr)

                ! --- qmc/psips group ---
                call h5gcreate_f(group_id, gpsips, subgroup_id, ierr)
                call h5gopen_f(group_id, gpsips, subgroup_id, ierr)

                ! Don't write out the entire array for storing particles but
                ! rather only the slots in use...
                call hdf5_write(subgroup_id, ddets, kinds, shape(walker_dets(:,:tot_walkers)), &
                                 walker_dets(:,:tot_walkers))

                call hdf5_write(subgroup_id, dpops, kinds, shape(walker_population(:,:tot_walkers)), &
                                 walker_population(:,:tot_walkers))

                call hdf5_write(subgroup_id, ddata, kinds, shape(walker_data(:,:tot_walkers)), &
                                 walker_data(:,:tot_walkers))
                if (non_blocking_comm) then
                    call hdf5_write(subgroup_id, dspawn, kinds, shape(received_list%sdata(:,:received_list%head(0,0))), &
                                    received_list%sdata(:,:received_list%head(0,0)))
                    call hdf5_write(subgroup_id, dnspawn, received_list%head(0,0))
                end if
                call hdf5_write(subgroup_id, dproc_map, kinds, shape(par_info%load%proc_map), par_info%load%proc_map)

                ! Can't use c_loc on a assumed shape array.  It's small, so just
                ! copy it.
                allocate(tmp_pop(size(total_population)))
                tmp_pop = total_population
                call hdf5_write(subgroup_id, dtot_pop, kinds, shape(tmp_pop), tmp_pop)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/state group ---
                call h5gcreate_f(group_id, gstate, subgroup_id, ierr)
                call h5gopen_f(group_id, gstate, subgroup_id, ierr)

                    call hdf5_write(subgroup_id, dncycles, ncycles)

                    call hdf5_write(subgroup_id, dshift, kinds, shape(shift), shift)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/reference group ---
                call h5gcreate_f(group_id, gref, subgroup_id, ierr)
                call h5gopen_f(group_id, gref, subgroup_id, ierr)

                    call hdf5_write(subgroup_id, dref, kinds, shape(f0), f0)

                    call hdf5_write(subgroup_id, dhsref, kinds, shape(hs_f0), hs_f0)

                    tmp = D0_population_cycle
                    call hdf5_write(subgroup_id, dref_pop, kinds, shape(tmp), tmp)

                call h5gclose_f(subgroup_id, ierr)

            call h5gclose_f(group_id, ierr)

            ! --- rng group ---
            call h5gcreate_f(file_id, grng, group_id, ierr)
            call h5gopen_f(file_id, grng, group_id, ierr)
            call h5gclose_f(group_id, ierr)

            ! And terminate HDF5.
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
#else
            call warning('dump_restart_hdf5', '# Not compiled with HDF5 support.  Cannot write out restart file.')
#endif

        end subroutine dump_restart_hdf5

        subroutine read_restart_hdf5(ri)

            ! Read QMC data from restart file.

            ! In:
            !    ri: restart information.  ri%restart_stem and ri%read_id are used.

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_read
#endif
            use errors, only: stop_all
            use const

            use fciqmc_data, only: walker_dets, walker_population, walker_data,  &
                                   shift, tot_nparticles, f0, hs_f0,             &
                                   D0_population, mc_cycles_done, tot_walkers,   &
                                   par_info, received_list
            use calc, only: calc_type, exact_diag, lanczos_diag, mc_hilbert_space, &
                            non_blocking_comm
            use parallel, only: nprocs

            type(restart_info_t), intent(in) :: ri

#ifndef DISABLE_HDF5
            ! HDF5 kinds
            type(hdf5_kinds_t) :: kinds
            ! HDF5 handles
            integer(hid_t) :: file_id, group_id, subgroup_id, dset_id, dspace_id

            character(255) :: restart_file
            integer :: restart_version_restart, calc_type_restart, nprocs_restart
            integer :: i0_length_restart
            type(c_ptr) :: ptr
            integer :: ierr
            real(p), target :: tmp(1)
            logical :: exists

            integer(HSIZE_T) :: dims(size(shape(walker_dets))), maxdims(size(shape(walker_dets)))


            ! Initialise HDF5 and open file.
            call h5open_f(ierr)
            call init_restart_hdf5(ri, .false., restart_file, kinds)
            call h5fopen_f(restart_file, H5F_ACC_RDONLY_F, file_id, ierr)
            if (ierr/=0) then
               call stop_all('read_restart_hdf5', "Unable to open restart file.")
            endif

            ! --- metadata group ---
            call h5gopen_f(file_id, gmetadata, group_id, ierr)

                call hdf5_read(group_id, dnprocs, nprocs_restart)

                call hdf5_read(group_id, drestart, restart_version_restart)

                call hdf5_read(group_id, di0_length, i0_length_restart)

                call hdf5_read(group_id, dcalc, calc_type_restart)

                ! [todo] - Allow restart files for one calculation types to be used to
                ! [todo] - restart a (suitably compatible) different calculation.
                ! AJWT (correctly) doesn't like the low-level handling of the
                ! calc_type bit string.  It's not very modular and doesn't
                ! really belong in the restart code.
                ! However, this will all change as the purity work progresses.
                ! [todo] - refactor calc_type handling into a reusable procedure.

                ! Different calc types are either not compatible or require
                ! hyperslabs (fewer particle types) or require copying (more
                ! particle types).
                ! Clear the flags for non-QMC calculations (which aren't
                ! restarted anyway and don't affect the QMC calculation).
                calc_type_restart = ieor(calc_type, calc_type_restart)
                calc_type_restart = iand(calc_type_restart, not(exact_diag))
                calc_type_restart = iand(calc_type_restart, not(lanczos_diag))
                calc_type_restart = iand(calc_type_restart, not(mc_hilbert_space))
                if (calc_type_restart /= 0) &
                    call stop_all('read_restart_hdf5', &
                                  'Restarting with different calculation types not supported.  Please implement.')
                ! Different restart versions require graceful handling of the
                ! additions/removals.
                if (restart_version /= restart_version_restart) &
                    call stop_all('read_restart_hdf5', &
                                  'Restarting from a different restart version not supported.  Please implement.')
                ! Different processor counts requires figuring out if
                ! a determinant should be on the processor or not (and reading
                ! in chunks).
                if (nprocs /= nprocs_restart) &
                    call stop_all('read_restart_hdf5', &
                                  'Restarting on a different number of processors not supported.  Please implement.')

                if (i0_length /= i0_length_restart) &
                    call stop_all('read_restart_hdf5', &
                                  'Restarting with a different (i0) bit string length not supported.  Please implement.')
            call h5gclose_f(group_id, ierr)

            ! --- qmc group ---
            call h5gopen_f(file_id, gqmc, group_id, ierr)

                ! --- qmc/psips group ---
                call h5gopen_f(group_id, gpsips, subgroup_id, ierr)

                ! Figure out how many determinants we wrote out...
                ! walker_dets has rank 2, so need not look that up!
                call h5dopen_f(subgroup_id, ddets, dset_id, ierr)
                call h5dget_space_f(dset_id, dspace_id, ierr)
                call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
                call h5dclose_f(dset_id, ierr)
                ! Number of determinants is the last index...
                tot_walkers = dims(size(dims))

                call hdf5_read(subgroup_id, ddets, kinds, shape(walker_dets), walker_dets)

                call hdf5_read(subgroup_id, dpops, kinds, shape(walker_population), walker_population)

                call hdf5_read(subgroup_id, ddata, kinds, shape(walker_data), walker_data)

                call hdf5_read(subgroup_id, dtot_pop, kinds, shape(tot_nparticles), tot_nparticles)

                if (non_blocking_comm) then

                    call hdf5_read(subgroup_id, dspawn, kinds, shape(received_list%sdata), received_list%sdata)

                    call hdf5_read(subgroup_id, dnspawn, received_list%head(0,0))

                end if

                call h5lexists_f(subgroup_id, dproc_map, exists, ierr)
                if (exists) call hdf5_read(subgroup_id, dproc_map, kinds, shape(par_info%load%proc_map), par_info%load%proc_map)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/state group ---
                call h5gopen_f(group_id, gstate, subgroup_id, ierr)

                    call hdf5_read(subgroup_id, dncycles, mc_cycles_done)

                    call hdf5_read(subgroup_id, dshift, kinds, shape(shift), shift)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/reference group ---
                call h5gopen_f(group_id, gref, subgroup_id, ierr)

                    call hdf5_read(subgroup_id, dref, kinds, shape(f0), f0)

                    call hdf5_read(subgroup_id, dhsref, kinds, shape(hs_f0), hs_f0)

                    call hdf5_read(subgroup_id, dref_pop, kinds, shape(tmp), tmp)
                    D0_population = tmp(1)

                call h5gclose_f(subgroup_id, ierr)

            call h5gclose_f(group_id, ierr)

            ! --- rng group ---
            call h5gopen_f(file_id, grng, group_id, ierr)
            call h5gclose_f(group_id, ierr)

            ! And terminate HDF5.
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
#else
            call stop_all('read_restart_hdf5', '# Not compiled with HDF5 support.  Cannot read in restart file.')
#endif

        end subroutine read_restart_hdf5

end module restart_hdf5
