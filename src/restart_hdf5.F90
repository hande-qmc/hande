Module restart_hdf5

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
! [review] - AJWT: adding a hyphen is the last resort to help understanding
    !            resort                # if present and true, the psip information must be re-sorted before use.
    !      state/
    !            shift                 # shift (energy offset/population control)
    !            ncycles               # number of Monte Carlo cycles performed
    !            hash seed             # hash seed passed to hash function to assign a state to a processor
    !            move frequency        # (log2 of the) frequency at which the processor location is modified in CCMC
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
    public :: dump_restart_hdf5, read_restart_hdf5, restart_info_t, init_restart_info_t, redistribute_restart_hdf5

    type restart_info_t
        ! If write_id is negative, then set Y=-ID-1 in parse_input, where Y is in the restart
        ! filename below.  If non-negative, generate Y such that the restart filename is unique.
        integer :: write_id ! ID number to write to.
        ! As for write_id but if negative find the highest possible value of Y  out of the
        ! existing restart files (assuming they have sequential values of Y).
        integer :: read_id  ! ID number to write to.
        ! Stem to use for creating restart filenames (of the format restart_stem.Y.pX.H5,
        ! where X is the processor rank and Y is a common positive integer related to
        ! write_id/read_id.
        character(8), private :: restart_stem = 'HANDE.RS'
    end type restart_info_t

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
                               dtot_pop = 'total population',       &
                               dproc_map = 'processor map',         &
                               dspawn = 'received_list',            &
                               dnspawn = 'nspawn',                  &
                               dresort = 'psip_resort',             &
                               dshift = 'shift',                    &
                               dncycles = 'ncycles',                &
                               dhash_seed = 'hash_seed',            &
                               dmove_freq = 'move_freq',            &
                               dref = 'reference determinant',      &
                               dref_pop = 'reference population @ t-1', &
                               dhsref = 'Hilbert space reference determinant'

    contains

        subroutine init_restart_info_t(ri, write_id, read_id)

            ! In:
            !    write_id: id used to write restart file out to.  Default: lowest non-existing file.
            !    read_id: id used to read restart file from.  Default: highest existing file.
            ! Out:
            !    ri: restart_info_t object with fields appropriately set.  (See above.)

            type(restart_info_t), intent(out) :: ri
            integer, intent(in), optional :: write_id, read_id

            ri = restart_info_t(0,0)
            ! ri%read_id or ri%write_id should be non-negative if input is huge(0) (i.e. unset)
            if (present(write_id)) then
                if (write_id == huge(0)) then
                    ri%write_id = 0
                else
                    ri%write_id = -write_id-1
                end if
            end if
            if (present(read_id)) then
                if (read_id == huge(0)) then
                    ri%read_id = 0
                else
                    ri%read_id = -read_id-1
                end if
            end if

        end subroutine init_restart_info_t

#ifndef DISABLE_HDF5
        subroutine init_restart_hdf5(ri, write_mode, filename, kinds, ip, verbose)

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
            !    ip (optional): processor index for the restart filename.  Default: iproc,
            !        the processor index assigned by MPI.
            !    verbose (optional): write output.  Default: true if on parent processor,
            !        false otherwise.
            ! Out:
            !    filename: name of the restart file.
            !    kinds (optional): derived type containing HDF5 types which correspond to the
            !        non-standard integer and real kinds used in HANDE.

            use hdf5_helper, only: hdf5_kinds_t, hdf5_kinds_init

            use parallel, only: nprocs, iproc, parent
            use utils, only: int_fmt, get_unique_filename

            type(restart_info_t), intent(in) :: ri
            logical, intent(in) :: write_mode
            character(*), intent(out) :: filename
            type(hdf5_kinds_t), intent(out), optional :: kinds
            integer, intent(in), optional :: ip
            logical, intent(in), optional :: verbose

            character(14) :: proc_suf
            integer :: id, ierr, ip_loc
            logical :: exists, verbose_loc

            verbose_loc = parent
            if (present(verbose)) verbose_loc = verbose
            ip_loc = iproc
            if (present(ip)) ip_loc = ip

            if (write_mode) then
                id = ri%write_id
            else
                id = ri%read_id
            end if

            ! Figure out filename: restart_stem.Y.pX.H5, where Y is related to id and X is the processor rank.
            write (proc_suf,'(".p",'//int_fmt(ip_loc,0)//',".H5")') ip_loc
            call get_unique_filename(trim(ri%restart_stem), trim(proc_suf), write_mode, min(id,0), filename)

            ! New HDF5 files have a '.H5' suffix. However, older HANDE restart
            ! files do not have this. Therefore, if the above file does not
            ! exist, try without '.H5'.
            inquire(file=filename, exist=exists)
            if ((.not. write_mode) .and. (.not. exists)) then
                write (proc_suf,'(".p",'//int_fmt(iproc,0)//')') ip_loc
                call get_unique_filename(trim(ri%restart_stem), trim(proc_suf), write_mode, min(id,0), filename)
            end if

            if (verbose_loc) then
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

            if (present(kinds)) call hdf5_kinds_init(kinds)

        end subroutine init_restart_hdf5
#endif

        subroutine dump_restart_hdf5(ri, qs, ncycles, total_population, nb_comm)

            ! Write out a restart file.

            ! In:
            !    ri: restart information.  ri%restart_stem and ri%write_id are used.
            !    qs: QMC state to write to restart file.
            !    ncycles: number of Monte Carlo cycles performed.
            !    total_population: the total population of each particle type.
            !    nb_comm: true if using non-blocking communications, in which case
            !       information from qs%spawn_store%spawn_recv is also written out
            !       to the restart file.

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_write
#endif
            use const
            use, intrinsic :: iso_c_binding
            use errors, only: stop_all
            use parallel, only: nprocs, iproc, parent, nthreads
            use utils, only: get_unique_filename, int_fmt

            use calc, only: calc_type, GLOBAL_META
            use errors, only: warning
            use qmc_data, only: qmc_state_t

            type(restart_info_t), intent(in) :: ri
            type(qmc_state_t), intent(in) :: qs
            integer, intent(in) :: ncycles
            real(p), intent(in) :: total_population(:)
            logical, intent(in) :: nb_comm
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
            real(p), allocatable, target :: tmp_pop(:)
            real(p), target :: tmp(1)


            ! Initialise HDF5 and open file.
            call h5open_f(ierr)
            call init_restart_hdf5(ri, .true., restart_file, kinds)
            ! NOTE: if file exists (ie user requested we re-use an existing file), then it is overwritten.
            call h5fcreate_f(restart_file, H5F_ACC_TRUNC_F, file_id, ierr)

            ! --- metadata group ---
            call h5gcreate_f(file_id, gmetadata, group_id, ierr)

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

                ! --- qmc/psips group ---
                call h5gcreate_f(group_id, gpsips, subgroup_id, ierr)

                ! Don't write out the entire array for storing particles but
                ! rather only the slots in use...
                call hdf5_write(subgroup_id, ddets, kinds, shape(qs%psip_list%states(:,:qs%psip_list%nstates)), &
                                 qs%psip_list%states(:,:qs%psip_list%nstates))

                call hdf5_write(subgroup_id, dpops, kinds, shape(qs%psip_list%pops(:,:qs%psip_list%nstates)), &
                                 qs%psip_list%pops(:,:qs%psip_list%nstates))

                call hdf5_write(subgroup_id, ddata, kinds, shape(qs%psip_list%dat(:,:qs%psip_list%nstates)), &
                                 qs%psip_list%dat(:,:qs%psip_list%nstates))
                if (nb_comm) then
                    associate(sdata=>qs%spawn_store%spawn_recv%sdata, head=>qs%spawn_store%spawn_recv%head)
                        call hdf5_write(subgroup_id, dspawn, kinds, shape(sdata(:,:head(0,0))), sdata(:,:head(0,0)))
                        call hdf5_write(subgroup_id, dnspawn, head(0,0))
                    end associate
                end if
                call hdf5_write(subgroup_id, dproc_map, kinds, shape(qs%par_info%load%proc_map%map), qs%par_info%load%proc_map%map)

                ! Can't use c_loc on a assumed shape array.  It's small, so just
                ! copy it.
                allocate(tmp_pop(size(total_population)))
                tmp_pop = total_population
                call hdf5_write(subgroup_id, dtot_pop, kinds, shape(tmp_pop), tmp_pop)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/state group ---
                call h5gcreate_f(group_id, gstate, subgroup_id, ierr)

                    call hdf5_write(subgroup_id, dncycles, ncycles)

                    call hdf5_write(subgroup_id, dshift, kinds, shape(qs%shift), qs%shift)

                    call hdf5_write(subgroup_id, dhash_seed, qs%spawn_store%spawn%hash_seed)

                    call hdf5_write(subgroup_id, dmove_freq, qs%spawn_store%spawn%move_freq)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/qs%ref group ---
                call h5gcreate_f(group_id, gref, subgroup_id, ierr)

                    call hdf5_write(subgroup_id, dref, kinds, shape(qs%ref%f0), qs%ref%f0)

                    call hdf5_write(subgroup_id, dhsref, kinds, shape(qs%ref%hs_f0), qs%ref%hs_f0)

                    tmp = qs%estimators%D0_population
                    call hdf5_write(subgroup_id, dref_pop, kinds, shape(tmp), tmp)

                call h5gclose_f(subgroup_id, ierr)

            call h5gclose_f(group_id, ierr)

            ! --- rng group ---
            call h5gcreate_f(file_id, grng, group_id, ierr)
            call h5gclose_f(group_id, ierr)

            ! And terminate HDF5.
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
#else
            if (parent) call warning('dump_restart_hdf5', '# Not compiled with HDF5 support.  Cannot write out restart file.')
#endif

        end subroutine dump_restart_hdf5

        subroutine read_restart_hdf5(ri, nb_comm, qs)

            ! Read QMC data from restart file.

            ! In:
            !    ri: restart information.  ri%restart_stem and ri%read_id are used.
            !    nb_comm: true if using non-blocking communications.
            ! In/Out:
            !    qs: qmc_state_t object.  Allocated on input, contains info
            !       from the restart file on exit.  If nb_comm is true, then
            !       information is also read into qs%spawn_store%spawn_recv.

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dtype_equal, dset_shape
#endif
            use errors, only: stop_all
            use const

            use calc, only: calc_type, exact_diag, lanczos_diag, mc_hilbert_space
            use parallel, only: nprocs
            use spawn_data, only: spawn_t
            use qmc_data, only: qmc_state_t
            use sort, only: qsort

            type(restart_info_t), intent(in) :: ri
            logical, intent(in) :: nb_comm
            type(qmc_state_t), intent(inout) :: qs

#ifndef DISABLE_HDF5
            ! HDF5 kinds
            type(hdf5_kinds_t) :: kinds
            ! HDF5 handles
            integer(hid_t) :: file_id, group_id, subgroup_id, dset_id, dspace_id

            character(255) :: restart_file
            integer :: restart_version_restart, calc_type_restart, nprocs_restart
            integer :: i0_length_restart
            type(c_ptr) :: ptr
            integer :: ierr, resort
            real(p), target :: tmp(1)
            logical :: exists

            integer(HSIZE_T) :: dims(size(shape(qs%psip_list%states))), maxdims(size(shape(qs%psip_list%states)))


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
                                  'Restarting with a different DET_SIZE is not supported.  Please implement.')
            call h5gclose_f(group_id, ierr)

            ! --- qmc group ---
            call h5gopen_f(file_id, gqmc, group_id, ierr)

                ! --- qmc/psips group ---
                call h5gopen_f(group_id, gpsips, subgroup_id, ierr)

                ! Figure out how many determinants we wrote out...
                ! qs%psip_list%states has rank 2, so need not look that up!
                call dset_shape(subgroup_id, ddets, dims)
                ! Number of determinants is the last index...
                qs%psip_list%nstates = dims(size(dims))

                call hdf5_read(subgroup_id, ddets, kinds, shape(qs%psip_list%states), qs%psip_list%states)

                if (.not. dtype_equal(subgroup_id, dpops, kinds%int_p)) &
                    call stop_all('read_restart_hdf5', &
                                  'Restarting with a different POP_SIZE is not supported.  Please implement.')
                call hdf5_read(subgroup_id, dpops, kinds, shape(qs%psip_list%pops), qs%psip_list%pops)

                call hdf5_read(subgroup_id, ddata, kinds, shape(qs%psip_list%dat), qs%psip_list%dat)

                call hdf5_read(subgroup_id, dtot_pop, kinds, shape(qs%psip_list%tot_nparticles), qs%psip_list%tot_nparticles)

                if (nb_comm) then
                    associate(recv=>qs%spawn_store%spawn_recv)
                        call hdf5_read(subgroup_id, dspawn, kinds, shape(recv%sdata), recv%sdata)
                        call hdf5_read(subgroup_id, dnspawn, recv%head(0,0))
                    end associate
                end if

                call h5lexists_f(subgroup_id, dproc_map, exists, ierr)
                if (exists) call hdf5_read(subgroup_id, dproc_map, kinds, shape(qs%par_info%load%proc_map%map), &
                                           qs%par_info%load%proc_map%map)

                call h5lexists_f(subgroup_id, dresort, exists, ierr)
                if (exists) then
                    call hdf5_read(subgroup_id, dresort, resort)
                    associate(pl=>qs%psip_list)
! [review] - AJWT: the docs above mention if true for resort, but does 1==true in hdf5 world (or indeed in FORTRAN)?
                        if (resort == 1) call qsort(pl%nstates, pl%states, pl%pops, pl%dat)
                    end associate
                end if

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/state group ---
                call h5gopen_f(group_id, gstate, subgroup_id, ierr)

                    call hdf5_read(subgroup_id, dncycles, qs%mc_cycles_done)

                    call hdf5_read(subgroup_id, dshift, kinds, shape(qs%shift), qs%shift)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/reference group ---
                call h5gopen_f(group_id, gref, subgroup_id, ierr)

                    call hdf5_read(subgroup_id, dref, kinds, shape(qs%ref%f0), qs%ref%f0)

                    call hdf5_read(subgroup_id, dhsref, kinds, shape(qs%ref%hs_f0), qs%ref%hs_f0)

                    call hdf5_read(subgroup_id, dref_pop, kinds, shape(tmp), tmp)
                    qs%estimators%D0_population = tmp(1)

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

! [review] - AJWT: This routine is dauntingly long and dense, but I shall persevere!

        subroutine redistribute_restart_hdf5(ri, nprocs_target)

            ! Create a new set of restart files for a different number of processors than they
            ! were created from.

            ! In/Out:
            !    ri: restart information.  ri%restart_stem, ri%read_id and ri%write_id are used.
            !        On output, ri%read_id is updated (if set) to point to the new id.

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_read, hdf5_write, dset_shape, dtype_equal
            use checking
            use errors, only: stop_all
            use const, only: i0, int_p, p
            use parallel, only: parent

            use calc, only: ccmc_calc, init_proc_map_t
            use qmc_data, only: ccmc_in_t
            use spawn_data, only: proc_map_t
            use spawning, only: assign_particle_processor
#endif

            type(restart_info_t), intent(in) :: ri
            integer, intent(in) :: nprocs_target

#ifndef DISABLE_HDF5

            ! Number of new restart files to work on at a time.
! [review] - AJWT: (2nd read) - I see what this parameter is now - can it have a name describing what it's the max of!
            integer, parameter :: nmax = 10
            ! Max determinants we'll assign to a processor in RAM before writing out to disk.
            ! (Note: memory usage is O(nmax*nchunk).
            integer, parameter :: nchunk = 100000

            integer(hid_t) :: orig_id, orig_group_id, orig_subgroup_id
            integer(hid_t) :: group_id, subgroup_id, dset_id
            integer(hid_t) :: new_id
            character(255) :: tmp_name
            character(255), allocatable :: orig_names(:), new_names(:)
            type(hdf5_kinds_t) :: kinds
            integer :: nprocs_read, ierr, i, inew, icurr, iproc_max, idet, ndets, ip, nmoved, calc_type_restart
            real(p), allocatable :: psip_data(:,:)
            integer(i0), allocatable :: psip_dets(:,:)
            integer(int_p), allocatable :: psip_pop(:,:)
            integer(hsize_t) :: dims(2)

            integer :: hash_shift, hash_seed, label_length, move_freq, slot_pos
            integer, allocatable :: ihead(:)
            integer(i0), allocatable :: psip_dets_new(:,:,:)
            real(p), allocatable :: psip_data_new(:,:,:)
            integer(int_p), allocatable :: psip_pop_new(:,:,:)
            logical :: exists
            type(ccmc_in_t) :: ccmc_in_defaults
            type(proc_map_t) :: pm_dummy

            ! Hard code 1 load-balancing slot per processor for simplicity.  If the user wishes to use multiple
            ! slots, we should allow this to change when reading in the redistributed restart files.
            call init_proc_map_t(1, pm_dummy)

            if (.not.parent) call stop_all('redistribute_restart_hdf5', &
                            'Restart redistribution must currently be performed in serial.  Please improve.')
! [review] - AJWT: While I think the code tells me that the filenames for write and read must be different, the error
! [review] - AJWT: message didn't convey that information.
            if (ri%write_id < 0 .and. ri%write_id == ri%read_id) &
                call stop_all('redistribute_restart_hdf5', 'Cannot redistribute restart information in place.')

            ! Find the number of processors used to produce the original set of files.
! [review] - AJWT: which is put into nprocs_read
            call h5open_f(ierr)
            call init_restart_hdf5(ri, .false., tmp_name, kinds, 0, .false.)
            call h5fopen_f(tmp_name, H5F_ACC_RDONLY_F, orig_id, ierr)
            if (ierr/=0) then
               call stop_all('redistribute_restart_hdf5', "Unable to open restart file.")
            endif
            call h5gopen_f(orig_id, gmetadata, orig_group_id, ierr)
            call hdf5_read(orig_group_id, dnprocs, nprocs_read)
            call h5gclose_f(orig_group_id, ierr)
            call h5fclose_f(orig_id, ierr)

            ! Create filenames and HDF5 IDs for all old and new files.
            allocate(orig_names(0:nprocs_read-1))
            do i = 0, nprocs_read-1
                call init_restart_hdf5(ri, .false., orig_names(i), ip=i, verbose=i==0)
            end do
            allocate(new_names(0:nprocs_target-1))
            do i = 0, nprocs_target-1
                call init_restart_hdf5(ri, .true., new_names(i), ip=i, verbose=i==0)
                call h5fcreate_f(new_names(i), H5F_ACC_TRUNC_F, new_id, ierr)
                call h5fclose_f(new_id, ierr)
            end do

! [review] - AJWT: Open the original files
            call h5fopen_f(orig_names(0), H5F_ACC_RDONLY_F, orig_id, ierr)
            call h5gopen_f(orig_id, gqmc, orig_group_id, ierr)
            call h5gopen_f(orig_group_id, gpsips, orig_subgroup_id, ierr)

            ! Get info relating to assigning states to processors.
! [review] - AJWT: Backwards compatability options to be removed at some point.
! [review] - AJWT: We should probably print a warning when these are used (i.e. not previously set in restart files).
            hash_seed = 7 ! hard-coded default at time of writing (so will work with past and future restart files)
            move_freq = 0 ! true unless doing CCMC.
! [review] - AJWT: There must be some nicer wrappers we can create - this is almost unreadable!  What about some functions?
            call h5gopen_f(orig_id, gmetadata, orig_group_id, ierr)
                call hdf5_read(orig_group_id, dcalc, calc_type_restart)
            call h5gclose_f(orig_group_id, ierr)
            call h5gopen_f(orig_id, gqmc, orig_group_id, ierr)
            call h5gopen_f(orig_group_id, gstate, orig_subgroup_id, ierr)
! [review] - AJWT: Shouldn't this be something like dhash_shift not dncycles?
            call hdf5_read(orig_subgroup_id, dncycles, hash_shift)
            call h5lexists_f(orig_subgroup_id, dhash_seed, exists, ierr)
            if (exists) call hdf5_read(orig_subgroup_id, dhash_seed, hash_seed)
            call h5lexists_f(orig_subgroup_id, dmove_freq, exists, ierr)
            if (exists) then
                call hdf5_read(orig_subgroup_id, dmove_freq, move_freq)
            else if (iand(calc_type_restart, ccmc_calc) /= 0) then
                ! Only relevant in CCMC.  Require user to set it in input file manually.
                move_freq = ccmc_in_defaults%move_freq
            end if
            call h5gclose_f(orig_subgroup_id, ierr)
            call h5gopen_f(orig_group_id, gpsips, orig_subgroup_id, ierr)

            ! Write out metadata to each new file.
            ! Can just copy it from the first old restart file as it is the same on all files...
            do i = 0, nprocs_target-1
                call h5fopen_f(new_names(i), H5F_ACC_RDWR_F, new_id, ierr)
                call h5ocopy_f(orig_id, gmetadata, new_id, gmetadata, ierr)
                call h5ocopy_f(orig_id, grng, new_id, grng, ierr)
                ! ...and non-psip-specific groups in the qmc/ group.
! [review] - AJWT: I see that the indentation is meant to help understand structure, though I'm not sure it does much
                call h5gcreate_f(new_id, gqmc, group_id, ierr)
                call h5gopen_f(new_id, gqmc, group_id, ierr)
                    call h5ocopy_f(orig_group_id, gstate, group_id, gstate, ierr)
                    call h5ocopy_f(orig_group_id, gref, group_id, gref, ierr)
                    ! ...and create (but don't fill) the psips group except the constant (new) processor map and total population.
                    call h5gcreate_f(group_id, gpsips, subgroup_id, ierr)
                    call h5gopen_f(group_id, gpsips, subgroup_id, ierr)
                        call hdf5_write(subgroup_id, dproc_map, kinds, shape(pm_dummy%map), pm_dummy%map)
                        call hdf5_write(subgroup_id, dresort, 1)
                        call h5ocopy_f(orig_subgroup_id, dtot_pop, subgroup_id, dtot_pop, ierr)
                    call h5gclose_f(subgroup_id, ierr)
                call h5gclose_f(group_id, ierr)
                ! Update info.
                ! NOTE: we don't modify the date/time/UUID...
                call h5gopen_f(new_id, gmetadata, group_id, ierr)
                    call h5dopen_f(group_id, dnprocs, dset_id, ierr)
                    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nprocs_target, [0_HSIZE_T,0_HSIZE_T], ierr)
                    call h5dclose_f(dset_id, ierr)
                call h5gclose_f(group_id, ierr)
                call h5fclose_f(new_id, ierr)
            end do

            call h5gclose_f(orig_subgroup_id, ierr)
            call h5gclose_f(orig_group_id, ierr)
            call h5fclose_f(orig_id, ierr)

            ! Read the old restart file for each processor in turn and place the psip
            ! information into the new restart file for the appropriate processor.
            ! NOTE: we do not do any load balancing here (and ignore any that was done).
! [review] - AJWT: (3rd read) A description of how there is caching and new files are only 
! [review] - AJWT:            written out a few chunks at a time in blocks on nmax would be helpful.
            do i = 0, nprocs_read-1
                call h5fopen_f(orig_names(i), H5F_ACC_RDONLY_F, orig_id, ierr)
                call h5gopen_f(orig_id, gqmc, orig_group_id, ierr)
                call h5gopen_f(orig_group_id, gpsips, orig_subgroup_id, ierr)

! [review] - AJWT: How can one be sure if all the data has been moved?
! [review] - AJWT: What if somebody adds a bit of data, but doens't change this section?
                    call dset_shape(orig_subgroup_id, ddets, dims)
                    allocate(psip_dets(dims(1), dims(2)), stat=ierr)
                    call check_allocate('psip_dets', size(psip_dets), ierr)
                    allocate(psip_dets_new(dims(1), min(nchunk,dims(2)), 0:min(nprocs_target, nmax)-1), stat=ierr)
                    call check_allocate('psip_dets_new', size(psip_dets), ierr)

                    call dset_shape(orig_subgroup_id, dpops, dims)
                    allocate(psip_pop(dims(1), dims(2)), stat=ierr)
                    call check_allocate('psip_pop', size(psip_pop), ierr)
                    allocate(psip_pop_new(dims(1), min(nchunk,dims(2)), 0:min(nprocs_target, nmax)-1), stat=ierr)
                    call check_allocate('psip_pop_new', size(psip_pop_new), ierr)

                    call dset_shape(orig_subgroup_id, ddata, dims)
                    allocate(psip_data(dims(1), dims(2)), stat=ierr)
                    call check_allocate('psip_data', size(psip_data), ierr)
                    allocate(psip_data_new(dims(1), min(nchunk,dims(2)), 0:min(nprocs_target, nmax)-1), stat=ierr)
                    call check_allocate('psip_data_new', size(psip_data_new), ierr)

                    ! Read.
                    if (.not. dtype_equal(orig_subgroup_id, ddets, kinds%i0)) &
                        call stop_all('redistribute_restart_hdf5', &
                                      'Restarting with a different DET_SIZE is not supported.  Please implement.')
                    if (.not. dtype_equal(orig_subgroup_id, dpops, kinds%int_p)) &
                        call stop_all('redistribute_restart_hdf5', &
                                      'Restarting with a different POP_SIZE is not supported.  Please implement.')

                    call hdf5_read(orig_subgroup_id, ddets, kinds, shape(psip_dets), psip_dets)
                    call hdf5_read(orig_subgroup_id, dpops, kinds, shape(psip_pop), psip_pop)
                    call hdf5_read(orig_subgroup_id, ddata, kinds, shape(psip_data), psip_data)

                    ! Distribute.
                    ! [todo] - non-blocking information.
                    ndets = dims(2)
                    label_length = size(psip_dets, dim=1)
! [review] - AJWT: (3rd read) The fact nmax is defined at the top of the routine and not used until now might indicate that this routine
! [review] - AJWT:              could be broken up
                    do inew = 0, nprocs_target-1, nmax
                        nmoved = 0
! [review] - AJWT: Call me a bear of little brain, but I don't know what this variable does, nor why inew+nmax might be relevant
! [review] - AJWT: (2nd read) naming nmax something sensible would definitely help this!

                        iproc_max = min(inew+nmax,nprocs_target)-1
! [review] - AJWT: Another cryptically named variable whose purpose is to be fathomed.  Might it 
                        allocate(ihead(inew:iproc_max))
                        ihead = 0
                        do idet = 1, ndets
                            ! Get processor index
! [review] - AJWT: Not really knowing much about the load balancing, I assume that slot_pos is irrelevant
                            call assign_particle_processor(psip_dets(:,idet), label_length, hash_seed, hash_shift, move_freq, &
                                                           nprocs_target, ip, slot_pos, pm_dummy%map, pm_dummy%nslots)
                            if (ip > iproc_max) then
                                ! Leave in cache
                                psip_dets(:,idet-nmoved) = psip_data(:,idet)
                                psip_pop(:,idet-nmoved) = psip_data(:,idet)
                                psip_data(:,idet-nmoved) = psip_data(:,idet)
                            else
                                nmoved = nmoved + 1
                                ihead(ip) = ihead(ip) + 1
                                psip_dets_new(:,ihead(ip),ip) = psip_dets(:,idet)
                                psip_pop_new(:,ihead(ip),ip) = psip_pop(:,idet)
                                psip_data_new(:,ihead(ip),ip) = psip_data(:,idet)
                                if (ihead(ip) == nchunk) then
                                    ! Dump out what we've found so far.
                                    call write_psip_info(new_names(ip), kinds, psip_dets_new(:,:,ip), psip_pop_new(:,:,ip), &
                                                         psip_data_new(:,:,ip))
                                    ihead(ip) = 0
                                end if
                            end if
                        end do
! [review] - AJWT: By everything else, I read the code to mean the psips for procs inew...iproc_max which haven't yet been dumped.
                        ! Dump out everything else...
                        do ip = inew, iproc_max
                            call write_psip_info(new_names(ip), kinds, psip_dets_new(:,:ihead(ip),ip), &
                                                 psip_pop_new(:,:ihead(ip),ip), psip_data_new(:,:ihead(ip),ip))
                        end do
                        ndets = ndets - nmoved
                        deallocate(ihead)
                    end do

                    deallocate(psip_dets, stat=ierr)
                    call check_deallocate('psip_dets', ierr)
                    deallocate(psip_dets_new, stat=ierr)
                    call check_deallocate('psip_dets_new', ierr)
                    deallocate(psip_pop, stat=ierr)
                    call check_deallocate('psip_pop', ierr)
                    deallocate(psip_pop_new, stat=ierr)
                    call check_deallocate('psip_pop_new', ierr)
                    deallocate(psip_data, stat=ierr)
                    call check_deallocate('psip_data', ierr)
                    deallocate(psip_data_new, stat=ierr)
                    call check_deallocate('psip_data_new', ierr)

                call h5gclose_f(orig_subgroup_id, ierr)
                call h5gclose_f(orig_group_id, ierr)
                call h5fclose_f(orig_id, ierr)
            end do

! [review] - AJWT: Might it be worth checking that ndets==0 here?

            if (parent) write (6,'()')

            contains

! [review] - AJWT: At least there's one subroutine.  The name is fairly intuitive, bit a comment about the inputs wouldn't go amiss.
                subroutine write_psip_info(fname, kinds, psip_dets, psip_pop, psip_data)

                    use const, only: i0, int_p, p

                    use hdf5
                    use hdf5_helper, only: hdf5_write, dset_shape, hdf5_kinds_t

                    character(*), intent(in) :: fname
                    type(hdf5_kinds_t), intent(in) :: kinds
                    integer(i0), intent(in) :: psip_dets(:,:)
                    integer(int_p), intent(in) :: psip_pop(:,:)
                    real(p), intent(in) :: psip_data(:,:)

                    integer(hid_t) :: file_id, group_id, subgroup_id
                    integer :: ierr
                    integer(hsize_t) :: dims(2)

                    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, ierr)
                    call h5gopen_f(file_id, gqmc, group_id, ierr)
                    call h5gopen_f(group_id, gpsips, subgroup_id, ierr)

                    call hdf5_write(subgroup_id, ddets, kinds, shape(psip_dets), psip_dets, .true.)
                    call hdf5_write(subgroup_id, dpops, kinds, shape(psip_pop), psip_pop, .true.)
                    call hdf5_write(subgroup_id, ddata, kinds, shape(psip_data), psip_data, .true.)

                    call h5gclose_f(subgroup_id, ierr)
                    call h5gclose_f(group_id, ierr)
                    call h5fclose_f(file_id, ierr)

            end subroutine write_psip_info

#else
            call stop_all('redistribute_restart_hdf5', '# Not compiled with HDF5 support.  Cannot manipulate restart files.')
#endif

        end subroutine redistribute_restart_hdf5

end module restart_hdf5
