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
    !            resort                # if present and true, the psip information must be re-sorted before use.
    !            scaling factor        # population scaling factor for real amplitudes.
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
    !      basis/
    !            nbasis                # Number of basis functions
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
    public :: dump_restart_hdf5, read_restart_hdf5, restart_info_t, init_restart_info_t, redistribute_restart_hdf5, &
              dump_restart_file_wrapper

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
                               grng = 'rng',            &
                               gbasis = 'basis'

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
                               dhsref = 'Hilbert space reference determinant', &
                               dscaling = 'population scale factor', &
                               dnbasis = 'nbasis'

    contains

        subroutine init_restart_info_t(ri, write_id, read_id)

            ! In (optional):
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
        subroutine init_restart_hdf5(ri, write_mode, filename, kinds, ip, verbose, fname_id)

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
            !    fname_id (optional): integer used to label the restart file, i.e. Y in
            !        restart_stem.Y.pX.H5.

            use hdf5_helper, only: hdf5_kinds_t, hdf5_kinds_init

            use parallel, only: nprocs, iproc, parent
            use utils, only: int_fmt, get_unique_filename

            type(restart_info_t), intent(in) :: ri
            logical, intent(in) :: write_mode
            character(*), intent(out) :: filename
            type(hdf5_kinds_t), intent(out), optional :: kinds
            integer, intent(in), optional :: ip
            logical, intent(in), optional :: verbose
            integer, intent(out), optional :: fname_id

            character(14) :: proc_suf
            integer :: id, ip_loc
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
            call get_unique_filename(trim(ri%restart_stem), trim(proc_suf), write_mode, min(id,0), filename, fname_id)

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

        subroutine dump_restart_hdf5(ri, qs, ncycles, total_population, nbasis, nb_comm)

            ! Write out a restart file.

            ! In:
            !    ri: restart information.  ri%restart_stem and ri%write_id are used.
            !    qs: QMC state to write to restart file.
            !    ncycles: number of Monte Carlo cycles performed.
            !    total_population: the total population of each particle type.
            !    nbasis: number of basis functions
            !    nb_comm: true if using non-blocking communications, in which case
            !       information from qs%spawn_store%spawn_recv is also written out
            !       to the restart file.

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_write
            use parallel, only: nprocs
            use calc, only: calc_type, GLOBAL_META
#else
            use parallel, only: parent
#endif
            use const
            use, intrinsic :: iso_c_binding
            use errors, only: stop_all
            use utils, only: get_unique_filename, int_fmt

            use errors, only: warning
            use qmc_data, only: qmc_state_t

            type(restart_info_t), intent(in) :: ri
            type(qmc_state_t), intent(in) :: qs
            integer, intent(in) :: ncycles
            real(p), intent(in) :: total_population(:)
            integer, intent(in) :: nbasis
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
            ! Used for array sizes
            integer :: ishape(2)
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

                ! It would be convenient to use array slices here, but unfortunately they might (depending on
                ! compiler) cause the entire psip array to be copied to a temporary.  Instead we pass the
                ! unsliced arrays and indicate the limit of the final dimension of the array to be printed out.
                ishape = shape(qs%psip_list%states)
                call hdf5_write(subgroup_id, ddets, kinds, ishape, qs%psip_list%states, qs%psip_list%nstates)

                ishape=shape(qs%psip_list%pops)
                call hdf5_write(subgroup_id, dpops, kinds, ishape, qs%psip_list%pops, qs%psip_list%nstates)

                ishape=shape(qs%psip_list%dat)
                call hdf5_write(subgroup_id, ddata, kinds, ishape, qs%psip_list%dat, qs%psip_list%nstates)

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

                ! Always write as int_64 for the sake of conversions
                call hdf5_write(subgroup_id, dscaling, kinds, [1], [int(qs%psip_list%pop_real_factor,int_64)])

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

            ! --- basis group ---
            call h5gcreate_f(file_id, gbasis, group_id, ierr)
            call hdf5_write(group_id, dnbasis, nbasis)
            call h5gclose_f(group_id, ierr)

            ! And terminate HDF5.
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
#else
            if (parent) call warning('dump_restart_hdf5', '# Not compiled with HDF5 support.  Cannot write out restart file.')
#endif

        end subroutine dump_restart_hdf5

        subroutine read_restart_hdf5(ri, nbasis, nb_comm, qs)

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
            use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dtype_equal, dset_shape, hdf5_path
            use restart_utils, only: convert_dets, convert_ref, convert_pops, change_pop_scaling, change_nbasis
            use calc, only: calc_type, exact_diag, lanczos_diag, mc_hilbert_space
            use parallel, only: nprocs
#endif
            use errors, only: stop_all, warning
            use const

            use spawn_data, only: spawn_t
            use qmc_data, only: qmc_state_t
            use sort, only: qsort
            use qmc_common, only: redistribute_particles

            type(restart_info_t), intent(in) :: ri
            logical, intent(in) :: nb_comm
            integer, intent(in) :: nbasis
            type(qmc_state_t), intent(inout) :: qs

#ifndef DISABLE_HDF5
            ! HDF5 kinds
            type(hdf5_kinds_t) :: kinds
            ! HDF5 handles
            integer(hid_t) :: file_id, group_id, subgroup_id

            character(255) :: restart_file
            integer :: restart_version_restart, calc_type_restart, nprocs_restart
            integer :: i0_length_restart, nbasis_restart
            integer :: ierr
            real(p), target :: tmp(1)
            logical :: exists, resort
            integer(int_64) :: restart_scale_factor(1)

            integer(HSIZE_T) :: dims(size(shape(qs%psip_list%states)))


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
                                  'Restarting on a different number of processors not supported.  &
                                  &Use the redistribute function.')

            call h5gclose_f(group_id, ierr)

            ! --- basis group ---
            call h5lexists_f(file_id, gbasis, exists, ierr)
            if (exists) then
                call hdf5_read(file_id, hdf5_path(gbasis, dnbasis), nbasis_restart)
                if (nbasis_restart > nbasis) &
                    call stop_all('read_restart_hdf5', &
                                  'Restarting with a smaller basis not supported.  Please implement.')
            else
                ! assume not changing basis
                nbasis_restart = nbasis
            end if

            ! --- qmc group ---
                call h5gopen_f(file_id, gqmc, group_id, ierr)
                ! --- qmc/psips group ---
                call h5gopen_f(group_id, gpsips, subgroup_id, ierr)

                ! Figure out how many determinants we wrote out...
                ! qs%psip_list%states has rank 2, so need not look that up!
                call dset_shape(subgroup_id, ddets, dims)
                ! Number of determinants is the last index...
                qs%psip_list%nstates = int(dims(size(dims)))

                if (i0_length == i0_length_restart) then
                    if (nbasis == nbasis_restart) then
                        call hdf5_read(subgroup_id, ddets, kinds, shape(qs%psip_list%states), qs%psip_list%states)
                    else
                        ! Change array bounds to restart with a larger basis
                        ! Assume that basis functions 1..nbasis_restart correspond to the original basis
                        call change_nbasis(subgroup_id, ddets, kinds, qs%psip_list%states)
                    end if
                else
                    if (nbasis /= nbasis_restart) &
                        call stop_all('read_restart_hdf5', &
                                      'Changing DET_SIZE and basis size simultaneously not supported.  Please implement.')
                    call convert_dets(subgroup_id, ddets, kinds, qs%psip_list%states)
                end if

                if (.not. dtype_equal(subgroup_id, dpops, kinds%int_p) .and. (bit_size(0_int_p) == 32)) &
                    call warning('read_restart_hdf5', &
                                  'Converting populations from 64 to 32 bit integers.  Overflow may occur. '// &
                                  'Compile HANDE with the CPPFLAG -DPOP_SIZE=64 to use 64-bit populations.')

                call h5lexists_f(subgroup_id, dscaling, exists, ierr)
                if (exists) then
                    call hdf5_read(subgroup_id, dscaling, kinds, shape(restart_scale_factor), restart_scale_factor)
                else
                    ! Assume the scaling factor is unchanged if absent
                    restart_scale_factor = qs%psip_list%pop_real_factor
                end if

                if (dtype_equal(subgroup_id, dpops, kinds%int_p)) then
                    call hdf5_read(subgroup_id, dpops, kinds, shape(qs%psip_list%pops), qs%psip_list%pops)
                else
                    call convert_pops(subgroup_id, dpops, kinds, qs%psip_list%pops, restart_scale_factor(1))
                end if

                if (restart_scale_factor(1) /= qs%psip_list%pop_real_factor) then
                    call change_pop_scaling(qs%psip_list%pops, restart_scale_factor(1), int(qs%psip_list%pop_real_factor,int_64))
                end if

                ! Need to redistribute across processors if int(nbasis/32) changed
                associate(spawn=>qs%spawn_store%spawn, pm=>qs%spawn_store%spawn%proc_map, pl=>qs%psip_list)
                    if (nbasis /= nbasis_restart .and. nprocs > 1) call redistribute_particles(pl%states, pl%pop_real_factor, &
                                                                            pl%pops, pl%nstates, pl%nparticles, spawn)
                end associate

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
                        if (resort) call qsort(pl%nstates, pl%states, pl%pops, pl%dat)
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

                    qs%ref%f0 = 0
                    qs%ref%hs_f0 = 0
                    if (i0_length == i0_length_restart) then
                        call hdf5_read(subgroup_id, dref, kinds, shape(qs%ref%f0), qs%ref%f0)
                        call hdf5_read(subgroup_id, dhsref, kinds, shape(qs%ref%hs_f0), qs%ref%hs_f0)
                    else
                        call convert_ref(subgroup_id, dref, kinds, qs%ref%f0)
                        call convert_ref(subgroup_id, dhsref, kinds, qs%ref%hs_f0)
                    end if

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

        subroutine redistribute_restart_hdf5(ri, nprocs_target, sys)

            ! Create a new set of restart files for a different number of processors than they
            ! were created from.

            ! In/Out:
            !    ri: restart information.  ri%restart_stem, ri%read_id and ri%write_id are used.
            !        On output, ri%read_id is updated (if set) to point to the new id.
            ! In:
            !    number of processors the restart files are to be split over (ie the number of
            !        processors the user wishes to restart the calculation on).
            !    sys (optional): a sys_t object, used to get the basis size.  Only necessary if
            !        changing DET_SIZE for an old restart file.

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_read, hdf5_write, dset_shape, dtype_equal, hdf5_path
            use checking
            use errors, only: warning, stop_all
            use parallel

            use calc, only: ccmc_calc, init_proc_map_t
            use qmc_data, only: ccmc_in_t, particle_t
            use spawn_data, only: proc_map_t
            use particle_t_utils, only: init_particle_t, dealloc_particle_t
            use spawning, only: assign_particle_processor
            use restart_utils, only: convert_dets, convert_ref, convert_pops
#else
            use errors, only: stop_all
#endif
            use system, only: sys_t

            type(restart_info_t), intent(in) :: ri
            integer, intent(in) :: nprocs_target
            type(sys_t), intent(in), optional :: sys

#ifndef DISABLE_HDF5

            ! Number of new restart files to work on at a time.
            integer, parameter :: nmax_files = 10
            ! Max determinants we'll assign to a processor in RAM before writing out to disk.
            ! (Note: memory usage is O(nmax_files*nchunk).
            integer, parameter :: nchunk = 100000

            integer(hid_t) :: orig_id, orig_group_id, orig_subgroup_id
            integer(hid_t) :: group_id, subgroup_id, dset_id
            integer(hid_t) :: new_id
            character(255) :: tmp_name
            character(255), allocatable :: orig_names(:), new_names(:)
            type(hdf5_kinds_t) :: kinds
            integer :: nprocs_read, ierr, i, iproc_min, iproc_max, idet, ndets, ip, nmoved, calc_type_restart
            integer(hsize_t) :: dims(2)

            integer :: hash_shift, hash_seed, move_freq, slot_pos, storage_type, nlinks, max_corder, write_id
            integer :: max_nstates, tensor_label_len, i0_length_restart, nbasis
            integer :: iproc_target_start, iproc_target_end, string_len
            integer, allocatable :: istate_proc(:)
            type(particle_t) :: psip_read, psip_new(0:nmax_files-1)
            logical :: exists
            type(ccmc_in_t) :: ccmc_in_defaults
            type(proc_map_t) :: pm_dummy
            type(restart_info_t) :: ri_write
            integer(i0), allocatable :: f0(:)
            integer(int_64) :: restart_scale_factor(1)

            ! Each processor reads from every restart file but only writes to
            ! a (unique) subset, [iproc_target_start,iproc_target_end].
            ! Thus parallelisation reduces time spent writing out and a little
            ! time on processing the list read in, but each processor still
            ! analyses the each restart file in full (but goes over it fewer
            ! times, perhaps).  Thus, one should not expect this naive algorithm
            ! to scale well...
            iproc_target_end = - 1
            do i = 0, iproc
                iproc_target_start = iproc_target_end + 1
                iproc_target_end = iproc_target_start + nprocs_target/nprocs - 1
                if (i < mod(nprocs_target,nprocs)) iproc_target_end = iproc_target_end + 1
            end do

            ! Hard code 1 load-balancing slot per processor for simplicity.  If the user wishes to use multiple
            ! slots, we should allow this to change when reading in the redistributed restart files.
            call init_proc_map_t(1, pm_dummy, nprocs_target)

            if (ri%write_id < 0 .and. ri%write_id == ri%read_id) &
                call stop_all('redistribute_restart_hdf5', 'Cannot write redistributed restart information to the file(s) from &
                                                           &which the information is read.')

            ! Find the number of processors used to produce the original set of files.
            call h5open_f(ierr)
            call init_restart_hdf5(ri, .false., tmp_name, kinds, 0, .false.)
            call h5fopen_f(tmp_name, H5F_ACC_RDONLY_F, orig_id, ierr)
            if (ierr/=0) call stop_all('redistribute_restart_hdf5', "Unable to open restart file.")
            call hdf5_read(orig_id, hdf5_path(gmetadata,dnprocs), nprocs_read)
            call h5fclose_f(orig_id, ierr)

            ! Create filenames and HDF5 IDs for all old and new files.
            allocate(orig_names(0:nprocs_read-1))
            do i = 0, nprocs_read-1
                call init_restart_hdf5(ri, .false., orig_names(i), ip=i, verbose=i==0.and.parent)
            end do
            allocate(new_names(iproc_target_start:iproc_target_end))
            ! Base new filestem on the one assigned to processor 0 (guaranteed to exist for all sets of restart files.)
            ri_write = ri
            if (parent) then
                call init_restart_hdf5(ri_write, .true., new_names(iproc_target_start), ip=0, verbose=parent, fname_id=write_id)
                ri_write%write_id = -write_id-1
            end if
#ifdef PARALLEL
            call MPI_Bcast(ri_write%write_id, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
#endif
            do i = iproc_target_start, iproc_target_end
                call init_restart_hdf5(ri_write, .true., new_names(i), ip=i, verbose=.false.)
                call h5fcreate_f(new_names(i), H5F_ACC_TRUNC_F, new_id, ierr)
                call h5fclose_f(new_id, ierr)
            end do

            ! Open the original file from proc=0 to get required metadata.  (The choice is arbitrary as
            ! all restart files contain the same metadata.)
            call h5fopen_f(orig_names(0), H5F_ACC_RDONLY_F, orig_id, ierr)

            call h5lexists_f(orig_id, hdf5_path(gqmc, gpsips, dspawn), exists, ierr)
            if (exists) then
                call stop_all('redistribute_restart_hdf5', 'Redistribution from non-blocking calculations &
                                                           &not currently implemented.  Please fix.')
            end if

            ! Get info relating to assigning states to processors.
            hash_seed = 7 ! hard-coded default at time of writing (so will work with past and future restart files)
            move_freq = 0 ! true unless doing CCMC.
            call hdf5_read(orig_id, hdf5_path(gmetadata, dcalc), calc_type_restart)

            ! CARE: as an implementation detail, the CCMC code uses the number of
            ! Monte Carlo cycles (stored in dncycles) as the hash shift.
            call hdf5_read(orig_id, hdf5_path(gqmc, gstate, dncycles), hash_shift)

            call h5lexists_f(orig_id, hdf5_path(gqmc, gstate, dhash_seed), exists, ierr)
            if (exists) then
                call hdf5_read(orig_id, hdf5_path(gqmc, gstate, dhash_seed), hash_seed)
            else if (iand(calc_type_restart, ccmc_calc) /= 0 .and. parent) then
                call warning('redistribute_restart_hdf5', &
                             'hash_seed not found in the restart file.  Using a default hard-coded value.')
            end if

            call h5lexists_f(orig_id, hdf5_path(gqmc, gstate, dmove_freq), exists, ierr)
            if (exists) then
                call hdf5_read(orig_id, hdf5_path(gqmc, gstate, dmove_freq), move_freq)
            else if (iand(calc_type_restart, ccmc_calc) /= 0) then
                ! Only relevant in CCMC.  Require user to set it in input file manually.
                move_freq = ccmc_in_defaults%move_freq
                if (parent) call warning('redistribute_restart_hdf5', &
                                         'move_freq not found in the restart file.  Using a default hard-coded value.')
            end if

            call h5gopen_f(orig_id, gqmc, orig_group_id, ierr)

            ! Check whether the integer type used for determinants is the same as the current
            ! settings.
            call hdf5_read(orig_id, hdf5_path(gmetadata, di0_length), i0_length_restart)
            if (i0_length /= i0_length_restart) then
                ! Try to get determinant string length (needed to convert DET_SIZE) from file,
                ! otherwise the user must supply a system object.
                call h5lexists_f(orig_id, hdf5_path(gqmc, gbasis, dnbasis), exists, ierr)
                if (exists) then
                    call hdf5_read(orig_id, hdf5_path(gqmc, gbasis, dnbasis), nbasis)
                    string_len = ceiling(real(nbasis)/i0_length)
                else if (present(sys)) then
                    string_len = sys%basis%string_len
                else
                    call stop_all('redistribute_restart_hdf5','A system object must be supplied to change DET_SIZE.')
                end if
                allocate(f0(string_len))
            end if

            ! Write out metadata to each new file.
            ! Can just copy it from the first old restart file as it is the same on all files...
            do i = iproc_target_start, iproc_target_end
                call h5fopen_f(new_names(i), H5F_ACC_RDWR_F, new_id, ierr)
                ! /metadata, /basis and /rng
                call h5ocopy_f(orig_id, gmetadata, new_id, gmetadata, ierr)
                call h5lexists_f(orig_id, gbasis, exists, ierr)
                if (exists) call h5ocopy_f(orig_id, gbasis, new_id, gbasis, ierr)
                ! Update determinant integer kind if necessary.
                if (i0_length /= i0_length_restart) then
                    call h5dopen_f(new_id, hdf5_path(gmetadata, di0_length), dset_id, ierr)
                    ! Can't use hdf5_write to replace an existing dataset
                    ! [todo] - can combine this into a single h5dwrite_f?
                    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, i0_length, [0_HSIZE_T,0_HSIZE_T], ierr)
                    call h5dclose_f(dset_id, ierr)
                end if
                call h5ocopy_f(orig_id, grng, new_id, grng, ierr)
                ! ...and non-psip-specific groups in the /qmc group.
                call h5gcreate_f(new_id, gqmc, group_id, ierr)
                call h5gopen_f(new_id, gqmc, group_id, ierr)
                    ! /qmc/state and /qmc/reference
                    call h5ocopy_f(orig_group_id, gstate, group_id, gstate, ierr)
                    if (i0_length == i0_length_restart) then
                        call h5ocopy_f(orig_group_id, gref, group_id, gref, ierr)
                    else
                        ! Need to convert to a different datatype
                        call h5gopen_f(orig_group_id, gref, orig_subgroup_id, ierr)
                        call h5gcreate_f(group_id, gref, subgroup_id, ierr)

                        call convert_ref(orig_subgroup_id, dref, kinds, f0)
                        call hdf5_write(subgroup_id, dref, kinds, shape(f0), f0)

                        call convert_ref(orig_subgroup_id, dhsref, kinds, f0)
                        call hdf5_write(subgroup_id, dhsref, kinds, shape(f0), f0)

                        call h5ocopy_f(orig_group_id, hdf5_path(gref, dref_pop), group_id, hdf5_path(gref, dref_pop), ierr)

                        call h5gclose_f(subgroup_id, ierr)
                        call h5gclose_f(orig_subgroup_id, ierr)
                    end if

                    ! ...and create the /qmc/psips group and fill in the constant (new) processor map and total population.
                    ! (scaling factor updated below)
                    call h5gcreate_f(group_id, gpsips, subgroup_id, ierr)
                    call hdf5_write(group_id, hdf5_path(gpsips, dproc_map), kinds, shape(pm_dummy%map), pm_dummy%map)
                    call h5ocopy_f(orig_group_id, hdf5_path(gpsips, dtot_pop), group_id, hdf5_path(gpsips, dtot_pop), ierr)
                    call hdf5_write(group_id, hdf5_path(gpsips, dresort), .true.)
                call h5gclose_f(group_id, ierr)

                ! Update info.  NOTE: we don't modify the date/time/UUID, just processor count...
                call h5dopen_f(new_id, hdf5_path(gmetadata, dnprocs), dset_id, ierr)
                call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nprocs_target, [0_HSIZE_T,0_HSIZE_T], ierr)
                call h5dclose_f(dset_id, ierr)

                call h5fclose_f(new_id, ierr)
            end do

            call h5fclose_f(orig_id, ierr)

            ! Read the old restart file for each processor in turn and place the psip
            ! information into the new restart file for the appropriate processor.
            ! NOTE: we do not do any load balancing here (and ignore any that was done).
            do i = 0, nprocs_read-1

                ! For each old restart file, we read the entire particle info into RAM.  We then consider a small (nmax_files)
                ! number of processors in nprocs_target at a time and find the particles assigned to that subset of processors.  To
                ! avoid needing needing nmax_files times the amount of memory for holding the particle info, we store up to nchunk
                ! psips on each target processor and periodically write out to disk when we have found that many particles.

                call h5fopen_f(orig_names(i), H5F_ACC_RDONLY_F, orig_id, ierr)
                call h5gopen_f(orig_id, gqmc, orig_group_id, ierr)
                call h5gopen_f(orig_group_id, gpsips, orig_subgroup_id, ierr)

                    ! Check no-one's added to the psips group without (at least) modifying the
                    ! following to handle it.
                    call h5gget_info_f(orig_subgroup_id, storage_type, nlinks, max_corder, ierr)
                    if (nlinks < 4 .or. nlinks > 8) then
                        ! Current datasets in psips group: dtot_pop, dresort, dproc_map, (all handled above), ddets, dpops, ddata (all handled below)
                        ! and dspawn (not currently handled (see above).
                        call stop_all('redistribute_restart_hdf5', &
                                      'psips group in restart file contains an unexpected number of datasets.  Please investigate!')
                    end if

                    call dset_shape(orig_subgroup_id, ddets, dims)
                    if (i0_length == i0_length_restart) then
                        tensor_label_len = int(dims(1))
                    else
                        tensor_label_len = string_len
                    end if
                    max_nstates = int(dims(2))
                    call dset_shape(orig_subgroup_id, dpops, dims)
                    psip_read%nspaces = int(dims(1))
                    call dset_shape(orig_subgroup_id, ddata, dims)
                    psip_read%info_size = int(dims(1)) - psip_read%nspaces

                    psip_new%nspaces = psip_read%nspaces
                    psip_new%info_size = psip_read%info_size

                    call init_particle_t(max_nstates, 0, tensor_label_len, .false., .false., psip_read, .false.)
                    do iproc_min = 0, min(nmax_files-1, iproc_target_end-iproc_target_start)
                        call init_particle_t(max_nstates, 0, tensor_label_len, .false., .false., psip_new(iproc_min), .false.)
                    end do

                    ! Read.
                    if (.not. dtype_equal(orig_subgroup_id, dpops, kinds%int_p) .and. (bit_size(0_int_p) == 32)) &
                        call warning('redistribute_restart_hdf5', &
                                      'Converting populations from 64 to 32 bit integers.  Overflow may occur. '// &
                                      'Compile HANDE with the CPPFLAG -DPOP_SIZE=64 to use 64-bit populations.')

                    if (i0_length == i0_length_restart) then
                        call hdf5_read(orig_subgroup_id, ddets, kinds, shape(psip_read%states), psip_read%states)
                    else
                        call convert_dets(orig_subgroup_id, ddets, kinds, psip_read%states)
                    end if

                    call h5lexists_f(orig_subgroup_id, dscaling, exists, ierr)
                    if (exists) then
                        call hdf5_read(orig_subgroup_id, dscaling, kinds, shape(restart_scale_factor), restart_scale_factor)
                    else
                        restart_scale_factor = 1
                    end if

                    if (dtype_equal(orig_subgroup_id, dpops, kinds%int_p)) then
                        call hdf5_read(orig_subgroup_id, dpops, kinds, shape(psip_read%pops), psip_read%pops)
                    else
                        call convert_pops(orig_subgroup_id, dpops, kinds, psip_read%pops, restart_scale_factor(1))
                    end if
                    call hdf5_read(orig_subgroup_id, ddata, kinds, shape(psip_read%dat), psip_read%dat)

                    ! Distribute.
                    ! [todo] - non-blocking information.
                    ndets = int(dims(2))
                    do iproc_min = iproc_target_start, iproc_target_end, nmax_files
                        nmoved = 0
                        iproc_max = min(iproc_min+nmax_files-1,iproc_target_end)
                        ! istate_proc(i) is the running total number of states found in the current original restart file that
                        ! belong on processor i in the target set of processors.
                        allocate(istate_proc(iproc_min:iproc_max))
                        istate_proc = 0
                        do idet = 1, ndets
                            ! Get processor index (slot_pos is not relevant here as not redoing any load balancing).
                            call assign_particle_processor(psip_read%states(:,idet), tensor_label_len, hash_seed, hash_shift, &
                                                           move_freq, nprocs_target, ip, slot_pos, pm_dummy%map, pm_dummy%nslots)
                            if (ip < iproc_target_start .or. ip > iproc_target_end) then
                                ! Being handled by another processor.  Safely ignore.
                                nmoved = nmoved + 1
                            else if (ip > iproc_max) then
                                ! Leave in cache
                                psip_read%states(:,idet-nmoved) = psip_read%states(:,idet)
                                psip_read%pops(:,idet-nmoved) = psip_read%pops(:,idet)
                                psip_read%dat(:,idet-nmoved) = psip_read%dat(:,idet)
                            else
                                nmoved = nmoved + 1
                                istate_proc(ip) = istate_proc(ip) + 1
                                ! Index for ip relative to iproc_min due to bounds used in psip_new.
                                associate(pl_new=>psip_new(ip-iproc_min))
                                    pl_new%states(:,istate_proc(ip)) = psip_read%states(:,idet)
                                    pl_new%pops(:,istate_proc(ip)) = psip_read%pops(:,idet)
                                    pl_new%dat(:,istate_proc(ip)) = psip_read%dat(:,idet)
                                    if (istate_proc(ip) == nchunk) then
                                        ! Dump out what we've found so far for target processor ip.
                                        call write_psip_info(new_names(ip), kinds, pl_new%states, pl_new%pops, pl_new%dat, &
                                            restart_scale_factor)
                                        istate_proc(ip) = 0
                                    end if
                                end associate
                            end if
                        end do
                        ! Dump out the particles for the current set of processors (iproc_min..iproc_max) that are still in the cache.
                        do ip = iproc_min, iproc_max
                            associate(pl_new=>psip_new(ip-iproc_min), istate=>istate_proc(ip))
                                call write_psip_info(new_names(ip), kinds, pl_new%states(:,:istate), pl_new%pops(:,:istate), &
                                                     pl_new%dat(:,:istate), restart_scale_factor)
                            end associate
                        end do
                        ndets = ndets - nmoved
                        deallocate(istate_proc)
                    end do

                    call dealloc_particle_t(psip_read)
                    do iproc_min = 0, min(nmax_files-1, iproc_target_end-iproc_target_start)
                        call dealloc_particle_t(psip_new(iproc_min))
                    end do

                call h5gclose_f(orig_subgroup_id, ierr)
                call h5gclose_f(orig_group_id, ierr)
                call h5fclose_f(orig_id, ierr)

                if (ndets /= 0) call stop_all('redistribute_restart_hdf5', &
                                              'Failed to redistribute all determinants.  Something went seriously wrong!')

            end do

            if (parent) write (6,'()')

#ifdef PARALLEL
            ! Just in case we go on to use the restart files produced in the same HANDE
            ! run, make sure processing has finished, otherwise a processor may attempt to
            ! read in a file which another processor is writing out to.
            call mpi_barrier(mpi_comm_world, ierr)
#endif

            contains

                subroutine write_psip_info(fname, kinds, psip_dets, psip_pop, psip_data, pop_scaling_factor)

                    ! Write out particle information (state label, population, associated data) to a restart file.

                    ! In:
                    !    fname: HDF5 filename to write to.  Must exist *and* have the qmc/psips group structure already created.
                    !    kinds: hdf5_kinds_t object to convert between Fortran/HANDE kinds and HDF5 types.
                    !    psip_dets: representation of determinants (or similar) labelling the set of occupied states.
                    !    psip_pop: population (in potentially several spaces) on each determinant/state.
                    !    psip_data: system/calculation-specific data associated with each state.
                    !    restart_scale_factor: factor used to encode populations in fixed precision.

                    ! Note that the second dimension of the psip_* arrays is assumed to be identical and every element of the
                    ! arrays passed in is written out.

                    use const, only: i0, int_p, p

                    use hdf5
                    use hdf5_helper, only: hdf5_write, dset_shape, hdf5_kinds_t

                    character(*), intent(in) :: fname
                    type(hdf5_kinds_t), intent(in) :: kinds
                    integer(i0), intent(in) :: psip_dets(:,:)
                    integer(int_p), intent(in) :: psip_pop(:,:)
                    real(p), intent(in) :: psip_data(:,:)
                    integer(int_64), intent(in) :: pop_scaling_factor(1)

                    integer(hid_t) :: file_id, group_id, subgroup_id
                    integer :: ierr
                    logical :: exists

                    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, ierr)
                    call h5gopen_f(file_id, gqmc, group_id, ierr)
                    call h5gopen_f(group_id, gpsips, subgroup_id, ierr)

                    call h5lexists_f(subgroup_id, dscaling, exists, ierr)
                    if (.not.exists) call hdf5_write(subgroup_id, dscaling, kinds, shape(pop_scaling_factor), pop_scaling_factor)
                    call hdf5_write(subgroup_id, ddets, kinds, shape(psip_dets), psip_dets, append=.true.)
                    call hdf5_write(subgroup_id, dpops, kinds, shape(psip_pop), psip_pop, append=.true.)
                    call hdf5_write(subgroup_id, ddata, kinds, shape(psip_data), psip_data, append=.true.)

                    call h5gclose_f(subgroup_id, ierr)
                    call h5gclose_f(group_id, ierr)
                    call h5fclose_f(file_id, ierr)

            end subroutine write_psip_info

#else
            call stop_all('redistribute_restart_hdf5', '# Not compiled with HDF5 support.  Cannot manipulate restart files.')
#endif

        end subroutine redistribute_restart_hdf5

        subroutine dump_restart_file_wrapper(qs, dump_restart_shift, dump_freq, ntot_particles, ireport, ncycles, &
                                             nbasis, ri_freq, ri_shift, nb_comm)

            ! Check if a restart file needs to be written, and if so then do so.

            ! In:
            !     qs: qmc_state_t object.  Particle and related info written out (if desired).
            !     dump_freq: How often (in iterations) to write out a restart file.  Pass in
            !         huge(0) to (effectively) disable.
            !     ntot_particles: total number of particles in each space.
            !     ireport: index of current report loop.
            !     ncycles: the number of iterations per report loop.
            !     nbasis: the number of basis functions
            !     ri_freq: restart_info_t object for periodically writing out the restart file.
            !     ri_freq: restart_info_t object for writing out the restart file once the shift
            !         is turned on.
            !     nb_comm: true if using non-blocking communications.
            ! In/Out:
            !     dump_restart_shift: should we dump a restart file just before
            !         the shift turns on?  If true and a restart file is written out, then
            !         returned as false.

            use const, only: p
            use qmc_data, only: qmc_state_t

            type(qmc_state_t), intent(in) :: qs
            logical, intent(inout) :: dump_restart_shift
            real(p), intent(in) :: ntot_particles(qs%psip_list%nspaces)
            integer, intent(in) :: ireport, ncycles, dump_freq, nbasis
            type(restart_info_t), intent(in) :: ri_freq, ri_shift
            logical, intent(in) :: nb_comm

            if (dump_restart_shift .and. any(qs%vary_shift)) then
                dump_restart_shift = .false.
                call dump_restart_hdf5(ri_shift, qs, qs%mc_cycles_done+ncycles*ireport, &
                                       ntot_particles, nbasis, nb_comm)
            else if (mod(ireport*ncycles,dump_freq) == 0) then
                call dump_restart_hdf5(ri_freq, qs, qs%mc_cycles_done+ncycles*ireport, &
                                       ntot_particles, nbasis, nb_comm)
            end if

        end subroutine dump_restart_file_wrapper

end module restart_hdf5
