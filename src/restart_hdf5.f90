module restart_hdf5
! [review] -  AJWT: We should consider future usage at this point before the format is entrenched!
!               I'll add this to the dev list email.
    ! Restart functionality based on the HDF5 library.  Note: this is only
    ! for QMC (ie FCIQMC, DMQMC or CCMC) calculations).

    ! We save things we absolutely need to restart the calculation, some useful
    ! metadata (to make it possible to figure out where the restart file came
    ! from) and some small data items to make life easier and avoid recomputing
    ! them.

    ! WARNING: We use some of the Fortran 2003 interfaces so HDF5 must be
    ! compiled with them enabled (i.e. --enable-fortran --enable-fortran2003 in
    ! the configure line).

    ! See HDF5 documentation and tutorials (http://www.hdfgroup.org/HDF5/).
    ! It's a bit hard to get going (not all examples are
    ! correct/helpful/self-explanatory!) but fortunately we restrict ourselves
    ! to just a simple usage...

    ! The HDF5 structure we use is:

! [review] -  AJWT: Does order matter?
! [review] -  AJWT: How do we keep this specification consistent with what is actually done?
    ! /                                # ROOT/
    !
    !  metadata/
    !           restart version        # Version of restart module used to produce the restart file.
    !           hande version          # git sha1 hash.  For info only (not used).
! [review] -  AJWT: Probably best to say 'not currently used on read-in'
    !           date                   # For info only (not used).
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
    public :: dump_restart_hdf5, read_restart_hdf5, restart_info_global

    type restart_info_t
! [review] -  AJWT: The comments in parse_input are not helpful in this regard!  More please.
        ! See comments in parse_input regarding the read_id and write_id.
        integer :: write_id ! ID number to write to.
        integer :: read_id  ! ID number to write to.
        integer :: write_restart_freq
        character(255) :: restart_stem = 'HANDE.RS' ! Stem to use for creating restart filenames (of the format restart_stem.Y.pX, where X is the processor rank and Y is a common integer given by write_id or read_id.
    end type restart_info_t

    ! Global restart info store until we have a calc type which is passed
    ! around...
    type(restart_info_t) :: restart_info_global = restart_info_t(0,0,huge(0))

    ! Version id of the restart file *produced*.  Please increment if you add
    ! anything to dump_restart_hdf5!
! [review] -  AJWT: Although the git history will protect the meaning of the version number, is it advisable to document this elsewhere too?
! [review] -  AJWT: I presume we will deal with conflicts in this organically
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
                               dref = 'reference determinant',      &
                               dref_pop = 'reference population @ t-1', &
                               dhsref = 'Hilbert space reference determinant'

    contains

        subroutine dump_restart_hdf5(ncycles, total_population)

            ! Write out a restart file.

            ! In:
            !    ncycles: number of Monte Carlo cycles performed.
            !    total_population: the total population of each particle type.

            use hdf5
            use const
            use, intrinsic :: iso_c_binding
            use report, only: VCS_VERSION, GLOBAL_UUID
            use parallel, only: nprocs, iproc, parent
            use utils, only: get_unique_filename, int_fmt

            use fciqmc_data, only: walker_dets, walker_population, walker_data, &
                                   shift, f0, hs_f0, tot_walkers,               &
                                   D0_population_cycle
            use calc, only: calc_type

            integer, intent(in) :: ncycles
            integer(lint), intent(in) :: total_population(:)
! [review] -  AJWT: This 255 character limit seems a trifle out-dated!
            character(255) :: restart_file

            ! HDF5 kinds
            integer(hid_t) :: h5_i0, h5_p, h5_lint
            ! HDF5 handles
            integer(hid_t) :: file_id, group_id, subgroup_id

            integer :: date_time(8)
            character(19) :: date_str
            integer :: ierr
            type(c_ptr) :: ptr
            integer(HSIZE_T) :: dshape2(2)
            integer(lint), allocatable, target :: tmp_pop(:)
            character(10) :: proc_suf
            real(p), target :: tmp(1)
! [review] -  AJWT: By now I'm getting the sinking feeling from the variable list that this procedure is quite monolithic!

! [review] -  AJWT: A comment as to the format of the filename would be helpful.  I'd've written one, but I couldn't immediately figure it out
! [review] -  AJWT: Having looked back I saw:   the format is restart_stem.Y.pX, where X is the processor rank and Y is a common integer given by write_id or read_id.
            ! Figure out filename.
            write (proc_suf,'(".p",'//int_fmt(iproc,0)//')') iproc
! [review] -  AJWT: Might ri be an input parameter whose value defaults to restart_info_global?
            associate(ri => restart_info_global)
                if (ri%write_id < 0) then
                    call get_unique_filename(trim(ri%restart_stem), trim(proc_suf), .true., ri%write_id, restart_file)
                else
                    call get_unique_filename(trim(ri%restart_stem), trim(proc_suf), .true., 0, restart_file)
                end if
            end associate

            if (parent) then
                if (nprocs > 1) then
                    write (6,'(1X,"#",1X,"Writing restart file to",1X,a,1X,"family.")') trim(restart_file)
                else
                    write (6,'(1X,"#",1X,"Writing restart file to",1X,a)') trim(restart_file)//'.'
                end if
            end if

            ! Initialise HDF5 and open file.
            ! NOTE: if file exists, then it is overwritten.
! [review] -  AJWT: But the get_unique_filename above ensures this doesn't happen?
            call h5open_f(ierr)
            call h5fcreate_f(restart_file, H5F_ACC_TRUNC_F, file_id, ierr)

! [review] -  AJWT: I don't immediately see what these are used for.
! [review] -  AJWT: These kinds are used in the write_* functions.
            ! i0 kind?  What a nice interface!
            h5_i0 = h5kind_to_type(i0, H5_INTEGER_KIND)
            h5_p = h5kind_to_type(p, H5_REAL_KIND)
            h5_lint = h5kind_to_type(lint, H5_INTEGER_KIND)

! [review] -  AJWT: This is getting a bit cryptic (but ok if you bear with it)
            ! --- metadata group ---
            call h5gcreate_f(file_id, gmetadata, group_id, ierr)
            call h5gopen_f(file_id, gmetadata, group_id, ierr)

                call write_string(group_id, dhande, VCS_VERSION)

! [review] -  AJWT: This doesn't appear to agree with the comments at the top
                call write_string(group_id, duuid, GLOBAL_UUID)

                call date_and_time(values=date_time)
! [review] -  AJWT: What does this actually look like?
                write (date_str,'(2(i0.2,":"),i0.2,1X,2(i0.2,"/"),i4)') date_time(5:7), date_time(3:1:-1)
                call write_string(group_id, ddate, date_str)

                call write_integer(group_id, dnprocs, nprocs)

                call write_integer(group_id, di0_length, i0_length)

                call write_integer(group_id, drestart, restart_version)

                call write_integer(group_id, dcalc, calc_type)

            call h5gclose_f(group_id, ierr)

            ! --- qmc group ---
            call h5gcreate_f(file_id, gqmc, group_id, ierr)
            call h5gopen_f(file_id, gqmc, group_id, ierr)

                ! --- qmc/psips group ---
                call h5gcreate_f(group_id, gpsips, subgroup_id, ierr)
                call h5gopen_f(group_id, gpsips, subgroup_id, ierr)

                ! Don't write out the entire array for storing particles but
                ! rather only the slots in use...
                dshape2(1) = size(walker_dets, dim=1, kind=HSIZE_T)
                dshape2(2) = tot_walkers
                ptr = c_loc(walker_dets)
                call write_ptr(subgroup_id, ddets, h5_i0, size(shape(walker_dets)), dshape2, ptr)

                dshape2(1) = size(walker_population, dim=1, kind=HSIZE_T)
                ptr = c_loc(walker_population)
                call write_ptr(subgroup_id, dpops, H5T_NATIVE_INTEGER, size(shape(walker_population)), dshape2, ptr)

                dshape2(1) = size(walker_data, dim=1, kind=HSIZE_T)
                ptr = c_loc(walker_data)
                call write_ptr(subgroup_id, ddata, h5_p, size(shape(walker_data)), dshape2, ptr)

                ! Can't use c_loc on a assumed shape array.  It's small, so just
                ! copy it.
                allocate(tmp_pop(size(total_population)))
                tmp_pop = total_population
                ptr = c_loc(tmp_pop)
                call write_ptr(subgroup_id, dtot_pop, h5_lint, size(shape(tmp_pop)), shape(tmp_pop, HSIZE_T), ptr)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/state group ---
                call h5gcreate_f(group_id, gstate, subgroup_id, ierr)
                call h5gopen_f(group_id, gstate, subgroup_id, ierr)

                    call write_integer(subgroup_id, dncycles, ncycles)

                    ptr = c_loc(shift)
                    call write_ptr(subgroup_id, dshift, h5_p, size(shape(shift)), shape(shift,HSIZE_T), ptr)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/reference group ---
                call h5gcreate_f(group_id, gref, subgroup_id, ierr)
                call h5gopen_f(group_id, gref, subgroup_id, ierr)

                    ptr = c_loc(f0)
                    call write_ptr(subgroup_id, dref, h5_i0, size(shape(f0)), shape(f0,HSIZE_T), ptr)

                    ptr = c_loc(hs_f0)
                    call write_ptr(subgroup_id, dhsref, h5_i0, size(shape(hs_f0)), shape(hs_f0,HSIZE_T), ptr)
                    tmp = D0_population_cycle
                    ptr = c_loc(tmp)
                    call write_ptr(subgroup_id, dref_pop, h5_p, 1, [1_HSIZE_T], ptr)

                call h5gclose_f(subgroup_id, ierr)

            call h5gclose_f(group_id, ierr)

            ! --- rng group ---
            call h5gcreate_f(file_id, grng, group_id, ierr)
            call h5gopen_f(file_id, grng, group_id, ierr)
            call h5gclose_f(group_id, ierr)

            ! And terminate HDF5.
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)

        end subroutine dump_restart_hdf5

        subroutine read_restart_hdf5()

            ! Read QMC data from restart file.

            use hdf5
            use errors, only: stop_all
            use const
            use parallel, only: nprocs, iproc, parent
            use utils, only: int_fmt, get_unique_filename

            use fciqmc_data, only: walker_dets, walker_population, walker_data,  &
                                   shift, tot_nparticles, f0, hs_f0,             &
                                   D0_population, mc_cycles_done, tot_walkers
            use calc, only: calc_type, exact_diag, lanczos_diag, mc_hilbert_space

            ! HDF5 kinds
            integer(hid_t) :: h5_i0, h5_p, h5_lint
            ! HDF5 handles
            integer(hid_t) :: file_id, group_id, subgroup_id, dset_id, dspace_id

            character(255) :: restart_file
            integer :: restart_version_restart, calc_type_restart, nprocs_restart
            integer :: i0_length_restart
            type(c_ptr) :: ptr
            integer :: ierr
            character(10) :: proc_suf
            real(p), target :: tmp(1)

            integer(HSIZE_T) :: dims(size(shape(walker_dets))), maxdims(size(shape(walker_dets)))

! [review] -  AJWT: This seems like needless duplication of what happens in dump_restart_hdf5 which could be put in a procedure
            ! Figure out filename.
            write (proc_suf,'(".p",'//int_fmt(iproc,0)//')') iproc
            associate(ri => restart_info_global)
                if (ri%read_id < 0) then
                    call get_unique_filename(trim(ri%restart_stem), trim(proc_suf), .false., ri%read_id, restart_file)
                else
                    call get_unique_filename(trim(ri%restart_stem), trim(proc_suf), .false., 0, restart_file)
                end if
            end associate

            if (parent) then
                if (nprocs > 1) then
                    write (6,'(1X,"Reading restart file from the",1X,a,1X,"family."/)') trim(restart_file)
                else
                    write (6,'(1X,"Reading restart file from",1X,a,/)') trim(restart_file)//'.'
                end if
            end if

            ! Initialise HDF5 and open file.
            call h5open_f(ierr)
            call h5fopen_f(restart_file, H5F_ACC_RDONLY_F, file_id, ierr)

            ! i0 kind?  What a nice interface!
            h5_i0 = h5kind_to_type(i0, H5_INTEGER_KIND)
            h5_p = h5kind_to_type(p, H5_REAL_KIND)
            h5_lint = h5kind_to_type(lint, H5_INTEGER_KIND)

! [review] -  AJWT: Here endeth the duplication
            ! --- metadata group ---
            call h5gopen_f(file_id, gmetadata, group_id, ierr)

                call read_integer(group_id, dnprocs, nprocs_restart)

                call read_integer(group_id, drestart, restart_version_restart)

                call read_integer(group_id, di0_length, i0_length_restart)

                call read_integer(group_id, dcalc, calc_type_restart)


! [review] -  AJWT: While bit strings are nice, I think the code below lacks modularity.
!              I foresee a time when the calc_type format will change, and the code
!              below will be a pain.  Perhaps some sort of interface for this.
!             In particular, I don't think the restart_read should be dealing with this sort of thing.
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

                ptr = c_loc(walker_dets)
                call read_ptr(subgroup_id, ddets, h5_i0, ptr)

                ptr = c_loc(walker_population)
                call read_ptr(subgroup_id, dpops, H5T_NATIVE_INTEGER, ptr)

                ptr = c_loc(walker_data)
                call read_ptr(subgroup_id, ddata, h5_p, ptr)

                ptr = c_loc(tot_nparticles)
                call read_ptr(subgroup_id, dtot_pop, h5_lint, ptr)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/state group ---
                call h5gopen_f(group_id, gstate, subgroup_id, ierr)

                    call read_integer(subgroup_id, dncycles, mc_cycles_done)

                    ptr = c_loc(shift)
                    call read_ptr(subgroup_id, dshift, h5_p, ptr)

                call h5gclose_f(subgroup_id, ierr)

                ! --- qmc/reference group ---
                call h5gopen_f(group_id, gref, subgroup_id, ierr)

                    ptr = c_loc(f0)
                    call read_ptr(subgroup_id, dref, h5_i0, ptr)

                    ptr = c_loc(hs_f0)
                    call read_ptr(subgroup_id, dhsref, h5_i0, ptr)

                    ptr = c_loc(tmp)
                    call read_ptr(subgroup_id, dref_pop, h5_p, ptr)
                    D0_population = tmp(1)

                call h5gclose_f(subgroup_id, ierr)

            call h5gclose_f(group_id, ierr)

            ! --- rng group ---
            call h5gopen_f(file_id, grng, group_id, ierr)
            call h5gclose_f(group_id, ierr)

            ! And terminate HDF5.
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)

        end subroutine read_restart_hdf5

        ! === Helper procedures: writing ===

! [review] -  AJWT: The write_* lends it self to being overloaded to a single function write_hdf5.  Similarly read_*
        subroutine write_string(id, dset, string)

            ! Write a string to an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    string: string to write out.

            use hdf5

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset, string

            integer(hid_t) :: type_id, dspace_id, dset_id
            integer :: ierr
            character(len(string)) :: sarr(1)

! [review] -  AJWT: This indicates perhaps a misunderstanding of the format which might not be good news
!               for compatability.  Any chance this can be looked at again?
            ! Can't figure out how to write one string.  Just copy into an array
            ! (which I can get to work, bizarrely!).  Only have a few strings to
            ! write out, so performance overhead is essentially 0.
            sarr = string

            ! Set up fortran string type...
            call h5tcopy_f(H5T_STRING, type_id, ierr)
            call h5tset_strpad_f(type_id, H5T_STR_NULLPAD_F, ierr)

            ! Create space and write string.
            call h5screate_simple_f(1, [1_HSIZE_T], dspace_id, ierr)
            call h5dcreate_f(id, dset, type_id, dspace_id, dset_id, ierr)
            call h5dwrite_vl_f(dset_id, type_id, sarr, [len(sarr(1),HSIZE_T),1_HSIZE_T], &
                               [len(sarr(1),SIZE_T)], ierr, dspace_id)
            call h5sclose_f(dspace_id, ierr)
            call h5dclose_f(dset_id, ierr)

            ! Release fortran string type.
            call h5tclose_f(type_id, ierr)

        end subroutine write_string

        subroutine write_integer(id, dset, val)

            ! Write an integer to an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    val: integer to write out.

            use hdf5

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer, intent(in) :: val

            integer(hid_t) :: dspace_id, dset_id
            integer :: ierr

            call h5screate_f(H5S_SCALAR_F, dspace_id, ierr)
            call h5dcreate_f(id, dset, H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)

            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, val, [0_HSIZE_T,0_HSIZE_T], ierr)

            call h5dclose_f(dset_id, ierr)
            call h5sclose_f(dspace_id, ierr)

        end subroutine write_integer

! [review] -  AJWT: This looks potentially dangerous - if dtype differs from the actual type of arr_ptr
!              Might an overloaded interface not be better?
        subroutine write_ptr(id, dset, dtype, arr_rank, arr_dim, arr_ptr)

            ! Write an array to an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    dtype: HDF5 data type of array.
            !    arr_rank: rank of array.
            !    arr_dim: size of array along each dimension.
            !    arr_ptr: C pointer to first element in array to be written out.

            ! NOTE: get dtype from h5kind_to_type if not using a native HDF5
            ! Fortran type.

            use hdf5
            use, intrinsic :: iso_c_binding

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer(hid_t), intent(in) :: dtype
            integer, intent(in) :: arr_rank
            integer(hsize_t), intent(in) :: arr_dim(:)
            type(c_ptr), intent(in) :: arr_ptr

            integer :: ierr
            integer(hid_t) :: dspace_id, dset_id

            call h5screate_simple_f(arr_rank, arr_dim, dspace_id, ierr)
            call h5dcreate_f(id, dset, dtype, dspace_id, dset_id, ierr)

            call h5dwrite_f(dset_id, dtype, arr_ptr, ierr)

            call h5dclose_f(dset_id, ierr)
            call h5sclose_f(dspace_id, ierr)

        end subroutine write_ptr

        ! === Helper procedures: reading ===

        subroutine read_integer(id, dset, val)

            ! Read an integer from an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            ! Out:
            !    val: integer read from HDF5 file.

            use hdf5

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer, intent(out) :: val

            integer(hid_t) :: dset_id
            integer :: ierr

            call h5dopen_f(id, dset, dset_id, ierr)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, val, [0_HSIZE_T,0_HSIZE_T], ierr)
            call h5dclose_f(dset_id, ierr)

        end subroutine read_integer

! [review] -  AJWT: I feel the frisson of data overflow and horrible bugs for the future in this code.
!               Perhaps at least a length check?
        subroutine read_ptr(id, dset, dtype, arr_ptr)

            ! Read an array from an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    dtype: HDF5 data type of array.
            ! In/Out:
            !    arr_ptr: C pointer to first element in array to read.  On
            !        output, the dataset is store in the array pointed to by
            !        arr_ptr.

            ! NOTE: get dtype from h5kind_to_type if not using a native HDF5
            ! Fortran type.

            use hdf5
            use, intrinsic :: iso_c_binding

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer(hid_t), intent(in) :: dtype
            type(c_ptr), intent(inout) :: arr_ptr

            integer :: ierr
            integer(hid_t) :: dset_id

            call h5dopen_f(id, dset, dset_id, ierr)
            call h5dread_f(dset_id, dtype, arr_ptr, ierr)
            call h5dclose_f(dset_id, ierr)

        end subroutine read_ptr

end module restart_hdf5
