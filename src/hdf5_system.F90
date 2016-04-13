module hdf5_system

    ! System information dumping based around HDF5 library (shamelessly based
    ! on approach used in restart_hdf5). Can be used with any HANDE calculation
    ! utilising an FCIDUMP to initialise a system (system type read_in).
    !
    ! Format has attempted to strike a balance between outputting everything
    ! possible and regenerating the entire setup. Where existing routines can
    ! easily reproduce the information required the information has not been
    ! included in the output file. This is by no means set in stone and can
    ! be easily adjusted.

    ! For warnings and tutorials on using HDF5 functionality please see the start of
    ! restart_hdf5.F90.

    ! HDF5 struncture used is:

    ! /                                 # ROOT/
    !
    ! metadata/
    !       hande version               # git sha1 hash.  For info only (not currently used on read-in).
    !       uuid                        # UUID of calculation.  For info only (not currently used on read-in).
    !       date                        # For info only (not currently used on read-in).
    !       sysdump version             # Version of hdf5_system module used to produce file.
    !
    ! sys/                              # see system.f90 for descriptions of same-named parameters
    !   system
    !   symmetry
    !   nsym
    !   sym0
    !   sym_max
    !   nsym_tot
    !   sym0_tot
    !   sym_max_tot
    !   CAS
    !
    !   basis/
    !       nbasis
    !                                   # Arrays containing information for
    !                                   # ith basis function at ith position.
    !       basis_spatial_index
    !       basis_symmetry
    !       basis_symmetry_index
    !       basis_symmetry_spin_index
    !       basis_ms
    !       basis_lz
    !       basis_sp_eigv
    !
    !   read_in/
    !       fcidump     (string)
    !       uhf         (bool)
    !       ecore       (real)
    !       uselz       (bool)
    !       comp        (bool)
    !       pg_mask     (int)

    !       integrals/                  # See write_1body_integrals/
    !                                   # write_coulomb_integrals for storage
    !                                   # structure information.
    !           one_body
    !           coulomb_ints
    !           one_body_im
    !           coulomb_ints_im

    ! where XXX/ indicates a group called XXX, YYY indicates a dataset called
    ! YYY and a nested structure indicates group membership and # is used to
    ! denote a comment.

    implicit none

    private
    public :: read_system_hdf5, dump_system_hdf5, get_filename

    ! Version ID of sysdump, just for compatability later on. Increment for any
    ! changes to functionality.
    integer, parameter :: sysdump_version = 0

    ! Group names.
    character(*), parameter ::  gmetadata = 'metadata',     &
                                gsys = 'system',            &
                                gbasis = 'basis',           &
                                gread_in = 'read_in',       &
                                gintegrals = 'integrals'

    ! Dataspace names.
    character(*), parameter ::  dsysdump = 'sysdump version',   &
                                dhande = 'hande version',       &
                                ddate = 'date',                 &
                                duuid = 'uuid',                 &

                                dsystem = 'system',             &
                                dnel = 'nelectrons',            &
                                dms = 'Ms',                     &
                                dsymmetry = 'symmetry',         &
                                dnsym = 'nsym',                 &
                                dsym0 = 'sym0',                 &
                                dsym_max = 'sym_max',           &
                                dnsym_tot = 'nsym_tot',         &
                                dsym0_tot = 'sym0_tot',         &
                                dsym_max_tot = 'sym_max_tot',   &
                                dcas = 'CAS',                   &

                                dfcidump = 'fcidump',           &
                                duhf = 'uhf',                   &
                                decore = 'ecore',               &
                                duselz = 'uselz',               &
                                dintegrals = 'integrals',       &
                                dcomp = 'comp',                 &

                                done_body = 'one_body',         &
                                done_body_im = 'one_body_im',   &
                                dcoulomb_ints = 'coulomb_ints', &
                                dcoulomb_ints_im =              &
                                              'coulomb_ints_im',&
                                dnbasis = 'nbasis',             &
                                dbasis_spat_ind =               &
                                        'basis_spatial_index',  &
                                dbasis_sym = 'basis_symmetry',  &
                                dbasis_sym_index =                &
                                        'basis_symmetry_index', &
                                dbasis_sym_spin_index =           &
                                        'basis_symmetry_spin_index',&
                                dbasis_ms = 'basis_ms',         &
                                dbasis_lz = 'basis_lz',         &
                                dbasis_sp_eigv = 'basis_sp_eigv',&
                                dpg_mask = 'pg_mask'

    contains

#ifndef DISABLE_HDF5
        subroutine init_system_hdf5(write_mode, sys, filename, kinds, verbose)

            ! Initialise HDF5 system functionality:
            ! * set filename;
            ! * print information line;
            ! * create HDF5 types

            ! NOTE: HDF5 library must be opened (h5open_f) before init_system_hdf5 is
            ! called and not closed between calling init_restart_hdf5 and operating on
            ! the restart file to ensure the HDF5 types match those calculated here.

            ! In:
            !   write_mode: true if writing out a system file, false for reading one in.
            !   sys: sys_t type. If writing out contains all required information. If
            !       reading in contains minimal read in info (nel, fcidump).
            !   verbose (optional): write output. Default: true.
            ! In/Out:
            !   filename: name of HDF5 system file.  Used if set to a non-empty string,
            !       otherwise set to a file based upon the original integral file and
            !       (if relevant) the CAS information plus a .H5 suffix.
            ! Out:
            !   kinds: derived tpe containing HDF5 types which correspond to the
            !       non-standard integer and real kinds used in HANDE.

            use hdf5_helper, only: hdf5_kinds_t, hdf5_kinds_init
            use system, only: sys_t

            logical, intent(in) :: write_mode
            type(sys_t), intent(in) :: sys

            type(hdf5_kinds_t), intent(out), optional :: kinds
            logical, intent(in), optional ::  verbose
            character(255), intent(inout) :: filename

            logical :: verbose_loc

            verbose_loc = .true.
            if (present(verbose)) verbose_loc = verbose

            call get_filename(write_mode, sys, filename)

            if (verbose_loc) then
                if (write_mode) then
                    write (6,'(1X,"Writing HDF5 system file to",1X,a,/)') trim(filename)
                else
                    write (6,'(1X,"Reading HDF5 system file from",1X,a,/)') trim(filename)
                end if
            end if

            if (present(kinds)) call hdf5_kinds_init(kinds)

        end subroutine init_system_hdf5

#endif
        subroutine get_filename(write_mode, sys, filename)

            ! Sets filename to write to when writing hdf5 system file. Need to separate
            ! to enable setting of filename to write to in all threads (not just parent)
            ! if want to be able topass out autogenerated filename.
            ! In:
            !   write_mode: true if writing out a system file, false for reading one in.
            !   sys: sys_t type. If writing out contains all required information. If
            !       reading in contains minimal read in info (nel, fcidump).
            ! In/Out:
            !   filename: name of HDF5 system file.  Used if set to a non-empty string,
            !       otherwise set to a file based upon the original integral file and
            !       (if relevant) the CAS information plus a .H5 suffix.

            use system, only: sys_t

            logical, intent(in) :: write_mode
            type(sys_t), intent(in) :: sys

            character(255), intent(inout) :: filename
#ifndef DISABLE_HDF5

            if (write_mode) then
                if (filename == '') then
                    ! If we used a CAS give clear indication in filename.
                    if (sys%CAS(1) == -1 .and. sys%CAS(2) == -1) then
                        write (filename, "(a,a)") trim(sys%read_in%fcidump), ".H5"
                    else
                        write (filename, '(a,"-CAS",i0,",",i0,".H5")') trim(sys%read_in%fcidump), sys%CAS
                    end if
                end if
            else
                ! If we're already here the filename ends in .H5 so no changes needed.
                filename = trim(sys%read_in%fcidump)
            end if
#else
            filename = ''
#endif

        end subroutine get_filename

        subroutine dump_system_hdf5(sys, filename)

            ! Produces HDF5 file containing system information generated from intdump.

            ! In:
            !   sys: derived type containing all information required to set up system
            !       again in binary form.
            ! In/Out:
            !   filename (optional): name of HDF5 to write system information to.  The input
            !       value is used if not an empty string and a name is created based upon
            !       sys%read_in%fcidump and (if relevant) the CAS information otherwise.

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_write
            use calc, only: GLOBAL_META
#else
            use parallel, only: parent
#endif
            use const
            use errors, only: stop_all, warning
            use system, only: sys_t

            type(sys_t), intent(in) :: sys
            character(255), intent(inout) :: filename

#ifndef DISABLE_HDF5

            integer :: date_time(8)
            character(19) :: date_str
            integer :: ierr

            ! HDF5 kinds
            type(hdf5_kinds_t) :: kinds
            ! HDF5 handles
            integer(hid_t) :: file_id, group_id, subgroup_id, subsubgroup_id

            ! Initialise HDF5 and open file.
            call h5open_f(ierr)

            call init_system_hdf5(.true., sys, filename, kinds)

            ! NOTE: if file exists (ie user requested we re-use an existing file), then it is overwritten.
            call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)

            ! --- metadata group ---
            call h5gcreate_f(file_id, gmetadata, group_id, ierr)

                call hdf5_write(group_id, dhande, GLOBAL_META%git_sha1)

                call hdf5_write(group_id, duuid, GLOBAL_META%uuid)

                ! Print out current time and date as HH:MM:SS DD/MM/YYYY.
                call date_and_time(values=date_time)
                write (date_str,'(2(i0.2,":"),i0.2,1X,2(i0.2,"/"),i4)') date_time(5:7), date_time(3:1:-1)
                call hdf5_write(group_id, ddate, date_str)

                call hdf5_write(group_id, dsysdump, sysdump_version)

            call h5gclose_f(group_id, ierr)

            ! --- sys group ---
            call h5gcreate_f(file_id, gsys, group_id, ierr)

            call hdf5_write(group_id, dsystem, sys%system)

            ! Write either original number of electrons or modified number of electrons
            ! depending upon whether using a CAS.
            call hdf5_write(group_id, dnel, sys%nel)

            call hdf5_write(group_id, dms, sys%Ms)

            call hdf5_write(group_id, dcas, kinds, [2_int_64], sys%CAS)

                ! --- basis subgroup ---
                call h5gcreate_f(group_id, gbasis, subgroup_id, ierr)

                associate(nbasis=>sys%basis%nbasis)
                    call hdf5_write(subgroup_id, dnbasis, sys%basis%nbasis)
                    call hdf5_write(subgroup_id, dbasis_spat_ind, kinds, [int(nbasis, kind=int_64)],&
                                sys%basis%basis_fns(:)%spatial_index)
                    call hdf5_write(subgroup_id, dbasis_sym, kinds, [int(nbasis, kind=int_64)],&
                                sys%basis%basis_fns(:)%sym)
                    call hdf5_write(subgroup_id, dbasis_sym_index, kinds, [int(nbasis, kind=int_64)],&
                                sys%basis%basis_fns(:)%sym_index)
                    call hdf5_write(subgroup_id, dbasis_sym_spin_index, kinds,&
                                [int(nbasis, kind=int_64)],  sys%basis%basis_fns(:)%sym_spin_index)
                    call hdf5_write(subgroup_id, dbasis_ms, kinds, [int(nbasis, kind=int_64)],&
                                sys%basis%basis_fns(:)%ms)
                    call hdf5_write(subgroup_id, dbasis_lz, kinds, [int(nbasis, kind=int_64)],&
                                sys%basis%basis_fns(:)%lz)
                    call hdf5_write(subgroup_id, dbasis_sp_eigv, kinds, [int(nbasis, kind=int_64)],&
                                sys%basis%basis_fns(:)%sp_eigv)

                end associate

                call h5gclose_f(subgroup_id, ierr)

                ! --- read_in subgroup ---
                call h5gcreate_f(group_id, gread_in, subgroup_id, ierr)

                call hdf5_write(subgroup_id, dfcidump, sys%read_in%fcidump)
                call hdf5_write(subgroup_id, duhf, sys%read_in%uhf)
                call hdf5_write(subgroup_id, decore, kinds, [1_int_64], &
                                            [sys%read_in%Ecore])
                call hdf5_write(subgroup_id, duselz, sys%read_in%uselz)
                call hdf5_write(subgroup_id, dcomp, sys%read_in%comp)

                ! Need to pass this value to be able to reinitiate symmetry
                call hdf5_write(subgroup_id, dpg_mask, sys%read_in%pg_sym%pg_mask)

                    ! --- integrals subsubgroup ---
                    ! Write each spin/sym block separately, use interface functions for
                    ! clarity of code. See below for info on functionality.
                    call h5gcreate_f(subgroup_id, gintegrals, subsubgroup_id, ierr)

                    call write_1body_integrals(subsubgroup_id, done_body, kinds, &
                            sys%read_in%pg_sym%nbasis_sym_spin, &
                            sys%read_in%one_e_h_integrals%integrals)
                    call write_coulomb_integrals(subsubgroup_id, dcoulomb_ints, kinds, &
                            sys%read_in%coulomb_integrals%integrals)
                    if (sys%read_in%comp) then
                        call write_1body_integrals(subsubgroup_id, done_body_im, kinds, &
                                sys%read_in%pg_sym%nbasis_sym_spin, &
                                sys%read_in%one_e_h_integrals_imag%integrals)
                        call write_coulomb_integrals(subsubgroup_id, dcoulomb_ints_im, kinds, &
                                sys%read_in%coulomb_integrals_imag%integrals)
                    end if

                    call h5gclose_f(subsubgroup_id, ierr)

                call h5gclose_f(subgroup_id, ierr)

            call h5gclose_f(group_id, ierr)

            ! Terminate HDF5
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
#else
            filename = ''
            if (parent)  call warning('dump_system_hdf5', '# Not compiled with HDF5 support. Cannot write out &
                                    sysdump file.')
#endif
        end subroutine dump_system_hdf5

        subroutine read_system_hdf5(sys, verbose)

            ! Read system data from HDF5 file.

            ! In:
            !   sys: system object containing only common sys options
            !       (as in set_common_sys_options) and fcidump.
            !   verbose (optional): write output.

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_read, hdf5_path
#endif
            use parallel

            use const
            use errors, only: stop_all, warning
            use system, only: sys_t
            use point_group_symmetry, only: init_pg_symmetry, print_pg_symmetry_info
            use checking, only: check_allocate, check_deallocate
            use basis, only: write_basis_fn_header, write_basis_fn, write_basis_fn_title
            use basis_types, only: init_basis_strings, print_basis_metadata
            use determinants, only: init_determinants
            use excitations, only: init_excitations
            use read_in_system, only: read_in_one_body
            use molecular_integrals, only: init_one_body_t, init_two_body_t, broadcast_one_body_t, broadcast_two_body_t

            type(sys_t), intent(inout) :: sys
            logical, optional, intent(in) :: verbose

#ifndef DISABLE_HDF5
            character(255) :: filename
            integer :: ierr, sysdump_dump_version, cas(2)

            ! HDF5 kinds
            type(hdf5_kinds_t) :: kinds
            ! HDF5 handles
            integer(hid_t) :: file_id, group_id, subgroup_id, subsubgroup_id

            ! Needed to enable use of current hdf5_helper functions for
            ! reading reals.
            real(p) :: ecore(1)

            integer :: i, nel, ms
            logical :: verbose_t
            integer, allocatable :: sp_fcidump_rank(:)

            verbose_t = .true.
            if (present(verbose)) verbose_t = verbose

            ! A little care must be taken with reading in the file to accomodate
            ! (e.g.) shared memory allocators.  We first read in the data to be
            ! broadcast to all processors, then allocate the (potentially
            ! shared) integral arrays and then read and broadcast the integrals.
            ! All reading is done on the parent processor but all processors
            ! must take part in (e.g.) MPI-3 shared memory allocation together,
            ! hence this ordering.
            if (parent) then
                ! Initialise HDF5 and open file.
                call h5open_f(ierr)
                call init_system_hdf5(.false., sys, filename, kinds)

                call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)
                if (ierr /= 0) then
                    call stop_all('read_system_hdf5', "Unable to open system file.")
                endif

                ! --- metadata group ---

                call hdf5_read(file_id, hdf5_path(gmetadata, dsysdump), sysdump_dump_version)

                if (sysdump_dump_version /= sysdump_version) then
                    call stop_all('read_system_hdf5', "Reading between different sysdump versions not supported.")
                end if

                ! --- sys group ---
                call h5gopen_f(file_id, gsys, group_id, ierr)

                call hdf5_read(group_id, dsystem, sys%system)

                call hdf5_read(group_id, dnel, nel)
                call hdf5_read(group_id, dms, ms)

                if ((sys%nel == 0) .and. sys%Ms == huge(1)) then
                    sys%nel = nel
                    sys%Ms = ms
                else if (sys%nel == 0) then
                    call stop_all('read_system_hdf5', "Specified total spin but not number &
                            &electrons. Please specify neither to use values contained within &
                            & system dump or both.")
                else if (sys%Ms == huge(1)) then
                    call stop_all('read_system_hdf5', "Specified number electrons but not total &
                            &spin. Please specify neither to use values contained within &
                            &system dump or both.")

                else if (nel /= sys%nel .or. ms /= sys%Ms) then
                    call warning('read_system_hdf5', "# WARNING: Overwriting nel and Ms within &
                                        &hdf5 file.")
                    write (*,*)
                end if

                ! --- basis subgroup ---
                call h5gopen_f(group_id, gbasis, subgroup_id, ierr)

                    associate(nbasis=>sys%basis%nbasis)

                        call hdf5_read(subgroup_id, dnbasis, nbasis)

                        allocate(sys%basis%basis_fns(nbasis),  stat=ierr)
                        call check_allocate('sys%basis%basis_fns', nbasis, ierr)

                        call hdf5_read(subgroup_id, dbasis_spat_ind, kinds, [int(nbasis, kind=int_64)],&
                                    sys%basis%basis_fns(:)%spatial_index)
                        call hdf5_read(subgroup_id, dbasis_sym, kinds, [int(nbasis, kind=int_64)],&
                                    sys%basis%basis_fns(:)%sym)
                        call hdf5_read(subgroup_id, dbasis_sym_index, kinds, [int(nbasis, kind=int_64)],&
                                    sys%basis%basis_fns(:)%sym_index)
                        call hdf5_read(subgroup_id, dbasis_sym_spin_index, kinds,&
                                    [int(nbasis, kind=int_64)],  sys%basis%basis_fns(:)%sym_spin_index)
                        call hdf5_read(subgroup_id, dbasis_ms, kinds, [int(nbasis, kind=int_64)],&
                                    sys%basis%basis_fns(:)%ms)
                        call hdf5_read(subgroup_id, dbasis_lz, kinds, [int(nbasis, kind=int_64)],&
                                    sys%basis%basis_fns(:)%lz)
                        call hdf5_read(subgroup_id, dbasis_sp_eigv, kinds, [int(nbasis, kind=int_64)],&
                                    sys%basis%basis_fns(:)%sp_eigv)
                    end associate

                call h5gclose_f(subgroup_id, ierr)

                sys%nvirt = sys%basis%nbasis - sys%nel

                call hdf5_read(group_id, dcas, kinds, [2_int_64], cas)

                if (sys%CAS(1) /= -1 .or. sys%CAS(2) /= -1) then
                    call stop_all('read_system_hdf5', 'Attempted to start calculation &
                            &using a CAS and HDF5 file initiation. This is not currently &
                            &supported. Use original INTDUMP to generate new &
                            &HDF5 file with given CAS instead.')
                end if
                sys%CAS = cas

                ! --- read_in subgroup ---
                call h5gopen_f(group_id, gread_in, subgroup_id, ierr)

                    call hdf5_read(subgroup_id, duhf, sys%read_in%uhf)
                    ! Workaround reading real value from HDF5.
                    call hdf5_read(subgroup_id, decore, kinds, [1_int_64], &
                                                ecore)
                    sys%read_in%Ecore = ecore(1)
                    call hdf5_read(subgroup_id, duselz, sys%read_in%uselz)
                    call hdf5_read(subgroup_id, dcomp, sys%read_in%comp)
                    ! All the symmetry information we need to initiate.
                    call hdf5_read(subgroup_id, dpg_mask, sys%read_in%pg_sym%pg_mask)

                call h5gclose_f(subgroup_id, ierr)
                call h5gclose_f(group_id, ierr)
            end if

#ifdef PARALLEL
            ! Distribute values needed for initialisation on other processes.
            call MPI_BCast(sys%nel, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
            call MPI_BCast(sys%symmetry, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
            call MPI_BCast(sys%read_in%uhf, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
            call MPI_BCast(sys%basis%nbasis, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
            call MPI_BCast(sys%CAS, 2, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
            ! Broadcast symmetry values.
            call MPI_BCast(sys%read_in%pg_sym%pg_mask, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
            ! Allocate basis on non-parent processors.
            if (.not. parent) then
                    allocate(sys%basis%basis_fns(sys%basis%nbasis),  stat = ierr)
                    call check_allocate('sys%basis%basis_fns', sys%basis%nbasis, ierr)
                    sys%nvirt = sys%basis%nbasis - sys%nel
            end if
            ! Broadcast basis.
            associate(nbasis=>sys%basis%nbasis)
                call MPI_BCast(sys%basis%basis_fns(:)%spatial_index, nbasis, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
                call MPI_BCast(sys%basis%basis_fns(:)%sym, nbasis, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
                call MPI_BCast(sys%basis%basis_fns(:)%sym_index, nbasis, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
                call MPI_BCast(sys%basis%basis_fns(:)%sym_spin_index, nbasis, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
                call MPI_BCast(sys%basis%basis_fns(:)%ms, nbasis, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
                call MPI_BCast(sys%basis%basis_fns(:)%lz, nbasis, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
                call MPI_BCast(sys%basis%basis_fns(:)%sp_eigv, nbasis, MPI_PREAL, root, MPI_COMM_WORLD, ierr)
            end associate

            ! Broadcast read_in parameters.
            call MPI_BCast(sys%read_in%Ecore, 1, MPI_PREAL, root, MPI_COMM_WORLD, ierr)
            call MPI_BCast(sys%read_in%uselz, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
            call MPI_BCast(sys%read_in%comp, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
#endif
            ! Initialise various system parameters on all nodes simultaneously.
            ! Do system initialisation that hasn't been read in, hopefully
            ! in same order as in conventional initialisation. Also write
            ! out read_in info for easy checking to compare to original.
            call init_basis_strings(sys%basis)
            call init_determinants(sys, sys%nel)
            call init_excitations(sys%basis)
            call init_pg_symmetry(sys, .true.)

            call init_one_body_t(sys%read_in%uhf, sys%read_in%pg_sym%gamma_sym, sys%read_in%pg_sym%nbasis_sym_spin, &
                                .false., sys%read_in%one_e_h_integrals)
            call init_two_body_t(sys%read_in%uhf, sys%basis%nbasis, sys%read_in%pg_sym%gamma_sym, sys%read_in%comp, &
                                .false., sys%read_in%coulomb_integrals)
            if (sys%read_in%comp) then
                call init_one_body_t(sys%read_in%uhf, sys%read_in%pg_sym%gamma_sym, &
                                    sys%read_in%pg_sym%nbasis_sym_spin, .true., &
                                    sys%read_in%one_e_h_integrals_imag)
                call init_two_body_t(sys%read_in%uhf, sys%basis%nbasis, sys%read_in%pg_sym%gamma_sym,&
                                     sys%read_in%comp, .true., sys%read_in%coulomb_integrals_imag)
            end if

            if (parent) then
                ! --- sys/read_in/integrals subsubgroup ---
                call h5gopen_f(file_id, hdf5_path(gsys, gread_in, gintegrals), subsubgroup_id, ierr)

                    call read_1body_integrals(subsubgroup_id, done_body, kinds, &
                        sys%read_in%uhf, sys%read_in%pg_sym%nbasis_sym_spin, &
                        sys%read_in%one_e_h_integrals)

                    call read_coulomb_integrals(subsubgroup_id, dcoulomb_ints, &
                        kinds, sys%read_in%coulomb_integrals)

                    if (sys%read_in%comp) then
                        call read_1body_integrals(subsubgroup_id, done_body_im, &
                            kinds, sys%read_in%uhf, sys%read_in%pg_sym%nbasis_sym_spin, &
                            sys%read_in%one_e_h_integrals_imag)

                        call read_coulomb_integrals(subsubgroup_id, dcoulomb_ints_im, &
                            kinds, sys%read_in%coulomb_integrals_imag)
                    end if
                call h5gclose_f(subsubgroup_id, ierr)
                call h5fclose_f(file_id, ierr)
                call h5close_f(ierr)
            end if

            ! Broadcast integrals.
            call broadcast_one_body_t(sys%read_in%one_e_h_integrals, root)
            call broadcast_two_body_t(sys%read_in%coulomb_integrals, root)
            if (sys%read_in%comp) then
                call broadcast_one_body_t(sys%read_in%one_e_h_integrals_imag, root)
                call broadcast_two_body_t(sys%read_in%coulomb_integrals_imag, root)
            end if
            if (parent) then
                if (verbose_t) then
                    call write_basis_fn_header(sys)
                    do i = 1, sys%basis%nbasis
                        call write_basis_fn(sys, sys%basis%basis_fns(i), ind=i, &
                                                new_line=.true.)
                    end do
                    write (6,'(/,1X,a8,f18.12)') 'E_core =', sys%read_in%Ecore
                else
                    call write_basis_fn_title()
                end if

                call print_basis_metadata(sys%basis, sys%nel, .false.)
                call print_pg_symmetry_info(sys)
            end if

            if (sys%read_in%dipole_int_file /= '') then
                if (parent) then
                    if (sys%read_in%uhf) then
                        allocate(sp_fcidump_rank(0:sys%basis%nbasis), stat=ierr)
                        call check_allocate('sp_fcidump_rank', sys%basis%nbasis+1, ierr)
                    else
                        allocate(sp_fcidump_rank(0:sys%basis%nbasis/2), stat=ierr)
                        call check_allocate('sp_fcidump_rank', sys%basis%nbasis/2 + 1, ierr)
                    end if
                    do i = lbound(sp_fcidump_rank, dim=1), ubound(sp_fcidump_rank, dim=1)
                        sp_fcidump_rank(i) = i
                    end do
                end if
                call read_in_one_body(sys%read_in%dipole_int_file, sys%basis%nbasis, sys%basis%basis_fns, &
                                      sys%read_in%pg_sym, sys%read_in%uhf, sp_fcidump_rank, sys%nel - sys%cas(1), &
                                      sys%read_in%one_body_op_integrals, sys%read_in%dipole_core)
                deallocate(sp_fcidump_rank, stat=ierr)
                call check_deallocate('sp_fcidump_rank', ierr)
            end if
#else
            if (parent)  call stop_all('read_system_hdf5', '# Not compiled with HDF5 support. Cannot read in &
                                    &sysdump file.')
#endif
        end subroutine read_system_hdf5

! --- utility functions to aid reading out of specific data structures ---

#ifndef DISABLE_HDF5
        subroutine write_1body_integrals(id, dname, kinds, nbasis_sym_spin, integs)

            ! Writes one-body integral values out in spin & symmetry chunks, as stored.
            ! In:
            !   id: hdf5 group id to write in.
            !   dname: name of dataset values will belong to.
            !   kinds: derived tpe containing HDF5 types which correspond to the
            !       non-standard integer and real kinds used in HANDE.
            !   nbasis_sym_spin: (i,j) gives no. basis functions with spin i, sym j.
            !   integs: integral storage from one_body_t

            use hdf5
            use hdf5_helper, only: hdf5_write, hdf5_kinds_t
            use base_types, only: alloc_rp1d
            use const, only: int_64

            character(*), intent(in) :: dname
            type(alloc_rp1d) :: integs(:,:)
            type(hdf5_kinds_t), intent(in) :: kinds
            integer(hid_t), intent(in) :: id
            integer, allocatable, intent(in) :: nbasis_sym_spin(:,:)

            integer :: shpe(2), ispin, isym

            character(155) :: dentr_name

            shpe = shape(integs)
            do ispin = 1, shpe(1)
                do isym = lbound(nbasis_sym_spin, dim=2), ubound(nbasis_sym_spin, dim=2)
                    call get_onebody_name(dname, ispin, isym, dentr_name)
                    call hdf5_write(id, dentr_name, kinds, shape(integs(ispin, isym + 1)%v, kind=int_64), integs(ispin, isym + 1)%v)
                end do
            end do

        end subroutine write_1body_integrals

        subroutine write_coulomb_integrals(id, dname, kinds, integs)

            ! Writes Coulomb integral value out in spin chunks, as stored.
            ! In:
            !   id: hdf5 group id to write in.
            !   dname: name of dataset values will belong to.
            !   kinds: derived tpe containing HDF5 types which correspond to the
            !       non-standard integer and real kinds used in HANDE.
            !   integs: integral storage from one_body_t.

            use hdf5
            use hdf5_helper, only: hdf5_write, hdf5_kinds_t
            use base_types, only: alloc_rp1d
            use const, only: int_64

            character(*), intent(in) :: dname

            type(alloc_rp1d) :: integs(:)
            type(hdf5_kinds_t), intent(in) :: kinds
            integer(hid_t), intent(in) :: id
            integer :: shpe(1), ispin
            character(155) :: dentr_name

            shpe = shape(integs)
            do ispin = 1, shpe(1)
                call get_coulomb_name(dname, ispin, dentr_name)
                call hdf5_write(id, dentr_name, kinds, shape(integs(ispin)%v, kind=int_64), integs(ispin)%v)
            end do

        end subroutine write_coulomb_integrals

        subroutine read_1body_integrals(id, dname, kinds, uhf, nbasis_sym_spin, store)

            ! Reads one body integrals from hdf5 previously output by hande.

            ! In:
            !   id: hdf5 group id to write in.
            !   dname: name of dataset values will belong to.
            !   kinds: derived tpe containing HDF5 types which correspond to the
            !       non-standard integer and real kinds used in HANDE.
            !   uhf: bool, sets whether integrals results from uhf calculation.
            !   nbasis_sym_spin: (i,j) gives no. basis functions with spin i, sym j.
            ! In/Out:
            !   store: Fully allocated one_body_t passed in, returned with all
            !       appropriate integral values stored.

            ! NB. must be called after allocating integral arrays within store.

            use hdf5
            use hdf5_helper, only: hdf5_read, hdf5_kinds_t
            use molecular_integral_types, only: one_body_t
            use checking, only: check_allocate
            use const, only: int_64

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dname
            type(hdf5_kinds_t), intent(in) :: kinds
            logical, intent(in) :: uhf
            integer, allocatable, intent(in) :: nbasis_sym_spin(:,:)
            type(one_body_t), intent(inout) :: store

            integer :: ispin, isym, nspin, dummy(2)
            integer(int_64) :: s(1)
            character(155) :: dentr_name

            if (uhf) then
                nspin = 2
            else
                nspin = 1
            end if

            do ispin = 1, nspin
                ! Workaround behaviour of lbound/ubound on allocatable arrays
                ! (sometimes 0-indexed, sometimes not)
                dummy = shape(nbasis_sym_spin)
                do isym = 0, dummy(2) - 1
                    call get_onebody_name(dname, ispin, isym, dentr_name)
                    s(1) = int(nbasis_sym_spin(ispin, isym), kind=int_64) * &
                           (int(nbasis_sym_spin(ispin, isym), kind=int_64) + 1_int_64) / 2_int_64
                    call hdf5_read(id, dentr_name, kinds, s, store%integrals(ispin, isym)%v)
                end do
            end do

        end subroutine read_1body_integrals

        subroutine read_coulomb_integrals(id, dname, kinds, store)
            ! Reads coulomb integrals from hdf5 previously output by hande.

            ! In:
            !   id: hdf5 group id to write in.
            !   dname: name of dataset values will belong to.
            !   kinds: derived tpe containing HDF5 types which correspond to the
            !       non-standard integer and real kinds used in HANDE.
            ! In/Out:
            !   store: Fully allocated two_body_t passed in, returned with all
            !       appropriate integral values stored.

            ! NB. must be called after allocating integral arrays within store.

            use hdf5
            use hdf5_helper, only: hdf5_read, hdf5_kinds_t
            use molecular_integral_types, only: two_body_t
            use const, only: int_64

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dname
            type(hdf5_kinds_t), intent(in) :: kinds
            type(two_body_t), intent(inout) :: store

            integer :: ispin, shpe(1)
            character(155) :: dentr_name

            shpe = shape(store%integrals)

            do ispin = 1, shpe(1)
                call get_coulomb_name(dname, ispin, dentr_name)
                call hdf5_read(id, dentr_name, kinds, shape(store%integrals(ispin)%v, kind=int_64), &
                                store%integrals(ispin)%v)
            end do

        end subroutine read_coulomb_integrals

        subroutine get_onebody_name(dname, ispin, isym, entry_name)

            ! Generates name for set parameters when storing one-body integral chunks.
            ! Name follows pattern "dname_ispin[ispin]_isym[isym]"
            ! In:
            !   dname: base name for parameter group (differs between real/imaginary components).
            !   ispin: index of spin chunk.
            !   isym: index of symmetry chunk.
            ! Out:
            !   entry_name: name for storage of information under.

            character(*), intent(in) :: dname
            integer, intent(in) :: ispin, isym
            character(155), intent(out) :: entry_name

            write (entry_name, "(a,a,i2.2,a,i2.2)") trim(dname), "_ispin", ispin, "_isym", isym

        end subroutine get_onebody_name

        subroutine get_coulomb_name(dname, ispin, entry_name)

            ! Generates name for set parameters when storing coulomb integral chunks.
            ! Name follows pattern "dname_ispin[ispin]"
            ! In:
            !   dname: base name for parameter group (differs between real/imaginary components).
            !   ispin: index of spin chunk.
            ! Out:
            !   entry_name: name for storage of information under.

            character(*), intent(in) :: dname
            integer, intent(in) :: ispin
            character(155), intent(out) :: entry_name

            write (entry_name, "(a,a,i2.2)") trim(dname), "_ispin", ispin

        end subroutine get_coulomb_name

#endif
end module hdf5_system
