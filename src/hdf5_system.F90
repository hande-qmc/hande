module hdf5_system

    implicit none



    integer, parameter :: sysdump_version = 0

    ! Group names.
    character(*), parameter ::  gmetadata = 'metadata',     &
                                gsys = 'system',            &
                                gbasis = 'basis',           &
                                gread_in = 'read_in',       &
                                gintegrals = 'integrals',   &
                                gpg_sym = 'pg_sym'


    ! Dataspace names.
    character(*), parameter ::  dsysdump = 'sysdump version',   &
                                dhande = 'hande version',       &
                                ddate = 'date',                 &
                                duuid = 'uuid',                 &

                                dsystem = 'system',             &
                                dsymmetry = 'symmetry',         &
                                dnsym = 'nsym',                 &
                                dsym0 = 'sym0',                 &
                                dsym_max = 'sym_max',           &
                                dnsym_tot = 'nsym_tot',         &
                                dsym0_tot = 'sym0_tot',         &
                                dsym_max_tot = 'sym_max_tot',   &
                                dcas = 'CAS',                   &
                                dmax_number_excitations =       &
                                       'max_number_excitations',&

                                dfcidump = 'fcidump',           &
                                duhf = 'uhf',                   &
                                decore = 'ecore',               &
                                duselz = 'uselz',               &
                                dintegrals = 'integrals',       &
                                dop_sym = 'op_sym',             &
                                dimag = 'imag',                 &
                                dcomp = 'comp',                 &

                                done_body = 'one_body',         &
                                done_body_im = 'one_body_im',   &
                                dcoulomb_ints = 'coulomb_ints', &
                                dcoulomb_ints_im =              &
                                              'coulomb_ints_im',&

                                dgamma_sym = 'gamma_sym',       &
                                dnbasis_sym = 'nbasis_sym',     &
                                dnbasis_sym_spin =              &
                                              'nbasis_sym_spin',&
                                dsym_spin_basis_fns =           &
                                           'sym_spin_basis_fns',&
                                dpg_mask = 'pg_mask',           &
                                dlz_mask = 'Lz_mask',           &
                                dlz_offset = 'lz_offset',       &
                                dlz_divisor = 'lz_divisor',     &

                                dnbasis = 'nbasis',             &
                                dstring_len = 'string_len',     &
                                dtensor_label_len =             &
                                             'tensor_label_len',&
                                dalpha_mask = 'alpha_mask',     &
                                dbeta_mask = 'beta_mask',       &
                                dbit_lookup = 'bit_lookup',     &
                                dbasis_lookup = 'basis_lookup', &
                                dexcit_mask = 'excit_mask',     &

                                dbasis_spat_ind =               &
                                        'basis_spatial_index',  &
                                dbasis_sym = 'basis_symmetry',  &
                                dbasis_sym_index =                &
                                        'basis_symmetry_index', &
                                dbasis_sym_spin_index =           &
                                        'basis_symmetry_spin_index',&
                                dbasis_ms = 'basis_ms',         &
                                dbasis_lz = 'basis_lz',         &
                                dbasis_sp_eigv = 'basis_sp_eigv'

    contains

#ifndef DISABLE_HDF5
        subroutine init_system_hdf5(write_mode, sys, filename, kinds, verbose)


            use hdf5_helper, only: hdf5_kinds_t, hdf5_kinds_init
            use utils, only: get_unique_filename
            use system, only: sys_t

            logical, intent(in) :: write_mode
            type(sys_t), intent(in) :: sys

            type(hdf5_kinds_t), intent(out), optional :: kinds
            logical, intent(in), optional ::  verbose

            character(255), intent(out) :: filename
            integer :: fname_id

            logical :: verbose_loc, exists


            if (present(verbose)) verbose_loc = verbose


            if (write_mode) then
                if (sys%CAS(1) == -1 .and. sys%CAS(2) == -1) then
                    write (filename, "(a,a)") trim(sys%read_in%fcidump), ".H5"
                else
                    write (filename, "(a,a,i2,a,i2,a)") trim(sys%read_in%fcidump), "-", sys%CAS(1), ",", sys%CAS(2), "CAS.H5"
                end if
            else
                filename = trim(sys%read_in%fcidump)
            end if



            if (present(kinds)) call hdf5_kinds_init(kinds)

        end subroutine init_system_hdf5

#endif


        subroutine dump_system_hdf5(sys)

            ! Produces HDF5 file containing system information generated from intdump.

            ! In:



            ! Out:

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

#ifndef DISABLE_HDF5
            character(255) :: integral_filename


            integer :: date_time(8)
            character(19) :: date_str
            integer :: ierr, nbasis

            ! HDF5 kinds
            type(hdf5_kinds_t) :: kinds
            ! HDF5 handles
            integer(hid_t) :: file_id, group_id, subgroup_id, subsubgroup_id

            ! Initialise HDF5 and open file.
            call h5open_f(ierr)
            call init_system_hdf5(.true., sys, integral_filename, kinds)

            call h5fcreate_f(integral_filename, H5F_ACC_TRUNC_F, file_id, ierr)

            ! --- metadata group ---
            call h5gcreate_f(file_id, gmetadata, group_id, ierr)

                call hdf5_write(group_id, dhande, GLOBAL_META%git_sha1)

                call hdf5_write(group_id, duuid, GLOBAL_META%uuid)

                ! Print out current time and date as HH:MM:SS DD/MM/YYYY.
                !write (date_str,'(2(i0.2,":"),i0.2,1X,2(i0.2,"/"),i4)') date_time(5:7), date_time(3:1:-1)
                !call hdf5_write(group_id, ddate, date_str)

                call hdf5_write(group_id, dsysdump, sysdump_version)

            call h5gclose_f(group_id, ierr)

            ! --- sys group ---
            call h5gcreate_f(file_id, gsys, group_id, ierr)

            call hdf5_write(group_id, dsystem, sys%system)

            call hdf5_write(group_id, dsymmetry, sys%symmetry)
            call hdf5_write(group_id, dnsym, sys%nsym)
            call hdf5_write(group_id, dsym0, sys%sym0)
            call hdf5_write(group_id, dsym_max, sys%sym_max)

            call hdf5_write(group_id, dnsym_tot, sys%nsym_tot)
            call hdf5_write(group_id, dsym0_tot, sys%sym0_tot)
            call hdf5_write(group_id, dsym_max_tot, sys%sym_max_tot)


            call hdf5_write(group_id, dcas, kinds, [2], sys%CAS)

            call hdf5_write(group_id, dmax_number_excitations,&
                                            sys%max_number_excitations)



                ! --- basis subgroup ---
                call h5gcreate_f(group_id, gbasis, subgroup_id, ierr)

                nbasis = sys%basis%nbasis

                call hdf5_write(subgroup_id, dnbasis, sys%basis%nbasis)
                call hdf5_write(subgroup_id, dbasis_spat_ind, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%spatial_index)
                call hdf5_write(subgroup_id, dbasis_sym, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%sym)
                call hdf5_write(subgroup_id, dbasis_sym_index, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%sym_index)
                call hdf5_write(subgroup_id, dbasis_sym_spin_index, kinds,&
                            [nbasis],  sys%basis%basis_fns(:)%sym_spin_index)
                call hdf5_write(subgroup_id, dbasis_ms, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%ms)
                call hdf5_write(subgroup_id, dbasis_lz, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%lz)
                call hdf5_write(subgroup_id, dbasis_sp_eigv, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%sp_eigv)


                call h5gclose_f(subgroup_id, ierr)

                ! --- read_in subgroup ---
                call h5gcreate_f(group_id, gread_in, subgroup_id, ierr)

                call hdf5_write(subgroup_id, dfcidump, sys%read_in%fcidump)
                call hdf5_write(subgroup_id, duhf, sys%read_in%uhf)
                call hdf5_write(subgroup_id, decore, kinds, [1], &
                                            [sys%read_in%Ecore])
                call hdf5_write(subgroup_id, duselz, sys%read_in%uselz)
                call hdf5_write(subgroup_id, dcomp, sys%read_in%comp)

                    ! --- integrals subsubgroup ---
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

                    ! Currently removed; can regenerate fairly painlessly with basis info.
                    ! --- pg sym subsubgroup ---
                    !call h5gcreate_f(subgroup_id, gpg_sym, subsubgroup_id, ierr)
                    !call hdf5_write(subsubgroup_id, dgamma_sym, &
                    !                sys%read_in%pg_sym%gamma_sym)
                    !call hdf5_write(subsubgroup_id, dnbasis_sym, kinds, &
                    !                shape(sys%read_in%pg_sym%nbasis_sym), &
                    !                sys%read_in%pg_sym%nbasis_sym)
                    !call hdf5_write(subsubgroup_id, dnbasis_sym_spin, kinds, &
                    !                shape(sys%read_in%pg_sym%nbasis_sym_spin), &
                    !                sys%read_in%pg_sym%nbasis_sym_spin)

                    !call hdf5_write(subsubgroup_id, dsym_spin_basis_fns, kinds, &
                    !                shape(sys%read_in%pg_sym%sym_spin_basis_fns), &
                    !                sys%read_in%pg_sym%sym_spin_basis_fns)
                    !call hdf5_write(subsubgroup_id, dpg_mask, &
                    !                sys%read_in%pg_sym%pg_mask)
                    !call hdf5_write(subsubgroup_id, dlz_mask, &
                    !                sys%read_in%pg_sym%lz_mask)
                    !call hdf5_write(subsubgroup_id, dlz_offset, &
                    !                sys%read_in%pg_sym%lz_offset)
                    !call hdf5_write(subsubgroup_id, dlz_divisor, &
                    !                sys%read_in%pg_sym%lz_divisor)
                    !call h5gclose_f(subsubgroup_id, ierr)
                call h5gclose_f(subgroup_id, ierr)

            call h5gclose_f(group_id, ierr)

            ! Terminate HDF5
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)

#else
            if (parent)  call warning('dump_system_hdf5', '# Not compiled with HDF5 support. Cannot write out &
                                    sysdump file.')
#endif
        end subroutine dump_system_hdf5

        subroutine read_system_hdf5(filename, sys, verbose)

#ifndef DISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_read
#endif
            use parallel, only: parent

            use const
            use errors, only: stop_all, warning
            use system, only: sys_t
            use point_group_symmetry, only: init_pg_symmetry, print_pg_symmetry_info
            use checking, only: check_allocate
            use basis, only: write_basis_fn_header, write_basis_fn, write_basis_fn_title
            use basis_types, only: init_basis_strings, print_basis_metadata
            use determinants, only: init_determinants
            use excitations, only: init_excitations


            type(sys_t), intent(inout) :: sys
            logical, optional, intent(in) :: verbose

#ifndef DISABLE_HDF5

            character(255), intent(out) :: filename

            integer :: ierr, nbasis, sysdump_dump_version, cas(2)
            type(hdf5_kinds_t) :: kinds
            integer(hid_t) :: file_id, group_id, subgroup_id, subsubgroup_id
            real(p) :: ecore(1)

            integer :: i
            logical :: exists, verbose_t

            verbose_t = .true.
            if (present(verbose)) verbose_t = verbose


            ! Initialise HDF5 and open file.
            call h5open_f(ierr)
            call init_system_hdf5(.false., sys, filename, kinds)

            call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)
            if (ierr /= 0) then
                call stop_all('read_system_hdf5', "Unable to open restart file.")
            endif

            ! --- metadata group ---
            call h5gopen_f(file_id, gmetadata, group_id, ierr)

                call hdf5_read(group_id, dsysdump, sysdump_dump_version)

                if (sysdump_dump_version /= sysdump_version) then
                    call stop_all('read_system_hdf5', "Reading between different &
                        &sysdump versions not supported.")
                end if
            call h5gclose_f(group_id, ierr)

            ! --- sys group ---
            call h5gopen_f(file_id, gsys, group_id, ierr)

            call hdf5_read(group_id, dsystem, sys%system)

            call hdf5_read(group_id, dsymmetry, sys%symmetry)
            call hdf5_read(group_id, dnsym, sys%nsym)
            call hdf5_read(group_id, dsym0, sys%sym0)
            call hdf5_read(group_id, dsym_max, sys%sym_max)

            call hdf5_read(group_id, dnsym_tot, sys%nsym_tot)
            call hdf5_read(group_id, dsym0_tot, sys%sym0_tot)
            call hdf5_read(group_id, dsym_max_tot, sys%sym_max_tot)

            call hdf5_read(group_id, dmax_number_excitations,&
                                            sys%max_number_excitations)

                ! --- basis subgroup ---
                call h5gopen_f(group_id, gbasis, subgroup_id, ierr)

                call hdf5_read(subgroup_id, dnbasis, nbasis)
                sys%basis%nbasis = nbasis

                allocate(sys%basis%basis_fns(nbasis),  stat = ierr)
                call check_allocate('sys%basis%basis_fns', nbasis, ierr)

                call hdf5_read(subgroup_id, dbasis_spat_ind, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%spatial_index)
                call hdf5_read(subgroup_id, dbasis_sym, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%sym)
                call hdf5_read(subgroup_id, dbasis_sym_index, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%sym_index)
                call hdf5_read(subgroup_id, dbasis_sym_spin_index, kinds,&
                            [nbasis],  sys%basis%basis_fns(:)%sym_spin_index)
                call hdf5_read(subgroup_id, dbasis_ms, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%ms)
                call hdf5_read(subgroup_id, dbasis_lz, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%lz)
                call hdf5_read(subgroup_id, dbasis_sp_eigv, kinds, [nbasis],&
                            sys%basis%basis_fns(:)%sp_eigv)

                call h5gclose_f(subgroup_id, ierr)

            sys%nvirt = sys%basis%nbasis - sys%nel

            call hdf5_read(group_id, dcas, kinds, [2], cas)

            if ((sys%CAS(1) /= -1 .or. sys%CAS(2) /= -1) .and. (sys%CAS(1) /= cas(1) &
                                            .or. sys%CAS(2) /= cas(2))) then
                call stop_all('read_system_hdf5', 'attempting to start calculation &
                        &with different CAS to that used in HDF5 file creation; not &
                        &currently supported. Use original INTDUMP to generate new &
                        &HDF5 file.')
            end if
            sys%CAS = cas


                ! --- read in subgroup ---
                call h5gopen_f(group_id, gread_in, subgroup_id, ierr)

                call hdf5_read(subgroup_id, duhf, sys%read_in%uhf)
                call hdf5_read(subgroup_id, decore, kinds, [1], &
                                            ecore)
                sys%read_in%Ecore = ecore(1)
                call hdf5_read(subgroup_id, duselz, sys%read_in%uselz)
                call hdf5_read(subgroup_id, dcomp, sys%read_in%comp)

                call init_basis_strings(sys%basis)
                call init_determinants(sys, sys%nel)
                call init_excitations(sys%basis)


                call init_pg_symmetry(sys)

                if (parent .and. verbose_t) then
                    call write_basis_fn_header(sys)
                    do i = 1, sys%basis%nbasis
                        call write_basis_fn(sys, sys%basis%basis_fns(i), ind=i, &
                                                new_line=.true.)
                    end do
                    write (6,'(/,1X,a8,f18.12)') 'E_core =', sys%read_in%Ecore
                else if (parent) then
                    call write_basis_fn_title()
                end if

                call print_basis_metadata(sys%basis, sys%nel, .false.)
                call print_pg_symmetry_info(sys)

                    ! ---integrals subsubgroup ---
                    call h5gopen_f(subgroup_id, gintegrals, subsubgroup_id, ierr)

                    call read_1body_integrals(subsubgroup_id, done_body, kinds, &
                        sys%read_in%uhf, sys%read_in%pg_sym%gamma_sym, &
                        sys%read_in%pg_sym%nbasis_sym_spin, .false., &
                        sys%read_in%one_e_h_integrals)

                    call read_coulomb_integrals(subsubgroup_id, dcoulomb_ints, &
                        kinds, sys%read_in%uhf, sys%basis%nbasis, &
                        sys%read_in%pg_sym%gamma_sym, sys%read_in%comp, .false., &
                        sys%read_in%coulomb_integrals)

                    if (sys%read_in%comp) then
                        call read_1body_integrals(subsubgroup_id, done_body_im, &
                            kinds, sys%read_in%uhf, sys%read_in%pg_sym%gamma_sym, &
                            sys%read_in%pg_sym%nbasis_sym_spin, .true., &
                            sys%read_in%one_e_h_integrals_imag)

                        call read_coulomb_integrals(subsubgroup_id, dcoulomb_ints_im, &
                            kinds, sys%read_in%uhf, sys%basis%nbasis, &
                            sys%read_in%pg_sym%gamma_sym, sys%read_in%comp, .true., &
                            sys%read_in%coulomb_integrals_imag)
                    end if
                    call h5gclose_f(subsubgroup_id, ierr)

                call h5gclose_f(subgroup_id, ierr)

            call h5gclose_f(group_id, ierr)

            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)

#else

            if (parent)  call warning('dump_system_hdf5', '# Not compiled with HDF5 support. Cannot write out &
                                    sysdump file.')
#endif
        end subroutine read_system_hdf5

! --- utility functions to aid reading out of specific data structures ---

#ifndef DISABLE_HDF5
        subroutine write_1body_integrals(id, dname, kinds, nbasis_sym_spin, integs)
            use hdf5
            use hdf5_helper, only: hdf5_write, hdf5_kinds_t
            use base_types, only: alloc_rp1d
            character(*), intent(in) :: dname

            type(alloc_rp1d) :: integs(:,:)
            type(hdf5_kinds_t), intent(in) :: kinds
            integer(hid_t), intent(in) :: id
            integer, allocatable, intent(in) :: nbasis_sym_spin(:,:)

            integer :: shpe(2), i, j

            character(155) :: dentr_name

            shpe = shape(integs)
            do i = 1, shpe(1)
                do j = lbound(nbasis_sym_spin, dim=2), ubound(nbasis_sym_spin, dim=2)
                    write (dentr_name, "(a,a,i2.2,a,i2.2)") trim(dname), "_ispin", i, "_isym", j
                    call hdf5_write(id, dentr_name, kinds, &
                                shape(integs(i,j + 1)%v), integs(i,j + 1)%v)
                end do
            end do
        end subroutine write_1body_integrals

        subroutine write_coulomb_integrals(id, dname, kinds, integs)
            use hdf5
            use hdf5_helper, only: hdf5_write, hdf5_kinds_t
            use base_types, only: alloc_rp1d
            character(*), intent(in) :: dname

            type(alloc_rp1d) :: integs(:)
            type(hdf5_kinds_t), intent(in) :: kinds
            integer(hid_t), intent(in) :: id

            integer :: shpe(1), i, j

            character(155) :: dentr_name

            shpe = shape(integs)
            do i = 1, shpe(1)
                write (dentr_name, "(a,a,i2.2)") dname, "_ispin", i
                call hdf5_write(id, dentr_name, kinds, &
                            shape(integs(i)%v), integs(i)%v)
            end do

        end subroutine write_coulomb_integrals

        subroutine read_1body_integrals(id, dname, kinds, uhf, op_sym, nbasis_sym_spin, imag, store)
            ! Reads one body integrals from hdf5 previously output by hande.

            ! In:
            !   id:
            !   dname:
            !   kinds:
            !   uhf:
            !   op_sym:
            !   nbasis_sym_spin:
            !   imag:
            ! Out:
            !   store:

            ! NB. should be called after reading all other system information in,
            ! as requires other information (uhf, nbasis_sym_spin) to size info
            ! to be read in.
            use hdf5
            use hdf5_helper, only: hdf5_read, hdf5_kinds_t, dset_shape
            use molecular_integrals, only: init_one_body_t
            use molecular_integral_types, only: one_body_t
            use checking, only: check_allocate

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dname
            type(hdf5_kinds_t), intent(in) :: kinds
            logical, intent(in) :: uhf, imag
            integer, intent(in) :: op_sym
            integer, allocatable, intent(in) :: nbasis_sym_spin(:,:)
            type(one_body_t), intent(out) :: store

            integer :: ispin, isym, s(1), nspin, ierr, dummy(2)
            integer(hsize_t) :: dum(1)
            character(155) :: dentr_name

            if (uhf) then
                nspin = 2
            else
                nspin = 1
            end if

            call init_one_body_t(uhf, op_sym, nbasis_sym_spin, imag, store)

            do ispin = 1, nspin
                dummy = shape(nbasis_sym_spin)
                do isym = 0, dummy(2) - 1
                    write (dentr_name, "(a,a,i2.2,a,i2.2)") trim(dname), "_ispin", ispin, "_isym", isym
                    s(1) = nbasis_sym_spin(ispin, isym) * (nbasis_sym_spin(ispin, isym) + 1) / 2
                    call hdf5_read(id, dentr_name, &
                                kinds, s,&
                                store%integrals(ispin, isym)%v)
                end do
            end do
        end subroutine read_1body_integrals

        subroutine read_coulomb_integrals(id, dname, kinds, uhf, nbasis, &
                                    op_sym, comp, imag, store)

            use hdf5
            use hdf5_helper, only: hdf5_read, hdf5_kinds_t, dset_shape
            use molecular_integrals, only: init_two_body_t
            use molecular_integral_types, only: two_body_t

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dname
            type(hdf5_kinds_t), intent(in) :: kinds
            logical, intent(in) :: uhf, imag, comp
            integer, intent(in) :: op_sym, nbasis
            type(two_body_t), intent(out) :: store

            integer :: ispin, shpe(1)
            integer(hsize_t) :: s(1)
            character(155) :: dentr_name

            call init_two_body_t(uhf, nbasis, op_sym, comp, imag, store)

            shpe = shape(store%integrals)

            do ispin = 1, shpe(1)

                write (dentr_name, "(a,a,i2.2)") trim(dname), "_ispin", ispin
                call dset_shape(id, dentr_name, s)
                call hdf5_read(id, dentr_name, kinds, &
                        shape(store%integrals(ispin)%v), &
                        store%integrals(ispin)%v)
            end do

        end subroutine read_coulomb_integrals
#endif
end module hdf5_system
