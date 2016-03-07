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
                                dnel = 'nel',                   &
                                dnvirt = 'nvirt',               &
                                dMs = 'Ms',                     &
                                dnalpha = 'nalpha',             &
                                dnbeta = 'nbeta',               &
                                dnvirt_alpha = 'nvirtalpha',    &
                                dnvirt_beta = 'nnvirtbeta',     &
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

#ifndef IDISABLE_HDF5
        subroutine init_system_hdf5(write_mode, filename, kinds, verbose)


            use hdf5_helper, only: hdf5_kinds_t, hdf5_kinds_init
            use utils, only: get_unique_filename


            logical, intent(in) :: write_mode

            type(hdf5_kinds_t), intent(out), optional :: kinds
            logical, intent(in), optional ::  verbose

            character(*), intent(out) :: filename
            integer :: fname_id

            logical :: verbose_loc, exists
            integer :: id
            character(4) :: stem
            character(6) :: suffix

            stem = 'STEM'
            suffix = 'SUFFIX'
            id = 1


            if (present(verbose)) verbose_loc = verbose

            call get_unique_filename(stem, suffix, write_mode, id, filename, fname_id)


            if (present(kinds)) call hdf5_kinds_init(kinds)

        end subroutine init_system_hdf5

#endif


        subroutine dump_system_info(sys)

            ! Produces HDF5 file containing system information generated from intdump.

            ! In:



            ! Out:

#ifndef IDISABLE_HDF5
            use hdf5
            use hdf5_helper, only: hdf5_kinds_t, hdf5_write

            use calc, only: GLOBAL_META
#else
            use parallel, only: parent
#endif
            use const
            use errors, only: stop_all, warning
            use checking, only: check_allocate, check_deallocate

            use system, only: sys_t

            type(sys_t), intent(in) :: sys

#ifndef IDISABLE_HDF5
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
            call init_system_hdf5(.true., integral_filename, kinds)

            call h5fcreate_f(integral_filename, H5F_ACC_TRUNC_F, file_id, ierr)

            ! --- metadata group ---
            call h5gcreate_f(file_id, gmetadata, group_id, ierr)

                call hdf5_write(group_id, dhande, GLOBAL_META%git_sha1)

                call hdf5_write(group_id, duuid, GLOBAL_META%uuid)

                ! Print out current time and date as HH:MM:SS DD/MM/YYYY.
                write (date_str,'(2(i0.2,":"),i0.2,1X,2(i0.2,"/"),i4)') date_time(5:7), date_time(3:1:-1)
                call hdf5_write(group_id, ddate, date_str)

                call hdf5_write(group_id, dsysdump, sysdump_version)

            call h5gclose_f(group_id, ierr)

            ! --- sys group ---
            call h5gcreate_f(file_id, gsys, group_id, ierr)

            call hdf5_write(group_id, dsystem, sys%system)

            call hdf5_write(group_id, dnel, sys%nel)

            call hdf5_write(group_id, dnvirt, sys%nvirt)

            call hdf5_write(group_id, dms, sys%Ms)

            call hdf5_write(group_id, dnalpha, sys%nalpha)
            call hdf5_write(group_id, dnbeta, sys%nbeta)

            call hdf5_write(group_id, dnvirt_alpha, sys%nvirt_alpha)
            call hdf5_write(group_id, dnvirt_beta, sys%nvirt_beta)

            call hdf5_write(group_id, dsymmetry, sys%symmetry)

            call hdf5_write(group_id, dnsym, sys%nsym)
            call hdf5_write(group_id, dsym0, sys%sym0)
            call hdf5_write(group_id, dsym_max, sys%sym_max)

            call hdf5_write(group_id, dnsym_tot, sys%nsym_tot)
            call hdf5_write(group_id, dsym0_tot, sys%sym0_tot)
            call hdf5_write(group_id, dsym_max_tot, sys%sym_max_tot)


            call hdf5_write(group_id, dcas, kinds, [nbasis], sys%CAS)

            call hdf5_write(group_id, dmax_number_excitations,&
                                            sys%max_number_excitations)



                ! --- basis subgroup ---
                call h5gcreate_f(group_id, gbasis, subgroup_id, ierr)

                nbasis = sys%basis%nbasis


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
                            sys%read_in%one_e_h_integrals%integrals)
                    call write_coulomb_integrals(subsubgroup_id, dcoulomb_ints, kinds, &
                            sys%read_in%coulomb_integrals%integrals)
                    if (sys%read_in%comp) then
                        call write_1body_integrals(subsubgroup_id, done_body_im, kinds, &
                                sys%read_in%one_e_h_integrals_imag%integrals)
                        call write_coulomb_integrals(subsubgroup_id, dcoulomb_ints_im, kinds, &
                                sys%read_in%coulomb_integrals_imag%integrals)
                    end if

                    call h5gclose_f(subsubgroup_id, ierr)
                    ! --- pg sym subsubgroup ---
                    call h5gcreate_f(subgroup_id, gpg_sym, subsubgroup_id, ierr)
                    call hdf5_write(subsubgroup_id, dgamma_sym, &
                                    sys%read_in%pg_sym%gamma_sym)
                    call hdf5_write(subsubgroup_id, dnbasis_sym, kinds, &
                                    shape(sys%read_in%pg_sym%nbasis_sym), &
                                    sys%read_in%pg_sym%nbasis_sym)
                    call hdf5_write(subsubgroup_id, dnbasis_sym_spin, kinds, &
                                    shape(sys%read_in%pg_sym%nbasis_sym_spin), &
                                    sys%read_in%pg_sym%nbasis_sym_spin)

                    call hdf5_write(subsubgroup_id, dsym_spin_basis_fns, kinds, &
                                    shape(sys%read_in%pg_sym%sym_spin_basis_fns), &
                                    sys%read_in%pg_sym%sym_spin_basis_fns)
                    call hdf5_write(subsubgroup_id, dpg_mask, &
                                    sys%read_in%pg_sym%pg_mask)
                    call hdf5_write(subsubgroup_id, dlz_mask, &
                                    sys%read_in%pg_sym%lz_mask)
                    call hdf5_write(subsubgroup_id, dlz_offset, &
                                    sys%read_in%pg_sym%lz_offset)
                    call hdf5_write(subsubgroup_id, dlz_divisor, &
                                    sys%read_in%pg_sym%lz_divisor)
                    call h5gclose_f(subsubgroup_id, ierr)
                call h5gclose_f(subgroup_id, ierr)

            call h5gclose_f(group_id, ierr)

            ! Terminate HDF5
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)

#else
            if (parent)  call warning('dump_system_hdf5', '# Not compiled with HDF5 support. Cannot write out &
                                    sysdump file.')
#endif
        end subroutine dump_system_info

        subroutine write_1body_integrals(id, dname, kinds, integs)
            use hdf5
            use hdf5_helper, only: hdf5_write, hdf5_kinds_t
            use base_types, only: alloc_rp1d
            character(*), intent(in) :: dname

            type(alloc_rp1d) :: integs(:,:)
            type(hdf5_kinds_t), intent(in) :: kinds
            integer(hid_t), intent(in) :: id

            integer :: shpe(2), i, j

            character(155) :: dentr_name

            shpe = shape(integs)
            do i = 1, shpe(1)
                do j = 1, (shpe(2))
                    write (dentr_name, *) dname, "ispin", i, "isym", j
                    call hdf5_write(id, dentr_name, kinds, &
                                shape(integs(i,j)%v), integs(i,j)%v)
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
                write (dentr_name, *) dname, "ispin", i
                call hdf5_write(id, dentr_name, kinds, &
                            shape(integs(i)%v), integs(i)%v)
            end do

        end subroutine write_coulomb_integrals

end module hdf5_system
