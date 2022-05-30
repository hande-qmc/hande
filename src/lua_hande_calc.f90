module lua_hande_calc

! Lua wrappers to for calculation procedures.  See top-level comments in
! lua_hande for more details about working with the Lua API.

implicit none

contains

    ! --- Calculation wrappers ---

    function lua_fci(L) result(nresult) bind(c)

        ! Run an FCI calculation.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    fci {
        !           sys = system,        -- required
        !           fci = { ... },
        !           davidson = { ... },
        !           reference = { ... },
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_exists, aot_table_close

        use calc, only: calc_type, exact_diag, davidson_diag
        use fci_lapack, only: do_fci_lapack
        use fci_davidson, only: do_fci_davidson
        use fci_utils, only: fci_in_t
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args, register_timing
        use reference_determinant, only: reference_t
        use system, only: sys_t

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state
        integer :: opts
        real :: t1, t2
        type(sys_t), pointer :: sys
        type(fci_in_t) :: fci_in
        type(reference_t) :: ref
        character(12), parameter :: keys(4) = [character(12) :: 'sys', 'fci', 'davidson', 'reference']

        call cpu_time(t1)

        lua_state = flu_copyptr(l)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)
        call read_fci_in(lua_state, opts, sys%basis, fci_in)
        call read_reference_t(lua_state, opts, ref, sys)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        if (fci_in%using_davidson) then
            calc_type = davidson_diag
            call do_fci_davidson(sys, fci_in, ref)
        else
            calc_type = exact_diag
            call do_fci_lapack(sys, fci_in, ref)
        end if

        nresult = 0

        call cpu_time(t2)
        call register_timing(lua_state, "FCI calculation", t2-t1)

    end function lua_fci

    function lua_hilbert_space(L) result(nresult) bind(c)

        ! Run a Monte Carlo calculation to estimate the size of the Hilbert
        ! space.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    hilbert_space {
        !       sys = system,      -- required
        !       hilbert = { ... }, -- required
        !       output = { ... },
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer

        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_exists, aot_table_open, aot_table_close, aot_table_set_val

        use parallel, only: parent
        use errors, only: stop_all
        use const, only: dp
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args, register_timing
        use lua_hande_calc_utils, only: init_output_unit, end_output_unit
        use system, only: sys_t
        use qmc_data, only: output_in_t

        use calc, only: calc_type, mc_hilbert_space
        use hilbert_space, only: estimate_hilbert_space

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state

        type(sys_t), pointer :: sys
        integer :: truncation_level, nattempts, ncycles, rng_seed
        type(output_in_t) :: output_in
        integer, allocatable :: ref_det(:)
        integer :: opts, io_unit, ret_vals
        real :: t1, t2
        real(dp) :: hs_mean, hs_se
        character(12), parameter :: keys(3) = [character(12) :: 'sys', 'hilbert', 'output']

        call cpu_time(t1)

        lua_state = flu_copyptr(l)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)
        call read_hilbert_args(lua_state, opts, sys%nel, nattempts, ncycles, truncation_level, ref_det, rng_seed)
        call read_output_in_t(lua_state, opts, output_in)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        ! AOTUS returns a vector of size 0 to denote a non-existent vector.
        if (size(ref_det) == 0) deallocate(ref_det)
        if (allocated(ref_det)) then
            if (size(ref_det) /= sys%nel .and. parent) call stop_all('lua_hilbert_space', &
                            'Reference determinant does not match the number of electrons in system.')
        end if

        calc_type = mc_hilbert_space
        call init_output_unit(output_in, sys, io_unit)
        call estimate_hilbert_space(sys, truncation_level, nattempts, ncycles, ref_det, rng_seed, io_unit, hs_mean, hs_se)
        call end_output_unit(output_in%out_filename, io_unit)

        nresult = 1
        call aot_table_open(lua_state, thandle=ret_vals)
        call aot_table_set_val(hs_mean, lua_state, ret_vals, 'mean')
        call aot_table_set_val(hs_se, lua_state, ret_vals, 'std. err.')

        call cpu_time(t2)
        call register_timing(lua_state, "Hilbert space estimation", t2-t1)

    end function lua_hilbert_space

    function lua_canonical_estimates(L) result(nresult) bind(c)

        ! Run a Monte Carlo calculation to estimate various mean-field properties of a
        ! system in the canonical ensemble.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    canonical_estimates {
        !       sys = system,      -- required
        !       canonical_estimates = { ... }, -- required
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use const, only: p

        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use errors, only: stop_all
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args, register_timing
        use system, only: sys_t

        use calc, only: calc_type, mc_canonical_estimates
        use canonical_estimates, only: estimate_canonical

        integer(c_int) :: nresult
        type(c_ptr), value :: l

        type(flu_state) :: lua_state

        type(sys_t), pointer :: sys
        integer :: opts, rng_seed, ncycles, nattempts
        real :: t1, t2
        logical :: fermi_temperature, all_spin_sectors
        real(p) :: beta
        character(19), parameter :: keys(2) = [character(19) :: 'sys', 'canonical_estimates']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)
        call read_canonical_args(lua_state, opts, fermi_temperature, beta, nattempts, ncycles, rng_seed, all_spin_sectors)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = mc_canonical_estimates
        call estimate_canonical(sys, fermi_temperature, beta, nattempts, ncycles, all_spin_sectors, rng_seed)

        ! [todo] - return estimate of various canonical mean-field energies and error to lua.
        nresult = 0

        call cpu_time(t2)
        call register_timing(lua_state, "Canonical energy estimation", t2-t1)

    end function lua_canonical_estimates

    function lua_simple_fciqmc(L) result(nresult) bind(c)

        ! Run an FCIQMC calculation using the simple algorithm (slow, serial,
        ! memory hungry).

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    simple_fciqmc {
        !       sys = system,         -- required
        !       qmc = { ... },        -- required
        !       restart = { ... },
        !       reference = { ... },
        !       sparse = true/false,
        !       qmc_state = qmc_state,
        !    }

        ! Returns a qmc_state.

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close, aot_get_val

        use simple_fciqmc, only: do_simple_fciqmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args, register_timing
        use qmc_data, only: qmc_in_t, restart_in_t, qmc_state_t
        use reference_determinant, only: reference_t
        use system, only: sys_t

        use calc, only: calc_type, simple_fciqmc_calc, fciqmc_calc
        use calc_system_init, only: set_spin_polarisation, set_symmetry_aufbau

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(restart_in_t) :: restart_in
        type(reference_t) :: reference
        type(qmc_state_t), pointer :: qmc_state_restart, qmc_state_out
        logical :: use_sparse_hamil, have_qmc_state

        integer :: opts, err
        real :: t1, t2
        character(12), parameter :: keys(6) = [character(12) :: 'sys', 'qmc', 'restart', 'reference', 'sparse', 'qmc_state']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, sys)
        ! If using Aufbau determined symmetry need to do after setting spin polarisation.
        if (sys%aufbau_sym) call set_symmetry_aufbau(sys, 6)

        ! Get main table.
        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in, .true.)
        call read_restart_in(lua_state, opts, restart_in)
        call read_reference_t(lua_state, opts, reference, sys)
        call aot_get_val(use_sparse_hamil, err, lua_state, opts, 'sparse', default=.true.)
        call warn_unused_args(lua_state, keys, opts)
        call get_qmc_state(lua_state, have_qmc_state, qmc_state_restart)
        call aot_table_close(lua_state, opts)

        calc_type = simple_fciqmc_calc + fciqmc_calc
        allocate(qmc_state_out)
        if (have_qmc_state) then
            call do_simple_fciqmc(sys, qmc_in, restart_in, reference, use_sparse_hamil, qmc_state_out, qmc_state_restart)
        else
            call do_simple_fciqmc(sys, qmc_in, restart_in, reference, use_sparse_hamil, qmc_state_out)
        end if

        call push_qmc_state(lua_state, qmc_state_out)
        nresult = 1

        call cpu_time(t2)
        call register_timing(lua_state, "Simple FCIQMC calculation", t2-t1)

    end function lua_simple_fciqmc

    function lua_fciqmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    fciqmc {
        !       sys = system,          -- required
        !       qmc = { ... },         -- required
        !       fciqmc = { ... },
        !       semi_stoch = { ... },
        !       restart = { ... },
        !       load_bal = { ... },
        !       reference = { ... },
        !       logging = { ... },
        !       output = { ... },
        !       blocking = { ... },
        !       qmc_state = qmc_state,
        !    }

        ! Returns a qmc_state.

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use dmqmc_data, only: dmqmc_in_t
        use fciqmc, only: do_fciqmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args, register_timing
        use lua_hande_calc_utils, only: init_output_unit, end_output_unit
        use qmc_data, only: qmc_in_t, fciqmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t, &
                            qmc_state_t, output_in_t, blocking_in_t
        use logging, only: logging_in_t
        use reference_determinant, only: reference_t
        use system, only: sys_t
        use state_histograms, only: state_histogram_in_t

        use calc, only: calc_type, fciqmc_calc
        use calc_system_init, only: set_spin_polarisation, set_symmetry_aufbau

        integer(c_int)  :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(fciqmc_in_t) :: fciqmc_in
        type(semi_stoch_in_t) :: semi_stoch_in
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference
        type(qmc_state_t), pointer :: qmc_state_restart, qmc_state_out
        type(logging_in_t) :: logging_in
        type(output_in_t) :: output_in
        type(blocking_in_t) :: blocking_in
        type(state_histogram_in_t) :: state_hist_in

        logical :: have_restart_state

        integer :: opts, io_unit
        real :: t1, t2
        character(15), parameter :: keys(12) = [character(15) :: 'sys', 'qmc', 'fciqmc', 'semi_stoch', 'restart', &
                                                                'load_bal', 'reference', 'qmc_state', 'logging', &
                                                                'output', 'blocking', 'state_histogram']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)

        call read_qmc_in(lua_state, opts, qmc_in)
        call read_fciqmc_in(lua_state, opts, fciqmc_in)
        call read_blocking_in(lua_state, opts, blocking_in)
        call read_semi_stoch_in(lua_state, opts, qmc_in, semi_stoch_in)
        call read_restart_in(lua_state, opts, restart_in)
        call read_load_bal_in(lua_state, opts, load_bal_in)
        call read_reference_t(lua_state, opts, reference, sys)
        call read_logging_in_t(lua_state, opts, logging_in)
        call read_output_in_t(lua_state, opts, output_in)
        call read_state_histogram_in(lua_state, opts, state_hist_in)

        call get_qmc_state(lua_state, have_restart_state, qmc_state_restart)

        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = fciqmc_calc
        allocate(qmc_state_out)

        call init_output_unit(output_in, sys, io_unit)

        ! If using Aufbau determined symmetry need to do after setting spin polarisation.
        if (sys%aufbau_sym) call set_symmetry_aufbau(sys, io_unit)

        if (have_restart_state) then
            call do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, load_bal_in, io_unit, reference, logging_in, &
                           blocking_in, qmc_state_out, state_hist_in, qmc_state_restart)
        else
            call do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, load_bal_in, io_unit, reference, logging_in, &
                           blocking_in, qmc_state_out, state_hist_in)
        end if

        call end_output_unit(output_in%out_filename, io_unit)

        ! Return qmc_state to the user.
        call push_qmc_state(lua_state, qmc_state_out)
        nresult = 1

        call cpu_time(t2)
        call register_timing(lua_state, "FCIQMC calculation", t2-t1)

    end function lua_fciqmc

    function lua_ccmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    ccmc{
        !       sys = system,         --  required
        !       qmc = { ... },        --  required
        !       ccmc = { ... }
        !       restart = { ... },
        !       reference = { ... },
        !       logging = { ... },
        !       output = { ... },
        !       blocking = { ... },
        !       qmc_state = qmc_state,
        !       psip_list = psip_list,
        !    }

        ! Returns a qmc_state.

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use dmqmc_data, only: dmqmc_in_t
        use ccmc, only: do_ccmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args, register_timing
        use lua_hande_calc_utils, only: init_output_unit, end_output_unit
        use qmc_data, only: qmc_in_t, ccmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t, &
                            qmc_state_t, output_in_t, blocking_in_t, particle_t
        use logging, only: logging_in_t
        use reference_determinant, only: reference_t
        use system, only: sys_t

        use calc, only: calc_type, ccmc_calc
        use calc_system_init, only: set_spin_polarisation, set_symmetry_aufbau
        use errors, only: warning
        use excitations, only: end_excitations, init_excitations

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(ccmc_in_t) :: ccmc_in
        type(semi_stoch_in_t) :: semi_stoch_in
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference
        type(qmc_state_t), pointer :: qmc_state_restart, qmc_state_out
        type(logging_in_t) :: logging_in
        type(output_in_t) :: output_in
        type(blocking_in_t) :: blocking_in
        type(particle_t), pointer :: psip_list_in

        logical :: have_restart_state, have_psip_list
        integer :: opts, io_unit
        real :: t1, t2
        character(10), parameter :: keys(10) = [character(10) :: 'sys', 'qmc', 'ccmc', 'restart', 'reference', 'qmc_state', &
                                                                'logging', 'output', 'blocking', 'psip_list']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in)
        call read_ccmc_in(lua_state, opts, ccmc_in, sys)

        ! Need to extend bit strings for additional excitation level information if needed.
        if (ccmc_in%even_selection) then
            sys%basis%info_string_len = 1
            sys%basis%tot_string_len = sys%basis%bit_string_len + sys%basis%info_string_len
            sys%basis%tensor_label_len = sys%basis%tot_string_len
            call end_excitations(sys%basis%excit_mask)
            call init_excitations(sys%basis)
        end if

        call read_blocking_in(lua_state, opts, blocking_in)
        ! note that semi-stochastic is not (yet) available in CCMC.
        !call read_semi_stoch_in(lua_state, opts, qmc_in, semi_stoch_in)
        call read_restart_in(lua_state, opts, restart_in)
        ! load balancing is not available in CCMC; must use default settings.
        call read_reference_t(lua_state, opts, reference, sys)

        call read_logging_in_t(lua_state, opts, logging_in)
        call read_output_in_t(lua_state, opts, output_in)


        call get_qmc_state(lua_state, have_restart_state, qmc_state_restart)
        call get_psip_list(lua_state, have_psip_list, psip_list_in)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = ccmc_calc
        allocate(qmc_state_out)

        call init_output_unit(output_in, sys, io_unit)
        ! If using Aufbau determined symmetry need to do after setting spin polarisation.
        if (sys%aufbau_sym) call set_symmetry_aufbau(sys, io_unit)

        if (have_restart_state) then
            call do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference, logging_in, blocking_in, &
                            io_unit, qmc_state_out, qmc_state_restart=qmc_state_restart)
        else if (have_psip_list) then
            call do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference, logging_in, blocking_in, &
                            io_unit, qmc_state_out, psip_list_in=psip_list_in)
        else
            call do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference, logging_in, blocking_in, &
                            io_unit, qmc_state_out)
        end if

        ! Reset bit string length for reuse.
        if (ccmc_in%even_selection) then
            sys%basis%info_string_len = 0
            sys%basis%tot_string_len = sys%basis%bit_string_len + sys%basis%info_string_len
            sys%basis%tensor_label_len = sys%basis%tot_string_len
            call end_excitations(sys%basis%excit_mask)
            call init_excitations(sys%basis)
        end if

        call end_output_unit(output_in%out_filename, io_unit)

        call push_qmc_state(lua_state, qmc_state_out)
        nresult = 1

        call cpu_time(t2)
        call register_timing(lua_state, "CCMC calculation", t2-t1)

    end function lua_ccmc

    function lua_uccmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    uccmc{
        !       sys = system,         --  required
        !       qmc = { ... },        --  required
        !       ccmc = { ... },
        !       uccmc = { ... }
        !       restart = { ... },
        !       reference = { ... },
        !       logging = { ... },
        !       output = { ... },
        !       blocking = { ... },
        !       qmc_state = qmc_state,
        !    }

        ! Returns a qmc_state.

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use dmqmc_data, only: dmqmc_in_t
        use uccmc, only: do_uccmc
        use trotterized_uccmc, only: do_trot_uccmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args, register_timing
        use lua_hande_calc_utils, only: init_output_unit, end_output_unit
        use qmc_data, only: qmc_in_t, ccmc_in_t, uccmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t, &
                            qmc_state_t, output_in_t, blocking_in_t
        use logging, only: logging_in_t
        use reference_determinant, only: reference_t
        use system, only: sys_t

        use calc, only: calc_type, uccmc_calc, trot_uccmc_calc
        use calc_system_init, only: set_spin_polarisation, set_symmetry_aufbau
        use excitations, only: end_excitations, init_excitations
        use errors, only: stop_all

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(ccmc_in_t) :: ccmc_in
        type(uccmc_in_t) :: uccmc_in
        type(semi_stoch_in_t) :: semi_stoch_in
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference
        type(qmc_state_t), pointer :: qmc_state_restart, qmc_state_out
        type(logging_in_t) :: logging_in
        type(output_in_t) :: output_in
        type(blocking_in_t) :: blocking_in

        logical :: have_restart_state
        integer :: opts, io_unit
        real :: t1, t2
        character(10), parameter :: keys(10) = [character(10) :: 'sys', 'qmc', 'ccmc', 'uccmc', 'restart', &
                                                'reference', 'qmc_state', 'logging', 'output', 'blocking']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)

        ![todo] can this be merged with lua_ccmc?
        call set_spin_polarisation(sys%basis%nbasis, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in)
        call read_ccmc_in(lua_state, opts, ccmc_in, sys)
        call read_uccmc_in(lua_state, opts, uccmc_in, sys)
        if (ccmc_in%even_selection) &
            call stop_all('lua_uccmc', 'Even selection is not compatible with UCCMC.')

        if (uccmc_in%trot) then
            sys%basis%info_string_len = 2
            sys%basis%tot_string_len = sys%basis%bit_string_len + sys%basis%info_string_len
            sys%basis%tensor_label_len = sys%basis%tot_string_len
            call end_excitations(sys%basis%excit_mask)
            call init_excitations(sys%basis)
        end if

        call read_restart_in(lua_state, opts, restart_in)
        call read_reference_t(lua_state, opts, reference, sys)

        call read_logging_in_t(lua_state, opts, logging_in)
        call read_output_in_t(lua_state, opts, output_in)

        call get_qmc_state(lua_state, have_restart_state, qmc_state_restart)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        if (uccmc_in%trot) then
            calc_type = trot_uccmc_calc
        else
            calc_type = uccmc_calc
        end if

        allocate(qmc_state_out)

        call init_output_unit(output_in, sys, io_unit)
        ! If using Aufbau determined symmetry need to do after setting spin polarisation.
        if (sys%aufbau_sym) call set_symmetry_aufbau(sys, io_unit)

        if (have_restart_state) then
            if (uccmc_in%trot) then
                call do_trot_uccmc(sys, qmc_in, ccmc_in, uccmc_in, restart_in, load_bal_in, reference, logging_in, &
                            io_unit, qmc_state_out, qmc_state_restart)
            else
                call do_uccmc(sys, qmc_in, ccmc_in, uccmc_in, restart_in, load_bal_in, reference, logging_in, &
                            io_unit, qmc_state_out, qmc_state_restart)
            end if
        else
            if (uccmc_in%trot) then
                call do_trot_uccmc(sys, qmc_in, ccmc_in, uccmc_in, restart_in, load_bal_in, reference, logging_in, &
                            io_unit, qmc_state_out)
            else
                call do_uccmc(sys, qmc_in, ccmc_in, uccmc_in,  restart_in, load_bal_in, reference, logging_in, &
                            io_unit, qmc_state_out)
            end if
        end if

        if (uccmc_in%trot) then
            sys%basis%info_string_len = 0
            sys%basis%tot_string_len = sys%basis%bit_string_len + sys%basis%info_string_len
            sys%basis%tensor_label_len = sys%basis%tot_string_len
            call end_excitations(sys%basis%excit_mask)
            call init_excitations(sys%basis)
        end if

        call end_output_unit(output_in%out_filename, io_unit)

        call push_qmc_state(lua_state, qmc_state_out)
        nresult = 1

        call cpu_time(t2)
        call register_timing(lua_state, "UCCMC calculation", t2-t1)

    end function lua_uccmc

    function lua_dmqmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    dmqmc{
        !       sys = system,         -- required
        !       qmc = { ... },        -- required
        !       dmqmc = { ... },
        !       ipdmqmc = { ... },
        !       operators = { ... },
        !       rdm = { ... },
        !       restart = { ... },
        !       reference = { ... },
        !       qmc_state = qmc_state,
        !    }

        ! Returns a qmc_state.

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close, aot_table_from_1Darray

        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t
        use dmqmc, only: do_dmqmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args, register_timing
        use qmc_data, only: qmc_in_t, restart_in_t, load_bal_in_t, qmc_state_t
        use reference_determinant, only: reference_t
        use system, only: sys_t, heisenberg, read_in, ueg
        use state_histograms, only: state_histogram_in_t

        use const, only: p
        use calc, only: calc_type, dmqmc_calc, doing_dmqmc_calc, dmqmc_energy_squared
        use checking, only: check_allocate
        use real_lattice, only: create_next_nearest_orbs
        use calc_system_init, only: set_spin_polarisation, set_symmetry_aufbau

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(dmqmc_in_t) :: dmqmc_in
        type(dmqmc_estimates_t) :: dmqmc_estimates
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference
        type(qmc_state_t), pointer :: qmc_state_out, qmc_state_restart
        type(state_histogram_in_t) :: state_hist_in

        logical :: have_restart_state
        integer :: opts, ierr
        real :: t1, t2
        real(p), allocatable :: sampling_probs(:)
        character(15), parameter :: keys(10) = [character(15) :: 'sys', 'qmc', 'dmqmc', 'ipdmqmc', 'operators', 'rdm', &
                                                                'restart', 'reference', 'qmc_state', 'state_histogram']

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in)
        call read_dmqmc_in(lua_state, sys%basis%nbasis, opts, sys%system, dmqmc_in, dmqmc_estimates%subsys_info)

        ! We are required to handle the tensor length ourselves.
        sys%basis%tensor_label_len = 2*sys%basis%tot_string_len

        if (doing_dmqmc_calc(dmqmc_energy_squared)) then
            ! Create info no longer set in init_real_space.
            allocate(sys%real_lattice%next_nearest_orbs(sys%basis%nbasis,sys%basis%nbasis), stat=ierr)
            call check_allocate('sys%real_lattice%next_nearest_orbs',sys%basis%nbasis*sys%basis%nbasis,ierr)
            call create_next_nearest_orbs(sys%basis, sys%real_lattice)
        end if

        ! [todo] - do spin polarisation in system setup.
        ! I'm sure we can handle this more gracefully than hack it in here...!
        if (dmqmc_in%all_spin_sectors) then
            if (sys%system /= read_in .and. sys%system /= ueg) sys%Ms = sys%lattice%nsites
            call set_spin_polarisation(sys%basis%nbasis, sys)
            if (sys%system == heisenberg) sys%max_number_excitations = sys%lattice%nsites/2
        else
            call set_spin_polarisation(sys%basis%nbasis, sys)
        end if
        ! If using Aufbau determined symmetry need to do after setting spin polarisation.
        if (sys%aufbau_sym) call set_symmetry_aufbau(sys, 6)

        ! Now system initialisation is complete (boo), act on the other options.
        call read_restart_in(lua_state, opts, restart_in)
        call read_reference_t(lua_state, opts, reference, sys)
        call get_qmc_state(lua_state, have_restart_state, qmc_state_restart)
        call read_state_histogram_in(lua_state, opts, state_hist_in)
        ! load balancing not implemented in DMQMC--just use defaults.
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = dmqmc_calc
        allocate(qmc_state_out)
        if (have_restart_state) then
            call do_dmqmc(sys, qmc_in, dmqmc_in, dmqmc_estimates, restart_in, load_bal_in, reference, qmc_state_out, &
                          sampling_probs, state_hist_in, qmc_state_restart)
        else
            call do_dmqmc(sys, qmc_in, dmqmc_in, dmqmc_estimates, restart_in, load_bal_in, reference, qmc_state_out, &
                          sampling_probs, state_hist_in)
        end if

        sys%basis%tensor_label_len = sys%basis%tot_string_len

        nresult = 1
        call push_qmc_state(lua_state, qmc_state_out)
        if (dmqmc_in%find_weights) then
            nresult = nresult + 1
            call aot_table_from_1Darray(lua_state, opts, sampling_probs)
        end if

        call cpu_time(t2)
        call register_timing(lua_state, "DMQMC calculation", t2-t1)

    end function lua_dmqmc

    function lua_mp1_mc(L) result(nresult) bind(c)

        ! Run a deterministic MP1 wavefunction calculation and produce the 
        ! exact MP2 energy, and return a stochastically coarse-grained psip_list object.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        ! psip_list = mp1_mc {
        !     sys = sys,
        !     mp1 = {
        !           D0_population = nD0,              -- required
        !           state_size = S,                   -- required
        !           real_amplitudes = true/false,
        !           spawn_cutoff = cutoff,
        !           rng_seed = seed,
        !           even_selection = true/false,
        !     },
        !     qmc = {...},
        !     ccmc = {...},
        ! }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close, aot_exists, aot_table_open, aot_get_val
        use aot_vector_module, only: aot_get_val

        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args, register_timing

        use system, only: sys_t
        use qmc_data, only: reference_t, qmc_state_t, particle_t, qmc_in_t, ccmc_in_t
        use mp1, only: mp1_in_t, sample_mp1_wfn
        use excitations, only: init_excitations, end_excitations

        use errors, only: warning

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys

        type(mp1_in_t) :: mp1_in
        type(reference_t) :: ref

        integer :: opts
        logical :: have_seed

        real :: t1, t2

        type(particle_t), pointer :: psip_list

        type(qmc_in_t) :: qmc_in
        type(ccmc_in_t) :: ccmc_in
        logical :: qmc_exists, ccmc_exists

        call cpu_time(t1)

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)

        qmc_exists = aot_exists(lua_state, opts, 'qmc')
        ccmc_exists = aot_exists(lua_state, opts, 'ccmc')

        if (qmc_exists) call read_qmc_in(lua_state, opts, qmc_in)
        if (ccmc_exists) call read_ccmc_in(lua_state, opts, ccmc_in, sys)
        call read_mp1_in(lua_state, opts, mp1_in, have_seed, qmc_exists)

        ref%ex_level = 2

        if (qmc_exists) then
            call warning('lua_mp1_mc', &
                'qmc table provided, overriding D0_normlisation, state_size, real_amplitudes, and rng_seed!')
            mp1_in%D0_norm = qmc_in%D0_population
            mp1_in%state_size = qmc_in%walker_length
            mp1_in%real_amplitudes = qmc_in%real_amplitudes
            ! Calling read_qmc_in fetches the seed provided in lua, or else generates a random seed
            mp1_in%seed = qmc_in%seed
            have_seed = .true.
            mp1_in%spawn_cutoff = qmc_in%spawn_cutoff
        end if

        if (ccmc_exists) then
            mp1_in%even_selection = ccmc_in%even_selection
        end if

        call aot_table_close(lua_state, opts)

        if (mp1_in%even_selection) then
            sys%basis%info_string_len = 1
            sys%basis%tot_string_len = sys%basis%bit_string_len + sys%basis%info_string_len
            sys%basis%tensor_label_len = sys%basis%tot_string_len
            call end_excitations(sys%basis%excit_mask)
            call init_excitations(sys%basis)
        end if

        allocate(psip_list)

        if (have_seed) then
            call sample_mp1_wfn(sys, mp1_in, ref, psip_list, mp1_in%seed)
        else
            call sample_mp1_wfn(sys, mp1_in, ref, psip_list)
        end if

        call push_psip_list(lua_state, psip_list)
        nresult = 1

        call cpu_time(t2)
        call register_timing(lua_state, "MP1 initialisation", t2-t1)

    end function lua_mp1_mc

    ! --- table-derived type wrappers ---

    subroutine read_fci_in(lua_state, opts, basis, fci_in)

        ! Read in the fci and (optionally) davidson tables to a fci_in_t object.

        ! fci = {
        !     write_hamiltonian = true/false,
        !     hamiltonian_file = filename,
        !     write_determinants = true/false,
        !     determinant_file = filename,
        !     write_nwfns = M,
        !     wfn_file = filename,
        !     nanalyse = N,
        !     blacs_block_size = block_size,
        !     hamiltonian_diagonal_only = true/false,
        !     rdm = { ... }, -- L-d vector containing the sites to include in subsystem A.
        ! }
        ! davidson = {
        !     ndavidson_eigv = int,
        !     ndavidson_trialvec = int,
        !     davidson_maxsize = int,
        !     davidson_tol = real,
        !     davidson_maxiter = int,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the fci and (optionally) davidson tables.
        !    basis: information about the single-particle basis set of the system.
        ! Out:
        !    fci_in: fci_in_t object containing generic fci/davidson input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use aot_vector_module, only: aot_get_val

        use lua_hande_utils, only: warn_unused_args
        use basis_types, only: basis_t
        use checking, only: check_allocate
        use fci_utils, only: fci_in_t
        use errors, only: stop_all

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(basis_t), intent(in) :: basis
        type(fci_in_t), intent(inout) :: fci_in

        integer :: fci_table, err, fci_nrdms
        integer, allocatable :: err_arr(:)

        character(25), parameter :: fci_keys(10) = [character(25) :: 'write_hamiltonian', 'hamiltonian_file', &
                                                                    'write_determinants', 'determinant_file', 'write_nwfns', &
                                                                    'wfn_file', 'nanalyse', 'blacs_block_size', 'rdm', &
                                                                    'hamiltonian_diagonal_only']
        character(18), parameter :: davidson_keys(6) = [character(18) :: 'using_davidson', 'ndavidson_eigv', 'ndavidson_trialvec', &
                                                                         'davidson_maxsize', 'davidson_tol', 'davidson_maxiter']

        if (aot_exists(lua_state, opts, 'fci')) then
            call aot_table_open(lua_state, opts, fci_table, 'fci')

            ! Optional arguments
            call aot_get_val(fci_in%write_hamiltonian, err, lua_state, fci_table, 'write_hamiltonian')
            call aot_get_val(fci_in%hamiltonian_file, err, lua_state, fci_table, 'hamiltonian_file')
            call aot_get_val(fci_in%write_determinants, err, lua_state, fci_table, 'write_determinants')
            call aot_get_val(fci_in%determinant_file, err, lua_state, fci_table, 'determinant_file')
            call aot_get_val(fci_in%print_fci_wfn, err, lua_state, fci_table, 'write_nwfns')
            call aot_get_val(fci_in%print_fci_wfn_file, err, lua_state, fci_table, 'wfn_file')
            call aot_get_val(fci_in%analyse_fci_wfn, err, lua_state, fci_table, 'nanalyse')
            call aot_get_val(fci_in%block_size, err, lua_state, fci_table, 'blacs_block_size')
            call aot_get_val(fci_in%hamiltonian_diagonal_only, err, lua_state, fci_table, 'hamiltonian_diagonal_only')

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, fci_table, 'rdm')) then
                ! Currently restricted to one RDM in a single FCI calculation.
                fci_nrdms = 1
                allocate(fci_in%subsys_info(fci_nrdms), stat=err)
                call check_allocate('fci_in%subsys_info', fci_nrdms, err)
                call aot_get_val(fci_in%subsys_info(fci_nrdms)%subsystem_A, err_arr, basis%nbasis, lua_state, fci_table, 'rdm')
                fci_in%subsys_info(fci_nrdms)%A_nsites = size(fci_in%subsys_info(fci_nrdms)%subsystem_A)
            end if

            call warn_unused_args(lua_state, fci_keys, fci_table)

            call aot_table_close(lua_state, fci_table)

        end if

        ! Davidson table: optional and indicates doing a Davidson calculation.
        if (aot_exists(lua_state, opts, 'davidson')) then
            fci_in%using_davidson = .true.
            call aot_table_open(lua_state, opts, fci_table, 'davidson')
            call aot_get_val(fci_in%ndavidson_eigv, err, lua_state, fci_table, 'ndavidson_eigv')
            call aot_get_val(fci_in%ndavidson_trialvec, err, lua_state, fci_table, 'ndavidson_trialvec')
            if (fci_in%ndavidson_trialvec <= fci_in%ndavidson_eigv) then 
                call stop_all('read_fci_in', 'Number of trial vectors smaller than number of eigenvalues solving for')
            end if
            call aot_get_val(fci_in%davidson_maxsize, err, lua_state, fci_table, 'davidson_maxsize')
            if (fci_in%davidson_maxsize <= fci_in%ndavidson_trialvec*2) then 
                call stop_all('read_fci_in', 'Max Davidson basis size smaller than twice the number of trial vectors')
            end if
            call aot_get_val(fci_in%davidson_tol, err, lua_state, fci_table, 'davidson_tol')
            call aot_get_val(fci_in%davidson_maxiter, err, lua_state, fci_table, 'davidson_maxiter')
            call warn_unused_args(lua_state, davidson_keys, fci_table)
            call aot_table_close(lua_state, fci_table)
        end if

    end subroutine read_fci_in

    subroutine read_hilbert_args(lua_state, opts, nel, nattempts, ncycles, ex_level, ref_det, rng_seed)

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle to the table which is input to the Lua hilbert_space
        !        routine.
        !    nel: The number of electrons in the system.
        ! Out:
        !    nattempts: number of samples, i.e. number of random determinants, to use each cycle.
        !    ncycles: number of Monte Carlo cycles to perform.
        !    ex_level: maximum excitation level relative to the reference
        !        determinant to include in the Hilbert space.
        !    ref_det: reference determinant.  If not supplied by the user then
        !        this will be deallocated on output.
        !    rng_seed: seed to initialise the random number generator.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use aot_vector_module, only: aot_get_val

        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args, get_rng_seed

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts, nel
        integer, intent(out) :: nattempts, ncycles, ex_level, rng_seed
        integer, allocatable :: ref_det(:)

        integer :: hilbert_table, err
        integer, allocatable :: err_arr(:)
        character(10), parameter :: keys(6) = [character(10) :: 'sys', 'nattempts', 'ncycles', 'ex_level', 'reference', 'rng_seed']

        if (.not. aot_exists(lua_state, opts, 'hilbert') .and. parent) &
            call stop_all('read_hilbert_args','"hilbert" table not present.')
        call aot_table_open(lua_state, opts, hilbert_table, 'hilbert')

        call aot_get_val(nattempts, err, lua_state, hilbert_table, 'nattempts')
        if (err /= 0 .and. parent) call stop_all('read_hilbert_args', 'Number of attempts not supplied.')
        call aot_get_val(ncycles, err, lua_state, hilbert_table, 'ncycles', default=20)
        call aot_get_val(ex_level, err, lua_state, hilbert_table, 'ex_level', default=-1)
        call aot_get_val(ref_det, err_arr, nel, lua_state, hilbert_table, key='reference')
        call get_rng_seed(lua_state, hilbert_table, rng_seed)

        call warn_unused_args(lua_state, keys, hilbert_table)
        call aot_table_close(lua_state, hilbert_table)

    end subroutine read_hilbert_args

    subroutine read_canonical_args(lua_state, opts, fermi_temperature, beta, nattempts, ncycles, rng_seed, all_spin_sectors)

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle to the table which is input to the Lua canonical_estimates
        !        routine.
        ! Out:
        !    fermi_temperature: if true, rescale beta as the inverse reduced
        !        temperature.
        !    beta: target temperature.
        !    nattempts: number of samples to use each cycle.
        !    ncycles: number of Monte Carlo cycles to perform.
        !    rng_seed: seed to initialise the random number generator.
        !    all_spin_sectors: True if averaging over spin.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use const, only: p
        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args, get_rng_seed

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        logical, intent(out) :: fermi_temperature, all_spin_sectors
        integer, intent(out) :: nattempts, ncycles, rng_seed
        real(p), intent(out) :: beta

        integer :: canonical_estimates_table, err
        character(17), parameter :: keys(6) = [character(17) :: 'nattempts', 'ncycles', 'beta', &
                                                                'fermi_temperature', 'rng_seed', 'all_spin_sectors']

        if (.not. aot_exists(lua_state, opts, 'canonical_estimates') .and. parent) &
            call stop_all('read_canonical_args','"canonical_estimates" table not present.')
        call aot_table_open(lua_state, opts, canonical_estimates_table, 'canonical_estimates')

        call aot_get_val(nattempts, err, lua_state, canonical_estimates_table, 'nattempts')
        if (err /= 0 .and. parent) call stop_all('read_canonical_args', 'nattempts: Number of attempts/cycle not supplied.')
        call aot_get_val(ncycles, err, lua_state, canonical_estimates_table, 'ncycles')
        if (err /= 0 .and. parent) call stop_all('read_canonical_args', 'ncycles: Number of cycles not supplied.')

        call aot_get_val(beta, err, lua_state, canonical_estimates_table, 'beta')
        if (err /= 0 .and. parent) call stop_all('read_canonical_args', 'beta: target temperature not supplied.')
        call aot_get_val(fermi_temperature, err, lua_state, canonical_estimates_table, 'fermi_temperature', default=.false.)

        call get_rng_seed(lua_state, canonical_estimates_table, rng_seed)
        call aot_get_val(all_spin_sectors, err, lua_state, canonical_estimates_table, 'all_spin_sectors', default=.false.)

        call warn_unused_args(lua_state, keys, canonical_estimates_table)
        call aot_table_close(lua_state, canonical_estimates_table)

    end subroutine read_canonical_args

    subroutine read_mp1_in(lua_state, opts, mp1_in, have_seed, qmc_exists)

        ! mp1 = {
        !     D0_population = nD0,              -- required
        !     state_size = S,                   -- required
        !     real_amplitudes = true/false,
        !     spawn_cutoff = cutoff,
        !     rng_seed = seed,
        !     even_selection = true/false,
        ! }
        ! In:
        !   opts: handle to the table which is input to the Lua canonical_estimates
        !        routine.
        !   qmc_exists: whether a qmc table is passed into mp1_mc. If yes, some options will be overwritten with those
        !        in the qmc table.
        ! In/Out:
        !   lua_state: flu/Lua state to which the HANDE API is added.
        ! Out:
        !   mp1_in: mp1_in_t object containing mp1 options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use mp1, only: mp1_in_t
        use parallel, only: parent
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(mp1_in_t), intent(out) :: mp1_in
        logical, intent(out) :: have_seed
        logical, intent(in) :: qmc_exists

        integer :: mp1_table, err
        character(15), parameter :: keys(6) = [character(15) :: 'D0_population', 'state_size', &
                                                                'real_amplitudes', 'spawn_cutoff', 'rng_seed','even_selection']

        if (.not. aot_exists(lua_state, opts, 'mp1') .and. .not. qmc_exists) then
            if (parent) call stop_all('read_mp1_args', '"mp1" table not present.')
        end if

        if (aot_exists(lua_state, opts, 'mp1')) then
            call aot_table_open(lua_state, opts, mp1_table, 'mp1')

            call aot_get_val(mp1_in%D0_norm, err, lua_state, mp1_table, 'D0_population')
            if (err /= 0 .and. parent .and. .not. qmc_exists) &
                call stop_all('read_mp1_args', 'D0_population: Internal normalisation not supplied.')
            call aot_get_val(mp1_in%state_size, err, lua_state, mp1_table, 'state_size')
            if (err /= 0 .and. parent .and. .not. qmc_exists) call stop_all('read_mp1_args', 'state_size not set.')

            call aot_get_val(mp1_in%real_amplitudes, err, lua_state, mp1_table, 'real_amplitudes')
            call aot_get_val(mp1_in%spawn_cutoff, err, lua_state, mp1_table, 'spawn_cutoff')
            call aot_get_val(mp1_in%even_selection, err, lua_state, mp1_table, 'even_selection')

            have_seed = aot_exists(lua_state, mp1_table, 'rng_seed')
            if (have_seed) call aot_get_val(mp1_in%seed, err, lua_state, mp1_table, 'rng_seed')

            call warn_unused_args(lua_state, keys, mp1_table)
            call aot_table_close(lua_state, mp1_table)
        end if

    end subroutine read_mp1_in

    subroutine read_qmc_in(lua_state, opts, qmc_in, short)

        ! Read in a qmc table to a qmc_in_t object.

        ! qmc = {
        !     tau = tau,                                  -- required
        !     init_pop = N,                               -- required unless restarting
        !     mc_cycles = ncycles,                        -- required
        !     nreports = nreport,                         -- required
        !     state_size = walker_length,                 -- required unless restarting
        !     spawned_state_size = spawned_walker_length, -- required unless restarting
        !     rng_seed = seed,
        !     target_population = pop,
        !     real_amplitudes = true/false,
        !     real_amplitude_force_32 = true/false,
        !     spawn_cutoff = cutoff,
        !     excit_gen =
        !     'renorm'/'renorm_spin'/'no_renorm'/'no_renorm_spin'/'power_pitzer'/'uniform_power_pitzer'('power_pitzer_orderM')/
        !     'heat_bath_power_pitzer'('power_pitzer_orderM_ij')/'heat_bath_power_pitzer_ref'('power_pitzer_orderN')/
        !     'uniform_cauchy_schwarz'('cauchy_schwarz_orderM')/'heat_bath_cauchy_schwarz'('cauchy_schwarz_orderM_ij')/
        !     'heat_bath'/'heat_bath_uniform_singles'('heat_bath_uniform')/'heat_bath_exact_singles'('heat_bath_single')
        !     power_pitzer_min_weight = power_pitzer_min_weight,
        !     tau_search = true/false,
        !     pattempt_single = prob,
        !     pattempt_double = prob,
        !     pattempt_update = true/false,
        !     pattempt_zero_accum_data = true/false,
        !     pattempt_parallel = prob,
        !     initial_shift = shift,
        !     shift_damping = damp_factor,
        !     shift_harmonic_forcing = harmonic_force_factor,
        !     shift_harmonic_crit_damp = true/false,
        !     initiator = true/false,
        !     initiator_threshold = pop,
        !     use_mpi_barriers = true/false,
        !     vary_shift_from = shift or "proje",
        !     vary_shift = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the qmc table.
        !    short (optional): if true, then don't read in otherwise critical
        !       options not required for simple_fciqmc.  Default: false.
        ! Out:
        !    qmc_in: qmc_in_t object containing generic QMC input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: qmc_in_t, excit_gen_renorm, excit_gen_no_renorm, excit_gen_renorm_spin, excit_gen_no_renorm_spin
        use qmc_data, only: excit_gen_power_pitzer, excit_gen_cauchy_schwarz_occ, excit_gen_cauchy_schwarz_occ_ij
        use qmc_data, only: excit_gen_power_pitzer_occ, excit_gen_power_pitzer_occ_ij, excit_gen_power_pitzer_orderN
        use qmc_data, only: excit_gen_heat_bath, excit_gen_heat_bath_uniform, excit_gen_heat_bath_single
        use lua_hande_utils, only: warn_unused_args, get_rng_seed
        use parallel, only: parent
        use errors, only: stop_all, warning

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(qmc_in_t), intent(out) :: qmc_in
        logical, intent(in), optional :: short

        integer :: qmc_table, err
        character(len=30) :: str
        logical :: skip, no_renorm

        character(32), parameter :: keys(35) = [character(32) :: 'tau', 'init_pop', 'mc_cycles', 'nreports', 'state_size', &
                                                                 'spawned_state_size', 'rng_seed', 'target_population', &
                                                                 'real_amplitudes', 'spawn_cutoff', 'no_renorm', 'tau_search', &
                                                                 'real_amplitude_force_32', &
                                                                 'pattempt_single', 'pattempt_double', 'pattempt_update', &
                                                                 'pattempt_zero_accum_data', &
                                                                 'pattempt_parallel', 'initial_shift', 'shift_damping', &
                                                                 'shift_harmonic_forcing_two_stage', &
                                                                 'shift_harmonic_forcing', 'shift_harmonic_crit_damp', &
                                                                 'initiator', 'initiator_threshold', 'use_mpi_barriers', &
                                                                 'vary_shift_from', 'excit_gen', 'power_pitzer_min_weight', &
                                                                 'reference_target', 'vary_shift', 'quasi_newton', &
                                                                 'quasi_newton_threshold', 'quasi_newton_value', &
                                                                 'quasi_newton_pop_control']

        if (present(short)) then
            skip = short
        else
            skip = .false.
        end if

        if (.not. aot_exists(lua_state, opts, 'qmc') .and. parent) call stop_all('read_qmc_in','"qmc" table not present.')

        call aot_table_open(lua_state, opts, qmc_table, 'qmc')

        ! Required arguments
        call aot_get_val(qmc_in%tau, err, lua_state, qmc_table, 'tau')
        if (err /= 0 .and. parent) call stop_all('read_qmc_in', 'tau not set.')
        call aot_get_val(qmc_in%ncycles, err, lua_state, qmc_table, 'mc_cycles')
        if (err /= 0 .and. parent) call stop_all('read_qmc_in', 'mc_cycles not set.')
        call aot_get_val(qmc_in%nreport, err, lua_state, qmc_table, 'nreports')
        if (err /= 0 .and. parent) call stop_all('read_qmc_in', 'nreports not set.')

        ! Usually required, but not always. (Checked later if needed.)
        call aot_get_val(qmc_in%D0_population, err, lua_state, qmc_table, 'init_pop')
        if (.not. skip) then
            call aot_get_val(qmc_in%walker_length, err, lua_state, qmc_table, 'state_size')
            call aot_get_val(qmc_in%spawned_walker_length, err, lua_state, qmc_table, 'spawned_state_size')
        end if

        ! Optional arguments (defaults set in derived type).
        call aot_get_val(qmc_in%real_amplitudes, err, lua_state, qmc_table, 'real_amplitudes')
        call aot_get_val(qmc_in%real_amplitude_force_32, err, lua_state, qmc_table, 'real_amplitude_force_32')
        call aot_get_val(qmc_in%spawn_cutoff, err, lua_state, qmc_table, 'spawn_cutoff')
        call aot_get_val(qmc_in%pattempt_single, err, lua_state, qmc_table, 'pattempt_single')
        call aot_get_val(qmc_in%pattempt_double, err, lua_state, qmc_table, 'pattempt_double')
        call aot_get_val(qmc_in%pattempt_update, err, lua_state, qmc_table, 'pattempt_update')
        call aot_get_val(qmc_in%pattempt_zero_accum_data, err, lua_state, qmc_table, 'pattempt_zero_accum_data')
        call aot_get_val(qmc_in%pattempt_parallel, err, lua_state, qmc_table, 'pattempt_parallel')
        call aot_get_val(qmc_in%tau_search, err, lua_state, qmc_table, 'tau_search')
        call aot_get_val(qmc_in%initial_shift, err, lua_state, qmc_table, 'initial_shift')
        call aot_get_val(qmc_in%shift_damping, err, lua_state, qmc_table, 'shift_damping')
        call aot_get_val(qmc_in%shift_harmonic_forcing_two_stage, err, lua_state, qmc_table, 'shift_harmonic_forcing_two_stage')
        call aot_get_val(qmc_in%shift_harmonic_forcing, err, lua_state, qmc_table, 'shift_harmonic_forcing')
        call aot_get_val(qmc_in%shift_harmonic_crit_damp, err, lua_state, qmc_table, 'shift_harmonic_crit_damp')
        call aot_get_val(qmc_in%target_particles, err, lua_state, qmc_table, 'target_population')
        call aot_get_val(qmc_in%initiator_approx, err, lua_state, qmc_table, 'initiator')
        call aot_get_val(qmc_in%initiator_pop, err, lua_state, qmc_table, 'initiator_threshold')
        call aot_get_val(qmc_in%use_mpi_barriers, err, lua_state, qmc_table, 'use_mpi_barriers')
        call aot_get_val(qmc_in%power_pitzer_min_weight, err, lua_state, qmc_table, 'power_pitzer_min_weight')
        call aot_get_val(qmc_in%quasi_newton, err, lua_state, qmc_table, 'quasi_newton')
        call aot_get_val(qmc_in%quasi_newton_threshold, err, lua_state, qmc_table, 'quasi_newton_threshold')
        call aot_get_val(qmc_in%quasi_newton_value, err, lua_state, qmc_table, 'quasi_newton_value')
        call aot_get_val(qmc_in%quasi_newton_pop_control, err, lua_state, qmc_table, 'quasi_newton_pop_control')


        if (aot_exists(lua_state, qmc_table, 'reference_target')) then
            qmc_in%target_reference = .true.
            call aot_get_val(qmc_in%target_particles, err, lua_state, qmc_table, 'reference_target')
            if (aot_exists(lua_state, qmc_table, 'target_population') .and. parent) call stop_all('read_qmc_in', &
                'Cannot provide both target_population and reference_target')
        end if

        qmc_in%vary_shift_present = aot_exists(lua_state, qmc_table, 'vary_shift')
        if (qmc_in%vary_shift_present) then
            call aot_get_val(qmc_in%vary_shift, err, lua_state, qmc_table, 'vary_shift')
        end if

        if (aot_exists(lua_state, qmc_table, 'no_renorm')) then
            call aot_get_val(no_renorm, err, lua_state, qmc_table, 'no_renorm')
            if (parent) call warning('read_qmc_in', 'no_renorm is deprecated.  Please use excit_gen instead.')
            if (no_renorm) then
                qmc_in%excit_gen = excit_gen_no_renorm
            else
                qmc_in%excit_gen = excit_gen_renorm
            end if
        end if

        if (aot_exists(lua_state, qmc_table, 'excit_gen')) then
            call aot_get_val(str, err, lua_state, qmc_table, 'excit_gen')
            select case(str)
            case('renorm')
                qmc_in%excit_gen = excit_gen_renorm
            case('renorm_spin')
                qmc_in%excit_gen = excit_gen_renorm_spin
            case('no_renorm')
                qmc_in%excit_gen = excit_gen_no_renorm
            case('no_renorm_spin')
                qmc_in%excit_gen = excit_gen_no_renorm_spin
            case('uniform_power_pitzer')
                qmc_in%excit_gen = excit_gen_power_pitzer_occ
            case('power_pitzer_orderM')
                qmc_in%excit_gen = excit_gen_power_pitzer_occ
            case('heat_bath_power_pitzer')
                qmc_in%excit_gen = excit_gen_power_pitzer_occ_ij
            case('power_pitzer_orderM_ij')
                qmc_in%excit_gen = excit_gen_power_pitzer_occ_ij
            case('power_pitzer')
                qmc_in%excit_gen = excit_gen_power_pitzer
            case('heat_bath_power_pitzer_ref')
                qmc_in%excit_gen = excit_gen_power_pitzer_orderN
            case('power_pitzer_orderN')
                qmc_in%excit_gen = excit_gen_power_pitzer_orderN
            case('uniform_cauchy_schwarz')
                qmc_in%excit_gen = excit_gen_cauchy_schwarz_occ
            case('cauchy_schwarz_orderM')
                qmc_in%excit_gen = excit_gen_cauchy_schwarz_occ
            case('heat_bath_cauchy_schwarz')
                qmc_in%excit_gen = excit_gen_cauchy_schwarz_occ_ij
            case('cauchy_schwarz_orderM_ij')
                qmc_in%excit_gen = excit_gen_cauchy_schwarz_occ_ij
            case('heat_bath')
                qmc_in%excit_gen = excit_gen_heat_bath
            case('heat_bath_uniform_singles')
                qmc_in%excit_gen = excit_gen_heat_bath_uniform
            case('heat_bath_uniform')
                qmc_in%excit_gen = excit_gen_heat_bath_uniform
            case('heat_bath_exact_singles')
                qmc_in%excit_gen = excit_gen_heat_bath_single
            case('heat_bath_single')
                qmc_in%excit_gen = excit_gen_heat_bath_single
            case default
                call stop_all('read_qmc_in', 'Invalid excit_gen setting: '//trim(str))
            end select
        end if

        ! If user sets initial shift and vary_shift_from, assume they know what
        ! they're doing.  Otherwise, vary the shift from the initial shift
        ! value.
        qmc_in%vary_shift_from = qmc_in%initial_shift

        ! Optional arguments requiring special care.
        call get_rng_seed(lua_state, qmc_table, qmc_in%seed)
        if (aot_exists(lua_state, qmc_table, 'vary_shift_from')) then
            call aot_get_val(qmc_in%vary_shift_from, err, lua_state, qmc_table, 'vary_shift_from')
            if (err /= 0) then
                ! Perhaps passed in a string?
                call aot_get_val(str, err, lua_state, qmc_table, 'vary_shift_from')
                qmc_in% vary_shift_from_proje = trim(str) == 'proje'
            end if
        end if

        call warn_unused_args(lua_state, keys, qmc_table)

        call aot_table_close(lua_state, qmc_table)

    end subroutine read_qmc_in

    subroutine read_fciqmc_in(lua_state, opts, fciqmc_in)

        ! Read in an fciqmc table (if it exists) to an fciqmc_in object.

        ! fciqmc = {
        !     non_blocking_comm = true/false,
        !     load_balancing = true/false,
        !     init_spin_inverse_reference_det = true/false,
        !     guiding_function = 'none'/'neel_singlet',
        !     trial_function = 'single_basis'/'neel_singlet',
        !     select_reference_det = true, -- OR
        !     select_reference_det = {
        !        update_every = mc_cycles,
        !        pop_factor = factor,
        !     }
        !     quadrature_initiator = true/false,
        !     replica_tricks = true/false,
        !     density_matrices = true/false,
        !     density_matrix_file = filename,
        !     density_matrix_report = report,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the fciqmc table.
        ! Out:
        !    fciqmc_in: fciqmc_in_t object containing FCIQMC-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: fciqmc_in_t, neel_singlet, neel_singlet_guiding
        use lua_hande_utils, only: warn_unused_args
        use parallel, only: parent
        use errors, only: stop_all

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(fciqmc_in_t), intent(out) :: fciqmc_in

        integer :: fciqmc_table, ref_det, err
        character(len=12) :: str
        logical :: ref_det_flag
        character(31), parameter :: keys(11) = [character(31) :: 'non_blocking_comm', 'load_balancing', 'guiding_function', &
                                                                 'init_spin_inverse_reference_det', 'trial_function', &
                                                                 'select_reference_det', 'quadrature_initiator', &
                                                                 'replica_tricks', 'density_matrices', 'density_matrix_file', &
                                                                 'density_matrix_report']

        if (aot_exists(lua_state, opts, 'fciqmc')) then

            call aot_table_open(lua_state, opts, fciqmc_table, 'fciqmc')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(fciqmc_in%non_blocking_comm, err, lua_state, fciqmc_table, 'non_blocking_comm')
            call aot_get_val(fciqmc_in%doing_load_balancing, err, lua_state, fciqmc_table, 'load_balancing')
            call aot_get_val(fciqmc_in%init_spin_inv_D0, err, lua_state, fciqmc_table, 'init_spin_inverse_reference_det')
            call aot_get_val(fciqmc_in%quadrature_initiator, err, lua_state, fciqmc_table, 'quadrature_initiator')
            call aot_get_val(fciqmc_in%replica_tricks, err, lua_state, fciqmc_table, 'replica_tricks')
            call aot_get_val(fciqmc_in%density_matrices, err, lua_state, fciqmc_table, 'density_matrices')
            call aot_get_val(fciqmc_in%density_matrix_file, err, lua_state, fciqmc_table, 'density_matrix_file')
            call aot_get_val(fciqmc_in%density_matrix_report, err, lua_state, fciqmc_table, 'density_matrix_report')

            ! If the user has asked to calculate an RDM, then replica_tricks
            ! *must* be on.
            if (fciqmc_in%density_matrices) then
                fciqmc_in%replica_tricks = .true.
            end if

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, fciqmc_table, 'select_reference_det')) then
                call aot_table_open(lua_state, fciqmc_table, ref_det, 'select_reference_det')
                if (ref_det == 0) then
                    ! Just passed a boolean (hopefully).
                    call aot_get_val(ref_det_flag, err, lua_state, fciqmc_table, 'select_reference_det', default=.false.)
                    if (ref_det_flag) fciqmc_in%select_ref_det_every_nreports = 20
                else
                    call aot_get_val(fciqmc_in%select_ref_det_every_nreports, err, lua_state, ref_det, 'update_every', default=20)
                    call aot_get_val(fciqmc_in%ref_det_factor, err, lua_state, ref_det, 'pop_factor')
                    call aot_table_close(lua_state, ref_det)
                end if
            end if
            if (aot_exists(lua_state, fciqmc_table, 'guiding_function')) then
                call aot_get_val(str, err, lua_state, fciqmc_table, 'guiding_function')
                select case(trim(str))
                case('none')
                    ! use default
                case('neel_singlet')
                    fciqmc_in%guiding_function = neel_singlet_guiding
                    fciqmc_in%trial_function = neel_singlet
                case default
                    if (parent) call stop_all('read_fciqmc_in', 'Unknown guiding function')
                end select
            end if
            if (aot_exists(lua_state, fciqmc_table, 'trial_function')) then
                call aot_get_val(str, err, lua_state, fciqmc_table, 'trial_function')
                select case(trim(str))
                case('none')
                    ! use default
                case('neel_singlet')
                    fciqmc_in%trial_function = neel_singlet
                case default
                    if (parent) call stop_all('read_fciqmc_in', 'Unknown trial wavefunction')
                end select
            end if

            call warn_unused_args(lua_state, keys, fciqmc_table)

            call aot_table_close(lua_state, fciqmc_table)

        end if

    end subroutine read_fciqmc_in

    subroutine read_semi_stoch_in(lua_state, opts, qmc_in, semi_stoch_in)

        ! Read in an semi_stoch table (if it exists) to an semi_stoch_in object.

        ! semi_stoch = {
        !     space = 'read'/'high'/'ci',       -- required if semi_stoch table exists.
        !     size = target_pop,                -- required if space == 'high'
        !     start_iteration = mc_cycles,
        !     shift_start_iteration = mc_cycle, -- overrides start_iteration
        !     separate_annihilation = true/false,
        !     write = true/false/id,
        !     read = id,
        !     ci_space = { ... },               -- required if space == 'ci'.  Reference table.
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        !    qmc_in: qmc_in_t object containing qmc input options.  If semi-stochastic
        !        options are present, then qmc_in%real_amplitudes is forced to be true.
        ! In:
        !    opts: handle for the table containing the semi_stoch table.
        ! Out:
        !    semi_stoch_in: semi_stoch_in_t object containing semi-stochastic-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use parallel, only: parent
        use errors, only: stop_all, warning

        use qmc_data, only: qmc_in_t, semi_stoch_in_t, high_pop_determ_space, read_determ_space, ci_determ_space, &
                            semi_stoch_combined_annihilation, semi_stoch_separate_annihilation
        use lua_hande_utils, only: warn_unused_args, get_flag_and_id

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(qmc_in_t), intent(inout) :: qmc_in
        type(semi_stoch_in_t), intent(out) :: semi_stoch_in

        integer :: semi_stoch_table, err
        character(len=10) :: str
        logical :: separate_annihilation
        character(21), parameter :: keys(9) = [character(21) :: 'space', 'size', 'start_iteration', 'separate_annihilation', &
                                                                'shift_start_iteration', 'write_determ_space', 'write', 'read', &
                                                                'ci_space']

        if (aot_exists(lua_state, opts, 'semi_stoch')) then

            if (parent .and. .not. qmc_in%real_amplitudes) &
                call warning('read_semi_stoch_in', 'Enabling real_amplitudes as required by semi_stoch.')
            qmc_in%real_amplitudes = .true.

            call aot_table_open(lua_state, opts, semi_stoch_table, 'semi_stoch')

            ! Required arguments.
            if (aot_exists(lua_state, semi_stoch_table, 'space')) then
                call aot_get_val(str, err, lua_state, semi_stoch_table, 'space')
                select case(str)
                case('read')
                    semi_stoch_in%space_type = read_determ_space
                    semi_stoch_in%target_size = -1
                case('high')
                    semi_stoch_in%space_type = high_pop_determ_space
                    call aot_get_val(semi_stoch_in%target_size, err, lua_state, semi_stoch_table, 'size')
                    if (err /= 0 .and. parent) call stop_all('read_semi_stoch_in', 'Target space size not given.')
                case('ci')
                    semi_stoch_in%space_type = ci_determ_space
                    call read_reference_t(lua_state, semi_stoch_table, semi_stoch_in%ci_space, ref_table_name='ci_space')
                    if (semi_stoch_in%ci_space%ex_level < 0) &
                        call stop_all('read_semi_stoch_in', 'Max excitation level not given for ci_space')
                case default
                    if (parent) call stop_all('read_semi_stoch_in', 'Unknown semi stochastic space type')
                end select
            else
                if (parent) call stop_all('read_semi_stoch_in', 'Deterministic space not specified.')
            end if

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(semi_stoch_in%start_iter, err, lua_state, semi_stoch_table, 'start_iteration')
            call aot_get_val(separate_annihilation, err, lua_state, semi_stoch_table, 'separate_annihilation', default=.true.)
            call aot_get_val(semi_stoch_in%read_id, err, lua_state, semi_stoch_table, 'read')
            call get_flag_and_id(lua_state, semi_stoch_table, semi_stoch_in%write_determ_space, semi_stoch_in%write_id, 'write')
            if (separate_annihilation) then
                semi_stoch_in%projection_mode = semi_stoch_separate_annihilation
            else
                semi_stoch_in%projection_mode = semi_stoch_combined_annihilation
            end if

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, semi_stoch_table, 'write_determ_space')) then
                call aot_get_val(semi_stoch_in%write_determ_space, err, lua_state, semi_stoch_table, 'write_determ_space')
                if (parent) call warning('read_semi_stoch_in', 'write_determ_space is deprecated.  Please use write instead.')
            end if
            if (aot_exists(lua_state, semi_stoch_table, 'shift_start_iteration')) then
                semi_stoch_in%start_iter = huge(0) ! fixed up once the shift comes on...
                call aot_get_val(semi_stoch_in%shift_iter, err, lua_state, semi_stoch_table, 'shift_start_iteration')
            end if

            call warn_unused_args(lua_state, keys, semi_stoch_table)

            call aot_table_close(lua_state, semi_stoch_table)

        end if

    end subroutine read_semi_stoch_in

    subroutine read_ccmc_in(lua_state, opts, ccmc_in, sys)

        ! Read in an ccmc table (if it exists) to an ccmc_in object.

        ! ccmc = {
        !     move_frequency = frequency,
        !     cluster_multispawn_threshold = threshold,
        !     full_non_composite = true/false,
        !     linked = true/false,
        !     vary_shift_reference = true/false,
        !     density_matrices = true/false,
        !     density_matrix_file = filename,
        !     even_selection = true/false,
        !     multiref = true/false,
        !     n_secondary_ref = number of additional references,
        !     secondary_ref1,...,secondary_ref999 ={...},
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the ccmc table.
        !    sys:  sys_t object containing information of current system. 
        ! Out:
        !    ccmc_in: ccmc_in_t object containing ccmc-specific input options.

        use, intrinsic :: iso_fortran_env, only: iostat_end
        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: ccmc_in_t
        use lua_hande_utils, only: warn_unused_args
        use system, only: sys_t
        use errors, only: stop_all

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type (sys_t), intent(in) :: sys
        type(ccmc_in_t), intent(out) :: ccmc_in

        integer :: ccmc_table, err, i, ir, ios
        character(28), parameter :: keys(15) = [character(28) :: 'move_frequency', 'cluster_multispawn_threshold', &
                                                                'full_non_composite', 'linked', 'vary_shift_reference', &
                                                                'density_matrices', 'density_matrix_file', 'even_selection', &
                                                                'multiref', 'n_secondary_ref', 'mr_acceptance_search', &
                                                                'mr_excit_lvl','mr_secref_file','mr_n_frozen','mr_read_in']
        character(23) :: string ! 32 bit integer has 10 digits, should be more than enough
        ! secondary_refX keywords are not hardcoded in, so we dynamically add them into the
        ! array of allowed keys 
        character(23), dimension(:), allocatable :: secondary_ref_keys
        character(28), dimension(:), allocatable :: keys_concat
        character(10) :: str
        logical :: secref_exist

        if (aot_exists(lua_state, opts, 'ccmc')) then

            call aot_table_open(lua_state, opts, ccmc_table, 'ccmc')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(ccmc_in%move_freq, err, lua_state, ccmc_table, 'move_frequency')
            call aot_get_val(ccmc_in%cluster_multispawn_threshold, err, lua_state, ccmc_table, 'cluster_multispawn_threshold')
            call aot_get_val(ccmc_in%full_nc, err, lua_state, ccmc_table, 'full_non_composite')
            call aot_get_val(ccmc_in%linked, err, lua_state, ccmc_table, 'linked')
            call aot_get_val(ccmc_in%vary_shift_reference, err, lua_state, ccmc_table, 'vary_shift_reference')
            call aot_get_val(ccmc_in%density_matrices, err, lua_state, ccmc_table, 'density_matrices')
            call aot_get_val(ccmc_in%density_matrix_file, err, lua_state, ccmc_table, 'density_matrix_file')
            call aot_get_val(ccmc_in%even_selection, err, lua_state, ccmc_table, 'even_selection')
            call aot_get_val(ccmc_in%multiref, err, lua_state, ccmc_table, 'multiref')
            if (ccmc_in%multiref) then
                call aot_get_val(ccmc_in%mr_read_in, err, lua_state, ccmc_table, 'mr_read_in')

                call aot_get_val(ccmc_in%n_secondary_ref, err, lua_state, ccmc_table, 'n_secondary_ref')
                if (ccmc_in%n_secondary_ref == 0 .and. .not. ccmc_in%mr_read_in) then
                    call stop_all('read_ccmc_in', 'Number of secondary references unspecified.') 
                end if

                if (ccmc_in%mr_read_in) then
                    call aot_get_val(ccmc_in%mr_secref_file, err, lua_state, ccmc_table, 'mr_secref_file')
                    if (err.ne.0) then
                        call stop_all('read_ccmc_in','mr_read_in set but mr_secref_file unset.')
                    end if
                    inquire(file=ccmc_in%mr_secref_file, exist=secref_exist)
                    if (.not. secref_exist) call stop_all('read_ccmc_in','mr_secref_file not found!')

                    ! Find out how many secondary references are in the file
                    open(newunit=ir, file=ccmc_in%mr_secref_file, status='old', form='formatted', action='read')
                    ccmc_in%n_secondary_ref = 0
                    do
                        read(ir, *, iostat=ios)
                        if (ios==iostat_end) exit
                        ccmc_in%n_secondary_ref = ccmc_in%n_secondary_ref + 1
                    end do
                    close(ir)

                    call aot_get_val(ccmc_in%mr_n_frozen, err, lua_state, ccmc_table, 'mr_n_frozen')
                    call aot_get_val(ccmc_in%mr_excit_lvl, err, lua_state, ccmc_table, 'mr_excit_lvl')
                    if (ccmc_in%mr_excit_lvl.eq.-1) then
                        call stop_all('read_ccmc_in','mr_read_in set but mr_excit_lvl unset.')
                    endif
                    keys_concat = keys
                end if

                allocate(secondary_ref_keys(ccmc_in%n_secondary_ref))

                if (.not. ccmc_in%mr_read_in) then
                    allocate(ccmc_in%secondary_refs(ccmc_in%n_secondary_ref))
                    do i = 1, ccmc_in%n_secondary_ref
                        ! I0 makes sure there are no whitespaces around the number string
                        write (string, '(A13,I0)') 'secondary_ref', i ! up to 2.15E9 secondary references can be provided      
                        ! trim makes sure there are no trailing whitespaces 
                        call read_reference_t(lua_state, ccmc_table, ccmc_in%secondary_refs(i), sys, trim(string))
                             if (.not. allocated(ccmc_in%secondary_refs(i)%occ_list0)) then
                                 call stop_all('read_ccmc_in', 'Uninitialised secondary reference determinant.') 
                             end if
                             if (ccmc_in%secondary_refs(i)%ex_level == -1 .or. ccmc_in%secondary_refs(i)%ex_level == sys%nel) then
                                 call stop_all('read_ccmc_in', 'Uninitialised secondary reference excitation level.')
                             end if
                             secondary_ref_keys(i) = trim(string)
                    end do
                    keys_concat = [keys,secondary_ref_keys]
                end if

                call aot_get_val(str, err, lua_state, ccmc_table, 'mr_acceptance_search')
                select case (str)
                    case ('bk_tree')
                        ccmc_in%mr_acceptance_search = 1
                    case default
                        ccmc_in%mr_acceptance_search = 0
                end select
            else
                keys_concat = keys
            end if
            call warn_unused_args(lua_state, keys_concat, ccmc_table)

            call aot_table_close(lua_state, ccmc_table)

        end if

    end subroutine read_ccmc_in

    subroutine read_uccmc_in(lua_state, opts, uccmc_in, sys)

        ! Read in an ccmc table (if it exists) to an ccmc_in object.

        ! uccmc = {
        !   pow_trunc = pow_trunc,
        !   variational_energy = true/false,
        !   average_wfn = true/false,
        !   trotterized = true/false,
        !   threshold = threshold,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the ccmc table.
        !    sys:  sys_t object containing information of current system. 
        ! Out:
        !    ccmc_in: ccmc_in_t object containing ccmc-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: uccmc_in_t
        use lua_hande_utils, only: warn_unused_args
        use system, only: sys_t
        use errors, only: stop_all

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type (sys_t), intent(in) :: sys
        type(uccmc_in_t), intent(out) :: uccmc_in

        integer :: uccmc_table, err, i
        character(28), parameter :: keys(5) = [character(28) :: 'pow_trunc', 'variational_energy', &
                                                           'average_wfn', 'trotterized', 'threshold']

        if (aot_exists(lua_state, opts, 'uccmc')) then

            call aot_table_open(lua_state, opts, uccmc_table, 'uccmc')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(uccmc_in%pow_trunc, err, lua_state, uccmc_table, 'pow_trunc')
            call aot_get_val(uccmc_in%variational_energy, err, lua_state, uccmc_table, 'variational_energy')
            call aot_get_val(uccmc_in%average_wfn, err, lua_state, uccmc_table, 'average_wfn')
            call aot_get_val(uccmc_in%trot, err, lua_state, uccmc_table, 'trotterized')
            call aot_get_val(uccmc_in%threshold, err, lua_state, uccmc_table, 'threshold')
            call warn_unused_args(lua_state, keys, uccmc_table)

            call aot_table_close(lua_state, uccmc_table)

        end if

    end subroutine read_uccmc_in

    subroutine read_dmqmc_in(lua_state, nbasis, opts, system_name, dmqmc_in, subsys_info)

        ! Read in DMQMC-related input options from various tables to a dmqmc_in_t object.

        ! dmqmc = {
        !     replica_tricks = true/false,
        !     fermi_temperature = true/false,
        !     all_sym_sectors = true/false,
        !     all_spin_sectors = true/false,
        !     beta_loops = Nb,
        !     final_beta = bf,
        !     sampling_weights = { ... },
        !     vary_weights = N,
        !     find_weights = true/false,
        !     find_weights_start = iteration,
        !     symmetrize = true/false,
        !     initiator_level = ilevel,
        !     symmetric_bloch = true/false,
        !     walker_scale_factor = factor,
        ! }
        ! ipdmqmc = { -- sets ipdmqmc to true
        !     target_beta = b,
        !     initial_matrix = 'free_electron'/'hartree_fock',
        !     grand_canonical_initialisation = true/false,
        !     metropolis_attempts = nattempts,
        !     symmetric_interaction_picture = true/false,
        !     piecewise_shift = ps,
        !     count_diagonal_occupations = true/false,
        !     skip_gci_reference_check = true/false,
        ! }
        ! operators = {
        !     renyi2 = true/false,
        !     energy = true/false,
        !     energy2 = true/false,
        !     staggered_magnetisation = true/false,
        !     correlation = { ... }
        !     excit_dist = true/false,
        !     excit_dist_start = iteration,
        !     kinetic_energy = true/false,
        !     H0_energy = true/false,
        !     potential_energy = true/false,
        !     HI_energy = true/false,
        !     mom_dist = kmax,
        !     structure_factor = kmax,
        !     ref_projected_energy = true/false,
        ! }
        ! rdm = {
        !     spawned_state_size = X,
        !     rdms = { { ... } , { ... } , { ... } ... }
        !     ground_state = true/false,
        !     ground_state_start = iteration,
        !     instantaneous = true/false,
        !     write = true/false,
        !     concurrence = true/false,
        !     von_neumann = true/false,
        !     renyi2 = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the above tables.
        !    nbasis: number of single-particle basis functions in the basis set of the system.
        !    system_name: system being studied.
        ! Out:
        !    dmqmc_in: dmqmc_in_t object containing DMQMC input options.
        !    subsys_info: array of subsys_t objects containing the subsystem
        !       for each RDM of interest. Unallocated if no RDMs are specified.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_exists, aot_table_open, aot_table_close, aot_get_val, aot_table_length
        use aot_vector_module, only: aot_get_val

        use calc
        use dmqmc_data, only: dmqmc_in_t, subsys_t, hartree_fock_dm, free_electron_dm
        use checking, only: check_allocate
        use parallel, only: parent
        use errors, only: stop_all, warning
        use lua_hande_utils, only: warn_unused_args
        use system, only: heisenberg

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: nbasis, opts, system_name
        type(dmqmc_in_t), intent(out) :: dmqmc_in
        type(subsys_t), intent(out), allocatable :: subsys_info(:)

        integer :: table, subtable, err, i
        integer, allocatable :: err_arr(:)
        logical :: op, bloch_symmetry_unset
        character(len=13) :: str

        character(30), parameter :: dmqmc_keys(15) = [character(30) :: 'replica_tricks', 'fermi_temperature', 'all_sym_sectors',  &
                                                                       'all_spin_sectors', 'beta_loops', 'sampling_weights',      &
                                                                       'find_weights', 'find_weights_start', 'symmetrize',        &
                                                                       'vary_weights', 'initiator_level', 'symmetric_bloch',      &
                                                                       'walker_scale_factor', 'final_beta', 'piecewise_shift']
        character(30), parameter :: ip_keys(9)     = [character(30) :: 'target_beta', 'initial_beta', 'initial_matrix',           &
                                                                       'grand_canonical_initialisation', 'metropolis_attempts',   &
                                                                       'symmetric_interaction_picture', 'piecewise_beta',         &
                                                                       'skip_gci_reference_check', 'count_diagonal_occupations']
        character(30), parameter :: op_keys(13)    = [character(30) :: 'renyi2', 'energy', 'energy2', 'staggered_magnetisation',  &
                                                                       'correlation', 'excit_dist', 'kinetic_energy',             &
                                                                       'H0_energy', 'potential_energy', 'HI_energy', 'mom_dist',  &
                                                                       'structure_factor', 'ref_projected_energy']
        character(30), parameter :: rdm_keys(9)    = [character(30) :: 'spawned_state_size', 'rdms', 'ground_state',              &
                                                                       'ground_state_start', 'instantaneous', 'write',            &
                                                                       'concurrence', 'von_neumann', 'renyi2']

        dmqmc_calc_type = 0
        bloch_symmetry_unset = .false.

        ! All optional and straightfoward except the vector quantities.

        if (aot_exists(lua_state, opts, 'dmqmc')) then
            call aot_table_open(lua_state, opts, table, 'dmqmc')
            call aot_get_val(dmqmc_in%replica_tricks, err, lua_state, table, 'replica_tricks')
            call aot_get_val(dmqmc_in%fermi_temperature, err, lua_state, table, 'fermi_temperature')
            call aot_get_val(dmqmc_in%all_sym_sectors, err, lua_state, table, 'all_sym_sectors')
            call aot_get_val(dmqmc_in%all_spin_sectors, err, lua_state, table, 'all_spin_sectors')
            call aot_get_val(dmqmc_in%beta_loops, err, lua_state, table, 'beta_loops')
            call aot_get_val(dmqmc_in%final_beta, err, lua_state, table, 'final_beta')
            call aot_get_val(dmqmc_in%piecewise_shift, err, lua_state, table, 'piecewise_shift')
            call aot_get_val(dmqmc_in%walker_scale_factor, err, lua_state, table, 'walker_scale_factor')
            call aot_get_val(dmqmc_in%find_weights, err, lua_state, table, 'find_weights')
            call aot_get_val(dmqmc_in%find_weights_start, err, lua_state, table, 'find_weights_start')
            call aot_get_val(dmqmc_in%half_density_matrix, err, lua_state, table, 'symmetrize')
            if (aot_exists(lua_state, table, 'symmetric_bloch')) then
                call aot_get_val(dmqmc_in%symmetric_bloch, err, lua_state, table, 'symmetric_bloch', default=.true.)
            else
                bloch_symmetry_unset = .true.
            end if
            dmqmc_in%symmetric = dmqmc_in%symmetric_bloch
            if (aot_exists(lua_state, table, 'sampling_weights')) then
                dmqmc_in%weighted_sampling = .true.
                ! Certainly can't have more excitation levels than basis functions, so that's a handy upper-limit.
                call aot_get_val(dmqmc_in%sampling_probs, err_arr, nbasis, lua_state, table, 'sampling_weights')
            end if
            if (aot_exists(lua_state, table, 'initiator_level')) then
                call aot_get_val(dmqmc_in%initiator_level, err, lua_state, table, 'initiator_level')
            end if
            dmqmc_in%vary_weights = aot_exists(lua_state, table, 'vary_weights')
            call aot_get_val(dmqmc_in%finish_varying_weights, err, lua_state, table, 'vary_weights')
            call warn_unused_args(lua_state, dmqmc_keys, table)
            call aot_table_close(lua_state, table)
        end if

        ! If using the 'find_weights' option then weighted sampling should
        ! always be used, even if not asked for via a separate input.
        if (dmqmc_in%find_weights) dmqmc_in%weighted_sampling = .true.

        if (aot_exists(lua_state, opts, 'ipdmqmc')) then
            dmqmc_in%ipdmqmc = .true.
            call aot_table_open(lua_state, opts, table, 'ipdmqmc')
            if (aot_exists(lua_state, table, 'initial_beta')) then
                call aot_get_val(dmqmc_in%target_beta, err, lua_state, table, 'initial_beta')
                if (parent) call warning('read_dmqmc_in', 'initial_beta is deprecated.  Please use target_beta instead.')
            end if
            call aot_get_val(dmqmc_in%target_beta, err, lua_state, table, 'target_beta')
            if (aot_exists(lua_state, table, 'initial_matrix')) then
                call aot_get_val(str, err, lua_state, table, 'initial_matrix')
                select case(trim(str))
                case('free_electron')
                    dmqmc_in%initial_matrix = free_electron_dm
                case('hartree_fock')
                    dmqmc_in%initial_matrix = hartree_fock_dm
                case default
                    if (parent) call stop_all('read_dmqmc_in', 'Unknown  inital density matrix')
                end select
            end if
            call aot_get_val(dmqmc_in%symmetric_interaction_picture, err, lua_state, table, &
                             'symmetric_interaction_picture', default=.true.)
            if (bloch_symmetry_unset) dmqmc_in%symmetric_bloch = dmqmc_in%symmetric_interaction_picture
            dmqmc_in%symmetric = dmqmc_in%symmetric_interaction_picture
            call aot_get_val(dmqmc_in%grand_canonical_initialisation, err, lua_state, table, 'grand_canonical_initialisation')
            call aot_get_val(dmqmc_in%metropolis_attempts, err, lua_state, table, 'metropolis_attempts')
            call aot_get_val(dmqmc_in%count_diagonal_occupations, err, lua_state, table, 'count_diagonal_occupations')
            call aot_get_val(dmqmc_in%skip_gci_reference_check, err, lua_state, table, 'skip_gci_reference_check')
            call aot_get_val(dmqmc_in%piecewise_shift, err, lua_state, table, 'piecewise_shift')
            call aot_get_val(dmqmc_in%walker_scale_factor, err, lua_state, table, 'walker_scale_factor')
            call warn_unused_args(lua_state, ip_keys, table)
            call aot_table_close(lua_state, table)
        end if

        if (aot_exists(lua_state, opts, 'operators')) then
            call aot_table_open(lua_state, opts, table, 'operators')
            call aot_get_val(op, err, lua_state, table, 'renyi2', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_full_r2
            call aot_get_val(op, err, lua_state, table, 'energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_energy
            call aot_get_val(op, err, lua_state, table, 'energy2', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_energy_squared
            call aot_get_val(op, err, lua_state, table, 'staggered_magnetisation', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_staggered_magnetisation
            call aot_get_val(op, err, lua_state, table, 'kinetic_energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_kinetic_energy
            call aot_get_val(op, err, lua_state, table, 'H0_energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_H0_energy
            call aot_get_val(op, err, lua_state, table, 'potential_energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_potential_energy
            call aot_get_val(op, err, lua_state, table, 'HI_energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_HI_energy
            call aot_get_val(op, err, lua_state, table, 'ref_projected_energy', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_ref_proj_energy
            if (aot_exists(lua_state, table, 'correlation')) then
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_correlation
                call aot_get_val(dmqmc_in%correlation_sites, err_arr, nbasis, lua_state, table, 'correlation')
            end if
            if (aot_exists(lua_state, table, 'mom_dist')) then
                call aot_get_val(dmqmc_in%mom_dist_kmax, err, lua_state, table, 'mom_dist')
                dmqmc_in%calc_mom_dist = .true.
            end if
            if (aot_exists(lua_state, table, 'structure_factor')) then
                call aot_get_val(dmqmc_in%struc_fac_qmax, err, lua_state, table, 'structure_factor')
                dmqmc_in%calc_struc_fac = .true.
            end if
            call aot_get_val(dmqmc_in%calc_excit_dist, err, lua_state, table, 'excit_dist')
            call warn_unused_args(lua_state, op_keys, table)
            call aot_table_close(lua_state, table)
        end if

        if (aot_exists(lua_state, opts, 'rdm')) then
            call aot_table_open(lua_state, opts, table, 'rdm')

            ! instantaneous is optional but if present, spawned_state_size must
            ! be.
            call aot_get_val(dmqmc_in%rdm%calc_inst_rdm, err, lua_state, table, 'instantaneous')

            ! If doing rdms, must specify them and the spawned state size.
            call aot_get_val(dmqmc_in%rdm%spawned_length, err, lua_state, table, 'spawned_state_size')
            if (dmqmc_in%rdm%calc_inst_rdm .and. err /= 0) then
                if (parent) call stop_all('read_dmqmc_in', 'spawned_state_size not specified in rdm table.')
            end if
            if (aot_exists(lua_state, table, 'rdms')) then
                dmqmc_in%rdm%doing_rdm = .true.
                call aot_table_open(lua_state, table, subtable, 'rdms')
                dmqmc_in%rdm%nrdms = aot_table_length(lua_state, subtable)
                allocate(subsys_info(dmqmc_in%rdm%nrdms), stat=err)
                call check_allocate('subsys_info', dmqmc_in%rdm%nrdms, err)
                do i = 1, dmqmc_in%rdm%nrdms
                    call aot_get_val(subsys_info(i)%subsystem_A, err_arr, nbasis, lua_state, subtable, pos=i)
                    subsys_info(i)%A_nsites = size(subsys_info(i)%subsystem_A)
                end do
            else
                if (parent) call stop_all('read_dmqmc_in', 'rdms not specified in rdm table.')
            end if

            ! Optional arguments.
            call aot_get_val(dmqmc_in%rdm%calc_ground_rdm, err, lua_state, table, 'ground_state')
            call aot_get_val(dmqmc_in%start_av_rdm, err, lua_state, table, 'ground_state_start')
            call aot_get_val(dmqmc_in%rdm%output_rdm, err, lua_state, table, 'write')
            call aot_get_val(dmqmc_in%rdm%doing_concurrence, err, lua_state, table, 'concurrence')
            call aot_get_val(dmqmc_in%rdm%doing_vn_entropy, err, lua_state, table, 'von_neumann')
            call aot_get_val(op, err, lua_state, table, 'renyi2', default=.false.)
            if (op) dmqmc_calc_type = dmqmc_calc_type + dmqmc_rdm_r2
            call warn_unused_args(lua_state, rdm_keys, table)
            call aot_table_close(lua_state, table)
        else
            dmqmc_in%rdm%nrdms = 0
        end if
        ! Can't average over spin sectors alone in normal dmqmc calculation.
        if (dmqmc_in%all_spin_sectors .and. .not. dmqmc_in%all_sym_sectors .and. system_name /= heisenberg .and. parent) &
                                            & call stop_all('read_dmqmc_in', 'specified symmetry averaging not imlemented')

    end subroutine read_dmqmc_in

    subroutine read_load_bal_in(lua_state, opts, load_bal_in)

        ! Read in an load_bal table (if it exists) to an load_bal_in object.

        ! load_bal = {
        !     nslots = N,
        !     min_pop = pop,
        !     target = percent,
        !     max_attempts = max,
        !     write = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the load_bal table.
        ! Out:
        !    load_bal_in: load_bal_in_t object containing load-balancing-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: load_bal_in_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(load_bal_in_t), intent(out) :: load_bal_in

        integer :: load_bal_table, err
        character(12), parameter :: keys(5) = [character(12) :: 'nslots', 'min_pop', 'target', 'max_attempts', 'write']

        if (aot_exists(lua_state, opts, 'load_bal')) then

            call aot_table_open(lua_state, opts, load_bal_table, 'load_bal')

            ! Optional arguments (defaults set in derived type).
            load_bal_in%nslots = 20 ! Default for when load balancing in use
            call aot_get_val(load_bal_in%nslots, err, lua_state, load_bal_table, 'nslots')
            call aot_get_val(load_bal_in%pop, err, lua_state, load_bal_table, 'min_pop')
            call aot_get_val(load_bal_in%percent, err, lua_state, load_bal_table, 'target')
            call aot_get_val(load_bal_in%max_attempts, err, lua_state, load_bal_table, 'max_attempts')
            call aot_get_val(load_bal_in%write_info, err, lua_state, load_bal_table, 'write')

            call warn_unused_args(lua_state, keys, load_bal_table)

            call aot_table_close(lua_state, load_bal_table)

        end if

    end subroutine read_load_bal_in

    subroutine read_reference_t(lua_state, opts, ref, sys, ref_table_name)

        ! Read in a reference table (if it exists) to a reference_t object.

        ! reference = {
        !     det = { ... }, -- N-electron vector
        !     hilbert_space_det = { ... }, -- N-electron vector
        !     ex_level = truncation_level
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the reference table.
        !    sys (optional): sys_t object describing the system of interest.  Used to default to
        !        full space if present.
        !    ref_table_name (optional): name of table holding reference info.  Default: reference.
        ! Out:
        !    ref: reference_t object containing the input options describing the
        !        reference.  Note that ex_level is set to the number of electrons
        !        if not provided (incl. if the reference table is not present in
        !        opts) and sys is present and -1 otherwise.

        use flu_binding, only: flu_State
        use aot_vector_module, only: aot_get_val
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use errors, only: stop_all

        use lua_hande_utils, only: warn_unused_args
        use reference_determinant, only: reference_t
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(reference_t), intent(out) :: ref
        type(sys_t), intent(in), optional :: sys
        character(*), intent(in), optional :: ref_table_name

        integer :: ref_table, err
        integer, allocatable :: err_arr(:)
        character(17), parameter :: keys(3) = [character(17) :: 'det', 'hilbert_space_det', 'ex_level']
        ! Avoid using allocatable strings here for old compiler support.
        character(100) :: ref_name

        ! Set to full space/a problematic value by default.
        ref%ex_level = -1
        if (present(sys)) ref%ex_level = sys%nel

        ref_name = 'reference'
        if (present(ref_table_name)) then
            if (len(ref_table_name) > len(ref_name)) call stop_all('read_reference_t', 'ref_name size too large.')
            ref_name = ref_table_name
        end if

        if (aot_exists(lua_state, opts, ref_name)) then

            call aot_table_open(lua_state, opts, ref_table, ref_name)

            ! Optional arguments requiring special handling.

            ! Check if deeterminants exist first to avoid aotus allocating components unnecessarily.
            if (aot_exists(lua_state, ref_table, 'det')) then
                call aot_get_val(ref%occ_list0, err_arr, sys%nel, lua_state, ref_table, 'det')
            end if
            if (aot_exists(lua_state, ref_table, 'hilbert_space_det')) then
                call aot_get_val(ref%hs_occ_list0, err_arr, sys%nel, lua_state, ref_table, 'hilbert_space_det')
            end if

            call aot_get_val(ref%ex_level, err, lua_state, ref_table, 'ex_level')
           
            call warn_unused_args(lua_state, keys, ref_table)

            call aot_table_close(lua_state, ref_table)
        end if

    end subroutine read_reference_t

    subroutine read_restart_in(lua_state, opts, restart_in)

        ! Read in a restart table (if it exists) to a restart_in_t object.

        ! restart = {
        !    read = true/false/id,
        !    write = true/false/id,
        !    write_shift = true/false/id,
        !    write_frequency = nreports,
        !    rng = true/false,
        ! }

        ! where id is an integer and sets the id for the restart file of the relevant type.
        ! If an id is provided, then the corresponding flag is also set to true.

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the restart table.
        ! Out:
        !    restart_in: restart_in_t object contianing the input options describing restart files.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_exists, aot_table_open, aot_table_close, aot_get_val

        use lua_hande_utils, only: warn_unused_args, get_flag_and_id
        use qmc_data, only: restart_in_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(restart_in_t), intent(out) :: restart_in

        integer :: err, restart_table
        character(15), parameter :: keys(5) = [character(15) :: 'read', 'write', 'write_shift', 'write_frequency', 'rng']

        if (aot_exists(lua_state, opts, 'restart')) then

            call aot_table_open(lua_state, opts, restart_table, 'restart')

            call aot_get_val(restart_in%write_freq, err, lua_state, restart_table, 'write_frequency')
            call aot_get_val(restart_in%restart_rng, err, lua_state, restart_table, 'rng')

            associate(r_in=>restart_in)
                call get_flag_and_id(lua_state, restart_table, r_in%read_restart, r_in%read_id, 'read')
                call get_flag_and_id(lua_state, restart_table, r_in%write_restart, r_in%write_id, 'write')
                call get_flag_and_id(lua_state, restart_table, r_in%write_restart_shift, r_in%write_shift_id, 'write_shift')
            end associate

            call warn_unused_args(lua_state, keys, restart_table)

            call aot_table_close(lua_state, restart_table)

        end if

    end subroutine read_restart_in

    subroutine read_output_in_t(lua_state, opts, output_in)

        ! Read in output table (if it exists) and set up.

        ! output = {
        !    filename = filename,
        !    reprint_sys = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the output table.
        ! Out:
        !    output_in: output_in_t object contianing the input options
        !       describing output files.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_exists, aot_table_open, aot_table_close, aot_get_val

        use lua_hande_utils, only: warn_unused_args, get_flag_and_id
        use qmc_data, only: output_in_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(output_in_t), intent(out) :: output_in

        integer :: err, output_table
        character(11), parameter :: keys(2) = [character(11) :: 'filename', 'reprint_sys']

        if (aot_exists(lua_state, opts, 'output')) then

            call aot_table_open(lua_state, opts, output_table, 'output')

            call aot_get_val(output_in%out_filename, err, lua_state, output_table, 'filename')
            call aot_get_val(output_in%reprint_sys_info, err, lua_state, output_table, 'reprint_sys')

            call warn_unused_args(lua_state, keys, output_table)

            call aot_table_close(lua_state, output_table)

        end if

    end subroutine read_output_in_t

    subroutine read_blocking_in(lua_state, opts, blocking_in)

        ! Read in a blocking table (if it exists) to a blocking_in_t object.

        ! blocking = {
        !    blocking_on_the_fly = true/false,
        !    start_save_frequency = log_2(nreports),
        !    start_point_number = N,
        !    filename = filename,
        !    start_point = niterations,
        !    error_limit = error_limit,
        !    blocks_used = blocks_used,
        !    min_blocks_used = min_blocks_used,
        !    auto_shift_damping = true/false,
        !    shift_damping_precision = ratio,
        !    force_shift_damping_opt = true/false,
        ! }

        ! If start_save_frequency, start_point_number, error_limit, min ratio and start_point
        ! is not specified default value is used.

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle of the table containing the blocking table.
        ! Out:
        !    blocking_in: blocking_in_t object containing the input options for
        !    blocking on the fly.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_exists, aot_table_open, aot_table_close
        use aot_table_module, only: aot_get_val

        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: blocking_in_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(blocking_in_t), intent(out) :: blocking_in

        integer :: err, blocking_table
        character(24), parameter :: keys(11) = [character(24) ::  'blocking_on_the_fly', 'start_save_frequency',   &
                                                                'start_point_number', 'filename', 'start_point', &
                                                                'error_limit', 'blocks_used', 'min_blocks_used', &
                                                                'auto_shift_damping', 'shift_damping_precision', &
                                                                'force_shift_damping_opt']

        if (aot_exists(lua_state, opts, 'blocking')) then

            call aot_table_open(lua_state, opts, blocking_table, 'blocking')

            call aot_get_val(blocking_in%blocking_on_the_fly, err, lua_state, blocking_table, 'blocking_on_the_fly')

            call aot_get_val(blocking_in%start_save_frequency, err, lua_state, blocking_table, 'start_save_frequency')

            call aot_get_val(blocking_in%start_point_number, err, lua_state, blocking_table, 'start_point_number')

            call aot_get_val(blocking_in%filename, err, lua_state, blocking_table, 'filename')

            call aot_get_val(blocking_in%start_point, err, lua_state, blocking_table, 'start_point')

            call aot_get_val(blocking_in%error_limit, err, lua_state, blocking_table, 'error_limit')

            call aot_get_val(blocking_in%blocks_used, err, lua_state, blocking_table, 'blocks_used')

            call aot_get_val(blocking_in%min_blocks_used, err, lua_state, blocking_table, 'min_blocks_used')

            call aot_get_val(blocking_in%auto_shift_damping, err, lua_state, blocking_table, 'auto_shift_damping')

            call aot_get_val(blocking_in%shift_damping_precision, err, lua_state, blocking_table, 'shift_damping_precision')

            call aot_get_val(blocking_in%force_shift_damping_opt, err, lua_state, blocking_table, 'force_shift_damping_opt')

            call warn_unused_args(lua_state, keys, blocking_table)

            call aot_table_close(lua_state, blocking_table)

        end if

    end subroutine read_blocking_in

    subroutine read_state_histogram_in(lua_state, opts, state_hist_in)

        ! Read in a state_histogram table to a state_histogram_in_t object.

        ! state_histogram = {
        !     report_frequency = nreports,
        !     decades = decades_to_bin,
        !     nbins = bins_per_decade,
        !     skip_memory_check = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the qmc table.
        ! Out:
        !    state_hist_in: Derived type object containing read in
        !        state histogram input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use state_histograms, only: state_histogram_in_t
        use lua_hande_utils, only: warn_unused_args
        use errors, only: stop_all, warning

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(state_histogram_in_t), intent(out) :: state_hist_in

        integer :: state_hist_table, err

        character(24), parameter :: keys(4) = [character(24) :: 'report_frequency', 'decades', 'nbins', &
                                                                'skip_memory_check']

        if (aot_exists(lua_state, opts, 'state_histogram')) then
            state_hist_in%state_histograms = .true.

            call aot_table_open(lua_state, opts, state_hist_table, 'state_histogram')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(state_hist_in%report_frequency, err, lua_state, state_hist_table, 'report_frequency')
            call aot_get_val(state_hist_in%decades, err, lua_state, state_hist_table, 'decades')
            call aot_get_val(state_hist_in%nbins, err, lua_state, state_hist_table, 'nbins')
            call aot_get_val(state_hist_in%skip_memory_check, err, lua_state, state_hist_table, 'skip_memory_check')

            call warn_unused_args(lua_state, keys, state_hist_table)

            call aot_table_close(lua_state, state_hist_table)
        end if

    end subroutine read_state_histogram_in

    subroutine get_qmc_state(lua_state, have_qmc_state, qmc_state)

        ! Get (if present) a qmc_state_t object passed in to resume a calculation.

        ! In/Out:
        !   lua_state: flu/Lua state to which the HANDE API is added.
        ! Out:
        !   have_qmc_state: true if qmc_state was passed in from lua.
        !   qmc_state: the qmc_state_t object passed in if have_qmc_state.

        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_table_top, aot_exists, aot_get_val
        use aot_table_ops_module, only: aot_table_open, aot_table_close

        use qmc_data, only: qmc_state_t
        use lua_hande_utils, only: get_userdata

        type(flu_State), intent(inout) :: lua_state
        logical, intent(out) :: have_qmc_state
        type(qmc_state_t), pointer, intent(out) :: qmc_state

        integer :: opts, qs_table
        type(c_ptr) :: qs_ptr

        opts = aot_table_top(lua_state)

        have_qmc_state = aot_exists(lua_state, opts, 'qmc_state')
        if (have_qmc_state) then
            ! Get qmc_state object
            call aot_table_open(lua_state, opts, qs_table, key='qmc_state')
            call get_userdata(lua_state, qs_table, "qmc_state", qs_ptr)
            call aot_table_close(lua_state, qs_table)
            call c_f_pointer(qs_ptr, qmc_state)
        end if

    end subroutine get_qmc_state

    subroutine get_psip_list(lua_state, have_psip_list, psip_list)

        ! Get (if present) a particle_t object passed in from an MP1 calculation to initialise a FCIQMC/CCMC calculation.

        ! In/Out:
        !   lua_state: flu/Lua state to which the HANDE API is added.
        ! Out:
        !   have_psip_list: true if psip_list was passed in from lua.
        !   psip_list: the particle_t object passed in if have_psip_list.

        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_table_top, aot_exists, aot_get_val
        use aot_table_ops_module, only: aot_table_open, aot_table_close

        use qmc_data, only: particle_t
        use lua_hande_utils, only: get_userdata

        type(flu_State), intent(inout) :: lua_state
        logical, intent(out) :: have_psip_list
        type(particle_t), pointer, intent(out) :: psip_list

        integer :: opts, pl_table
        type(c_ptr) :: pl_ptr

        opts = aot_table_top(lua_state)

        have_psip_list = aot_exists(lua_state, opts, 'psip_list')
        if (have_psip_list) then
            ! Get psip_list object
            call aot_table_open(lua_state, opts, pl_table, key='psip_list')
            call get_userdata(lua_state, pl_table, "psip_list", pl_ptr)
            call aot_table_close(lua_state, pl_table)
            call c_f_pointer(pl_ptr, psip_list)
        end if

    end subroutine get_psip_list

    subroutine push_qmc_state(lua_state, qmc_state)

        ! Add a table containing the passed qmc_state to the lua stack for returning

        ! In/Out:
        !   lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !   qmc_state: the qmc_state object to return to lua

        use, intrinsic :: iso_c_binding, only: c_loc
        use flu_binding, only: flu_State, flu_pushlightuserdata, flu_pushstring, flu_settable, flu_pushcclosure, fluL_setmetatable
        use aot_table_ops_module, only: aot_table_open, aot_table_close

        use qmc_data, only: qmc_state_t

        type(qmc_state_t), pointer, intent(in) :: qmc_state
        type(flu_state), intent(inout) :: lua_state

        integer :: table

        ! Create table to become qmc_state object
        call aot_table_open(lua_state, thandle=table)
        
        ! Add qmc_state pointer as t.qmc_state
        call flu_pushstring(lua_state, "qmc_state")
        call flu_pushlightuserdata(lua_state, c_loc(qmc_state))
        call flu_settable(lua_state, table)

        ! Add deallocation function as t:free()
        call flu_pushstring(lua_state, "free")
        call flu_pushcclosure(lua_state, lua_dealloc_qmc_state, 0)
        call flu_settable(lua_state, table)

        ! Set metatable to mark for finalisation.  Note metatable is created in register_lua_hande_api.
        call fluL_setmetatable(lua_state, "qmc_state")

    end subroutine push_qmc_state

    subroutine push_psip_list(lua_state, psip_list)

        ! Add a table containing the passed psip_list to the lua stack for returning

        ! In/Out:
        !   lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !   psip_list: the particle_t object to return to lua

        use, intrinsic :: iso_c_binding, only: c_loc
        use flu_binding, only: flu_State, flu_pushlightuserdata, flu_pushstring, flu_settable, flu_pushcclosure, fluL_setmetatable
        use aot_table_ops_module, only: aot_table_open, aot_table_close

        use qmc_data, only: particle_t

        type(particle_t), pointer, intent(in) :: psip_list
        type(flu_state), intent(inout) :: lua_state

        integer :: table

        ! Create table to become psip_list object
        call aot_table_open(lua_state, thandle=table)
        
        ! Add psip_list pointer as t.psip_list
        call flu_pushstring(lua_state, "psip_list")
        call flu_pushlightuserdata(lua_state, c_loc(psip_list))
        call flu_settable(lua_state, table)

        ! Add deallocation function as t:free()
        call flu_pushstring(lua_state, "free")
        call flu_pushcclosure(lua_state, lua_dealloc_psip_list, 0)
        call flu_settable(lua_state, table)

        ! Set metatable to mark for finalisation.  Note metatable is created in register_lua_hande_api.
        call fluL_setmetatable(lua_state, "psip_list")

    end subroutine push_psip_list
    
    function lua_dealloc_qmc_state(L) result(nresult) bind(c)

        ! Deallocate a qmc_state object.  Expects to be called from lua with a single argument --
        ! the qmc_state object to be deallocated.

        ! In/Out:
        !   L: lua state (bare C pointer).

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer, c_null_ptr
        use flu_binding, only: flu_State, flu_copyptr, flu_pushstring, flu_pushlightuserdata, flu_settable, &
                               flu_getmetatable, flu_pop
        use aot_table_ops_module, only: aot_table_top
        use aot_table_module, only: aot_get_val, aot_table_close
        use aot_top_module, only: aot_top_get_val

        use qmc_data, only: qmc_state_t
        use dealloc, only: dealloc_qmc_state_t
        use errors, only: stop_all
        use lua_hande_utils, only: get_userdata

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state
        integer :: qs_table
        type(c_ptr) :: qs_ptr
        type(qmc_state_t), pointer :: qs

        lua_state = flu_copyptr(L)

        qs_table = aot_table_top(lua_state)
        call get_userdata(lua_state, qs_table, "qmc_state", qs_ptr)
        call c_f_pointer(qs_ptr, qs)

        if (associated(qs)) then
            call dealloc_qmc_state_t(qs)
            deallocate(qs)
        end if

        ! Update table with deallocated pointer.
        call flu_pushstring(lua_state, "qmc_state")
        call flu_pushlightuserdata(lua_state, c_null_ptr)
        call flu_settable(lua_state, qs_table)

        call aot_table_close(lua_state, qs_table)

        nresult = 0

    end function lua_dealloc_qmc_state

    function lua_dealloc_psip_list(L) result(nresult) bind(c)

        ! Deallocate a particle_t object.  Expects to be called from lua with a single argument --
        ! the particle_t object to be deallocated.

        ! In/Out:
        !   L: lua state (bare C pointer).

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer, c_null_ptr
        use flu_binding, only: flu_State, flu_copyptr, flu_pushstring, flu_pushlightuserdata, flu_settable, &
                               flu_getmetatable, flu_pop
        use aot_table_ops_module, only: aot_table_top
        use aot_table_module, only: aot_get_val, aot_table_close
        use aot_top_module, only: aot_top_get_val

        use qmc_data, only: particle_t
        use particle_t_utils, only: dealloc_particle_t
        use errors, only: stop_all
        use lua_hande_utils, only: get_userdata

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state
        integer :: pl_table
        type(c_ptr) :: pl_ptr
        type(particle_t), pointer :: pl

        lua_state = flu_copyptr(L)

        pl_table = aot_table_top(lua_state)
        call get_userdata(lua_state, pl_table, "psip_list", pl_ptr)
        call c_f_pointer(pl_ptr, pl)

        if (associated(pl)) then
            call dealloc_particle_t(pl)
            deallocate(pl)
        end if

        ! Update table with deallocated pointer.
        call flu_pushstring(lua_state, "psip_list")
        call flu_pushlightuserdata(lua_state, c_null_ptr)
        call flu_settable(lua_state, pl_table)

        call aot_table_close(lua_state, pl_table)

        nresult = 0

    end function lua_dealloc_psip_list

    subroutine read_logging_in_t(lua_state, opts, logging_in)

        ! Read in options associated with the logging table (only for debug builds).

        ! logging = {
        !   calc = verbosity_level,
        !   spawn = verbosity_level,
        !   death = verbosity_level,
        !   }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the semi_stoch table.
        ! Out:
        !    logging_in: logging_in_t object containing logging-specific input options.


        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use lua_hande_utils, only: warn_unused_args

        use logging, only: logging_in_t

        use errors, only: stop_all
        use const, only: debug
        use parallel, only: parent

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts

        type(logging_in_t), intent(out) :: logging_in
        integer :: logging_table, err

        character(15), parameter :: keys(8) = [character(15) :: 'calc', 'spawn', 'death', &
                                                'stoch_selection', 'select', 'start', 'finish', 'write']

        if (aot_exists(lua_state, opts, 'logging')) then

            if (debug) then

                call aot_table_open(lua_state, opts, logging_table, 'logging')

                call aot_get_val(logging_in%calc, err, lua_state, logging_table, 'calc')

                call aot_get_val(logging_in%spawn, err, lua_state, logging_table, 'spawn')

                call aot_get_val(logging_in%death, err, lua_state, logging_table, 'death')

                call aot_get_val(logging_in%stoch_selection, err, lua_state, logging_table, 'stoch_selection')

                call aot_get_val(logging_in%selection, err, lua_state, logging_table, 'select')

                call aot_get_val(logging_in%start_iter, err, lua_state, logging_table, 'start')

                call aot_get_val(logging_in%end_iter, err, lua_state, logging_table, 'finish')

                call aot_get_val(logging_in%write_to, err, lua_state, logging_table, 'write')

                call warn_unused_args(lua_state, keys, logging_table)
                call aot_table_close(lua_state, logging_table)

            else
                if (parent) call stop_all('read_logging_in_t', &
                        'Tried to pass logging options to a non-debug build.')
            end if

        end if
    end subroutine read_logging_in_t

end module lua_hande_calc
