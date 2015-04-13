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
        !           sys = system,
        !           fci = { ... },
        !           lanczos = { ... },
        !           reference = { ... },
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_exists, aot_table_close

        use calc, only: calc_type, exact_diag, lanczos_diag
        use fci_lanczos, only: do_fci_lanczos
        use fci_lapack, only: do_fci_lapack
        use fci_utils, only: fci_in_t
        use lua_hande_system, only: get_sys_t
        use qmc_data, only: reference_t
        use system, only: sys_t

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state
        integer :: opts
        type(sys_t), pointer :: sys
        type(fci_in_t) :: fci_in
        type(reference_t) :: ref
        logical :: lanczos

        lua_state = flu_copyptr(l)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)
        lanczos = aot_exists(lua_state, opts, 'lanczos')
        call read_fci_in(lua_state, opts, sys%basis, fci_in)
        call read_reference_t(lua_state, opts, sys, ref)
        call aot_table_close(lua_state, opts)

        if (lanczos) then
            calc_type = lanczos_diag
            call do_fci_lanczos(sys, fci_in, ref, .false.)
        else
            calc_type = exact_diag
            call do_fci_lapack(sys, fci_in, ref)
        end if

    end function lua_fci

    function lua_hilbert_space(L) result(nresult) bind(c)

        ! Run a Monte Carlo calculation to estimate the size of the Hilbert
        ! space.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !     hilbert_args = {
        !         ncycles = n,
        !         ex_level = level,
        !         rng_seed = seed,
        !         reference = { ... },  -- nel-dimensional vector.
        !    }
        !
        !    hilbert_space {sys = sys_t, hilbert = hilbert_args}

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer

        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_exists, aot_table_close

        use errors, only: stop_all
        use lua_hande_system, only: get_sys_t
        use system, only: sys_t

        use calc, only: calc_type, mc_hilbert_space
        use hilbert_space, only: estimate_hilbert_space

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state

        type(c_ptr) :: sys_ptr
        type(sys_t), pointer :: sys
        integer :: truncation_level, ncycles, rng_seed
        integer, allocatable :: ref_det(:)
        integer :: opts, err
        logical :: have_seed

        lua_state = flu_copyptr(l)
        call get_sys_t(lua_state, sys)

        opts = aot_table_top(lua_state)
        call read_hilbert_args(lua_state, opts, sys%nel, ncycles, truncation_level, ref_det, rng_seed, have_seed)
        call aot_table_close(lua_state, opts)

        ! AOTUS returns a vector of size 0 to denote a non-existent vector.
        if (size(ref_det) == 0) deallocate(ref_det)
        if (allocated(ref_det)) then
            if (size(ref_det) /= sys%nel) call stop_all('lua_hilbert_space', &
                            'Reference determinant does not match the number of electrons in system.')
        end if

        calc_type = mc_hilbert_space
        if (have_seed) then
            call estimate_hilbert_space(sys, truncation_level, ncycles, ref_det, rng_seed)
        else
            call estimate_hilbert_space(sys, truncation_level, ncycles, ref_det)
        end if

        ! [todo] - return estimate of space and error to lua.
        nresult = 0

    end function lua_hilbert_space

    function lua_kinetic_energy(L) result(nresult) bind(c)

        ! Run a Monte Carlo calculation to estimate the kinetic energy of a
        ! system in the canonical ensemble.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    kin_args = {
        !        nattempts = npop,
        !        ncycles = n,
        !        beta = beta,
        !        fermi_temperature = true/false,
        !        rng_seed = seed
        !    }
        !
        !    kinetic_energy {sys = sys_t, kinetic = kin_args}

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use const, only: p

        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use calc, only: ms_in
        use errors, only: stop_all
        use lua_hande_system, only: get_sys_t
        use system, only: sys_t

        use calc, only: calc_type, mc_canonical_kinetic_energy
        use canonical_kinetic_energy, only: estimate_kinetic_energy

        integer(c_int) :: nresult
        type(c_ptr), value :: l

        type(flu_state) :: lua_state

        type(c_ptr) :: sys_ptr
        type(sys_t), pointer :: sys
        integer :: opts, err, rng_seed, ncycles, nattempts
        logical :: fermi_temperature, have_seed
        real(p) :: beta

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)

        if (sys%ueg%chem_pot == huge(1.0_p)) call stop_all('lua_kinetic_energy', &
                                                           'chem_pot: chemical potential not supplied.')

        opts = aot_table_top(lua_state)
        call read_kinetic_args(lua_state, opts, fermi_temperature, beta, nattempts, ncycles, rng_seed, have_seed)
        call aot_table_close(lua_state, opts)

        calc_type = mc_canonical_kinetic_energy
        if (have_seed) then
            call estimate_kinetic_energy(sys, fermi_temperature, beta, nattempts, ncycles, rng_seed)
        else
            call estimate_kinetic_energy(sys, fermi_temperature, beta, nattempts, ncycles)
        end if

        ! [todo] - return estimate of kinetic energy and error to lua.
        nresult = 0

    end function lua_kinetic_energy

    function lua_simple_fciqmc(L) result(nresult) bind(c)

        ! Run an FCIQMC calculation using the simple algorithm (slow, serial,
        ! memory hungry).

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    simple_fciqmc {
        !       sys = system,
        !       qmc = { ... },
        !       restart = { ... },
        !       reference = { ... },
        !    }

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use simple_fciqmc, only: do_simple_fciqmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: qmc_in_t, restart_in_t, reference_t
        use qmc, only: init_proc_pointers
        use system, only: sys_t

        use calc, only: ms_in
        use calc, only: calc_type, simple_fciqmc_calc, fciqmc_calc
        use system, only: set_spin_polarisation

        integer :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(restart_in_t) :: restart_in
        type(reference_t) :: reference

        integer :: opts
        character(10), parameter :: keys(4) = [character(10) :: 'sys', 'qmc', 'restart', 'reference']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in, .true.)
        call read_restart_in(lua_state, opts, restart_in)
        call read_reference_t(lua_state, opts, sys, reference)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = simple_fciqmc_calc + fciqmc_calc
        ! NOTE: sparse_hamil is currently disabled.
        call do_simple_fciqmc(sys, qmc_in, restart_in, reference, .false.)

        nresult = 0

    end function lua_simple_fciqmc

    function lua_fciqmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    fciqmc {
        !       sys = system,
        !       qmc = { ... },
        !       fciqmc = { ... },
        !       semi_stoch = { ... },
        !       restart = { ... },
        !       load_bal = { ... },
        !       reference = { ... },
        !    }

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use dmqmc_data, only: dmqmc_in_t
        use fciqmc, only: do_fciqmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: qmc_in_t, fciqmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t, reference_t
        use qmc, only: init_proc_pointers
        use system, only: sys_t

        use calc, only: ms_in
        use calc, only: calc_type, fciqmc_calc
        use system, only: set_spin_polarisation

        integer :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(fciqmc_in_t) :: fciqmc_in
        type(semi_stoch_in_t) :: semi_stoch_in
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference

        integer :: opts
        character(10), parameter :: keys(7) = [character(10) :: 'sys', 'qmc', 'fciqmc', 'semi_stoch', 'restart', &
                                                                'load_bal', 'reference']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)

        call read_qmc_in(lua_state, opts, qmc_in)
        call read_fciqmc_in(lua_state, opts, fciqmc_in)
        call read_semi_stoch_in(lua_state, opts, qmc_in, semi_stoch_in)
        call read_restart_in(lua_state, opts, restart_in)
        call read_load_bal_in(lua_state, opts, load_bal_in)
        call read_reference_t(lua_state, opts, sys, reference)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = fciqmc_calc
        call init_proc_pointers(sys, qmc_in, reference, fciqmc_in=fciqmc_in)
        call do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, load_bal_in, reference)

        nresult = 0

    end function lua_fciqmc

    function lua_ccmc(L) result(nresult) bind(c)

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    ccmc{
        !       sys = system,
        !       qmc = { ... },
        !       ccmc = { ... }
        !       restart = { ... },
        !       reference = { ... },
        !    }

        ! See interface documentation for the relevant read_TYPE procedure to
        ! understand the options available within a given subtable.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
        use flu_binding, only: flu_State, flu_copyptr
        use aot_table_module, only: aot_table_top, aot_table_close

        use dmqmc_data, only: dmqmc_in_t
        use ccmc, only: do_ccmc
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args
        use qmc_data, only: qmc_in_t, ccmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t, reference_t
        use qmc, only: init_proc_pointers
        use system, only: sys_t

        use calc, only: ms_in
        use calc, only: calc_type, ccmc_calc
        use system, only: set_spin_polarisation

        integer :: nresult
        type(c_ptr), value :: L

        type(flu_state) :: lua_state
        type(sys_t), pointer :: sys
        type(qmc_in_t) :: qmc_in
        type(ccmc_in_t) :: ccmc_in
        type(semi_stoch_in_t) :: semi_stoch_in
        type(restart_in_t) :: restart_in
        type(load_bal_in_t) :: load_bal_in
        type(reference_t) :: reference

        integer :: opts
        character(10), parameter :: keys(5) = [character(10) :: 'sys', 'qmc', 'ccmc', 'restart', 'reference']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        ! [todo] - do spin polarisation in system setup.
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        ! Get main table.
        opts = aot_table_top(lua_state)
        call read_qmc_in(lua_state, opts, qmc_in)
        call read_ccmc_in(lua_state, opts, ccmc_in)
        ! note that semi-stochastic is not (yet) available in CCMC.
        !call read_semi_stoch_in(lua_state, opts, qmc_in, semi_stoch_in)
        call read_restart_in(lua_state, opts, restart_in)
        ! load balancing is not available in CCMC; must use default settings.
        call read_reference_t(lua_state, opts, sys, reference)
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        calc_type = ccmc_calc
        call init_proc_pointers(sys, qmc_in, reference)
        call do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference)

        nresult = 0

    end function lua_ccmc

    ! --- table-derived type wrappers ---

    subroutine read_fci_in(lua_state, opts, basis, fci_in)

        ! Read in the fci and (optionally) lanczos tables to a fci_in_t object.

        ! fci = {
        !     write_hamiltonian = true/false,
        !     hamiltonian_file = filename,
        !     write_determinants = true/false,
        !     determinant_file = filename,
        !     nanalyse = N,
        !     blacs_block_size = block_size,
        !     rdm = { ... } -- L-d vector containing the sites to include in subsystem A.
        ! }
        ! lanczos = {
        !     neigv = N,
        !     nbasis = M,
        !     direct = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the fci and (optionally) the lanczos table(s).
        !    basis: information about the single-particle basis set of the system.
        ! Out:
        !    fci_in: fci_in_t object containing generic fci/lanczos input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use aot_vector_module, only: aot_get_val

        use basis_types, only: basis_t
        use checking, only: check_allocate
        use errors, only: stop_all
        use fci_utils, only: fci_in_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(basis_t), intent(in) :: basis
        type(fci_in_t), intent(inout) :: fci_in

        integer :: fci_table, err, fci_nrdms
        integer, allocatable :: err_arr(:)

        if (aot_exists(lua_state, opts, 'fci')) then
            call aot_table_open(lua_state, opts, fci_table, 'fci')

            ! Optional arguments
            call aot_get_val(fci_in%write_hamiltonian, err, lua_state, fci_table, 'write_hamiltonian')
            call aot_get_val(fci_in%hamiltonian_file, err, lua_state, fci_table, 'hamiltonian_file')
            call aot_get_val(fci_in%write_determinants, err, lua_state, fci_table, 'write_determinants')
            call aot_get_val(fci_in%determinant_file, err, lua_state, fci_table, 'determinant_file')
            call aot_get_val(fci_in%analyse_fci_wfn, err, lua_state, fci_table, 'nanalyse')
            call aot_get_val(fci_in%block_size, err, lua_state, fci_table, 'blacs_block_size')

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, fci_table, 'rdm')) then
                ! Currently restricted to one RDM in a single FCI calculation.
                fci_nrdms = 1
                allocate(fci_in%rdm_info(fci_nrdms), stat=err)
                call check_allocate('fci_in%rdm_info', fci_nrdms, err)
                call aot_get_val(fci_in%rdm_info(fci_nrdms)%subsystem_A, err_arr, basis%nbasis, lua_state)
                fci_in%rdm_info(fci_nrdms)%A_nsites = size(fci_in%rdm_info(fci_nrdms)%subsystem_A)
            end if

            call aot_table_close(lua_state, fci_table)

        end if

        ! Lanczos table: optional and indicates doing a Lanczos calculation.
        if (aot_exists(lua_state, opts, 'lanczos')) then
            call aot_table_open(lua_state, opts, fci_table, 'lanczos')
            call aot_get_val(fci_in%nlanczos_eigv, err, lua_state, fci_table, 'neigv')
            call aot_get_val(fci_in%lanczos_string_len, err, lua_state, fci_table, 'nbasis')
            call aot_get_val(fci_in%direct_lanczos, err, lua_state, fci_table, 'direct')
            call aot_table_close(lua_state, fci_table)
        end if

    end subroutine read_fci_in

    subroutine read_hilbert_args(lua_state, opts, nel, ncycles, ex_level, ref_det, rng_seed, have_seed)

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle to the table which is input to the Lua hilbert_space
        !        routine.
        !    nel: The number of electrons in the system.
        ! Out:
        !    ncycles: number of cycles i.e. number of random determinants to
        !        generate.
        !    ex_level: maximum excitation level relative to the reference
        !        determinant to include in the Hilbert space.
        !    ref_det: reference determinant.  If not supplied by the user then
        !        this will be deallocated on output.
        !    rng_seed: seed to initialise the random number generator.
        !    have_seed: True is the user inputs an RNG seed, false otherwise.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close
        use aot_vector_module, only: aot_get_val

        use const, only: p
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts, nel
        integer, intent(out) :: ncycles, ex_level, rng_seed
        integer, allocatable :: ref_det(:)
        logical, intent(out) :: have_seed

        integer :: hilbert_table, err
        integer, allocatable :: err_arr(:)
        character(10), parameter :: keys(5) = [character(10) :: 'sys', 'ncycles', 'ex_level', 'reference', 'rng_seed']

        if (.not. aot_exists(lua_state, opts, 'hilbert')) call stop_all('read_hilbert_args','"hilbert" table not present.')
        call aot_table_open(lua_state, opts, hilbert_table, 'hilbert')

        call aot_get_val(ncycles, err, lua_state, hilbert_table, 'ncycles')
        if (err /= 0) call stop_all('read_hilbert_args', 'Number of cycles not supplied.')
        call aot_get_val(ex_level, err, lua_state, hilbert_table, 'ex_level', default=-1)
        call aot_get_val(ref_det, err_arr, nel, lua_state, hilbert_table, key='reference')
        have_seed = aot_exists(lua_state, hilbert_table, 'rng_seed')
        call aot_get_val(rng_seed, err, lua_state, hilbert_table, 'rng_seed')

        call warn_unused_args(lua_state, keys, hilbert_table)
        call aot_table_close(lua_state, hilbert_table)

    end subroutine read_hilbert_args

    subroutine read_kinetic_args(lua_state, opts, fermi_temperature, beta, nattempts, ncycles, rng_seed, have_seed)

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle to the table which is input to the Lua kinetic_energy
        !        routine.
        ! Out:
        !    fermi_temperature: if true, rescale beta as the inverse reduced
        !        temperature.
        !    beta: target temperature.
        !    nattempts: number of samples to use each cycle.
        !    ncycles: number of Monte Carlo cycles to perform.
        !    rng_seed: seed to initialise the random number generator.
        !    have_seed: True is the user inputs an RNG seed, false otherwise.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use const, only: p
        use errors, only: stop_all
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        logical, intent(out) :: fermi_temperature, have_seed
        integer, intent(out) :: nattempts, ncycles, rng_seed
        real(p), intent(out) :: beta

        integer :: kinetic_table, err
        character(17), parameter :: keys(5) = [character(17) :: 'nattempts', 'ncycles', 'beta', &
                                                                'fermi_temperature', 'rng_seed']

        if (.not. aot_exists(lua_state, opts, 'kinetic')) call stop_all('read_kinetic_args','"kinetic" table not present.')
        call aot_table_open(lua_state, opts, kinetic_table, 'kinetic')

        call aot_get_val(nattempts, err, lua_state, kinetic_table, 'nattempts')
        if (err /= 0) call stop_all('read_kinetic_args', 'nattempts: Number of attempts/cycle not supplied.')
        call aot_get_val(ncycles, err, lua_state, kinetic_table, 'ncycles')
        if (err /= 0) call stop_all('read_kinetic_args', 'ncycles: Number of cycles not supplied.')

        call aot_get_val(beta, err, lua_state, kinetic_table, 'beta')
        if (err /= 0) call stop_all('read_kinetic_args', 'beta: target temperature not supplied.')
        call aot_get_val(fermi_temperature, err, lua_state, kinetic_table, 'fermi_temperature', default=.false.)

        have_seed = aot_exists(lua_state, kinetic_table, 'rng_seed')
        call aot_get_val(rng_seed, err, lua_state, kinetic_table, 'rng_seed')

        call warn_unused_args(lua_state, keys, kinetic_table)
        call aot_table_close(lua_state, kinetic_table)

    end subroutine read_kinetic_args

    subroutine read_qmc_in(lua_state, opts, qmc_in, short)

        ! Read in a qmc table to a qmc_in_t object.

        ! qmc = {
        !     tau = tau,                                  -- required
        !     init_pop = N,                               -- required
        !     mc_cycles = ncycles,                        -- required
        !     nreports = nreport,                         -- required
        !     state_size = walker_length,                 -- required
        !     spawned_state_size = spawned_walker_length, -- required
        !     rng_seed = seed,
        !     target_population = pop,
        !     real_amplitudes = true/false,
        !     spawn_cutoff = cutoff,
        !     no_renorm = true/false,
        !     tau_search = true/false,
        !     pattempt_single = prob,
        !     pattempt_double = prob,
        !     initial_shift = shift,
        !     shift_damping = damp_factor,
        !     initiator = true/false,
        !     initiator_threshold = pop,
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

        use calc, only: GLOBAL_META, gen_seed
        use qmc_data, only: qmc_in_t
        use lua_hande_utils, only: warn_unused_args
        use errors, only: stop_all

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(qmc_in_t), intent(out) :: qmc_in
        logical, intent(in), optional :: short

        integer :: qmc_table, err
        character(len=10) :: str
        logical :: skip

        if (present(short)) then
            skip = short
        else
            skip = .false.
        end if

        call aot_table_open(lua_state, opts, qmc_table, 'qmc')

        ! Required arguments
        call aot_get_val(qmc_in%tau, err, lua_state, qmc_table, 'tau')
        if (err /= 0) call stop_all('read_qmc_in', 'tau not set.')
        call aot_get_val(qmc_in%D0_population, err, lua_state, qmc_table, 'init_pop')
        if (err /= 0) call stop_all('read_qmc_in', 'init_pop not set.')
        call aot_get_val(qmc_in%ncycles, err, lua_state, qmc_table, 'mc_cycles')
        if (err /= 0) call stop_all('read_qmc_in', 'mc_cycles not set.')
        call aot_get_val(qmc_in%nreport, err, lua_state, qmc_table, 'nreports')
        if (err /= 0) call stop_all('read_qmc_in', 'nreports not set.')
        if (.not. skip) then
            call aot_get_val(qmc_in%walker_length, err, lua_state, qmc_table, 'state_size')
            if (err /= 0) call stop_all('read_qmc_in', 'state_size not set.')
            call aot_get_val(qmc_in%spawned_walker_length, err, lua_state, qmc_table, 'spawned_state_size')
            if (err /= 0) call stop_all('read_qmc_in', 'spawned_state_size not set.')
        end if

        ! Optional arguments (defaults set in derived type).
        call aot_get_val(qmc_in%real_amplitudes, err, lua_state, qmc_table, 'real_amplitudes')
        call aot_get_val(qmc_in%spawn_cutoff, err, lua_state, qmc_table, 'spawn_cutoff')
        call aot_get_val(qmc_in%no_renorm, err, lua_state, qmc_table, 'no_renorm')
        call aot_get_val(qmc_in%pattempt_single, err, lua_state, qmc_table, 'pattempt_single')
        call aot_get_val(qmc_in%pattempt_double, err, lua_state, qmc_table, 'pattempt_double')
        call aot_get_val(qmc_in%tau_search, err, lua_state, qmc_table, 'tau_search')
        call aot_get_val(qmc_in%initial_shift, err, lua_state, qmc_table, 'initial_shift')
        call aot_get_val(qmc_in%shift_damping, err, lua_state, qmc_table, 'shift_damping')
        call aot_get_val(qmc_in%target_particles, err, lua_state, qmc_table, 'target_population')
        call aot_get_val(qmc_in%initiator_approx, err, lua_state, qmc_table, 'initiator')
        call aot_get_val(qmc_in%initiator_pop, err, lua_state, qmc_table, 'initiator_threshold')

        ! If user sets initial shift and vary_shift_from, assume they know what
        ! they're doing.  Otherwise, vary the shift from the initial shift
        ! value.
        qmc_in%vary_shift_from = qmc_in%initial_shift

        ! Optional arguments requiring special care.
        if (aot_exists(lua_state, qmc_table, 'rng_seed')) then
            call aot_get_val(qmc_in%seed, err, lua_state, qmc_table, 'rng_seed')
        else
            qmc_in%seed = gen_seed(GLOBAL_META%uuid)
        end if
        if (aot_exists(lua_state, qmc_table, 'vary_shift_from')) then
            call aot_get_val(qmc_in%vary_shift_from, err, lua_state, qmc_table, 'vary_shift_from')
            if (err /= 0) then
                ! Perhaps passed in a string?
                call aot_get_val(str, err, lua_state, qmc_table, 'vary_shift_from')
                qmc_in% vary_shift_from_proje = trim(str) == 'proje'
            end if
        end if

        ! [todo] - check unused args.
        call warn_unused_args(lua_state, ['tau'], qmc_table)

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
        use errors, only: stop_all

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(fciqmc_in_t), intent(out) :: fciqmc_in

        integer :: fciqmc_table, ref_det, err
        character(len=12) :: str
        logical :: ref_det_flag

        if (aot_exists(lua_state, opts, 'fciqmc')) then

            call aot_table_open(lua_state, opts, fciqmc_table, 'fciqmc')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(fciqmc_in%non_blocking_comm, err, lua_state, fciqmc_table, 'non_blocking_comm')
            call aot_get_val(fciqmc_in%doing_load_balancing, err, lua_state, fciqmc_table, 'load_balancing')
            call aot_get_val(fciqmc_in%init_spin_inv_D0, err, lua_state, fciqmc_table, 'init_spin_inverse_reference_det')

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, fciqmc_table, 'select_reference_det')) then
                call aot_table_open(lua_state, fciqmc_table, ref_det, 'select_reference_det')
                if (ref_det == 0) then
                    ! Just passed a boolean (hopefully).
                    call aot_get_val(ref_det_flag, err, lua_state, ref_det, 'select_reference_det', default=.false.)
                    if (ref_det_flag) fciqmc_in%select_ref_det_every_nreports = 20
                else
                    call aot_table_close(lua_state, ref_det)
                    call aot_get_val(fciqmc_in%select_ref_det_every_nreports, err, lua_state, ref_det, 'update_every', default=20)
                    call aot_get_val(fciqmc_in%ref_det_factor, err, lua_state, ref_det, 'pop_factor')
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
                    call stop_all('read_fciqmc_in', 'Unknown guiding function')
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
                    write (6,*) trim(str)
                    call stop_all('read_fciqmc_in', 'Unknown trial wavefunction')
                end select
            end if

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], fciqmc_table)

            call aot_table_close(lua_state, fciqmc_table)

        end if

    end subroutine read_fciqmc_in

    subroutine read_semi_stoch_in(lua_state, opts, qmc_in, semi_stoch_in)

        ! Read in an semi_stoch table (if it exists) to an semi_stoch_in object.

        ! semi_stoch = {
        !     space = 'read'/'high',            -- required if semi_stoch table exists.
        !     size = target_pop,                -- required if space == 'high'
        !     start_iteration = mc_cycles,
        !     shift_start_iteration = mc_cycle, -- overrides start_iteration
        !     separate_annihilation = true/false,
        !     write_determ_space = true/false,
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
        use errors, only: stop_all

        use qmc_data, only: qmc_in_t, semi_stoch_in_t, high_pop_determ_space, read_determ_space
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(qmc_in_t), intent(inout) :: qmc_in
        type(semi_stoch_in_t), intent(out) :: semi_stoch_in

        integer :: semi_stoch_table, ref_det, err
        character(len=10) :: str
        logical :: ref_det_flag

        if (aot_exists(lua_state, opts, 'semi_stoch')) then

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
                    if (err /= 0) call stop_all('read_semi_stoch_in', 'Target space size not given.')
                end select
            else
                call stop_all('read_semi_stoch_in', 'Deterministic space not specified.')
            end if

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(semi_stoch_in%start_iter, err, lua_state, semi_stoch_table, 'start_iteration')
            call aot_get_val(semi_stoch_in%write_determ_space, err, lua_state, semi_stoch_table, 'write_determ_space')
            call aot_get_val(semi_stoch_in%separate_annihil, err, lua_state, semi_stoch_table, 'separate_annihilation')

            ! Optional arguments requiring special care.
            if (aot_exists(lua_state, semi_stoch_table, 'shift_start_iteration')) then
                semi_stoch_in%start_iter = huge(0) ! fixed up once the shift comes on...
                call aot_get_val(semi_stoch_in%shift_iter, err, lua_state, semi_stoch_table, 'shift_start_iteration')
            end if

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], semi_stoch_table)

            call aot_table_close(lua_state, semi_stoch_table)

        end if

    end subroutine read_semi_stoch_in

    subroutine read_ccmc_in(lua_state, opts, ccmc_in)

        ! Read in an ccmc table (if it exists) to an ccmc_in object.

        ! ccmc = {
        !     move_freq = frequency,
        !     cluster_multispawn_threshold = threshold,
        !     full_non_composite = true/false,
        !     linked = true/false,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the ccmc table.
        ! Out:
        !    ccmc_in: ccmc_in_t object containing ccmc-specific input options.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: ccmc_in_t
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(ccmc_in_t), intent(out) :: ccmc_in

        integer :: ccmc_table, ref_det, err
        character(len=10) :: str
        logical :: ref_det_flag

        if (aot_exists(lua_state, opts, 'ccmc')) then

            call aot_table_open(lua_state, opts, ccmc_table, 'ccmc')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(ccmc_in%move_freq, err, lua_state, ccmc_table, 'move_frequency')
            call aot_get_val(ccmc_in%cluster_multispawn_threshold, err, lua_state, ccmc_table, 'cluster_multispawn_threshold')
            call aot_get_val(ccmc_in%full_nc, err, lua_state, ccmc_table, 'full_non_composite')
            call aot_get_val(ccmc_in%linked, err, lua_state, ccmc_table, 'linked')

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], ccmc_table)

            call aot_table_close(lua_state, ccmc_table)

        end if

    end subroutine read_ccmc_in

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

        use qmc_data, only: load_bal_in_t
        use lua_hande_utils, only: warn_unused_args

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(load_bal_in_t), intent(out) :: load_bal_in

        integer :: load_bal_table, ref_det, err
        character(len=10) :: str
        logical :: ref_det_flag

        if (aot_exists(lua_state, opts, 'load_bal')) then

            call aot_table_open(lua_state, opts, load_bal_table, 'load_bal')

            ! Optional arguments (defaults set in derived type).
            call aot_get_val(load_bal_in%nslots, err, lua_state, load_bal_table, 'nslots')
            call aot_get_val(load_bal_in%pop, err, lua_state, load_bal_table, 'min_pop')
            call aot_get_val(load_bal_in%percent, err, lua_state, load_bal_table, 'target')
            call aot_get_val(load_bal_in%max_attempts, err, lua_state, load_bal_table, 'max_attempts')
            call aot_get_val(load_bal_in%write_info, err, lua_state, load_bal_table, 'write')

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], load_bal_table)

            call aot_table_close(lua_state, load_bal_table)

        end if

    end subroutine read_load_bal_in

    subroutine read_reference_t(lua_state, opts, sys, ref)

        ! Read in a reference table (if it exists) to a reference_t object.

        ! reference = {
        !     det = { ... }, -- N-electron vector
        !     hilbert_space_det = { ... }, -- N-electron vector
        !     ex_level = truncation_level,
        ! }

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    opts: handle for the table containing the reference table.
        !    sys: sys_t object describing the system of interest.
        ! Out:
        !    ref: reference_t object contianing the input options describing the
        !        reference.  Note that ex_level is set to the number of electrons
        !        if not provided (incl. if the reference table is not present in
        !        opts).

        use flu_binding, only: flu_State
        use aot_vector_module, only: aot_get_val
        use aot_table_module, only: aot_get_val, aot_exists, aot_table_open, aot_table_close

        use qmc_data, only: reference_t
        use system, only: sys_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(sys_t), intent(in) :: sys
        type(reference_t), intent(out) :: ref

        integer :: ref_table, err
        integer, allocatable :: err_arr(:)

        ! Set to full space by default.
        ref%ex_level = sys%nel

        if (aot_exists(lua_state, opts, 'reference')) then

            call aot_table_open(lua_state, opts, ref_table, 'reference')

            ! Optional arguments requiring special handling.

            ! Check if deeterminants exist first to avoid aotus allocating components unnecessarily.
            if (aot_exists(lua_state, ref_table, 'det')) then
                call aot_get_val(ref%occ_list0, err_arr, sys%nel, lua_state, ref_table, 'det')
            end if
            if (aot_exists(lua_state, ref_table, 'hilbert_space_det')) then
                call aot_get_val(ref%hs_occ_list0, err_arr, sys%nel, lua_state, ref_table, 'hilbert_space_det')
            end if

            call aot_get_val(ref%ex_level, err, lua_state, ref_table, 'ex_level')

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], load_bal_table)

            call aot_table_close(lua_state, ref_table)

        end if

    end subroutine read_reference_t

    subroutine read_restart_in(lua_state, opts, restart_in)

        ! Read in a restart table (if it exists) to a restart_in_t object.

        ! restart = {
        !    read = true/false/id,
        !    write = true/false/id,
        !    write_shift = true/false/id,
        !    write_frequency = niterations,
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

        use qmc_data, only: restart_in_t

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: opts
        type(restart_in_t), intent(out) :: restart_in

        integer :: err, restart_table

        if (aot_exists(lua_state, opts, 'restart')) then

            call aot_table_open(lua_state, opts, restart_table, 'restart')

            call aot_get_val(restart_in%write_freq, err, lua_state, restart_table, 'write_frequency')

            write (6,*) 'restart', restart_in%write_restart, restart_in%read_restart
            associate(r_in=>restart_in)
                call get_flag_and_id(lua_state, restart_table, r_in%read_restart, r_in%read_id, 'read')
                call get_flag_and_id(lua_state, restart_table, r_in%write_restart, r_in%write_id, 'write')
                call get_flag_and_id(lua_state, restart_table, r_in%write_restart_shift, r_in%write_shift_id, 'write_shift')
            end associate
            write (6,*) 'restart', restart_in%write_restart, restart_in%read_restart

            ! [todo] - check unused args.
            !call warn_unused_args(lua_state, [], load_bal_table)

            call aot_table_close(lua_state, restart_table)

        end if

        contains

            subroutine get_flag_and_id(lua_state, restart_table, flag, id, key)

                type(flu_State), intent(inout) :: lua_state
                integer, intent(in) :: restart_table
                logical, intent(inout) :: flag
                integer, intent(inout) :: id
                character(*), intent(in) :: key

                integer :: err

                if (aot_exists(lua_state, restart_table, key)) then
                    call aot_get_val(flag, err, lua_state, restart_table, key)
                    if (err /= 0) then
                        ! Passed an id instead.
                        flag = .true.
                        call aot_get_val(id, err, lua_state, restart_table, key)
                    end if
                end if

            end subroutine get_flag_and_id

    end subroutine read_restart_in

end module lua_hande_calc
