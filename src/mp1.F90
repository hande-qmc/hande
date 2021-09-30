module mp1

use const, only: p, int_64
use qmc_data, only: excit_gen_renorm

implicit none

! todo
! #. test/verify/fix
!      + MP2 energy
!      + convergence?
! #. fix todos
! #. output (pretty table, JSON)
! #. write restart file
! #. document

! [todo] - document
type mp1_in_t
    real(p) :: D0_norm
    integer(int_64) :: nattempts
    integer :: ncycles
    integer :: state_size
    logical :: real_amplitudes = .false.
    real(p) :: spawn_cutoff = 0.01_p

    ! Copied from qmc_data.f90::qmc_in_t
    ! BZ-[todo]: add actual input options for this for excit_gens other than default.
    integer :: excit_gen = excit_gen_renorm
    real(p) :: pattempt_single = -1, pattempt_double = -1
    logical :: pattempt_update = .false.
    logical :: pattempt_zero_accum_data = .false.
    real(p) :: pattempt_parallel = -1
end type mp1_in_t



contains

    subroutine sample_mp1_wfn(sys, mp1_in, ref_in, logging_in, tijab, rng_seed)

        ! Generate the mp1 wavefunction (1+T_2)|HF> exactly and use it as an initial guess for FCIQMC/CCMC
        ! In:
        !   sys: system being studied.
        !   mpi_in: input options relating to mp1.
        !   ref_in: current reference determinant, set in lua_hande_calc.
        !   logging_in: (currently untested) logging inputs.
        !   rng_seed:
        ! Out:
        !   tijab: particle_t object for use in FCIQMC/CCMC calculation

        use const, only: int_p, i0, depsilon

        use checking, only: check_allocate
        use errors, only: stop_all
        use json_out
        use dSFMT_interface, only: dSFMT_t, dSFMT_init, get_rand_close_open
        use parallel

        use system, only: sys_t, copy_sys_spin_info, sys_t_json
        use calc_system_init, only: set_spin_polarisation
        use qmc_data, only: particle_t, reference_t, annihilation_flags_t
        use qmc_data, only: load_bal_state_t, reference_t_json, qmc_in_t, qmc_state_t
        use qmc, only: init_proc_pointers, init_reference, init_excit_gen
        use ccmc_data, only: cluster_t, ex_lvl_dist_t
        use determinants, only: det_info_t, alloc_det_info_t
        use spawn_data, only: spawn_t, proc_map_t, alloc_spawn_t, dealloc_spawn_t

        use annihilation, only: direct_annihilation, remove_unoccupied_dets
        use calc, only: GLOBAL_META, gen_seed
        use load_balancing, only: init_proc_map_t
        use ccmc_utils, only: cumulative_population, find_D0
        use ccmc_selection, only: select_cluster
        use determinants, only: encode_det
        use determinant_decoders, only: decode_det_occ
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use hamiltonian, only: get_hmatel
        use hamiltonian_data
        use particle_t_utils, only: init_particle_t, dealloc_particle_t
        use proc_pointers, only: decoder_ptr
!        use qmc, only: init_sc0_ptr
        use reference_determinant, only: set_reference_det
        use search, only: binary_search
        use stoch_utils, only: stochastic_round_spawned_particle
        use spawning, only: create_spawned_particle, attempt_to_spawn, assign_particle_processor
        use logging, only: logging_t, logging_in_t
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(inout) :: sys
        type(mp1_in_t), intent(in) :: mp1_in
        type(reference_t), intent(in), optional :: ref_in
        type(logging_in_t), intent(in) :: logging_in
        type(particle_t), intent(out) :: tijab
        integer, intent(in), optional :: rng_seed

        type(spawn_t) :: spawn
        type(reference_t) :: ref
        type(det_info_t) :: cdet
        type(cluster_t) :: cluster
        type(excit_t), parameter :: null_excit = excit_t(0,[0,0],[0,0],.false.)
        type(excit_t) :: excit
        type(annihilation_flags_t) :: annihilation_flags
        type(proc_map_t) :: proc_map
        type(dSFMT_t) :: rng
        type(sys_t) :: sys_bak

        integer :: seed, max_nspawned_states, nspaces, state_size, psip_element_size, homo, lumo
        integer :: i, j, a, b, ia, ib, iocc_a, iocc_b, ierr, pos, istate, icycle
        integer(int_64) :: iattempt
        logical :: hit
        integer(int_p) :: nspawn, old_pop
        real(p) :: norm, spawn_cutoff, amplitude, denom
        type(hmatel_t) :: intgrl, ampl
        real(dp) :: emp2
        integer(i0) :: f(sys%basis%bit_string_len)
        integer :: excitor_sign, excitor_level
        integer :: max_cluster_size, slot
        real(p), allocatable :: cumulative_abs_real_pops(:)
        real(p) :: tot_abs_real_pop
        type(json_out_t) :: js
        integer :: io_unit
        type(ex_lvl_dist_t) :: ex_lvl_dist
        integer :: D0_pos, D0_proc, nD0_proc
        complex(p) :: D0_normalisation
        type(logging_t) :: logging_info
        type(excit_gen_data_t) :: excit_gen_data
        type(qmc_in_t) :: qmc_in_cast

        io_unit = 6
        if (parent) then
            write (io_unit,'(1X,"MP1 wavefunction MC")')
            write (io_unit,'(1X,"-------------------",/)')
        end if

        if (present(rng_seed)) then
            seed = rng_seed
        else
            seed = gen_seed(GLOBAL_META%uuid)
        end if
        call dSFMT_init(seed+iproc, 50000, rng)


        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, sys)

        if (debug) call init_logging(logging_in, logging_info, ref%ex_level)
        
        qmc_in_cast%excit_gen = mp1_in%excit_gen
        qmc_in_cast%pattempt_single = mp1_in%pattempt_single
        qmc_in_cast%pattempt_double = mp1_in%pattempt_double
        qmc_in_cast%pattempt_update = mp1_in%pattempt_update
        qmc_in_cast%pattempt_zero_accum_data = mp1_in%pattempt_zero_accum_data
        qmc_in_cast%pattempt_parallel = mp1_in%pattempt_parallel
        call init_proc_pointers(sys, qmc_in_cast, ref, io_unit)
        call init_reference(sys, ref_in, io_unit, ref)
        call init_excit_gen(sys, qmc_in_cast, ref, .false., excit_gen_data)
        
        ! Similarly for particle stores.  We use
        ! * tijab: exact amplitudes on reference and doubles.
        state_size = mp1_in%state_size
        if (state_size < 0) then
            ! Given in MB.  Convert to elements.
            ! bit string and population element.
            psip_element_size = (sys%basis%tensor_label_len*i0_length + int_p_length)/8
            ! data element
#ifdef SINGLE_PRECISION
            psip_element_size = psip_element_size + 4
#else
            psip_element_size = psip_element_size + 8
#endif
            state_size = -int((real(state_size,p)*10**6)/psip_element_size)
        end if

        tijab%nspaces = 1
        tijab%info_size = 0
        call init_particle_t(state_size, 0, sys%basis%tensor_label_len, .true., .false., tijab) ! Strictly speaking tijab need only hold the double amplitudes.

        ! Legacy global data (boo!) initialisation
!        call init_sc0_ptr(sys)

        ! annihilation_flags default to 'off'.  Only feature we might be using is real amplitudes.
        annihilation_flags%real_amplitudes = mp1_in%real_amplitudes
        ! Just need the simplest decoder as not doing anything based on the clusters produced beyond setting a bit string...
        decoder_ptr => decode_det_occ

        ! Initialise spawn store.  Use sane defaults for options irrelevant to sampling a wavefunction (i.e. disable initiator, only
        ! one particle space, default processor map, no MPI barriers).
        nspaces = 1
        spawn_cutoff = mp1_in%spawn_cutoff
        call init_proc_map_t(1, proc_map)
        if (.not. mp1_in%real_amplitudes) spawn_cutoff = 0.0_p
        ! Need at least one spawning slot per attempt so give ourselves a generous amount of breathing room for load inbalance.
        max_nspawned_states = 2*ceiling(real(mp1_in%state_size)/nprocs)*nprocs
        call alloc_spawn_t(sys%basis%tensor_label_len, sys%basis%nbasis, nspaces, .false., max_nspawned_states, spawn_cutoff, &
                           tijab%pop_real_factor, proc_map, 7, .false., spawn)
        ! spawn_t now holds all the processor map/load balancing info required so can safely deallocate proc_map.
        deallocate(proc_map%map)

        ! Start with the reference
        tijab%nstates = 1
        tijab%states(:,tijab%nstates) = ref%f0
        tijab%pops(1,tijab%nstates) = nint(mp1_in%D0_norm)*tijab%pop_real_factor
        ! We define the energy zero to be <HF|H|HF>.  The annihilation framework evaluates <D|H|D> for other determinants.
        tijab%dat(1,tijab%nstates) = 0.0_p
        tijab%nparticles(1) = real(sum(abs(tijab%pops(1,:tijab%nstates))),p)/tijab%pop_real_factor

        call alloc_det_info_t(sys, cdet)
        allocate(cluster%excitors(ref%ex_level+2))

        if (parent) then
            call json_object_init(js, tag=.true.)
            call sys_t_json(js, sys)
            call mp1_in_t_json(js, mp1_in)
            call reference_t_json(js, ref, terminal=.true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io,'()')
        end if

        homo = 1
        lumo = sys%basis%nbasis
        associate(bl=>sys%basis%bit_lookup, bfns=>sys%basis%basis_fns)
            do i = 2, sys%basis%nbasis
                if (btest(ref%f0(bl(2,i)), bl(1,i))) then
                    ! Occupied
                    if (bfns(i)%sp_eigv  > bfns(homo)%sp_eigv) homo = i
                else
                        if (bfns(i)%sp_eigv  < bfns(lumo)%sp_eigv) homo = i
                end if
            end do
            if (abs(bfns(lumo)%sp_eigv-bfns(homo)%sp_eigv) < depsilon) then
                call stop_all('sample_mp1_wfn', 'Cannot generate MP1 wavefunction for a system with a zero HOMO-LUMO gap.')
            end if
        end associate

        ! Start with generating N_0(1+T2) deterministically.
        excit%nexcit = 2
        iocc_a = 0
        iocc_b = 0
        emp2 = 0.0_dp
        do i = 1, sys%nel
            do j = i+1, sys%nel
                excit%from_orb = [ref%occ_list0(i), ref%occ_list0(j)]
                ! Trivial approach to parallelisation: split a over processors.
                ! (Particles are distributed by the annihilation framework anyway.)
                do ia = 1, sys%basis%nbasis-sys%nel, nprocs
                    ! Find the next unoccupied orbital.
                    do
                        a = ia + iocc_a
                        if (.not.btest(ref%f0(sys%basis%bit_lookup(2,a)), sys%basis%bit_lookup(1,a))) exit
                        iocc_a = iocc_a + 1
                    end do
                    do ib = ia+1, sys%basis%nbasis-sys%nel
                        ! Find the next unoccupied orbital greater than a.
                        ! We could optimise this by restricting the choice of b to conserve spin and spatial symmetries...
                        do
                            b = ib + iocc_b
                            if (.not.btest(ref%f0(sys%basis%bit_lookup(2,b)), sys%basis%bit_lookup(1,b))) exit
                            iocc_b = iocc_b + 1
                        end do
                        excit%to_orb = [a,b]
                        call create_excited_det(sys%basis, ref%f0, excit, f)
                        ! t_{ijab} = <ij||ab> / (e_i + e_j - e_a - e_b).
                        ! NOTE: sp_eigv is not the Hartree--Fock eigenvalues in all cases (e.g. Hubbard model in k-space, UEG!).
                        associate(bf => sys%basis%basis_fns)
                            intgrl = get_hmatel(sys, ref%f0, f)
                            denom = (1.0_p / (bf(i)%sp_eigv + bf(j)%sp_eigv - bf(a)%sp_eigv - bf(b)%sp_eigv))
                            ampl = intgrl * denom
                            if (sys%read_in%comp) then
                                emp2 = emp2 + real(conjg(intgrl%c)*ampl%c)
                            else
                                emp2 = emp2 + intgrl%r*ampl%r
                            end if
                            ampl = ampl * mp1_in%D0_norm
                        end associate
                        ! Note attempt_to_spawn creates particles of opposite sign, hence set parent_sign=-1 to undo this unwanted
                        ! sign change.
                        ! TODO - cope with complex valued spawns
                        nspawn = attempt_to_spawn(rng, 1.0_p, spawn%cutoff, tijab%pop_real_factor, ampl%r, 1.0_p, -1_int_p)
                        if (nspawn /= 0_int_p) call create_spawned_particle(sys%basis, ref, cdet, excit, nspawn, 1, spawn, f)
                    end do
                end do
            end do
        end do
        write (6,*) 'direct MP2 =', emp2

        call direct_annihilation(sys, rng, ref, annihilation_flags, tijab, spawn)
        if (parent) write (6,'()')

        call dealloc_spawn_t(spawn)

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine sample_mp1_wfn

    subroutine mp1_status(sys, ref, D0_norm, iteration, psip_list)

        use const, only: dp, int_p
        use qmc_data, only: reference_t, particle_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: ref
        real(p), intent(in) :: D0_norm
        integer, intent(in) :: iteration
        type(particle_t), intent(in) :: psip_list

        real(dp) :: emp2
        real(p) :: tot_nparticles
        integer :: tot_nstates

        tot_nstates = psip_list%nstates
        tot_nparticles = psip_list%nparticles(1)
        emp2 = calc_emp2(sys, ref%f0, D0_norm, psip_list)
#ifdef PARALLEL
        ! No attempt to optimise comms here.  This should not be done a huge number of times...
        ! [todo] - MPI
#endif

        write (6,'(1X,i6,2X,i12,2(2X,f16.8))') iteration, tot_nstates, tot_nparticles, emp2

    end subroutine mp1_status

    function calc_emp2(sys, f0, D0_norm, psip_list) result(emp2)

        ! Calculatates the projected energy of the system
        
        use const, only: dp, int_p, i0, depsilon
        use qmc_data, only: particle_t
        use system, only: sys_t
        use hamiltonian, only: get_hmatel
        use hamiltonian_data, only: hmatel_t
        use errors, only: stop_all

        real(dp) :: emp2
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: D0_norm
        type(particle_t), intent(in) :: psip_list

        integer :: istate
        type(hmatel_t) :: hmatel

        emp2 = 0.0_dp
        do istate = 1, psip_list%nstates
            if (any(psip_list%states(:,istate) /= f0)) then
                hmatel = get_hmatel(sys, f0, psip_list%states(:,istate))
                if (sys%read_in%comp) then
                    if (aimag(hmatel%c*psip_list%pops(2,istate)) > depsilon) then
                        ! Calling stop_all means we can't make this function pure, 
                        ! might be worth it for sanity and mp2 isn't the time consuming 
                        ! part of the calculation anyway
                        call stop_all('calc_emp2', 'Non-zero complex energy detected')
                    else
                        emp2 = emp2 + real(hmatel%c)*psip_list%pops(1,istate) + aimag(hmatel%c)*psip_list%pops(2,istate)
                    end if
                else
                    emp2 = emp2 + hmatel%r*psip_list%pops(1,istate)
                end if
            end if
        end do
        emp2 = emp2 / (psip_list%pop_real_factor*D0_norm)

    end function calc_emp2

    subroutine mp1_in_t_json(js, mp1_in, terminal)

        ! Serialise a mpi_in_t object in JSON format.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   mp1_in: mp1_in_t object containing MP1 sampling input values (including any defaults set).
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out

        type(json_out_t), intent(inout) :: js
        type(mp1_in_t), intent(in) :: mp1_in
        logical, intent(in), optional :: terminal

        call json_object_init(js, 'mp1')
        call json_write_key(js, 'D0_norm', mp1_in%D0_norm)
        call json_write_key(js, 'ncycles', mp1_in%ncycles)
        call json_write_key(js, 'nattempts', mp1_in%nattempts)
        call json_write_key(js, 'state_size', mp1_in%state_size)
        call json_write_key(js, 'real_amplitudes', mp1_in%real_amplitudes)
        call json_write_key(js, 'spawn_cutoff', mp1_in%spawn_cutoff, terminal=.true.)
        call json_object_end(js, terminal)

    end subroutine mp1_in_t_json

end module mp1
