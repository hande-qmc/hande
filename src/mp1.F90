module mp1

use const, only: p, int_64

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
end type mp1_in_t

contains

    subroutine sample_mp1_wfn(sys, mp1_in, ref_in, rng_seed)

        ! [todo] - document

        use const, only: int_p, i0, depsilon

        use checking, only: check_allocate
        use errors, only: stop_all
        use json_out
        use dSFMT_interface, only: dSFMT_t, dSFMT_init, get_rand_close_open
        use parallel

        use system, only: sys_t, copy_sys_spin_info, sys_t_json
        use calc_system_init, only: set_spin_polarisation
        use qmc_data, only: particle_t, reference_t, annihilation_flags_t, load_bal_state_t, reference_t_json
        use qmc, only: init_reference
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

        type(sys_t), intent(inout) :: sys
        type(mp1_in_t), intent(in) :: mp1_in
        type(reference_t), intent(in), optional :: ref_in
        integer, intent(in), optional :: rng_seed

        type(particle_t) :: tijab, psip_list
        type(spawn_t) :: spawn
        type(reference_t) :: ref
        type(det_info_t) :: cdet
        type(cluster_t) :: cluster
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

        io_unit = 6
        if (parent) then
            write (io_unit,'(1X,"MP1 wavefunction MC")')
            write (io_unit,'(1X,"-------------------",/)')
        end if

        if (debug) call init_logging(logging_in, logging_info, qs%ref%ex_level)

        if (present(rng_seed)) then
            seed = rng_seed
        else
            seed = gen_seed(GLOBAL_META%uuid)
        end if
        call dSFMT_init(seed+iproc, 50000, rng)


        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, sys)

        call init_reference(sys, ref_in, io_unit, ref)

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
        max_nspawned_states = 2*ceiling(real(mp1_in%nattempts)/nprocs)*nprocs
        call alloc_spawn_t(sys%basis%tensor_label_len, sys%basis%nbasis, nspaces, .false., max_nspawned_states, spawn_cutoff, &
                           psip_list%pop_real_factor, proc_map, 7, .false., spawn)
        ! spawn_t now holds all the processor map/load balancing info required so can safely deallocate proc_map.
        deallocate(proc_map%map)

        ! Similarly for particle stores.  We use
        ! * tijab: exact amplitudes on reference and doubles.
        ! * psip_list:  current set of amplitudes based on sampling e^{T2}|HF>.
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
            state_size = int((real(state_size,p)*10**6)/psip_element_size)
        end if
        tijab%nspaces = 1
        tijab%info_size = 0
        call init_particle_t(state_size, 0, sys%basis%tensor_label_len, .true., .false., tijab) ! Strictly speaking tijab need only hold the double amplitudes.
        psip_list%nspaces = 1
        psip_list%info_size = 0
        call init_particle_t(state_size, 0, sys%basis%tensor_label_len, mp1_in%real_amplitudes, .false., psip_list) ! [todo] - real_amplitude_force_32 option

        ! Start with the reference
        tijab%nstates = 1
        tijab%states(:,tijab%nstates) = ref%f0
        tijab%pops(1,tijab%nstates) = nint(mp1_in%D0_norm)*tijab%pop_real_factor
        ! We define the energy zero to be <HF|H|HF>.  The annihilation framework evaluates <D|H|D> for other determinants.
        tijab%dat(1,tijab%nstates) = 0.0_p
        tijab%nparticles(1) = real(sum(abs(tijab%pops(1,:tijab%nstates))),p)/tijab%pop_real_factor

        call alloc_det_info_t(sys, cdet)
        allocate(cluster%excitors(ref%ex_level+2))
        allocate(cumulative_abs_real_pops(size(psip_list%states,dim=2)), stat=ierr)
        call check_allocate('cumulative_abs_real_pops', size(cumulative_abs_real_pops), ierr)

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
                            emp2 = emp2 + real(conjg(intgrl%c)*ampl%c)
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
        spawn%head = spawn%head_start

        psip_list%nstates = tijab%nstates
        psip_list%states(:,1:psip_list%nstates) = tijab%states(:,1:psip_list%nstates)
        psip_list%pops(:,1:psip_list%nstates) = tijab%pops(:,1:psip_list%nstates)
        psip_list%dat(:,1:psip_list%nstates) = tijab%dat(:,1:psip_list%nstates)
        psip_list%nparticles = tijab%nparticles

        ! Exact MP1 wavefunction, e^T_2|HF> = 1 + T_2 + 1/2 T_2^2 + ...|HF>
        ! Sample this by randomly generating clusters from the (1+T_2)|HF> space with the appropriate amplitude.
        ! [todo] - describe in more detail.

        ! Set of amplitudes on reference and doubles are fixed so the the reference position and cumulative population information
        ! won't change for the set of excitors from which we sample clusters...
        call assign_particle_processor(ref%f0, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, &
                                       nprocs, D0_proc, slot, spawn%proc_map%map, spawn%proc_map%nslots)
        if (iproc == D0_proc) then
            D0_pos = -1
            nD0_proc = 1
            call find_D0(psip_list, ref%f0, D0_pos)
            D0_normalisation = cmplx(tijab%pops(1,D0_pos), tijab%pops(2,D0_pos), p)/tijab%pop_real_factor
        else
            D0_pos = -1
            nD0_proc = 0 ! No reference excitor on the processor.
        end if
#ifdef PARALLEL
        call mpi_bcast(D0_normalisation, 1, mpi_pcomplex, D0_proc, MPI_COMM_WORLD, ierr)
#endif
        call cumulative_population(tijab%pops, tijab%states(sys%basis%tot_string_len, :), tijab%nstates, D0_proc, D0_pos, &
             tijab%pop_real_factor, .false., sys%read_in%comp, cumulative_abs_real_pops, &
                                   tot_abs_real_pop, ex_lvl_dist)
        excit = excit_t(0, [0,0], [0,0], .false.)
        max_cluster_size = min(sys%nel, ref%ex_level+2, tijab%nstates-nD0_proc)

        call mp1_status(sys, ref, mp1_in%D0_norm, 0, psip_list)

        do icycle = 1, mp1_in%ncycles
            ! Note again minimal attempt at parallelisation.  This bit is meant to be fast compared to an FCIQMC/CCMC calculation!
            do iattempt = 1, mp1_in%nattempts, nprocs
                ! This is pretty crude: we generate clusters in the 1+T_2 space (despite already having that exactly) but also
                ! generate, with the appropriate weight, higher order terms as well.  It's convenient, despite the redundant work,
                ! to generate 1+T_2 again to avoid any normalisation issues in the cluster generation probabilities.
                call select_cluster(rng, sys, tijab, ref%f0, ref%ex_level, .false., mp1_in%nattempts, D0_normalisation, 0.0_p, &
                                    cumulative_abs_real_pops, tot_abs_real_pop, 0, max_cluster_size, cdet, cluster)
                if (cluster%excitation_level <= ref%ex_level) then
                    ! Create a particle with the amplitude of the cluster with the probability of selecting the cluster.
                    amplitude = cluster%amplitude*cluster%cluster_to_det_sign
                    nspawn = attempt_to_spawn(rng, 1.0_p, spawn%cutoff, psip_list%pop_real_factor, amplitude, &
                                              cluster%pselect, -1_int_p)
                    if (nspawn /= 0_int_p) call create_spawned_particle(sys%basis, ref, cdet, excit, nspawn, 1, spawn, cdet%f)
                end if
            end do

            call direct_annihilation(sys, rng, ref, annihilation_flags, psip_list, spawn)
            spawn%head = spawn%head_start

            ! Internal normalisation: reset N_0 to original value.
            if (iproc == D0_proc) then
                D0_pos = -1
                call find_D0(psip_list, ref%f0, D0_pos)
                norm = (real(psip_list%pops(1,D0_pos),p)/psip_list%pop_real_factor)/mp1_in%D0_norm
            end if
#ifdef PARALLEL
            call mpi_bcast(norm, 1, mpi_preal, D0_proc, MPI_COMM_WORLD, ierr)
#endif
            do istate = 1, psip_list%nstates
                ! Note: we remain in encoded format.
                amplitude = real(psip_list%pops(1,istate),p) / norm
                old_pop = psip_list%pops(1,istate)
                if (amplitude > 0.0_p) then
                    psip_list%pops(1,istate) = stochastic_round_spawned_particle(spawn%cutoff, amplitude, rng)
                else
                    psip_list%pops(1,istate) = -stochastic_round_spawned_particle(spawn%cutoff, -amplitude, rng)
                end if
                psip_list%nparticles(1) = psip_list%nparticles(1) &
                                              + real(abs(psip_list%pops(1,istate))-abs(old_pop),p)/psip_list%pop_real_factor
            end do
            call remove_unoccupied_dets(rng, psip_list, mp1_in%real_amplitudes)

            call mp1_status(sys, ref, mp1_in%D0_norm, icycle, psip_list)

        end do
        if (parent) write (6,'()')

        call dealloc_spawn_t(spawn)
        call dealloc_particle_t(tijab)
        call dealloc_particle_t(psip_list)

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

    pure function calc_emp2(sys, f0, D0_norm, psip_list) result(emp2)

        use const, only: dp, int_p, i0
        use qmc_data, only: particle_t
        use system, only: sys_t
        use hamiltonian, only: get_hmatel

        real(dp) :: emp2
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: D0_norm
        type(particle_t), intent(in) :: psip_list

        integer :: istate

        emp2 = 0.0_dp
        do istate = 1, psip_list%nstates
            if (any(psip_list%states(:,istate) /= f0)) then
                emp2 = emp2 + get_hmatel(sys, f0, psip_list%states(:,istate)) * psip_list%pops(1,istate)
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
