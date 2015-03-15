module qmc

! Launcher and initialisation routines for the various QMC algorithms.

use fciqmc_data

implicit none

contains

! --- QMC wrapper ---

    subroutine do_qmc(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, reference, load_bal_in, &
                      dmqmc_in)

        ! Initialise and run stochastic quantum chemistry procedures.

        ! In/Out:
        !    sys: system being studied.  This should(!) be returned unaltered on
        !         output from each procedure, but might be varied during the
        !         run if needed.
        !    qmc_in: input options relating to QMC methods.
        !    fciqmc_in: input options relating to FCIQMC.
        !    reference: the reference determinant.
        !    dmqmc_in: input options relating to DMQMC.
        ! In:
        !    ccmc_in: input options relating to CCMC.
        !    semi_stoch_in: Input options for the semi-stochastic adaptation.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.

        use calc

        use ccmc, only: do_ccmc
        use ct_fciqmc, only: do_ct_fciqmc
        use dmqmc, only: do_dmqmc
        use fciqmc, only: do_fciqmc
        use hellmann_feynman_sampling, only: do_hfs_fciqmc

        use qmc_data, only: qmc_in_t, fciqmc_in_t, ccmc_in_t, semi_stoch_in_t
        use qmc_data, only: restart_in_t, reference_t, load_bal_in_t
        use dmqmc_data, only: dmqmc_in_t
        use system, only: sys_t, copy_sys_spin_info, set_spin_polarisation
        use parallel, only: nprocs

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(fciqmc_in_t), intent(inout) :: fciqmc_in
        type(ccmc_in_t), intent(inout) :: ccmc_in
        type(semi_stoch_in_t), intent(in) :: semi_stoch_in
        type(restart_in_t), intent(in) :: restart_in
        type(reference_t), intent(inout) :: reference
        type(load_bal_in_t), intent(inout) :: load_bal_in
        type(dmqmc_in_t), intent(inout) :: dmqmc_in

        real(p) :: hub_matel
        type(sys_t) :: sys_bak

        ! Initialise procedure pointers.
        call init_proc_pointers(sys, qmc_in, dmqmc_in)

        ! Set spin variables.
        call copy_sys_spin_info(sys, sys_bak)
        call set_spin_polarisation(sys%basis%nbasis, ms_in, sys)

        ! Initialise data.
        call init_qmc(sys, qmc_in, fciqmc_in, restart_in, reference, load_bal_in, dmqmc_in)

        ! Calculation-specifc initialisation and then run QMC calculation.

        if (doing_calc(dmqmc_calc)) then
            call do_dmqmc(sys, qmc_in, dmqmc_in, restart_in, reference, load_bal_in)
        else if (doing_calc(ct_fciqmc_calc)) then
            call do_ct_fciqmc(sys, qmc_in, restart_in, reference, load_bal_in, hub_matel)
        else if (doing_calc(ccmc_calc)) then
            call do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, reference, load_bal_in)
        else
            ! Doing FCIQMC calculation (of some sort) using the original
            ! timestep algorithm.
            if (doing_calc(hfs_fciqmc_calc)) then
                call do_hfs_fciqmc(sys, qmc_in, restart_in, reference, load_bal_in)
            else
                call do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, load_bal_in, reference)
            end if
        end if

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine do_qmc

! --- Initialisation routines ---

    subroutine init_qmc(sys, qmc_in, fciqmc_in, restart_in, reference, load_bal_in, dmqmc_in)

        ! Initialisation for fciqmc calculations.
        ! Setup the spin polarisation for the system, initialise the RNG,
        ! allocate the required memory for the list of walkers and set the
        ! initial walker.

        ! In:
        !    sys: system being studied.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.
        ! In/Out:
        !    qmc_in: input options relating to QMC methods.
        !    fciqmc_in: input options relating to FCIQMC.
        !    reference: current reference determinant.
        !    dmqmc_in: input options relating to DMQMC.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use parallel
        use utils, only: int_fmt

        use basis, only: write_basis_fn
        use calc
        use dmqmc_procedures, only: init_dmqmc
        use determinants, only: decode_det, encode_det, write_det
        use energy_evaluation, only: calculate_hf_signed_pop, nparticles_start_ind
        use qmc_common, only: find_single_double_prob
        use reference_determinant, only: set_reference_det
        use hfs_data, only: O00, hf_signed_pop
        use proc_pointers, only: sc0_ptr, op0_ptr
        use spawn_data, only: alloc_spawn_t
        use spawning, only: assign_particle_processor
        use system
        use symmetry, only: symmetry_orb_list
        use momentum_symmetry, only: gamma_sym, sym_table
        use utils, only: factorial_combination_1
        use restart_hdf5, only: restart_info_global, read_restart_hdf5

        use qmc_data, only: qmc_in_t, fciqmc_in_t, restart_in_t, reference_t, load_bal_in_t
        use dmqmc_data, only: dmqmc_in_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(fciqmc_in_t), intent(inout) :: fciqmc_in
        type(restart_in_t), intent(in) :: restart_in
        type(reference_t), intent(inout) :: reference
        type(load_bal_in_t), intent(inout) :: load_bal_in
        type(dmqmc_in_t), intent(inout) :: dmqmc_in

        integer :: ierr
        integer :: i, j, D0_proc, D0_inv_proc, ipos, occ_list0_inv(sys%nel), slot
        integer :: step, size_main_walker, size_spawned_walker
        integer :: nwalker_int, nwalker_int_p, nwalker_real
        integer :: ref_sym ! the symmetry of the reference determinant
        integer(i0) :: f0_inv(sys%basis%string_len)
        integer(int_64) :: tmp_int_64

        if (parent) write (6,'(1X,a6,/,1X,6("-"),/)') 'FCIQMC'

        ! --- Array sizes depending upon QMC algorithms ---

        sampling_size = 1
        if (doing_calc(hfs_fciqmc_calc)) then
            sampling_size = sampling_size + 1
        else if (replica_tricks) then
            sampling_size = sampling_size + 1
        end if

        ! Each determinant occupies string_len kind=i0 integers,
        ! sampling_size kind=int_p integers, sampling_size kind=p reals and one
        ! integer. If the Neel singlet state is used as the reference state for
        ! the projected estimator, then a further 2 reals are used per
        ! determinant.
        if (trial_function == neel_singlet) then
            info_size = 2
        else
            info_size = 0
        end if
        nwalker_int = 1
        nwalker_int_p = sampling_size
        nwalker_real = sampling_size + info_size

        ! Thus the number of bits occupied by each determinant in the main
        ! walker list is given by string_len*i0_length+nwalker_int*32+
        ! nwalker_int_p*int_p_length+nwalker_real*32 (*64 if double precision).
        ! The number of bytes is simply 1/8 this.
#ifdef SINGLE_PRECISION
        size_main_walker = sys%basis%tensor_label_len*i0_length/8 + nwalker_int_p*int_p_length/8 + &
                           nwalker_int*4 + nwalker_real*4
#else
        size_main_walker = sys%basis%tensor_label_len*i0_length/8 + nwalker_int_p*int_p_length/8 + &
                           nwalker_int*4 + nwalker_real*8
#endif
        if (walker_length < 0) then
            ! Given in MB.  Convert.  Note: important to avoid overflow in the
            ! conversion!
            walker_length = int((-real(walker_length,p)*10**6)/size_main_walker)
        end if

        ! Each spawned_walker occupies spawned_size kind=int_s integers.
        if (qmc_in%initiator_approx) then
            size_spawned_walker = (sys%basis%tensor_label_len+sampling_size+1)*int_s_length/8
        else
            size_spawned_walker = (sys%basis%tensor_label_len+sampling_size)*int_s_length/8
        end if
        if (spawned_walker_length < 0) then
            ! Given in MB.  Convert.
            ! Note that we store 2 arrays.
            spawned_walker_length = int((-real(spawned_walker_length,p)*10**6)/(2*size_spawned_walker))
        end if

        if (parent) then
            write (6,'(1X,a53,f7.2)') &
                'Memory allocated per core for main walker list (MB): ', &
                size_main_walker*real(walker_length,p)/10**6
            write (6,'(1X,a57,f7.2)') &
                'Memory allocated per core for spawned walker lists (MB): ', &
                size_spawned_walker*real(2*spawned_walker_length,p)/10**6
            write (6,'(1X,a48,'//int_fmt(walker_length,1)//')') &
                'Number of elements per core in main walker list:', walker_length
            write (6,'(1X,a51,'//int_fmt(spawned_walker_length,1)//',/)') &
                'Number of elements per core in spawned walker list:', spawned_walker_length
        end if

        ! --- Memory allocation ---

        ! Allocate main walker lists.
        allocate(nparticles(sampling_size), stat=ierr)
        call check_allocate('nparticles', sampling_size, ierr)
        allocate(tot_nparticles(sampling_size), stat=ierr)
        call check_allocate('tot_nparticles', sampling_size, ierr)
        allocate(walker_dets(sys%basis%tensor_label_len,walker_length), stat=ierr)
        call check_allocate('walker_dets', sys%basis%string_len*walker_length, ierr)
        allocate(walker_population(sampling_size,walker_length), stat=ierr)
        call check_allocate('walker_population', sampling_size*walker_length, ierr)
        allocate(walker_data(sampling_size+info_size,walker_length), stat=ierr)
        call check_allocate('walker_data', size(walker_data), ierr)
        allocate(nparticles_proc(sampling_size, nprocs), stat=ierr)
        call check_allocate('nparticles_proc', nprocs*sampling_size, ierr)

        ! Allocate spawned walker lists and spawned walker times (ct_fciqmc only)
        if (mod(spawned_walker_length, nprocs) /= 0) then
            if (parent) write (6,'(1X,a68)') 'spawned_walker_length is not a multiple of the number of processors.'
            spawned_walker_length = ceiling(real(spawned_walker_length)/nprocs)*nprocs
            if (parent) write (6,'(1X,a35,'//int_fmt(spawned_walker_length,1)//',a1,/)') &
                                        'Increasing spawned_walker_length to',spawned_walker_length,'.'
        end if

        ! Allocate the shift.
        allocate(shift(sampling_size), stat=ierr)
        call check_allocate('shift', size(shift), ierr)
        shift = qmc_in%initial_shift

        allocate(vary_shift(sampling_size), stat=ierr)
        call check_allocate('vary_shift', size(vary_shift), ierr)

        if (doing_calc(ct_fciqmc_calc)) then
            allocate(spawn_times(spawned_walker_length),stat=ierr)
            call check_allocate('spawn_times',spawned_walker_length,ierr)
        end if

        ! Set the real encoding shift, depending on whether 32 or 64-bit integers
        ! are being used.
        if (qmc_in%real_amplitudes) then
            if (bit_size(0_int_p) == 64) then
                ! Allow a maximum population of 2^32, and a minimum fractional
                ! part of 2^-31.
                real_bit_shift = 31
            else if (bit_size(0_int_p) == 32) then
                ! Allow a maximum population of 2^20, and a minimum fractional
                ! part of 2^-11.
                real_bit_shift = 11
            end if
        else
            ! Allow no fractional part for walker populations.
            real_bit_shift = 0
        end if
        ! Store 2**real_bit_shift for ease.
        real_factor = 2_int_p**(int(real_bit_shift, int_p))

        ! If not using real amplitudes then we always want spawn_cutoff to be
        ! equal to 1.0, so overwrite the default.
        if (.not. qmc_in%real_amplitudes) spawn_cutoff = 0.0_p

        call alloc_spawn_t(sys%basis%tensor_label_len, sampling_size, qmc_in%initiator_approx, &
                         spawned_walker_length, spawn_cutoff, real_bit_shift, 7, use_mpi_barriers, qmc_spawn)
        if (fciqmc_in%non_blocking_comm) then
            call alloc_spawn_t(sys%basis%tensor_label_len, sampling_size, qmc_in%initiator_approx, &
                               spawned_walker_length, spawn_cutoff, real_bit_shift, 7, .false., received_list)
        end if

        if (nprocs == 1 .or. .not. fciqmc_in%doing_load_balancing) load_bal_in%nslots = 1
        call init_parallel_t(sampling_size, nparticles_start_ind-1, fciqmc_in%non_blocking_comm, par_info, &
                             load_bal_in%nslots)

        allocate(reference%f0(sys%basis%string_len), stat=ierr)
        call check_allocate('reference%f0',sys%basis%string_len,ierr)
        allocate(reference%hs_f0(sys%basis%string_len), stat=ierr)
        call check_allocate('reference%hs_f0', size(reference%hs_f0), ierr)

        ! --- Number of reports ---
        ! When using the propagate_to_beta option the number of iterations in imaginary
        ! time we want to do depends on what value of beta we are seeking. It's
        ! annoying to have to modify this in the input file, so just do it here.
        if (propagate_to_beta) qmc_in%nreport = int(ceiling(init_beta/(qmc_in%ncycles*qmc_in%tau)))

        ! --- Initial walker distributions ---
        ! Note occ_list could be set and allocated in the input.

        if (restart_in%read_restart) then
            if (.not.allocated(reference%occ_list0)) then
                allocate(reference%occ_list0(sys%nel), stat=ierr)
                call check_allocate('reference%occ_list0',sys%nel,ierr)
            end if
            call read_restart_hdf5(restart_info_global, fciqmc_in%non_blocking_comm, reference)
            ! Need to re-calculate the reference determinant data
            call decode_det(sys%basis, reference%f0, reference%occ_list0)
            if (trial_function == neel_singlet) then
                ! Set the Neel state data for the reference state, if it is being used.
                reference%H00 = 0.0_p
            else
                reference%H00 = sc0_ptr(sys, reference%f0)
            end if
            if (doing_calc(hfs_fciqmc_calc)) O00 = op0_ptr(sys, reference%f0)
        else

            ! Reference det

            ! Set the reference determinant to be the spin-orbitals with the lowest
            ! single-particle eigenvalues which satisfy the spin polarisation.
            ! Note: this is for testing only!  The symmetry input is currently
            ! ignored.
            if (sym_in < sys%sym_max) then
                call set_reference_det(sys, reference%occ_list0, .false., sym_in)
            else
                call set_reference_det(sys, reference%occ_list0, .false.)
            end if

            call encode_det(sys%basis, reference%occ_list0, reference%f0)

            if (allocated(reference%hs_occ_list0)) then
                call encode_det(sys%basis, reference%hs_occ_list0, reference%hs_f0)
            else
                allocate(reference%hs_occ_list0(sys%nel), stat=ierr)
                call check_allocate('reference%hs_occ_list0', size(reference%hs_occ_list0), ierr)
                reference%hs_occ_list0 = reference%occ_list0
                reference%hs_f0 = reference%f0
            end if

            ! Energy of reference determinant.
            reference%H00 = sc0_ptr(sys, reference%f0)
            if (doing_calc(hfs_fciqmc_calc)) O00 = op0_ptr(sys, reference%f0)

            ! In general FCIQMC, we start with psips only on the
            ! reference determinant, so set tot_walkers = 1 and
            ! initialise walker_population. For DMQMC, this is
            ! not required, as psips are spawned along the diagonal
            ! initially.
            if (.not.doing_calc(dmqmc_calc)) then
                tot_walkers = 1
                ! Zero all populations...
                walker_population(:,tot_walkers) = 0_int_p
                ! Set initial population of Hamiltonian walkers.
                walker_population(1,tot_walkers) = nint(qmc_in%D0_population)*real_factor
                ! Set the bitstring of this psip to be that of the
                ! reference state.
                walker_dets(:,tot_walkers) = reference%f0

                ! Determine and set properties for the reference state which we start on.
                ! By definition, when using a single determinant as a reference state:
                walker_data(1,tot_walkers) = 0.0_p
                ! Or if not using a single determinant:
                if (trial_function == neel_singlet) then
                    ! Set the Neel state data for the reference state, if it is being used.
                    walker_data(1,tot_walkers) = reference%H00
                    reference%H00 = 0.0_p

                    walker_data(sampling_size+1,tot_walkers) = sys%lattice%nsites/2
                    ! For a rectangular bipartite lattice, nbonds = ndim*nsites.
                    ! The Neel state cannot be used for non-bipartite lattices.
                    walker_data(sampling_size+2,tot_walkers) = sys%lattice%ndim*sys%lattice%nsites
                end if

                ! Finally, we need to check if the reference determinant actually
                ! belongs on this processor.
                ! If it doesn't, set the walkers array to be empty.
                call assign_particle_processor(reference%f0, sys%basis%string_len, qmc_spawn%hash_seed, &
                                               qmc_spawn%hash_shift, qmc_spawn%move_freq, nprocs, &
                                               D0_proc, slot, load_bal_in%nslots)
                if (D0_proc /= iproc) tot_walkers = 0
            end if

            ! For the Heisenberg model and open shell systems, it is often useful to
            ! have psips start on both the reference state and the spin-flipped version.
            if (fciqmc_in%init_spin_inv_D0) then

                ! Need to handle the Heisenberg model (consisting of spinors on
                ! lattice sites) and electron systems differently, as the
                ! Heisenberg model has no concept of unoccupied basis
                ! functions/holes.
                select case (sys%system)
                case (heisenberg)
                    ! Flip all spins in f0 to get f0_inv
                    f0_inv = not(reference%f0)
                    ! In general, the basis bit string has some padding at the
                    ! end which must be unset.  We need to clear this...
                    ! Loop over all bits after the last basis function.
                    i = sys%basis%bit_lookup(2,sys%basis%nbasis)
                    do ipos = sys%basis%bit_lookup(1,sys%basis%nbasis)+1, i0_end
                        f0_inv(i) = ibclr(f0_inv(i), ipos)
                    end do
                case default
                    ! Swap each basis function for its spin-inverse
                    ! This looks somewhat odd, but relies upon basis
                    ! functions alternating down (Ms=-1) and up (Ms=1).
                    do i = 1, sys%nel
                        if (mod(reference%occ_list0(i),2) == 1) then
                            occ_list0_inv(i) = reference%occ_list0(i) + 1
                        else
                            occ_list0_inv(i) = reference%occ_list0(i) - 1
                        end if
                    end do
                    call encode_det(sys%basis, occ_list0_inv, f0_inv)
                end select

                call assign_particle_processor(f0_inv, sys%basis%string_len, qmc_spawn%hash_seed, &
                                               qmc_spawn%hash_shift, qmc_spawn%move_freq, nprocs, &
                                               D0_inv_proc, slot, load_bal_in%nslots)

                ! Store if not identical to reference det.
                if (D0_inv_proc == iproc .and. any(reference%f0 /= f0_inv)) then
                    tot_walkers = tot_walkers + 1
                    ! Zero all populations for this determinant.
                    walker_population(:,tot_walkers) = 0_int_p
                    ! Set the population for this basis function.
                    walker_population(1,tot_walkers) = nint(qmc_in%D0_population)*real_factor
                    walker_data(1,tot_walkers) = sc0_ptr(sys, reference%f0) - reference%H00
                    select case(sys%system)
                    case(heisenberg)
                        if (trial_function /= single_basis) then
                            walker_data(1,tot_walkers) = 0
                        else
                            walker_data(1,tot_walkers) = sc0_ptr(sys, reference%f0) - reference%H00
                        end if
                    case default
                        walker_data(1,tot_walkers) = sc0_ptr(sys, reference%f0) - reference%H00
                    end select
                    walker_dets(:,tot_walkers) = f0_inv
                    ! If we are using the Neel state as a reference in the
                    ! Heisenberg model, then set the required data.
                    if (trial_function == neel_singlet) then
                        walker_data(sampling_size+1,tot_walkers) = 0
                        walker_data(sampling_size+2,tot_walkers) = 0
                    end if
                end if
            end if

        end if ! End of initialisation of reference state(s)/restarting from previous calculations

        ! Total number of particles on processor.
        ! Probably should be handled more simply by setting it to be either 0 or
        ! D0_population or obtaining it from the restart file, as appropriate.
        forall (i=1:sampling_size) nparticles(i) = sum(abs( real(walker_population(i,:tot_walkers),p)/real_factor))
        ! Should we already be in varyshift mode (e.g. restarting a calculation)?
#ifdef PARALLEL
        do i=1, sampling_size
            call mpi_allgather(nparticles(i), 1, MPI_PREAL, nparticles_proc(i,:), 1, MPI_PREAL, MPI_COMM_WORLD, ierr)
        end do
        ! When restarting a non-blocking calculation this sum will not equal
        ! tot_nparticles as some walkers have been communicated around the report
        ! loop. The correct total is in the restart file so get it from there.
        if (.not. restart_in%read_restart) forall(i=1:sampling_size) tot_nparticles(i) = sum(nparticles_proc(i,:))
#else
        tot_nparticles = nparticles
        nparticles_proc(:sampling_size,1) = nparticles(:sampling_size)
#endif

        ! Decide whether the shift should be turned on from the start.
        vary_shift = tot_nparticles >= qmc_in%target_particles

        if (doing_calc(hfs_fciqmc_calc)) then
#ifdef PARALLEL
            tmp_int_64 = calculate_hf_signed_pop()
            call mpi_allreduce(tmp_int_64, hf_signed_pop, sampling_size, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
            hf_signed_pop = calculate_hf_signed_pop()
#endif
        end if

        ! calculate the reference determinant symmetry
        ref_sym = symmetry_orb_list(sys, reference%occ_list0)

        ! If not set at input, set probability of selecting single or double
        ! excitations based upon the reference determinant and assume other
        ! determinants have a roughly similar ratio of single:double
        ! excitations.
        if (qmc_in%pattempt_single < 0 .or. qmc_in%pattempt_double < 0) then
            call find_single_double_prob(sys, reference%occ_list0, qmc_in%pattempt_single, qmc_in%pattempt_double)
        else
            ! renormalise just in case input wasn't
            qmc_in%pattempt_single = qmc_in%pattempt_single/(qmc_in%pattempt_single+qmc_in%pattempt_double)
            qmc_in%pattempt_double = 1.0_p - qmc_in%pattempt_single
        end if

        ! Calculate all the possible different amplitudes for the Neel singlet state
        ! and store them in an array
        if (trial_function == neel_singlet) then
            allocate(neel_singlet_amp(-1:(sys%lattice%nsites/2)+1), stat=ierr)
            call check_allocate('neel_singlet_amp',(sys%lattice%nsites/2)+1,ierr)

            neel_singlet_amp(-1) = 0
            neel_singlet_amp((sys%lattice%nsites/2)+1) = 0
            do i=0,(sys%lattice%nsites/2)
                neel_singlet_amp(i) = factorial_combination_1( (sys%lattice%nsites/2)-i , i )
                neel_singlet_amp(i) = -(2*mod(i,2)-1) * neel_singlet_amp(i)
            end do
        end if

        ! When doing a DMQMC calculation, call a routine to initialise all the required
        ! arrays, ie to store thermal quantities, and to initalise reduced density matrix
        ! quantities if necessary.
        if (doing_calc(dmqmc_calc)) then
            call init_dmqmc(sys, qmc_in, dmqmc_in)
        end if

        if (parent) then
            write (6,'(1X,a29,1X)',advance='no') 'Reference determinant, |D0> ='
            call write_det(sys%basis, sys%nel, reference%f0, new_line=.true.)
            write (6,'(1X,a16,f20.12)') 'E0 = <D0|H|D0> =',reference%H00
            if (doing_calc(hfs_fciqmc_calc)) write (6,'(1X,a17,f20.12)') 'O00 = <D0|O|D0> =',O00
            write(6,'(1X,a34)',advance='no') 'Symmetry of reference determinant:'
            select case(sys%system)
            case (hub_k)
                call write_basis_fn(sys, sys%basis%basis_fns(2*ref_sym), new_line=.true., print_full=.false.)
            case default
                write(6,'(i2)') ref_sym
            end select
            write (6,'(1X,a46,1X,f8.4)') 'Probability of attempting a single excitation:', qmc_in%pattempt_single
            write (6,'(1X,a46,1X,f8.4)') 'Probability of attempting a double excitation:', qmc_in%pattempt_double
            if (doing_calc(dmqmc_calc)) then
                write (6,'(1X,a54,'//int_fmt(int(qmc_in%D0_population,int_64),1)//')') &
                              'Initial population on the trace of the density matrix:', int(qmc_in%D0_population,int_64)
            else
                write (6,'(1X,a44,1X,f11.4,/)') &
                              'Initial population on reference determinant:',qmc_in%D0_population
                write (6,'(1X,a53,/)') 'Note that the correlation energy is relative to |D0>.'
            end if
            if (qmc_in%initiator_approx) then
                write (6,'(1X,a24)') 'Initiator method in use.'
                write (6,'(1X,a48,1X,f3.1,/)') &
                    'Population for a determinant to be an initiator:', qmc_in%initiator_pop
            end if
            write (6,'(1X,a46,/)') 'Information printed out every QMC report loop:'
            write (6,'(1X,a69)') 'Note that all particle populations are averaged over the report loop.'
            write (6,'(1X,a58,/)') 'Shift: the shift calculated at the end of the report loop.'
            if (.not. doing_calc(dmqmc_calc)) then
                write (6,'(1X,a46)') 'H_0j: <D_0|H|D_j>, Hamiltonian matrix element.'
                write (6,'(1X,a60)') 'N_j: population of Hamiltonian particles on determinant D_j.'
                if (doing_calc(hfs_fciqmc_calc)) then
                    write (6,'(1X,a43)') 'O_0j: <D_0|O|D_j>, operator matrix element.'
                    write (6,'(1X,a67)') "N'_j: population of Hellmann--Feynman particles on determinant D_j."
                end if
            else
                if (doing_dmqmc_calc(dmqmc_full_r2)) then
                    write (6, '(1X,a104)') 'Trace: The current total population on the diagonal elements of the &
                                         &first replica of the density matrix.'
                    write (6, '(1X,a107)') 'Trace 2: The current total population on the diagonal elements of the &
                                         &second replica of the density matrix.'
                else
                    write (6, '(1X,a83)') 'Trace: The current total population on the diagonal elements of the &
                                         &density matrix.'
                end if
                if (doing_dmqmc_calc(dmqmc_full_r2)) then
                    write (6, '(1X,a81)') 'Full S2: The numerator of the estimator for the Renyi entropy of the &
                                          &full system.'
                end if
                if (doing_dmqmc_calc(dmqmc_energy)) then
                    write (6, '(1X,a92)') '\sum\rho_{ij}H_{ji}: The numerator of the estimator for the expectation &
                                         &value of the energy.'
                end if
                if (doing_dmqmc_calc(dmqmc_energy_squared)) then
                    write (6, '(1X,a100)') '\sum\rho_{ij}H2{ji}: The numerator of the estimator for the expectation &
                                         &value of the energy squared.'
                end if
                if (doing_dmqmc_calc(dmqmc_correlation)) then
                    write (6, '(1X,a111)') '\sum\rho_{ij}S_{ji}: The numerator of the estimator for the expectation &
                                         &value of the spin correlation function.'
                end if
                if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
                    write (6, '(1X,a109)') '\sum\rho_{ij}M2{ji}: The numerator of the estimator for the expectation &
                                         &value of the staggered magnetisation.'
                end if
                if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
                    write (6, '(1x,a73)') 'RDM(n) S2: The numerator of the estimator for the Renyi entropy of RDM n.'
                end if
                if (calc_inst_rdm) then
                    write (6, '(1x,a83)') 'RDM(n) trace m: The current total population on the diagonal of replica m &
                                          &of RDM n.'
                end if
                if (calc_excit_dist) then
                    write (6, '(1x,a86)') 'Excit. level n: The fraction of particles on excitation level n of the &
                                          &density matrix.'
                end if
            end if
            write (6,'(1X,a61)') '# H psips: current total population of Hamiltonian particles.'
            if (doing_calc(hfs_fciqmc_calc)) then
                write (6,'(1X,a68)') '# HF psips: current total population of Hellmann--Feynman particles.'
            end if
            write (6,'(1X,"# states: number of many-particle states occupied.")')
            if (.not. fciqmc_in%non_blocking_comm) then
                write (6,'(1X,"# spawn_events: number of successful spawning events across all processors.")')
            end if
            write (6,'(1X,a56,/)') 'R_spawn: average rate of spawning across all processors.'
            write (6,'(1X,a41,/)') 'time: average time per Monte Carlo cycle.'
        end if

    end subroutine init_qmc

    subroutine init_proc_pointers(sys, qmc_in, dmqmc_in)

        ! Set function pointers for QMC calculations.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    dmqmc_in: input options relating to DMQMC.

        ! System and calculation data
        use calc
        use hfs_data
        use system
        use parallel, only: parent
        use qmc_data, only: qmc_in_t
        use dmqmc_data, only: dmqmc_in_t

        ! Procedures to be pointed to.
        use death, only: stochastic_death
        use determinants
        use dmqmc_estimators
        use dmqmc_procedures
        use energy_evaluation
        use excit_gen_mol
        use excit_gen_op_mol
        use excit_gen_hub_k
        use excit_gen_op_hub_k
        use excit_gen_real_lattice
        use excit_gen_ueg, only: gen_excit_ueg_no_renorm
        use hamiltonian_chung_landau, only: slater_condon0_chung_landau
        use hamiltonian_hub_k, only: slater_condon0_hub_k
        use hamiltonian_hub_real, only: slater_condon0_hub_real
        use hamiltonian_heisenberg, only: diagonal_element_heisenberg, diagonal_element_heisenberg_staggered
        use hamiltonian_molecular, only: slater_condon0_mol
        use hamiltonian_ueg, only: slater_condon0_ueg, kinetic_energy_ueg
        use heisenberg_estimators
        use ifciqmc, only: set_parent_flag, set_parent_flag_dummy
        use importance_sampling
        use operators
        use spawning

        ! Procedure pointers
        use proc_pointers

        ! Utilities
        use errors, only: stop_all

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in

        ! 0. In general, use the default spawning routine.
        spawner_ptr => spawn_standard

        ! 1. Set system-specific procedure pointers.
        !     * projected energy estimator
        !     * diagonal hamiltonian matrix element evaluation
        !     * spawning
        !     * excitation generators
        select case(sys%system)
        case(hub_k)

            decoder_ptr => decode_det_spinocc_spinunocc
            update_proj_energy_ptr => update_proj_energy_hub_k
            sc0_ptr => slater_condon0_hub_k
            spawner_ptr => spawn_lattice_split_gen
            if (qmc_in%no_renorm) then
                gen_excit_ptr%full => gen_excit_hub_k_no_renorm
                gen_excit_ptr%init => gen_excit_init_hub_k_no_renorm
                gen_excit_ptr%finalise => gen_excit_finalise_hub_k_no_renorm
            else
                gen_excit_ptr%full => gen_excit_hub_k
                gen_excit_ptr%init => gen_excit_init_hub_k
                gen_excit_ptr%finalise => gen_excit_finalise_hub_k
            end if

        case(hub_real, chung_landau)

            ! The Hubbard model in a local orbital basis and the Chung--Landau
            ! Hamiltonian have the same off-diagonal operator so we can use the
            ! same excitation generators and just a different function for the
            ! diagonal matrix elements.
            ! Note that the Chung--Landau model (as in Phys Rev B 85 (2012)
            ! 115115) is contains spinless fermions.
            decoder_ptr => decode_det_occ
            update_proj_energy_ptr => update_proj_energy_hub_real
            if (sys%system == hub_real) then
                sc0_ptr => slater_condon0_hub_real
            else
                sc0_ptr => slater_condon0_chung_landau
            end if

            if (qmc_in%no_renorm) then
                gen_excit_ptr%full => gen_excit_hub_real_no_renorm
            else
                gen_excit_ptr%full => gen_excit_hub_real
            end if

        case(heisenberg)

            ! Only need occupied orbitals list, as for the real Hubbard case.
            decoder_ptr => decode_det_occ
            ! Set which trial wavefunction to use for the energy estimator.
            select case(trial_function)
            case (single_basis)
                update_proj_energy_ptr => update_proj_energy_heisenberg_basic
            case (neel_singlet)
                update_proj_energy_ptr => update_proj_energy_heisenberg_neel_singlet
            end select

            ! Set whether the applied staggered magnetisation is non-zero.
            if (abs(sys%heisenberg%staggered_magnetic_field) > depsilon) then
                sc0_ptr => diagonal_element_heisenberg_staggered
            else
                sc0_ptr => diagonal_element_heisenberg
            end if

            ! Set which guiding wavefunction to use, if requested.
            if (qmc_in%no_renorm) then
                gen_excit_ptr%full => gen_excit_heisenberg_no_renorm
            else
                    gen_excit_ptr%full => gen_excit_heisenberg
            end if
            select case(guiding_function)
            case (neel_singlet_guiding)
                spawner_ptr => spawn_importance_sampling
                gen_excit_ptr%trial_fn => neel_trial_state
            end select

        case(read_in)

            update_proj_energy_ptr => update_proj_energy_mol
            sc0_ptr => slater_condon0_mol

            if (qmc_in%no_renorm) then
                gen_excit_ptr%full => gen_excit_mol_no_renorm
                decoder_ptr => decode_det_occ
            else
                gen_excit_ptr%full => gen_excit_mol
                decoder_ptr => decode_det_occ_symunocc
            end if

        case(ueg)

            update_proj_energy_ptr => update_proj_energy_ueg
            sc0_ptr => slater_condon0_ueg

            if (qmc_in%no_renorm) then
                gen_excit_ptr%full => gen_excit_ueg_no_renorm
                decoder_ptr => decode_det_occ
            else
                if (parent) then
                    write (6,'(1X,"WARNING: renormalised excitation generators not implemented.")')
                    write (6,'(1X,"WARNING: If this upsets you, please send patches (or bribe James with beer).",/)')
                end if
                gen_excit_ptr%full => gen_excit_ueg_no_renorm
                decoder_ptr => decode_det_occ
            end if

        case default
            call stop_all('init_proc_pointers','QMC not implemented for this system yet.')
        end select

        ! 2. Set calculation-specific procedure pointers

        ! 2: initiator-approximation
        if (qmc_in%initiator_approx) then
            set_parent_flag_ptr => set_parent_flag
            if (all(ras > 0)) then
                create_spawned_particle_ptr => create_spawned_particle_initiator_ras
            else if (truncate_space) then
                create_spawned_particle_ptr => create_spawned_particle_initiator_truncated
            else
                create_spawned_particle_ptr => create_spawned_particle_initiator
            end if
        else
            set_parent_flag_ptr => set_parent_flag_dummy
            if (all(ras > 0)) then
                create_spawned_particle_ptr => create_spawned_particle_ras
            else if (truncate_space) then
                create_spawned_particle_ptr => create_spawned_particle_truncated
            else
                create_spawned_particle_ptr => create_spawned_particle
            end if
        end if

        ! 2: density-matrix
        if (doing_calc(dmqmc_calc)) then

            ! Spawned particle creation. 
            if (half_density_matrix) then
                if (truncate_space) then
                    create_spawned_particle_dm_ptr => create_spawned_particle_truncated_half_density_matrix
                else
                    create_spawned_particle_dm_ptr => create_spawned_particle_half_density_matrix
                end if
            else
                if (truncate_space) then
                    create_spawned_particle_dm_ptr => create_spawned_particle_truncated_density_matrix
                else
                    create_spawned_particle_dm_ptr => create_spawned_particle_density_matrix
                end if
            end if

            ! Weighted importance sampling routines.
            if (dmqmc_in%weighted_sampling) then
                spawner_ptr => spawn_importance_sampling
                gen_excit_ptr%trial_fn => dmqmc_weighting_fn
            end if

            ! Expectation values.
            if (doing_dmqmc_calc(dmqmc_energy)) then
                if (propagate_to_beta) then
                    update_dmqmc_energy_and_trace_ptr => dmqmc_energy_and_trace_propagate
                else
                    update_dmqmc_energy_and_trace_ptr => dmqmc_energy_and_trace
                end if
            end if
            select case(sys%system)
            case(heisenberg)
                if (doing_dmqmc_calc(dmqmc_energy_squared)) &
                                         update_dmqmc_energy_squared_ptr => dmqmc_energy_squared_heisenberg
                if (doing_dmqmc_calc(dmqmc_correlation)) update_dmqmc_correlation_ptr => &
                                     dmqmc_correlation_function_heisenberg
                if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) &
                                         update_dmqmc_stag_mag_ptr => dmqmc_stag_mag_heisenberg
            case(ueg)
                if (propagate_to_beta) then
                    if (free_electron_trial) then
                        trial_dm_ptr => kinetic_energy_ueg
                    else
                        trial_dm_ptr => slater_condon0_ueg
                    end if
                end if
            case(hub_k)
                if (propagate_to_beta) then
                    if (free_electron_trial) then
                        trial_dm_ptr => kinetic0_hub_k
                    else
                        trial_dm_ptr => slater_condon0_hub_k
                    end if
                end if
            end select

        end if

        ! 2: Hellmann--Feynman operator sampling
        if (doing_calc(hfs_fciqmc_calc)) then
            select case(hf_operator)
            case(hamiltonian_operator)
                op0_ptr => sc0_ptr
                update_proj_hfs_ptr => update_proj_hfs_hamiltonian
                spawner_hfs_ptr => spawner_ptr
            case(kinetic_operator)
                update_proj_hfs_ptr => update_proj_hfs_diagonal
                spawner_hfs_ptr => spawn_null
                if (sys%system == hub_k) then
                    op0_ptr => kinetic0_hub_k
                else
                    call stop_all('init_proc_pointers','System not yet supported in HFS with operator given.')
                end if
            case(double_occ_operator)
                if (sys%system == hub_k) then
                    ! Shamelessly re-use the Hamiltonian excitation generators.
                    gen_excit_hfs_ptr%full => gen_excit_ptr%full
                    gen_excit_hfs_ptr%init => gen_excit_ptr%init
                    gen_excit_hfs_ptr%finalise => gen_excit_ptr%finalise
                    spawner_hfs_ptr => spawn_lattice_split_gen_importance_sampling
                    ! Scale the Hamiltonian matrix element to obtain the matrix
                    ! element of this operator.
                    gen_excit_hfs_ptr%trial_fn => gen_excit_double_occ_matel_hub_k
                    update_proj_hfs_ptr => update_proj_hfs_double_occ_hub_k
                    op0_ptr => double_occ0_hub_k
                else
                    call stop_all('init_proc_pointers','System not yet supported in HFS with operator given.')
                end if
            case(dipole_operator)
                if (sys%system == read_in) then
                    op0_ptr => one_body0_mol
                    update_proj_hfs_ptr => update_proj_hfs_one_body_mol
                    spawner_hfs_ptr => spawner_ptr
                    if (qmc_in%no_renorm) then
                        gen_excit_hfs_ptr%full => gen_excit_one_body_mol_no_renorm
                    else
                        gen_excit_hfs_ptr%full => gen_excit_one_body_mol
                    end if
                else
                    call stop_all('init_proc_pointers','System not yet supported in HFS with operator given.')
                end if
            case default
                call stop_all('init_proc_pointers','Operator given is not yet supported')
            end select
        end if

    end subroutine init_proc_pointers

end module qmc
