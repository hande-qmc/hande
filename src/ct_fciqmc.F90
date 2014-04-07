module ct_fciqmc

! Evolve the walker population using a continuous time algorithm (i.e. jump
! directly to the next spawning event without a timestep).

use fciqmc_data
use const, only: p, lint

implicit none

contains

    subroutine do_ct_fciqmc(sys, matel)

        ! In:
        !    sys: system being studied
        !    matel: off-diagonal Hamiltonian matrix element (ignoring sign due
        !       to permutations).  Either U (Bloch orbitals) or
        !       t (atomic/real-space orbitals).

        use annihilation, only: direct_annihilation
        use basis, only: basis_length
        use calc, only: seed, initiator_approximation
        use determinants, only: det_info, alloc_det_info
        use excitations, only: excit
        use qmc_common
        use proc_pointers
        use system, only: sys_t, hub_real, hub_k
        use interact

        use checking
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use parallel
        use utils, only: rng_init_info
        use restart_hdf5, only: restart_info_global, dump_restart_hdf5

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: matel ! either U or t, depending whether we are working in the real or k-space

        integer :: nspawned, ndeath, ireport, idet
        integer(lint) :: nattempts
        real(dp) :: nparticles_old(sampling_size)
        integer :: iparticle, tmp_pop, max_nexcitations, ierr, proc_id
        integer, allocatable :: current_pos(:,:) ! (nthreads, 0:max(1,nprocs-1))
        real(p) :: time, t_barrier, K_ii, R, sum_off_diag
        real :: t1, t2
        type(det_info) :: cdet
        type(excit) :: connection
        type(excit), allocatable :: connection_list(:)
        logical :: soft_exit
        real(p):: hmatel
        type(excit) :: D0_excit
        type(dSFMT_t) :: rng
        integer, parameter :: thread_id = 0
        integer :: spawned_pop

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

        ! index of spawning array which contains population
        spawned_pop = basis_length + 1

        if (sys%system == hub_k) then
            associate(sl=>sys%lattice)
                max_nexcitations = sys%nalpha*sys%nbeta*min(sl%nsites-sys%nalpha,sl%nsites-sys%nbeta)
            end associate
        else if (sys%system == hub_real) then
            max_nexcitations = 2*sys%lattice%ndim*sys%nel
        end if

        allocate(connection_list(max_nexcitations), stat=ierr)
        call check_allocate('connection_list', max_nexcitations, ierr)
        allocate(current_pos(nthreads,0:max(1,nprocs-1)), stat=ierr)
        call check_allocate('current_pos', size(current_pos), ierr)

        sum_off_diag = max_nexcitations*matel

        call alloc_det_info(sys, cdet)

        nparticles_old = tot_nparticles

        t_barrier = tau ! or we could just not bother with the t_barrier var...

        if (parent) call write_fciqmc_report_header()
        call initial_fciqmc_status(sys)

        ! time the report loop
        call cpu_time(t1)

        ! Main fciqmc loop
        do ireport = 1, nreport

            call init_report_loop()
            call init_mc_cycle(nattempts, ndeath)

            ! Loop over determinants in the walker list.
            do idet = 1, tot_walkers

                ! Get the determinant bitstring once so we do not need to keep
                ! doing it. Then find lists of orbitals.
                cdet%f => walker_dets(:,idet)
                cdet%data => walker_data(:,idet)
                call decoder_ptr(sys, cdet%f, cdet)

                tmp_pop = walker_population(1,idet)

                ! Evaluate the projected energy.
                call update_proj_energy_ptr(sys, f0, cdet, real(walker_population(1,idet),p), &
                                            D0_population_cycle, proj_energy, D0_excit, hmatel)

                ! Loop over each walker on the determinant.
                do iparticle = 1, abs(walker_population(1,idet))

                    ! Spawn until the next annihilation barrier.
                    time = 0.0_p
                    do

                        ! We pass R to the timestep generator. Luckily all R_ij,
                        ! i/=j are the same for the hubbard model (U or
                        ! t - stored in matel),  and there are nexcitations of them.
                        R = abs(walker_data(1,idet) - shift(1)) + sum_off_diag
                        time = time + timestep(rng, R)

                        if ( time > t_barrier ) exit

                        call ct_spawn(rng, sys, cdet, walker_data(1,idet), int(walker_population(1,idet)), &
                                      R, nspawned, connection)

                        if (nspawned /= 0) then

                            ! If the spawned walker and the parent (all the
                            ! walkers on a particular det. have the same sgn due
                            ! to annihilation) are of opposite sgn we get death.
                            ! If death then kill the walker immediately and move
                            ! onto the next one.
                            if (connection%nexcit == 0 .and. &
                                       walker_population(1,idet)*nspawned < 0.0_p) then
                                tmp_pop = tmp_pop + nspawned
                                ! abs(nspawned) guaranteed to be 1
                                nparticles(1) = nparticles(1) - 1
                                ! The walker is dead---no need to continue spawning to barrier.
                                ndeath = ndeath + 1
                                exit
                            end if

                            ! If there were some walkers spawned, append them to the
                            ! spawned array - maintaining processor blocks if going in
                            ! parallel. We now also have an extra "time" array giving
                            ! the birth time of the walker.
                            call create_spawned_particle_ct(cdet, connection, nspawned, 1, time)

                        end if

                    end do

                end do

                ! update the walker population from the death events
                walker_population(1,idet) = tmp_pop

            end do

            ! Now we advance all the spawned walkers to the barrier from their
            ! respective birth times. Any walkers spawned as a consequence of
            ! this  must be appened to the spawned array and themselves advanced
            ! to the barrier.

            ! Start the first element in each block in qmc_spawn%sdata.
            current_pos = qmc_spawn%head_start + 1
            do
                do proc_id = 0, nprocs-1

                    if (current_pos(thread_id,proc_id) <= qmc_spawn%head(0,proc_id)) then
                        ! Have spawned walkers in the block to be sent to
                        ! processor proc_id.  Need to advance them to the barrier.

                        ! decode the spawned walker bitstring
                        cdet%f = int(qmc_spawn%sdata(:basis_length,current_pos(thread_id,proc_id)), i0)
                        K_ii = sc0_ptr(sys, cdet%f) - H00
                        call decoder_ptr(sys, cdet%f,cdet)

                        ! Spawn from this walker & append to the spawned array until
                        ! we hit the barrier
                        time = spawn_times(current_pos(thread_id,proc_id))
                        do

                            R = abs(K_ii - shift(1)) + sum_off_diag
                            time = time + timestep(rng, R)

                            if ( time > t_barrier ) exit

                            call ct_spawn(rng, sys, cdet, K_ii, int(qmc_spawn%sdata(spawned_pop,current_pos(thread_id,proc_id))), &
                                          R, nspawned, connection)

                            if (nspawned /= 0) then

                                ! Handle walker death
                                if(connection%nexcit == 0 .and. &
                                        qmc_spawn%sdata(spawned_pop,current_pos(thread_id,proc_id))*nspawned < 0) then
                                    qmc_spawn%sdata(spawned_pop,current_pos(thread_id,proc_id)) = &
                                            qmc_spawn%sdata(spawned_pop,current_pos(thread_id,proc_id)) + nspawned
                                    ndeath = ndeath + 1
                                    exit ! The walker is dead - do not continue
                                end if

                                ! Add a walker to the end of the spawned walker list in the
                                ! appropriate block - this will increment the appropriate
                                ! spawning heads for the processors which were spawned on
                                call create_spawned_particle_ct(cdet, connection, nspawned, spawned_pop, time)

                            end if

                        end do

                        ! go on to the next element
                        current_pos(thread_id,proc_id) = current_pos(thread_id,proc_id) + 1

                    end if

                end do

                ! Spawned all children and future generations to the barrier?
                if (all(current_pos == qmc_spawn%head+1)) exit

            end do

            call end_mc_cycle(ndeath, nattempts)

            call direct_annihilation(sys, rng, initiator_approximation)

            call end_report_loop(ireport, nparticles_old, t1, soft_exit)

            if (soft_exit) exit

        end do

        if (parent) then
            call write_fciqmc_final(ireport)
            write(6,'()')
        end if

        call load_balancing_report()

        if (soft_exit) then
            mc_cycles_done = mc_cycles_done + ncycles*ireport
        else
            mc_cycles_done = mc_cycles_done + ncycles*nreport
        end if

        if (dump_restart_file) then
            call dump_restart_hdf5(restart_info_global, mc_cycles_done, nparticles_old)
            write (6,'()')
        end if

        deallocate(current_pos, stat=ierr)
        call check_deallocate('current_pos', ierr)
        deallocate(connection_list, stat=ierr)
        call check_deallocate('connection_list', ierr)

    end subroutine do_ct_fciqmc


    subroutine ct_spawn(rng, sys, cdet, K_ii, parent_sgn, R, nspawned, connection)

        ! Randomly select a (valid) excitation

        ! In:
        !    sys: system being studied.
        !    cdet: info on current determinant, |D>, that we will spawn from.
        !    K_ii: the diagonal matrix element for the determinant |D>,
        !        < D | H - E_HF - S | D >.
        !    parent_sgn: sgn on the parent determinant (i.e. +ve or -ve integer)
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    nspawned: +/- 1 as @ the end of each time "jump" we only spawn
        !        1 walker.
        !    connection: the excitation connection between the parent and child
        !        determinants

        use excitations, only: excit
        use determinants, only: det_info
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t, hub_real, hub_k
        use hamiltonian_hub_real, only: slater_condon1_hub_real_excit
        use hamiltonian_hub_k, only: slater_condon2_hub_k_excit
        use excit_gen_hub_k, only: choose_ij_hub_k, find_ab_hub_k

        type(sys_t), intent(in) :: sys
        type(det_info), intent(in) :: cdet
        real(p), intent(in) :: K_ii, R
        integer, intent(in) :: parent_sgn
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: nspawned
        type(excit), intent(out) :: connection

        real(p) :: rand, K_ij
        logical :: allowed_excitation
        integer :: i, j, a, b, ij_sym

        rand = get_rand_close_open(rng)*R

        if (rand < abs(K_ii - shift(1))) then
            connection%nexcit = 0 ! spawn onto the same determinant (death/cloning)
            K_ij = K_ii - shift(1)
        else
            ! Generate a random excitation and reject if it's forbidden (i.e.
            ! the orbitals are already occupied).
            if (sys%system == hub_k) then
                ! Choose a random (i,j) pair to excite from.
                call choose_ij_hub_k(rng, sys, cdet%occ_list_alpha, cdet%occ_list_beta, i ,j, ij_sym)
                ! Choose a random (a,b) pair to attempt to excite to.
                ! The symmetry of (a,b) is set by the symmetry of (i,j) and
                ! hence b is uniquely determined by the choice of i,j and a.
                ! We choose a to be an unoccupied alpha spin-orbital and then
                ! reject the spawning attempt if b is in fact occupied.
                call find_ab_hub_k(rng, sys, cdet%f, cdet%unocc_list_alpha, ij_sym, a, b, allowed_excitation)
                if (allowed_excitation) then
                    connection%nexcit = 2
                    connection%from_orb(1:2) = (/ i,j /)
                    connection%to_orb(1:2) = (/ a,b /)
                    call slater_condon2_hub_k_excit(sys, cdet%f, connection, K_ij)
                else
                    K_ij = 0.0_p
                end if
            else if (sys%system == hub_real) then
                connection%nexcit = 1
                call slater_condon1_hub_real_excit(sys, cdet%f, connection, K_ij)
            end if

        end if

        if (K_ij == 0.0_p) then
            nspawned = 0
        else if (K_ij < 0.0_p) then    ! child is same sign as parent
            nspawned = sign(1,parent_sgn)
        else
            nspawned = -sign(1,parent_sgn)
        end if

    end subroutine ct_spawn

    subroutine create_spawned_particle_ct(cdet, connection, nspawn, particle_type, spawn_time)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the type of particle created.  Must correspond to
        !        the desired element in the spawning array (i.e. be spawned_pop
        !        for Hamiltonian particles and spawned_hf_pop for
        !        Hellmann--Feynman particles).
        !    spawn_time: The amount of imaginary time which has elapsed since
        !        the previous annihilation barrier.

        use hashing
        use parallel, only: iproc, nprocs

        use basis, only: basis_length
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det
        use fciqmc_data, only: qmc_spawn

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type
        real(p), intent(in) :: spawn_time

        integer(i0) :: f_new(basis_length)
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

#ifdef PARALLEL
        ! (Extra credit for parallel calculations)
        ! Need to determine which processor the spawned walker should be sent
        ! to.  This communication is done during the annihilation process, after
        ! all spawning and death has occured..
        iproc_spawn = modulo(murmurhash_bit_string(f_new, basis_length, qmc_spawn%hash_seed), nprocs)
#endif

        ! Move to the next position in the spawning array.
        qmc_spawn%head(0,iproc_spawn) = qmc_spawn%head(0,iproc_spawn) + 1

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        qmc_spawn%sdata(:,qmc_spawn%head(0,iproc_spawn)) = 0_int_s
        qmc_spawn%sdata(:basis_length,qmc_spawn%head(0,iproc_spawn)) = int(f_new, int_s)
        qmc_spawn%sdata(particle_type,qmc_spawn%head(0,iproc_spawn)) = int(nspawn, int_s)
        spawn_times(qmc_spawn%head(0,iproc_spawn)) = spawn_time

    end subroutine create_spawned_particle_ct

    function timestep(rng, R) result(dt)

        ! In:
        !    R: \sum_i < D | H - E_0 - S | D_i >, the sum of all non-zero matrix
        !       elements connected to the current determinant, D.
        ! In/Out:
        !    rng: random number generator.
        ! Returns:
        !    dt: the (stochastic) time which elapses before the next spawning
        !        event.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        real(p) :: dt
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(in)  :: R

        dt = -(1.0_p/R)*log(get_rand_close_open(rng))

    end function timestep


end module ct_fciqmc
