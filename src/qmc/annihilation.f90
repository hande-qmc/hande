module annihilation

use const

implicit none

contains

    subroutine direct_annihilation(sys, rng, reference, annihilation_flags, psip_list, spawn, &
                                   nspawn_events, determ)

        ! Annihilation algorithm.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        ! In:
        !    sys: system being studied.
        !    reference: current reference determinant.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    rng: random number generator.
        !    psip_list: particle_t object containing psip information after the
        !       death step for the current iteration; combined with the spawned
        !       particles on exit.
        !    spawn: spawn_t object containing the set of spawned particles.
        !    determ (optional): Derived type containing information on the
        !       semi-stochastic part of the simulation.
        ! Out:
        !    nspawn_events (optional): number of successful spawning events on
        !       the processor.

        use parallel, only: iproc
        use spawn_data, only: spawn_t, annihilate_wrapper_spawn_t, calc_events_spawn_t, memcheck_spawn_t
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: semi_stoch_t, particle_t, annihilation_flags_t, semi_stoch_separate_annihilation
        use reference_determinant, only: reference_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        type(reference_t), intent(in) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer, optional, intent(out) :: nspawn_events
        type(semi_stoch_t), intent(inout), optional :: determ

        logical :: doing_semi_stoch

        doing_semi_stoch = .false.
        if (present(determ)) doing_semi_stoch = determ%doing_semi_stoch
        if (present(nspawn_events)) nspawn_events = calc_events_spawn_t(spawn)

        call memcheck_spawn_t(spawn, dont_warn=spawn%warned)

        ! If performing a semi-stochastic calculation then the annihilation
        ! process is slightly different, so call the correct routines depending
        ! on the situation.
        if (doing_semi_stoch) then
            if (determ%projection_mode == semi_stoch_separate_annihilation) then
                call deterministic_annihilation(rng, psip_list, determ)
                call annihilate_wrapper_spawn_t(spawn, annihilation_flags%initiator_approx)
            else
                call annihilate_wrapper_spawn_t(spawn, annihilation_flags%initiator_approx, determ%sizes(iproc))
            end if

            call annihilate_main_list_wrapper(sys, rng, reference, annihilation_flags, psip_list, spawn, determ_flags=determ%flags)
        else
            call annihilate_wrapper_spawn_t(spawn, annihilation_flags%initiator_approx)
            call annihilate_main_list_wrapper(sys, rng, reference, annihilation_flags, psip_list, spawn)
        end if

    end subroutine direct_annihilation

    subroutine direct_annihilation_received_list(sys, rng, reference, annihilation_flags, psip_list, spawn_recv)

        ! Annihilation algorithm for non-blocking communications.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! If doing load balancing as well the annihilation procedure changes a
        ! bit. If load balancing has been decided upon then proc_map will have been
        ! updated by now. This means any walkers spawned during the current iteration
        ! will have been added to the modified section of the spawned walker array.
        ! To take care of annihilation we first merge the received list into the
        ! main list and then redistribute sections of the main list according
        ! the the updated proc map, as in normal load balancing.
        ! These walkers will be added to the spawned list and then communicated,
        ! using non-blocking comms instead of the normal MPI_AlltoAll. So they
        ! will need to be evolved upon receipt, as is the case for nomal
        ! non-blocking communications.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        ! In:
        !    sys: system being studied.
        !    reference: current reference determinant.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    rng: random number generator.
        !    psip_list: particle_t object containing psip information after the
        !       death step for the current iteration; combined with the spawned
        !       particles on exit.
        !    spawn_recv: spawn_t object containing spawned particles received
        !        from other processors.

        use spawn_data, only: annihilate_wrapper_non_blocking_spawn, spawn_t
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: particle_t, annihilation_flags_t
        use reference_determinant, only: reference_t

        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(dSFMT_t), intent(inout) :: rng
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn_recv

        ! Perform annihilation inside received list. This involves annihilating
        ! walkers which were spawned onto this processor from other processors
        ! (not including the current processor) from  the previous iteration.
        ! They have since been evolved so they can be annihilated with the main list.
        ! First annihilate within spawn_recv.
        call annihilate_wrapper_non_blocking_spawn(spawn_recv, annihilation_flags%initiator_approx)
        ! Annihilate with main list.
        call annihilate_main_list_wrapper(sys, rng, reference, annihilation_flags, psip_list, spawn_recv)

    end subroutine direct_annihilation_received_list

    subroutine direct_annihilation_spawned_list(sys, rng, reference, annihilation_flags, psip_list, spawn, send_counts, &
                                                req_data_s, non_block_spawn, nspawn_events)

        ! Annihilation algorithm for non-blocking communications.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        ! In:
        !    sys: system being studied.
        !    reference: current reference determinant.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    rng: random number generator.
        !    psip_list: particle_t object containing psip information after the
        !       death step for the current iteration; combined with the spawned
        !       particles on exit.
        !    spawn: spawn_t object containing spawned particles.
        !    send_counts: array of messages sizes. Will be allocated in
        !       calculate_displacements and sent in non_blocking_send.
        !    req_data_s: array of requests for non-blocking send of walkers.
        ! Out:
        !    non_block_spawn: number of spawned particles on current processor
        !       during current MC cycle.
        ! Out (optional):
        !    nspawn_events (optional): number of successful spawning events on
        !       the processor.

        use parallel, only: nthreads, iproc
        use spawn_data, only: annihilate_wrapper_non_blocking_spawn, calculate_displacements, &
                              non_blocking_send, memcheck_spawn_t
        use sort, only: qsort
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: particle_t, annihilation_flags_t
        use reference_determinant, only: reference_t
        use spawn_data, only: spawn_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        type(reference_t), intent(in) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout) :: send_counts(0:)
        integer, intent(inout) :: req_data_s(0:)
        integer, intent(out) :: non_block_spawn(:)
        integer, optional, intent(out) :: nspawn_events

        integer, parameter :: thread_id = 0

        ! Need to calculate how many walkers we are going to send to all other
        ! processors. Need to do it now as spawn%head changes meaning upon annihilation.
        call calculate_displacements(spawn, send_counts, non_block_spawn)
        if (present(nspawn_events)) nspawn_events = non_block_spawn(1)

        call memcheck_spawn_t(spawn, dont_warn=spawn%warned)

        ! Perform annihilation within the spawned walker list.
        ! This involves locating, compressing and sorting the section of the spawned
        ! list which needs to be annihilated with the main list on this processor.
        call annihilate_wrapper_non_blocking_spawn(spawn, annihilation_flags%initiator_approx, iproc)
        ! Annihilate portion of spawned list with main list.
        call annihilate_main_list_wrapper(sys, rng, reference, annihilation_flags, psip_list, spawn, &
                                          spawn%head_start(thread_id, iproc)+nthreads)
        ! Communicate walkers spawned onto other processors during this
        ! evolution step to their new processors.
        call non_blocking_send(spawn, send_counts, req_data_s)

    end subroutine direct_annihilation_spawned_list

    subroutine annihilate_main_list_wrapper(sys, rng, reference, annihilation_flags, psip_list, spawn, &
                                            lower_bound, determ_flags)

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process during non-blocking
        ! communications.

        ! In:
        !    sys: system being studied.
        !    reference: current reference determinant.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    rng: random number generator.
        !    psip_list: particle_t object containing psip information after the
        !       death step for the current iteration; combined with the spawned
        !       particles on exit.
        !    spawn: spawn_t object containing spawned particles. For non-blocking
        !       communications a subsection of the spawned walker list will be annihilated
        !       with the main list, otherwise the entire list will be annihilated and merged.
        !    determ_flags (optional): A list of flags specifying whether determinants in
        !        the corresponding particle_t%states are deterministic or not.
        ! In (optional):
        !     lower_bound: starting point we annihiliate from in spawn_t object.

        use system, only: sys_t
        use spawn_data, only: spawn_t
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: particle_t, annihilation_flags_t
        use reference_determinant, only: reference_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        type(reference_t), intent(in) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        integer, optional, intent(in) :: lower_bound
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout), optional :: determ_flags(:)

        integer, parameter :: thread_id = 0
        integer :: spawn_start

        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if

        if (spawn%head(thread_id,0) >= spawn_start) then
            ! Have spawned walkers on this processor.

            if (annihilation_flags%initiator_approx) then
                call annihilate_main_list_initiator(psip_list, spawn, sys%basis%tensor_label_len, lower_bound)
            else
                call annihilate_main_list(psip_list, spawn, sys%basis%tensor_label_len, lower_bound)
            end if

            ! Remove determinants with zero walkers on them from the main
            ! walker list.
            call remove_unoccupied_dets(rng, psip_list, annihilation_flags%real_amplitudes, determ_flags)

            ! Remove low-population spawned walkers by stochastically
            ! rounding their population up to one or down to zero.
            if (annihilation_flags%real_amplitudes) then
                call round_low_population_spawns(rng, psip_list%pop_real_factor, spawn, lower_bound)
            end if

            ! Insert new walkers into main walker list.
            call insert_new_walkers(sys, psip_list, reference, annihilation_flags, spawn, determ_flags, lower_bound)

        else

            ! No spawned walkers so we only have to check to see if death has
            ! killed the entire population on a determinant.
            call remove_unoccupied_dets(rng, psip_list, annihilation_flags%real_amplitudes, determ_flags)

        end if

    end subroutine annihilate_main_list_wrapper

    subroutine annihilate_main_list(psip_list, spawn, tensor_label_len, lower_bound)

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        ! In:
        !    tensor_label_len: number of elements in the bit array describing the position
        !       of the particle in the space (i.e.  determinant label in vector/pair of
        !       determinants label in array).
        ! In/Out:
        !    psip_list: particle_t object containing psip information.
        !       On exit the particles spawned onto occupied sites have been
        !       annihilated with the existing populations.
        !    spawn: spawn_t obeject containing spawned particles to be annihilated with main
        !       list.  On exit contains particles spawned onto non-occupied
        !       sites.
        ! In (optional):
        !    lower_bound: starting point we annihiliate from in spawn_t object.
        !       Default: 1.

        use search, only: binary_search
        use spawn_data, only: spawn_t
        use qmc_data, only: particle_t

        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer, intent(in) :: tensor_label_len
        integer, intent(in), optional :: lower_bound

        integer :: i, pos, k, istart, iend, nannihilate, spawn_start
        integer(int_p) :: old_pop(psip_list%nspaces)
        integer(i0) :: f(tensor_label_len)

        logical :: hit
        integer, parameter :: thread_id = 0

        nannihilate = 0
        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if
        istart = 1
        iend = psip_list%nstates

        do i = spawn_start, spawn%head(thread_id,0)
            f = int(spawn%sdata(:tensor_label_len,i), i0)
            call binary_search(psip_list%states, f, istart, iend, hit, pos)
            if (hit) then
                ! Annihilate!
                old_pop = psip_list%pops(:,pos)
                psip_list%pops(:,pos) = psip_list%pops(:,pos) + &
                    int(spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes,i), int_p)
                nannihilate = nannihilate + 1
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                associate(nparticles=>psip_list%nparticles)
                    nparticles = nparticles + real(abs(psip_list%pops(:,pos)) - abs(old_pop),p)/psip_list%pop_real_factor
                end associate
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                k = i - nannihilate
                spawn%sdata(:,k) = spawn%sdata(:,i)
            end if
        end do

        spawn%head(thread_id,0) = spawn%head(thread_id,0) - nannihilate

    end subroutine annihilate_main_list

    subroutine annihilate_main_list_initiator(psip_list, spawn, tensor_label_len, lower_bound)

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list starting from lower bound in spawn.

        ! This version is for the initiator algorithm, where we also need to
        ! discard spawned walkers which are on previously unoccupied determinants
        ! and which are from non-initiator or non-sign-coherent events.

        ! In:
        !    tensor_label_len: number of elements in the bit array describing the position
        !       of the particle in the space (i.e.  determinant label in vector/pair of
        !       determinants label in array).
        ! In/Out:
        !    psip_list: particle_t object containing psip information.
        !       On exit the particles spawned onto occupied sites have been
        !       annihilated with the existing populations.
        !    spawn: spawn_t object we wish to annihilate with main list.
        !       On exit contains particles spawned onto non-occupied sites.
        ! In (Optional):
        !    lower_bound: starting point we annihiliate from in spawn_t object.

        use search, only: binary_search
        use qmc_data, only: particle_t
        use spawn_data, only: spawn_t

        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer, intent(in) :: tensor_label_len
        integer, intent(in), optional :: lower_bound

        integer :: i, ipart, pos, istart, iend, nannihilate, spawn_start
        integer(int_p) :: old_pop(psip_list%nspaces)
        integer(i0) :: f(tensor_label_len)
        logical :: hit, discard
        integer, parameter :: thread_id = 0

        nannihilate = 0
        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if
        istart = 1
        iend = psip_list%nstates
        do i = spawn_start, spawn%head(thread_id,0)
            f = int(spawn%sdata(:tensor_label_len,i), i0)
            call binary_search(psip_list%states, f, istart, iend, hit, pos)
            if (hit) then
                old_pop = psip_list%pops(:,pos)
                ! Need to take into account that the determinant might not have
                ! a non-zero population for all particle types.
                do ipart = 1, psip_list%nspaces
                    if (psip_list%pops(ipart,pos) /= 0_int_p) then
                        ! Annihilate!
                        psip_list%pops(ipart,pos) = psip_list%pops(ipart,pos) + &
                                                        int(spawn%sdata(ipart+spawn%bit_str_len,i), int_p)
                    else if (.not.btest(spawn%sdata(spawn%flag_indx,i),ipart-1)) then
                        ! Keep only if from a multiple spawning event or an
                        ! initiator.
                        ! If this is the case, then sdata(flag_indx,i)
                        ! does not have a bit set in corresponding to 2**(ipart-1)
                        ! where ipart+bit_str_len is the index of this walker type
                        ! in the spawn%sdata array.
                        psip_list%pops(ipart,pos) = int(spawn%sdata(ipart+spawn%bit_str_len,i), int_p)
                    end if
                end do
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                associate(nparticles=>psip_list%nparticles)
                    nparticles = nparticles + real(abs(psip_list%pops(:,pos)) - abs(old_pop),p)/psip_list%pop_real_factor
                end associate
                ! One more entry to be removed from the spawn%sdata array.
                nannihilate = nannihilate + 1
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                ! Keep only progeny spawned by initiator determinants
                ! or multiple sign-coherent events.  If neither of these
                ! conditions are met then the (j-1)-th bit of sdata(flag_indx,i) is set,
                ! where j is the particle index in spawn%sdata(bit_str_len+1:,:).
                discard = .true.
                do ipart = 1, spawn%ntypes
                    if (btest(spawn%sdata(spawn%flag_indx,i),ipart-1)) then
                        ! discard attempting spawnings from non-initiator walkers
                        ! onto unoccupied determinants.
                        ! note that the number of particles (psip_list%nparticles) was not
                        ! updated at the time of spawning, so doesn't change.
                        spawn%sdata(spawn%bit_str_len+ipart,i-nannihilate) = 0_int_s
                    else
                        ! keep!
                        spawn%sdata(spawn%bit_str_len+ipart,i-nannihilate) = spawn%sdata(spawn%bit_str_len+ipart,i)
                        discard = .false.
                    end if
                end do
                if (discard) then
                    ! Don't need to keep any particles from the current slot so can
                    ! just overwrite them...
                    nannihilate = nannihilate + 1
                else
                    ! Need to copy the bit string across...
                    spawn%sdata(:tensor_label_len,i-nannihilate) = spawn%sdata(:tensor_label_len,i)
                end if
            end if
        end do

        spawn%head(thread_id,0) = spawn%head(thread_id,0) - nannihilate

    end subroutine annihilate_main_list_initiator

    subroutine deterministic_annihilation(rng, psip_list, determ)

        ! Add in the deterministic spawnings to the main list.

        ! In:
        !    determ: Derived type containing information on the semi-stochastic
        !       part of the simulation.
        !    psip_list: particle_t object containing psip information.
        !       On exit the particles spawned onto occupied sites by
        !       the determinstic projection have been combined with the existing
        !       population.
        ! In/Out:
        !    rng: random number generator.

        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use qmc_data, only: semi_stoch_t, particle_t

        type(particle_t), intent(inout) :: psip_list
        type(dSFMT_t), intent(inout) :: rng
        type(semi_stoch_t), intent(in), optional :: determ

        integer :: i, ind
        integer(int_p) :: nspawn, spawn_sign, old_pop(psip_list%nspaces)
        real(p) :: scaled_amp

        do i = 1, size(determ%vector)
            ind = determ%indices(i)

            scaled_amp = determ%vector(i)*psip_list%pop_real_factor
            spawn_sign = 1
            if (scaled_amp < 0.0_p) spawn_sign = -1
            ! Stochastically round the scaled amplitude to the nearest integer
            ! in order to encode it.
            scaled_amp = abs(scaled_amp)
            nspawn = int(scaled_amp, int_p)
            scaled_amp = scaled_amp - nspawn
            if (scaled_amp > get_rand_close_open(rng)) nspawn = nspawn + 1_int_p

            ! Add in the now-encoded deterministic spawning amplitude.
            old_pop = psip_list%pops(:,ind)
            psip_list%pops(1,ind) = psip_list%pops(1,ind) + spawn_sign*nspawn
            psip_list%nparticles = psip_list%nparticles + real(abs(psip_list%pops(:,ind))-abs(old_pop),p)/psip_list%pop_real_factor
        end do

    end subroutine deterministic_annihilation

    subroutine remove_unoccupied_dets(rng, psip_list, real_amplitudes, determ_flags)

        ! Remove any determinants with 0 population.
        ! This can be done in a more efficient manner by doing it only when
        ! necessary...
        ! Also, before doing this, stochastically round up or down any
        ! populations which are less than one.
        ! The above steps are not performed for deterministic states which are
        ! kept in psip_list%states whatever their sign.

        ! In:
        !     real_amplitudes: true if using real particle amplitudes.
        ! In/Out:
        !     rng: random number generator.
        !     psip_list: psip information.  On exit unoccupied states in the
        !         stochastic space have been removed.
        !     determ_flags: A list of flags specifying whether determinants in
        !         psip_list%states are deterministic or not.

        use dSFMT_interface, only: dSFMT_t
        use stoch_utils, only: stochastic_round
        use qmc_data, only: particle_t

        type(dSFMT_t), intent(inout) :: rng
        type(particle_t), intent(inout) :: psip_list
        logical, intent(in) :: real_amplitudes
        integer, intent(inout), optional :: determ_flags(:)

        integer :: nzero, i, k
        integer(int_p) :: old_pop(psip_list%nspaces)
        logical :: determ_det

        nzero = 0
        do i = 1, psip_list%nstates

            determ_det = .false.
            if (present(determ_flags)) determ_det = determ_flags(i) == 0

            ! Stochastically round the walker populations up to real_factor
            ! (which is equal to 1 in the decoded representation) or down to
            ! zero. This is not done for deterministic states.
            if (real_amplitudes .and. (.not. determ_det)) then
                old_pop = psip_list%pops(:,i)
                call stochastic_round(rng, psip_list%pops(:,i), psip_list%pop_real_factor, psip_list%nspaces)
                associate(nparticles=>psip_list%nparticles)
                    nparticles = nparticles + real(abs(psip_list%pops(:,i)) - abs(old_pop),p)/psip_list%pop_real_factor
                end associate
            end if

            if (all(psip_list%pops(:,i) == 0_int_p) .and. (.not. determ_det)) then
                nzero = nzero + 1
            else if (nzero > 0) then
                k = i - nzero
                psip_list%states(:,k) = psip_list%states(:,i)
                psip_list%pops(:,k) = psip_list%pops(:,i)
                psip_list%dat(:,k) = psip_list%dat(:,i)
                if (present(determ_flags)) determ_flags(k) = determ_flags(i)
            end if
        end do
        psip_list%nstates = psip_list%nstates - nzero

    end subroutine remove_unoccupied_dets

    subroutine round_low_population_spawns(rng, pop_real_factor, spawn, lower_bound)

        ! Loop over all spawned walkers. For each walker with a population of
        ! less than one, round it up to one with a probability equal to its
        ! population, otherwise down to zero. This will ensure that the
        ! expectation value of the amplitudes on all determinants are the same
        ! as before, but will prevent many low-weight walkers remaining in the
        ! simulation.

        ! Note that all deterministic states are always kept in the main list.
        ! The walkers in the spawned list at this point only contain states
        ! not in the main list, so cannot contain any deterministic states.
        ! Therefore, as an optimisation, we don't need to check if determinants
        ! are deterministic or not before any rounding.

        ! In/Out:
        !    rng: random number generator.
        !    pop_real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        !    spawn: spawn_t object containing the set of spawned particles.
        ! In (optional):
        !    lower_bound: starting point we annihiliate from in spawn_t object.
        !       Default: 1.

        use dSFMT_interface, only: dSFMT_t
        use spawn_data, only: spawn_t
        use stoch_utils, only: stochastic_round

        type(dSFMT_t), intent(inout) :: rng
        integer(int_p), intent(in) :: pop_real_factor
        type(spawn_t), intent(inout) :: spawn
        integer, optional, intent(in) :: lower_bound

        integer :: i, k, nremoved, spawn_start
        integer(int_s) :: real_factor_s
        integer, parameter :: thread_id = 0

        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if

        real_factor_s = int(pop_real_factor, int_s)

        nremoved = 0
        ! [note] - It might be more efficient to combine this with insert_new_walkers.
        ! [note] - The number of particles to insert should be small by this point though...
        do i = spawn_start, spawn%head(thread_id,0)

            ! spawned_population holds the spawned population in its encoded
            ! form (see comments for psip_list%pops).
            associate(spawned_population => spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes, i))

                ! Stochastically round the walker populations up or down to
                ! real_factor (which is equal to 1 in the decoded representation).
                call stochastic_round(rng, spawned_population, real_factor_s, spawn%ntypes)

                ! If all the amplitudes for this determinant were zeroed then we
                ! don't want to add it to the main list.
                if (all(spawned_population == 0_int_s)) then
                    nremoved = nremoved + 1
                else
                    ! Shuffle this determinant down to fill in any newly opened
                    ! slots.
                    k = i - nremoved
                    spawn%sdata(:,k) = spawn%sdata(:,i)
                end if

            end associate

        end do

        spawn%head(thread_id,0) = spawn%head(thread_id,0) - nremoved

    end subroutine round_low_population_spawns

    subroutine insert_new_walkers(sys, psip_list, ref, annihilation_flags, spawn, determ_flags, lower_bound)

        ! Insert new walkers into the main walker list from the spawned list.
        ! This is done after all particles have been annihilated, so the spawned
        ! list contains only new walkers.

        ! In:
        !    sys: system being studied.
        !    ref: reference determinant --- the diagonal matrix elements are required.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    psip_list: psip information.
        !    spawn: spawn_t object containing list of particles spawned onto
        !        previously occupied sites.
        ! In (optional):
        !    determ_flags: A list of flags specifying whether determinants in
        !        psip_list%states are deterministic or not.
        !    lower_bound: starting point we annihiliate from in spawn_t object.

        use errors, only: stop_all
        use parallel, only: iproc
        use qmc_data, only: particle_t, annihilation_flags_t
        use reference_determinant, only: reference_t
        use search, only: binary_search
        use spawn_data, only: spawn_t
        use system, only: sys_t
        use utils, only: int_fmt

        use, intrinsic :: iso_fortran_env, only: error_unit

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(inout) :: psip_list
        type(reference_t), intent(in) :: ref
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout), optional :: determ_flags(:)
        integer, intent(in), optional :: lower_bound

        integer :: i, istart, iend, j, k, pos, spawn_start, disp
        real(p) :: real_population(psip_list%nspaces)

        logical :: hit
        integer, parameter :: thread_id = 0
        real :: fill_fraction

        ! Merge new walkers into the main list.

        ! Both the main list and the spawned list are sorted: the spawned list
        ! is sorted explicitly and the main list is sorted by construction via
        ! merging.

        ! 1. Find the position where the spawned walker should go.
        ! 2. Move all walkers above it to create the vacant slot for the new
        ! walker.  As we know how many elements we are inserting, we only need
        ! move a given walker at most once.
        ! 3. Insert the new walker at the bottom of the shifted block so it
        ! doesn't have to be moved again to accommodate other new walkers.

        ! We can make the search faster as we iterate through the spawned
        ! walkers in descending order, so once we know where one walker goes, we
        ! know that the next new walker has to go below it, allowing us to
        ! search through an ever-decreasing number of elements.

        if (present(lower_bound)) then
            spawn_start = lower_bound
            disp = lower_bound - 1
        else
            spawn_start = 1
            disp = 0
        end if

        ! Don't bother to perform these checks and print more error messages
        ! if we've run out of memory already.
        if (.not. psip_list%error) then
            fill_fraction = real(psip_list%nstates+(spawn%head(thread_id,0)-spawn_start+1))/size(psip_list%states,2)
            if (fill_fraction > 1.00) then
                write (error_unit,'(1X,"# Error: No space left in main particle array on processor",'//int_fmt(iproc,1)//',".")') &
                              iproc
                write (error_unit,'(1X,"# Error: HANDE will exit at the end of this report loop.")')
                write (error_unit,'(1X,"# Error: Note that spawning until the end of the report loop will be affected and&
                              & so results from this final loop may be slightly incorrect.")')
                write (error_unit,'(1X,"# Error: Some reconvergence time should be allowed if continuing from a&
                              & subsequent restart file.")')

                psip_list%error = .true.
            else if (fill_fraction > 0.95) then
                if (psip_list%warn) then
                    write (error_unit,'(1X,"# Warning: filled over 95% of main particle array on processor",'//int_fmt(iproc,1)&
                              //',".")') iproc
                    write (error_unit,'(1x,"This warning only prints once")')
                    psip_list%warn = .false.
                end if
                psip_list%warning_count = psip_list%warning_count + 1
            end if
        end if

        if (.not. psip_list%error) then
            istart = 1
            iend = psip_list%nstates
            do i = spawn%head(thread_id,0), spawn_start, -1

                ! spawned det is not in the main walker list.
                call binary_search(psip_list%states, int(spawn%sdata(:sys%basis%tensor_label_len,i), i0), &
                                   istart, iend, hit, pos)
                ! f should be in slot pos.  Move all determinants above it.
                do j = iend, pos, -1
                    ! i is the number of determinants that will be inserted below j.
                    k = j + i - disp
                    psip_list%states(:,k) = psip_list%states(:,j)
                    psip_list%pops(:,k) = psip_list%pops(:,j)
                    psip_list%dat(:,k) = psip_list%dat(:,j)
                    if (present(determ_flags)) determ_flags(k) = determ_flags(j)
                end do

                ! Insert new walker into pos and shift it to accommodate the number
                ! of elements that are still to be inserted below it.
                k = pos + i - 1 - disp

                ! The encoded spawned walker sign.
                associate(spawned_population => spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes, i), &
                        tbl=>sys%basis%tensor_label_len)
                    call insert_new_walker(sys, psip_list, annihilation_flags, k, int(spawn%sdata(:tbl,i), i0), &
                                           int(spawned_population, int_p), ref)
                    ! Extract the real sign from the encoded sign.
                    real_population = real(spawned_population,p)/psip_list%pop_real_factor
                    psip_list%nparticles = psip_list%nparticles + abs(real_population)
                end associate

                ! A deterministic state can never leave the main list so cannot be
                ! in the spawned list at this point. So set the flag to specify
                ! that this state is not deterministic.
                if (present(determ_flags)) determ_flags(k) = 1

                ! Next walker will be inserted below this one.
                iend = pos - 1
            end do

            ! Update psip_list%nstates.
            psip_list%nstates = psip_list%nstates + spawn%head(thread_id,0) - disp
        end if

    end subroutine insert_new_walkers

    subroutine insert_new_walker(sys, psip_list, annihilation_flags,  pos, det, population, ref)

        ! Insert a new determinant, det, at position pos in psip_list%states. Also
        ! insert a new population at position pos in psip_list%pops and
        ! calculate and insert all new values of psip_list%dat for the given
        ! walker. Note that all data at position pos in these arrays will be
        ! overwritten.

        ! In:
        !    sys: system being studied.
        !    psip_list: psip information.
        !    annihilation_flags: calculation specific annihilation flags.
        !    pos: The position in the walker arrays in which to insert the new
        !        data.
        !    det: The determinant to insert into psip_list%states.
        !    population: The population to insert into psip_list%pops.
        !    ref: reference determinant.

        use calc, only: doing_calc, hfs_fciqmc_calc, dmqmc_calc
        use heisenberg_estimators, only: neel_singlet_data
        use proc_pointers, only: sc0_ptr, op0_ptr, trial_dm_ptr
        use system, only: sys_t
        use qmc_data, only: particle_t, annihilation_flags_t, neel_singlet
        use reference_determinant, only: reference_t

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(inout) :: psip_list
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        integer, intent(in) :: pos
        integer(i0), intent(in) :: det(sys%basis%tensor_label_len)
        integer(int_p), intent(in) :: population(psip_list%nspaces)
        type(reference_t), intent(in) :: ref

        ! Insert the new determinant.
        psip_list%states(:,pos) = det
        ! Insert the new population.
        psip_list%pops(:,pos) = population
        ! Calculate and insert all new components of psip_list%dat.
        psip_list%dat(1,pos) = sc0_ptr(sys, det) - ref%H00
        associate(pl=>psip_list)
            if (annihilation_flags%trial_function == neel_singlet) &
                pl%dat(pl%nspaces+1:pl%nspaces+2,pos) = neel_singlet_data(sys, det)
        end associate
        if (doing_calc(hfs_fciqmc_calc)) then
            ! Set psip_list%dat(2:,k) = <D_i|O|D_i> - <D_0|O|D_0>.
            psip_list%dat(2,pos) = op0_ptr(sys, det) - ref%O00
        else if (doing_calc(dmqmc_calc)) then
            if (annihilation_flags%propagate_to_beta .and. .not. annihilation_flags%symmetric) then
                ! Store H^T_ii-H_jj so we can propagate with ~ 1 + \Delta\beta(H^T_ii - H_jj),
                ! where H^T_ii is the "trial" Hamiltonian.
                associate(bl=>sys%basis%string_len, pl=>psip_list)
                    pl%dat(1,pos) = -trial_dm_ptr(sys, pl%states((bl+1):(2*bl),pos)) + (pl%dat(1,pos) + ref%H00)
                end associate
            else if (annihilation_flags%propagate_to_beta .and. annihilation_flags%symmetric) then
                associate(bl=>sys%basis%string_len, pl=>psip_list)
                    pl%dat(1,pos) = -0.5_p*((trial_dm_ptr(sys, det)+trial_dm_ptr(sys, pl%states((bl+1):(2*bl),pos))) - &
                                     (pl%dat(1,pos) + ref%H00 + sc0_ptr(sys, pl%states((bl+1):(2*bl),pos))))
                end associate
            else
                ! Set the energy to be the average of the two induvidual energies.
                associate(bl=>sys%basis%string_len, pl=>psip_list)
                    pl%dat(1,pos) = (pl%dat(1,pos) + sc0_ptr(sys, pl%states((bl+1):(2*bl),pos)) - ref%H00)/2
                end associate
                if (annihilation_flags%replica_tricks) then
                    psip_list%dat(2:psip_list%nspaces,pos) = psip_list%dat(1,pos)
                end if
            end if
        end if

    end subroutine insert_new_walker

end module annihilation
