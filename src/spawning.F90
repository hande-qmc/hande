module spawning

! Module for procedures involved in the spawning step of the FCIQMC algorithm.

! We wish to spawn with probability
!   tau |H_ij|,
! where tau is the timestep of the simulation.
! The probability of doing something is the probability of selecting to attempt
! to do it multiplied by the probability of actually doing it, hence:
!   p_spawn = p_select tau*|H_ij|/p_gen
! and p_select = p_gen for normalised probabilities.
! p_select is included intrinsically in the algorithm, as we choose a random
! determinant, j, that is connected to the i-th determinant and attempt to spawn
! from a particle on i-th determinant onto the j-th determinant.
! Thus we compare to the probability tau*|H_ij|/p_gen in order to determine
! whether the spawning attempt is successful.

! This is just for top-level spawning routines and utility functions.  The
! actual work of generating a random excitation is done in the system-specific
! excit_gen_* modules.

! TODO: profile to discover how much time is spent obtaining the correct sign of
! a matrix element (i.e. find_permutation_* routines).  These can be avoided
! unless the spawning event is successful as we only need the correct sign of
! the matrix element if offspring are produced.

use const
implicit none

contains

!--- Spawning wrappers ---

    subroutine spawn(rng, sys, cdet, parent_sign, gen_excit_ptr, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant.

        ! This is just a thin wrapper around a system-specific excitation
        ! generator and a utility function.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !        gen_excit_ptr%full *must* be set to a procedure which generates
        !        a complete excitation.
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use excitations, only: excit
        use system, only: sys_t
        use proc_pointers, only: gen_excit_ptr_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel

        ! 1. Generate random excitation.
        call gen_excit_ptr%full(rng, sys, cdet, pgen, connection, hmatel)

        ! 2. Attempt spawning.
        nspawn = attempt_to_spawn(rng, hmatel, pgen, parent_sign)

    end subroutine spawn

    subroutine spawn_importance_sampling(rng, sys, cdet, parent_sign, gen_excit_ptr, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant.

        ! This subroutine applies a transformation to the Hamiltonian to achieve
        ! importance sampling of the stochastic process.

        ! This is just a thin wrapper around a system-specific excitation
        ! generator, trial-function specific transformation routine and
        ! a utility function.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !        gen_excit_ptr%full *must* be set to a procedure which generates
        !        a complete excitation.
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use system, only: sys_t
        use excitations, only: excit
        use proc_pointers, only: gen_excit_ptr_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel

        ! 1. Generate random excitation.
        call gen_excit_ptr%full(rng, sys, cdet, pgen, connection, hmatel)

        ! 2. Transform Hamiltonian matrix element by trial function.
        call gen_excit_ptr%trial_fn(sys, cdet, connection, hmatel)

        ! 3. Attempt spawning.
        nspawn = attempt_to_spawn(rng, hmatel, pgen, parent_sign)

    end subroutine spawn_importance_sampling

    subroutine spawn_lattice_split_gen(rng, sys, cdet, parent_sign, gen_excit_ptr, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant.

        ! This is just a thin wrapper around a system-specific excitation
        ! generator and a utility function.

        ! This subroutine applies a transformation to the Hamiltonian to achieve
        ! importance sampling of the stochastic process.

        ! For lattice models (principally the Hubbard model; the savings for the
        ! Heisenberg model are minimal) one can make a useful optimisation, as
        ! determining |H_ij| is fast (indeed, constant if non-zero) and does not
        ! require full knowledge of the excitation.  As only a small fraction of
        ! spawning events are successful, it is faster to not do any unnecessary
        ! work and test whether the spawning event is successful before
        ! finalising the excitation.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !        gen_excit_ptr%init and gen_excit_ptr%finalise *must* be set to
        !        a pair of procedures which generate a complete excitation.
        !        gen_excit_ptr%init must return (at least) the connecting matrix
        !        element and gen_excit_ptr%finalise must fill in the rest of the
        !        information about the excitation.
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use system, only: sys_t
        use excitations, only: excit
        use fciqmc_data, only: tau
        use proc_pointers, only: gen_excit_ptr_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, abs_hmatel, hmatel

        ! 1. Generate enough of a random excitation to determinant the
        ! generation probability and |H_ij|.
        call gen_excit_ptr%init(rng, sys, cdet, pgen, connection, abs_hmatel)

        ! 2. Attempt spawning.
        nspawn = nspawn_from_prob(rng, tau*abs_hmatel/pgen)

        if (nspawn /= 0) then

            ! 3. Complete excitation and find sign of connecting matrix element.
            call gen_excit_ptr%finalise(rng, sys, cdet, connection, hmatel)

            ! 4. Find sign of offspring.
            call set_child_sign(hmatel, parent_sign, nspawn)

        end if

    end subroutine spawn_lattice_split_gen

    subroutine spawn_lattice_split_gen_importance_sampling(rng, sys, cdet, parent_sign, gen_excit_ptr, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant.

        ! This is just a thin wrapper around a system-specific excitation
        ! generator and a utility function.

        ! For lattice models (principally the Hubbard model; the savings for the
        ! Heisenberg model are minimal) one can make a useful optimisation, as
        ! determining |H_ij| is fast (indeed, constant if non-zero) and does not
        ! require full knowledge of the excitation.  As only a small fraction of
        ! spawning events are successful, it is faster to not do any unnecessary
        ! work and test whether the spawning event is successful before
        ! finalising the excitation.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !        gen_excit_ptr%init and gen_excit_ptr%finalise *must* be set to
        !        a pair of procedures which generate a complete excitation.
        !        gen_excit_ptr%init must return (at least) the connecting matrix
        !        element and gen_excit_ptr%finalise must fill in the rest of the
        !        information about the excitation.
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use system, only: sys_t
        use excitations, only: excit
        use fciqmc_data, only: tau
        use proc_pointers, only: gen_excit_ptr_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, tilde_hmatel, hmatel

        ! 1. Generate enough of a random excitation to determinant the
        ! generation probability and |H_ij|.
        call gen_excit_ptr%init(rng, sys, cdet, pgen, connection, tilde_hmatel)

        ! 2. Transform Hamiltonian matrix element by trial function.
        call gen_excit_ptr%trial_fn(sys, cdet, connection, tilde_hmatel)

        ! 3. Attempt spawning.
        nspawn = nspawn_from_prob(rng, tau*abs(tilde_hmatel)/pgen)

        if (nspawn /= 0) then

            ! 4. Complete excitation and find sign of connecting matrix element.
            ! *NOTE*: this returns the original matrix element and *not* the
            ! matrix element after the trial function transformation.
            call gen_excit_ptr%finalise(rng, sys, cdet, connection, hmatel)

            ! 5. Find sign of offspring.
            ! Note that we don't care about the value of H_ij at this step, only
            ! the sign.
            call set_child_sign(tilde_hmatel*hmatel, parent_sign, nspawn)

        end if

    end subroutine spawn_lattice_split_gen_importance_sampling

    subroutine spawn_null(rng, sys, cdet, parent_sign, gen_excit_ptr, nspawn, connection)

        ! This is a null spawning routine for use with operators which are
        ! diagonal in the basis and hence only have a cloning step in the
        ! Hellmann-Feynman sampling.  It does *nothing*.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from (or not, in this case!).
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        !    gen_excit_ptr: procedure pointer to excitation generators.
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use system, only: sys_t
        use excitations, only: excit
        use proc_pointers, only: gen_excit_ptr_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        ! Just some null operations to avoid -Wall -Werror causing errors.
        connection%nexcit = huge(0)

        ! Return nspawn = 0 as we don't want to do any spawning.
        nspawn = 0

    end subroutine spawn_null

!--- Attempt spawning based upon random excitation ---

    function nspawn_from_prob(rng, probability) result(number_spawned)

        ! Generate the number spawned from a probability. If probability is greater than
        ! zero, then number spawned = int(probability) + stochastic{0,1}
        ! where the latter half of the RHS is a stochastic spawning from the remainder
        !
        ! In:
        !    probability: the spawning probability
        ! In/Out:
        !    rng: random number generator.
        !
        ! Returns:
        !    number_spawned: the number spawned from this probability

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        implicit none
        real(p), intent(in) :: probability
        type(dSFMT_t), intent(inout) :: rng
        integer             :: number_spawned
        real(p)             :: psuccess, pstochastic

        ! Generate random number
        psuccess = get_rand_close_open(rng)

        ! Multiple offspring
        number_spawned = int(probability)

        ! Stochastic offspring
        pstochastic = probability - number_spawned
        if (pstochastic > psuccess) number_spawned = number_spawned + 1

    end function nspawn_from_prob

    subroutine set_child_sign(hmatel, parent_sign, nspawn)

        ! Set the sign of the child walkers based upon the sign of the
        ! Hamiltonian matrix element connecting the parent determinant to the
        ! child determinant.
        !
        ! In:
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    nspawn: On input, the unsigned population of child walkers spawned
        !    from |D> to |D'>.  On output, the *signed* population of child
        !    walkers produced by this spawning attempt.

        real(p), intent(in) :: hmatel
        integer, intent(in) :: parent_sign
        integer, intent(inout) :: nspawn

        ! If H_ij is positive, then the spawned walker is of opposite
        ! sign to the parent, otherwise the spawned walkers if of the same
        ! sign as the parent.
        if (hmatel > 0.0_p) then
            nspawn = -sign(nspawn, parent_sign)
        else
            nspawn = sign(nspawn, parent_sign)
        end if

    end subroutine set_child_sign

    function attempt_to_spawn(rng, hmatel, pgen, parent_sign) result(nspawn)

        ! In:
        !    hmatel: Hamiltonian matrix element connecting a determinant and an
        !    excitation of that determinant.
        !    pgen: probability of generating the excitation.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    rng: random number generator.
        ! Returns:
        !    number of particles spawned.  0 indicates the spawning attempt was
        !    unsuccessful.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use fciqmc_data, only: tau

        integer :: nspawn

        real(p), intent(in) :: pgen, hmatel
        integer, intent(in) :: parent_sign
        type(dSFMT_t), intent(inout) :: rng

        real(p) :: pspawn

        ! Essentially combines nspawn_from_prob and set_child_sign into one
        ! convenient routine...

        ! 1. Calculate probability spawning is successful.
        pspawn = tau*abs(hmatel)/pgen

        ! 2. Attempt spawning.

        ! Need to take into account the possibilty of a spawning attempt
        ! producing multiple offspring...
        ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and
        ! then spawn a particle with probability pspawn-floor(pspawn).
        nspawn = int(pspawn)
        pspawn = pspawn - nspawn

        if (pspawn > get_rand_close_open(rng)) nspawn = nspawn + 1

        if (nspawn > 0) then

            ! 3. If H_ij is positive, then the spawned walker is of opposite
            ! sign to the parent, otherwise the spawned walkers if of the same
            ! sign as the parent.
            if (hmatel > 0.0_p) then
                nspawn = -sign(nspawn, parent_sign)
            else
                nspawn = sign(nspawn, parent_sign)
            end if

        end if

    end function attempt_to_spawn

!--- Assuming spawning is successful, create new particle appropriately ---

    function assign_particle_processor(particle_label, length, seed, shift, freq, np) result(particle_proc)

        ! In:
        !    particle_label: bit string which describes the location/basis
        !       function/etc of the particle (ie psip or excip).
        !    length: size of particle_label.
        !    seed: seed to pass to the hashing function.
        !    shift: value to add to the hash of the label before determining
        !       the processor to which the label is assigned.
        !    freq: frequency over which the result changes exactly once.
        !       See comments below.  Ignored if the shift is 0.  Must be smaller
        !       than 32.
        !    np: number of processors over which the particles are to be
        !       distributed.
        ! Returns:
        !    particle_proc: processor number (in range [0,np-1] on which the
        !       particle should reside.

        use hashing, only: murmurhash_bit_string

        integer :: particle_proc
        integer(i0), intent(in) :: particle_label(length)
        integer, intent(in) :: length, seed, shift, freq, np

        integer :: hash, offset
        integer(i0) :: mod_label(length)

        ! (Extra credit for parallel calculations)
        ! Hash the label to get a (hopefully uniform) distribution across all
        ! possible particle labels and then modulo it to assign each label in
        ! a (hopefully uniform) fashion.
        hash = murmurhash_bit_string(particle_label, length, seed)
        if (shift == 0) then
            ! p = hash(label) % np
            particle_proc = modulo(hash, np)
        else
            ! o = [ hash(label) + shift ] >> freq
            ! p = [ hash(label) + o ] % np
            ! Explanation:
            ! We wish to slowly vary the processor a label is assigned to.
            ! The shift is a fast(ish) varying value (e.g. the iteration
            ! number).
            ! [ hash(label) + shift ] >> freq changes exactly once in 2^freq
            ! consecutive values of the shift, i.e. when the freq lower bits of
            ! [ hash(label) + shift ] is greater than 2^freq.  We add this
            ! offset onto the label and rehash.  label+offset varies once every
            ! 2^freq values of the shift and hence the assigned processor
            ! changes at most once in this window.
            offset = ishft(hash+shift, -freq)
            mod_label = particle_label + offset
            hash = murmurhash_bit_string(mod_label, length, seed)
            particle_proc = modulo(hash, np)
        end if

    end function assign_particle_processor

    subroutine create_spawned_particle(cdet, connection, nspawn, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use hashing
        use parallel, only: iproc, nprocs, nthreads
        use omp_lib

        use basis, only: basis_length
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det
        use spawn_data, only: spawn_t

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        integer(i0) :: f_new(basis_length)
        integer :: iproc_spawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

        iproc_spawn = assign_particle_processor(f_new, basis_length, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs)

        ! Move to the next position in the spawning array.
        spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
        spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new
        spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn

    end subroutine create_spawned_particle

    subroutine create_spawned_particle_initiator(cdet, connection, nspawn, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use hashing
        use parallel, only: iproc, nprocs, nthreads
        use omp_lib

        use basis, only: basis_length
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det
        use spawn_data, only: spawn_t

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        integer(i0) :: f_new(basis_length)
        integer :: iproc_spawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

        iproc_spawn = assign_particle_processor(f_new, basis_length, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs)

        ! Move to the next position in the spawning array.
        spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
        spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new
        spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn
        ! initiator_flag: flag indicating the staturs of the parent determinant.
        !     initiator_flag = 0 indicates the parent is an initiator.
        !     initiator_flag = 1 indicates the parent is not an initiator.
        spawn%sdata(spawn%flag_indx,spawn%head(thread_id,iproc_spawn)) = cdet%initiator_flag

    end subroutine create_spawned_particle_initiator

    subroutine create_spawned_particle_truncated(cdet, connection, nspawn, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! TODO: comment

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use hashing
        use parallel, only: iproc, nprocs, nthreads
        use omp_lib

        use basis, only: basis_length
        use calc, only: truncation_level
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det, get_excitation_level
        use fciqmc_data, only: hs_f0
        use spawn_data, only: spawn_t

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        integer(i0) :: f_new(basis_length)
        integer :: iproc_spawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

        ! Only accept spawning if it's within the truncation level.
        if (get_excitation_level(hs_f0, f_new) <= truncation_level) then

            iproc_spawn = assign_particle_processor(f_new, basis_length, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs)

            ! Move to the next position in the spawning array.
            spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

            ! Set info in spawning array.
            ! Zero it as not all fields are set.
            spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
            spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new
            spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn

        end if

    end subroutine create_spawned_particle_truncated

    subroutine create_spawned_particle_initiator_truncated(cdet, connection, nspawn, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! TODO: comment

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use hashing
        use parallel, only: iproc, nprocs, nthreads
        use omp_lib

        use basis, only: basis_length
        use calc, only: truncation_level
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det, get_excitation_level
        use fciqmc_data, only: hs_f0
        use spawn_data, only: spawn_t

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        integer(i0) :: f_new(basis_length)
        integer :: iproc_spawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

        ! Only accept spawning if it's within the truncation level.
        if (get_excitation_level(hs_f0, f_new) <= truncation_level) then

            iproc_spawn = assign_particle_processor(f_new, basis_length, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs)

            ! Move to the next position in the spawning array.
            spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

            ! Set info in spawning array.
            ! Zero it as not all fields are set.
            spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
            spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new
            spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn
            ! initiator_flag: flag indicating the staturs of the parent determinant.
            !     initiator_flag = 0 indicates the parent is an initiator.
            !     initiator_flag = 1 indicates the parent is not an initiator.
            spawn%sdata(spawn%flag_indx,spawn%head(thread_id,iproc_spawn)) = cdet%initiator_flag

        end if

    end subroutine create_spawned_particle_initiator_truncated

    subroutine create_spawned_particle_ras(cdet, connection, nspawn, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! TODO: comment

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use hashing
        use parallel, only: iproc, nprocs, nthreads
        use omp_lib

        use basis, only: basis_length
        use bit_utils, only: count_set_bits
        use calc, only: truncation_level, ras1, ras3, ras1_min, ras3_max
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det, get_excitation_level, in_ras
        use spawn_data, only: spawn_t

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        integer(i0) :: f_new(basis_length)
        integer :: iproc_spawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

        ! Only accept spawning if it's within the RAS space.
        if (in_ras(ras1, ras3, ras1_min, ras3_max, f_new)) then

            iproc_spawn = assign_particle_processor(f_new, basis_length, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs)

            ! Move to the next position in the spawning array.
            spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

            ! Set info in spawning array.
            ! Zero it as not all fields are set.
            spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
            spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new
            spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn

        end if

    end subroutine create_spawned_particle_ras

    subroutine create_spawned_particle_initiator_ras(cdet, connection, nspawn, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! TODO: comment

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use hashing
        use parallel, only: iproc, nprocs, nthreads
        use omp_lib

        use basis, only: basis_length
        use bit_utils, only: count_set_bits
        use calc, only: truncation_level, ras1, ras3, ras1_min, ras3_max
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det, get_excitation_level, in_ras
        use spawn_data, only: spawn_t

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        integer(i0) :: f_new(basis_length)
        integer :: iproc_spawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

        ! Only accept spawning if it's within the RAS space.
        if (in_ras(ras1, ras3, ras1_min, ras3_max, f_new)) then

            iproc_spawn = assign_particle_processor(f_new, basis_length, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs)

            ! Move to the next position in the spawning array.
            spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

            ! Set info in spawning array.
            ! Zero it as not all fields are set.
            spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
            spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new
            spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn
            ! initiator_flag: flag indicating the staturs of the parent determinant.
            !     initiator_flag = 0 indicates the parent is an initiator.
            !     initiator_flag = 1 indicates the parent is not an initiator.
            spawn%sdata(spawn%flag_indx,spawn%head(thread_id,iproc_spawn)) = cdet%initiator_flag

        end if

    end subroutine create_spawned_particle_initiator_ras

    subroutine create_spawned_particle_density_matrix(f1, f2, connection, nspawn, spawning_end, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    f1: bitstring corresponding to the end which is currently
        !         being spawned from.
        !    f2: bitstring corresponding to the 'inactive' end, which
        !         was not involved in the current spawning step.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    spawning_end: Specifies which 'end' we are spawning from
        !        currently, ie, if the elements of the density matrix are
        !        \rho_{i,j}, are we spawning from the i end, 1, or the j end, 2.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use basis, only: basis_length, total_basis_length
        use errors, only: stop_all
        use excitations, only: excit, create_excited_det
        use hashing
        use parallel, only: iproc, nprocs, nthreads
        use spawn_data, only: spawn_t

        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer, intent(in) :: nspawn
        integer, intent(in) :: spawning_end
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        type(excit), intent(in) :: connection
        integer(i0) :: f_new(basis_length)
        integer(i0) :: f_new_tot(total_basis_length)

        integer :: iproc_spawn
        integer, parameter :: thread_id = 0

        ! Create bit string of new determinant. The entire two-ended
        ! bitstring is eventually stored in f_new_tot.
        call create_excited_det(f1, connection, f_new)
        f_new_tot = 0
        if (spawning_end==1) then
            f_new_tot(:basis_length) = f_new
            f_new_tot((basis_length+1):(total_basis_length)) = f2
        else
            f_new_tot(:basis_length) = f2
            f_new_tot((basis_length+1):(total_basis_length)) = f_new
        end if

        iproc_spawn = assign_particle_processor(f_new_tot, total_basis_length, spawn%hash_seed, &
                                                spawn%hash_shift, spawn%move_freq, nprocs)

        ! Move to the next position in the spawning array.
        spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

        ! spawn%head_start(0,1) holds the number of slots in the spawning array per processor.
        if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(0,iproc_spawn) >= spawn%head_start(0,1)) &
            call stop_all('create_spawned_particle_density_matrix',&
                           'There is no space left in the spawning array.')

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
        spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new_tot
        spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn

    end subroutine create_spawned_particle_density_matrix

    subroutine create_spawned_particle_half_density_matrix(f1, f2, connection, nspawn, spawning_end, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! If walker tries to spawn in the lower triangle of density matrix
        ! then reflect it to upper triangle by swapping bit strings for the
        ! two ends. The current position in the spawning array is updated.

        ! In:
        !    f1: bitstring corresponding to the end which is currently
        !         being spawned from.
        !    f2: bitstring corresponding to the 'inactive' end, which
        !         was not involved in the current spawning step.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    spawning_end: Specifies which 'end' we are spawning from
        !        currently, ie, if the elements of the density matrix are
        !        \rho_{i,j}, are we spawning from the i end, 1, or the j end, 2.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use determinants, only: det_compare
        use basis, only: basis_length, total_basis_length
        use calc, only: truncation_level
        use errors, only: stop_all
        use excitations, only: excit, create_excited_det
        use hashing
        use parallel, only: iproc, nprocs
        use spawn_data, only: spawn_t

        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer, intent(in) :: nspawn
        integer, intent(in) :: spawning_end
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        type(excit), intent(in) :: connection
        integer(i0) :: f_new(basis_length)
        integer(i0) :: f_new_tot(total_basis_length)

        integer :: iproc_spawn
        integer, parameter :: thread_id = 0

        ! Create bit string of new determinant. The entire two-ended
        ! bitstring is eventually stored in f_new_tot.
        call create_excited_det(f1, connection, f_new)
        f_new_tot = 0

        ! Test to see whether the new determinant resides in the upper
        ! triangle of the density matrix. If so keep bit string ends
        ! as they are. If not then swap bitstring ends so that the
        ! new psips are reflected into the upper triangle of the density
        ! matrix as they try to spawn.
        if (det_compare(f_new, f2, basis_length) == -1) then
            f_new_tot(:basis_length) = f_new
            f_new_tot((basis_length+1):(total_basis_length)) = f2
        else
            f_new_tot(:basis_length) = f2
            f_new_tot((basis_length+1):(total_basis_length)) = f_new
        end if

        iproc_spawn = assign_particle_processor(f_new_tot, total_basis_length, spawn%hash_seed, &
                                                spawn%hash_shift, spawn%move_freq, nprocs)

        ! Move to the next position in the spawning array.
        spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + 1

        ! spawn%head_start(0,1) holds the number of slots in the spawning array per processor.
        if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(0,iproc_spawn) >= spawn%head_start(0,1)) &
            call stop_all('create_spawned_particle_half_density_matrix',&
                           'There is no space left in the spawning array.')

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
        spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new_tot
        spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn

    end subroutine create_spawned_particle_half_density_matrix

    subroutine create_spawned_particle_truncated_half_density_matrix(f1, f2, connection, nspawn, &
                                                                          spawning_end, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! A spawned walker is only created on (f1', f2) if f1' and f2 do not differ by
        ! more than truncation_level basis functions, where f1' is obtained by
        ! applying the connection to f1.

        ! Note: This is the half density matrix version of create_spawned_particle_truncated_density_matrix
        ! It works in an ientical way apart from attempts to spawn on the lower triangle of the
        ! density matrix are reflected so that they are spawned on the upper triangle of the
        ! density matrix.

        ! In:
        !    f1: bitstring corresponding to the end which is currently
        !         being spawned from.
        !    f2: bitstring corresponding to the 'inactive' end, which
        !         was not involved in the current spawning step.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    spawning_end: Specifies which 'end' we are spawning from
        !        currently, ie, if the elements of the density matrix are
        !        \rho_{i,j}, are we spawning from the i end, 1, or the j end, 2.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use basis, only: basis_length, total_basis_length
        use determinants, only: det_compare
        use calc, only: truncation_level
        use errors, only: stop_all
        use excitations, only: excit, create_excited_det, get_excitation_level
        use hashing
        use parallel, only: iproc, nprocs
        use spawn_data, only: spawn_t

        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer, intent(in) :: nspawn
        integer, intent(in) :: spawning_end
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        type(excit), intent(in) :: connection
        integer(i0) :: f_new(basis_length)
        integer(i0) :: f_new_tot(total_basis_length)

        integer :: iproc_spawn
        integer, parameter :: thread_id = 0

        ! Create bit string of new determinant. The entire two-ended
        ! bitstring is eventually stored in f_new_tot.
        call create_excited_det(f1, connection, f_new)

        if (get_excitation_level(f2, f_new) <= truncation_level) then

            f_new_tot = 0
            ! Test to see whether the new determinant resides in the upper
            ! triangle of the density matrix. If so keep bit string ends
            ! as they are. If not then swap bitstring ends so that the
            ! new psips are reflected into the upper triangle of the density
            ! matrix as they try to spawn.
            if (det_compare(f_new, f2, basis_length) == -1) then
                f_new_tot(:basis_length) = f_new
                f_new_tot((basis_length+1):(total_basis_length)) = f2
            else
                f_new_tot(:basis_length) = f2
                f_new_tot((basis_length+1):(total_basis_length)) = f_new
            end if

            iproc_spawn = assign_particle_processor(f_new_tot, total_basis_length, spawn%hash_seed, &
                                                    spawn%hash_shift, spawn%move_freq, nprocs)

            ! Move to the next position in the spawning array.
            spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + 1

            ! spawn%head_start(0,1) holds the number of slots in the spawning array per processor.
            if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(0,iproc_spawn) >= spawn%head_start(0,1)) &
                call stop_all('create_spawned_particle_truncated_half_density_matrix',&
                               'There is no space left in the spawning array.')

            ! Set info in spawning array.
            ! Zero it as not all fields are set.
            spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
            spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new_tot
            spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn

        end if

    end subroutine create_spawned_particle_truncated_half_density_matrix

    subroutine create_spawned_particle_truncated_density_matrix(f1, f2, connection, nspawn, spawning_end, particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! A spawned walker is only created on (f1', f2) if f1' and f2 do not differ by
        ! more than truncation_level basis functions, where f1' is obtained by
        ! applying the connection to f1.

        ! In:
        !    f1: bitstring corresponding to the end which is currently
        !         being spawned from.
        !    f2: bitstring corresponding to the 'inactive' end, which
        !         was not involved in the current spawning step.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    spawning_end: Specifies which 'end' we are spawning from
        !        currently, ie, if the elements of the density matrix are
        !        \rho_{i,j}, are we spawning from the i end, 1, or the j end, 2.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use basis, only: basis_length, total_basis_length
        use calc, only: truncation_level
        use errors, only: stop_all
        use excitations, only: excit, create_excited_det, get_excitation_level
        use hashing
        use parallel, only: iproc, nprocs, nthreads
        use spawn_data, only: spawn_t

        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer, intent(in) :: nspawn
        integer, intent(in) :: spawning_end
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        type(excit), intent(in) :: connection
        integer(i0) :: f_new(basis_length)
        integer(i0) :: f_new_tot(total_basis_length)

        integer :: iproc_spawn
        integer, parameter :: thread_id = 0

        ! Create bit string of new determinant. The entire two-ended
        ! bitstring is eventually stored in f_new_tot.
        call create_excited_det(f1, connection, f_new)

        if (get_excitation_level(f2, f_new) <= truncation_level) then

            f_new_tot = 0
            if (spawning_end==1) then
                f_new_tot(:basis_length) = f_new
                f_new_tot((basis_length+1):(total_basis_length)) = f2
            else
                f_new_tot(:basis_length) = f2
                f_new_tot((basis_length+1):(total_basis_length)) = f_new
            end if

            iproc_spawn = assign_particle_processor(f_new_tot, total_basis_length, spawn%hash_seed, &
                                                    spawn%hash_shift, spawn%move_freq, nprocs)

            ! Move to the next position in the spawning array.
            spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

            ! spawn%head_start(0,1) holds the number of slots in the spawning array per processor.
            if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(0,iproc_spawn) >= spawn%head_start(0,1)) &
                call stop_all('create_spawned_particle_truncated_density_matrix', &
                               'There is no space left in the spawning array.')

            ! Set info in spawning array.
            ! Zero it as not all fields are set.
            spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0
            spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = f_new_tot
            spawn%sdata(total_basis_length+particle_type,spawn%head(thread_id,iproc_spawn)) = nspawn

        end if

    end subroutine create_spawned_particle_truncated_density_matrix

    subroutine create_spawned_particle_rdm(irdm, nspawn_in, particle_type, rdm_spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    irdm: the index of the rdm that this psip contributes to.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    rdm_spawn: rdm_spawn_t object to which the spanwed particle will be added.

        use bit_utils, only: operator(.bitstrgt.)
        use dmqmc_procedures, only: rdms
        use errors, only: stop_all
        use fciqmc_data, only: rdm_spawn_t
        use hashing
        use parallel, only: iproc, nprocs, nthreads
        use hash_table, only: hash_table_pos_t, lookup_hash_table_entry
        use hash_table, only: assign_hash_table_entry
        use utils, only: int_fmt
                              

        integer, intent(in) :: irdm
        integer, intent(in) :: nspawn_in
        integer, intent(in) :: particle_type
        type(rdm_spawn_t), intent(inout) :: rdm_spawn
        integer :: nspawn, rdm_bl

        integer(i0) :: f_new_tot(2*rdms(irdm)%rdm_basis_length)
        integer(i0) :: f_temp1(rdms(irdm)%rdm_basis_length), f_temp2(rdms(irdm)%rdm_basis_length)

        integer :: iproc_spawn
        ! WARNING!  The below algorithm is *not* suitable for conversion to
        ! thread-safety as each thread could be spawning onto the same RDM
        ! element, yet the hash table requires a given element to exist in one
        ! (and only one) location, which is not compatible with how the threaded
        ! spawning arrays work.
        integer, parameter :: thread_id = 0

        type(hash_table_pos_t) :: pos
        logical :: hit
        integer :: err_code

        rdm_bl = rdms(irdm)%rdm_basis_length
        ! nspawn will be doubles for diagonal elements.
        nspawn = nspawn_in

        ! To enforce that the rdm is symmetric, add this psip to both \rho_{ij} and
        ! \rho_{ji}. These will lead to an rdm with a trace twice what it would have
        ! been, but this is not a problem.
        f_new_tot = 0

        f_temp1 = rdms(irdm)%end1
        f_temp2 = rdms(irdm)%end2

        if (f_temp1 .bitstrgt. f_temp2) then
            ! If below the diagonal, swap the bitstrings so that the spawning occurs above it.
            f_new_tot(:rdm_bl) = f_temp2
            f_new_tot(rdm_bl+1:2*rdm_bl) = f_temp1
        else
            f_new_tot(:rdm_bl) = f_temp1
            f_new_tot(rdm_bl+1:2*rdm_bl) = f_temp2
            if (all(f_temp1 == f_temp2)) then
                ! Because off-diagonal elements have been doubled (elements above the diagonal taking
                ! contributions from both below and above it), we must double the diagonal elements too.
                nspawn = nspawn*2
            end if
        end if

        associate(spawn=>rdm_spawn%spawn, ht=>rdm_spawn%ht, bsl=>rdm_spawn%spawn%bit_str_len)

            iproc_spawn = assign_particle_processor(f_new_tot, 2*rdm_bl, spawn%hash_seed, &
                                                    spawn%hash_shift, spawn%move_freq, nprocs)

            call lookup_hash_table_entry(ht, f_new_tot, pos, hit)

            if (hit) then
                associate(indx => ht%table(pos%ientry,pos%islot))
                    ! Direct annihilation/cancellation in spawning array.
                    spawn%sdata(bsl+particle_type,indx) = spawn%sdata(bsl+particle_type,indx) + nspawn
                end associate
            else
                ! Move to the next position in the spawning array.
                call assign_hash_table_entry(ht, pos%islot, pos, err_code)
                if (err_code /= 0) then
                    write(6,'(1X,a9,'//int_fmt(err_code,1)//')') 'err_code:', err_code
                    call stop_all('create_spawned_particle_rdm','Error in assigning hash &
                                  &table entry.')
                end if
                ! Fix hash table to point to the head of the spawn data
                ! for this thread/processor.
                spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

                ! spawn%head_start(0,1) holds the number of slots in the spawning array per processor.
                if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(0,iproc_spawn) >= spawn%head_start(0,1)) &
                    call stop_all('create_spawned_particle_rdm','There is no space left in the RDM array.')

                ht%table(pos%ientry,pos%islot) = spawn%head(thread_id,iproc_spawn)
                associate(indx => ht%table(pos%ientry,pos%islot))
                    ! Set info in spawning array.
                    ! Zero it as not all fields are set.
                    spawn%sdata(:,indx) = 0
                    spawn%sdata(:bsl,indx) = f_new_tot
                    spawn%sdata(bsl+particle_type,indx) = nspawn
                end associate
            end if

        end associate

    end subroutine create_spawned_particle_rdm

end module spawning
