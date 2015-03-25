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

    subroutine spawn_standard(rng, sys, qmc_in, tau, spawn_cutoff, real_factor, cdet, parent_sign, &
                              gen_excit_ptr, weights, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant.

        ! This is just a thin wrapper around a system-specific excitation
        ! generator and a utility function.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    tau: timestep.
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !        gen_excit_ptr%full *must* be set to a procedure which generates
        !        a complete excitation.
        !    weights: importance sampling weights.
        ! Out:
        !    nspawn: number of particles spawned, in the encoded representation.
        !        0 indicates the spawning attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use qmc_data, only: qmc_in_t
        use system, only: sys_t
        use proc_pointers, only: gen_excit_ptr_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        real(p), intent(in) :: tau
        integer(int_p), intent(in) :: spawn_cutoff
        integer(int_p), intent(in) :: real_factor
        type(det_info_t), intent(in) :: cdet
        integer(int_p), intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        real(p), allocatable, intent(in) :: weights(:)
        integer(int_p), intent(out) :: nspawn
        type(excit_t), intent(out) :: connection

        real(p) :: pgen, hmatel

        ! 1. Generate random excitation.
        call gen_excit_ptr%full(rng, sys, qmc_in, cdet, pgen, connection, hmatel)

        ! 2. Attempt spawning.
        nspawn = attempt_to_spawn(rng, tau, spawn_cutoff, real_factor, hmatel, pgen, parent_sign)

    end subroutine spawn_standard

    subroutine spawn_importance_sampling(rng, sys, qmc_in, tau, spawn_cutoff, real_factor, cdet, parent_sign, &
                                         gen_excit_ptr, weights, nspawn, connection)

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
        !    qmc_in: input options relating to QMC methods.
        !    tau: timestep.
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !        gen_excit_ptr%full *must* be set to a procedure which generates
        !        a complete excitation.
        !    weights: importance sampling weights.
        ! Out:
        !    nspawn: number of particles spawned, in the encoded representation.
        !        0 indicates the spawning attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: gen_excit_ptr_t
        use qmc_data, only: qmc_in_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        real(p), intent(in) :: tau
        integer(int_p), intent(in) :: spawn_cutoff
        integer(int_p), intent(in) :: real_factor
        type(det_info_t), intent(in) :: cdet
        integer(int_p), intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        real(p), allocatable, intent(in) :: weights(:)
        integer(int_p), intent(out) :: nspawn
        type(excit_t), intent(out) :: connection

        real(p) :: pgen, hmatel

        ! 1. Generate random excitation.
        call gen_excit_ptr%full(rng, sys, qmc_in, cdet, pgen, connection, hmatel)

        ! 2. Transform Hamiltonian matrix element by trial function.
        call gen_excit_ptr%trial_fn(sys, cdet, connection, weights, hmatel)

        ! 3. Attempt spawning.
        nspawn = attempt_to_spawn(rng, tau, spawn_cutoff, real_factor, hmatel, pgen, parent_sign)

    end subroutine spawn_importance_sampling

    subroutine spawn_lattice_split_gen(rng, sys, qmc_in, tau, spawn_cutoff, real_factor, cdet, parent_sign, &
                                       gen_excit_ptr, weights, nspawn, connection)

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
        !    qmc_in: input options relating to QMC methods.
        !    tau: timestep.
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
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
        !    weights: importance sampling weights.
        ! Out:
        !    nspawn: number of particles spawned, in the encoded representation.
        !        0 indicates the spawning attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: gen_excit_ptr_t
        use qmc_data, only: qmc_in_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        real(p), intent(in) :: tau
        integer(int_p), intent(in) :: spawn_cutoff
        integer(int_p), intent(in) :: real_factor
        type(det_info_t), intent(in) :: cdet
        integer(int_p), intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        real(p), allocatable, intent(in) :: weights(:)
        integer(int_p), intent(out) :: nspawn
        type(excit_t), intent(out) :: connection

        real(p) :: pgen, abs_hmatel, hmatel

        ! 1. Generate enough of a random excitation to determinant the
        ! generation probability and |H_ij|.
        call gen_excit_ptr%init(rng, sys, qmc_in, cdet, pgen, connection, abs_hmatel)

        ! 2. Attempt spawning.
        nspawn = nspawn_from_prob(rng, spawn_cutoff, real_factor, tau*abs_hmatel/pgen)

        if (nspawn /= 0_int_p) then

            ! 3. Complete excitation and find sign of connecting matrix element.
            call gen_excit_ptr%finalise(rng, sys, qmc_in, cdet, connection, hmatel)

            ! 4. Find sign of offspring.
            call set_child_sign(hmatel, parent_sign, nspawn)

        end if

    end subroutine spawn_lattice_split_gen

    subroutine spawn_lattice_split_gen_importance_sampling(rng, sys, qmc_in, tau, spawn_cutoff, real_factor, cdet, parent_sign, &
                                                           gen_excit_ptr, weights, nspawn, connection)

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
        !    qmc_in: input options relating to QMC methods.
        !    tau: timestep.
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
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
        !    weights: importance sampling weights.
        ! Out:
        !    nspawn: number of particles spawned, in the encoded representation.
        !        0 indicates the spawning attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: gen_excit_ptr_t
        use qmc_data, only: qmc_in_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        real(p), intent(in) :: tau
        integer(int_p), intent(in) :: spawn_cutoff
        integer(int_p), intent(in) :: real_factor
        type(det_info_t), intent(in) :: cdet
        integer(int_p), intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        real(p), allocatable, intent(in) :: weights(:)
        integer(int_p), intent(out) :: nspawn
        type(excit_t), intent(out) :: connection

        real(p) :: pgen, tilde_hmatel, hmatel

        ! 1. Generate enough of a random excitation to determinant the
        ! generation probability and |H_ij|.
        call gen_excit_ptr%init(rng, sys, qmc_in, cdet, pgen, connection, tilde_hmatel)

        ! 2. Transform Hamiltonian matrix element by trial function.
        call gen_excit_ptr%trial_fn(sys, cdet, connection, weights, tilde_hmatel)

        ! 3. Attempt spawning.
        nspawn = nspawn_from_prob(rng, spawn_cutoff, real_factor, tau*abs(tilde_hmatel)/pgen)

        if (nspawn /= 0_int_p) then

            ! 4. Complete excitation and find sign of connecting matrix element.
            ! *NOTE*: this returns the original matrix element and *not* the
            ! matrix element after the trial function transformation.
            call gen_excit_ptr%finalise(rng, sys, qmc_in, cdet, connection, hmatel)

            ! 5. Find sign of offspring.
            ! Note that we don't care about the value of H_ij at this step, only
            ! the sign.
            call set_child_sign(tilde_hmatel*hmatel, parent_sign, nspawn)

        end if

    end subroutine spawn_lattice_split_gen_importance_sampling

    subroutine spawn_null(rng, sys, qmc_in, tau, spawn_cutoff, real_factor, cdet, parent_sign, gen_excit_ptr, weights, &
                          nspawn, connection)

        ! This is a null spawning routine for use with operators which are
        ! diagonal in the basis and hence only have a cloning step in the
        ! Hellmann-Feynman sampling.  It does *nothing*.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    tau: timestep.
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from (or not, in this case!).
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !    weights: importance sampling weights.
        ! Out:
        !    nspawn: number of particles spawned, in the encoded representation.
        !        0 indicates the spawning attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: gen_excit_ptr_t
        use qmc_data, only: qmc_in_t
        use dSFMT_interface, only: dSFMT_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        real(p), intent(in) :: tau
        integer(int_p), intent(in) :: spawn_cutoff
        integer(int_p), intent(in) :: real_factor
        type(det_info_t), intent(in) :: cdet
        integer(int_p), intent(in) :: parent_sign
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        real(p), allocatable, intent(in) :: weights(:)
        integer(int_p), intent(out) :: nspawn
        type(excit_t), intent(out) :: connection

        ! Just some null operations to avoid -Wall -Werror causing errors.
        connection%nexcit = huge(0)

        ! Return nspawn = 0 as we don't want to do any spawning.
        nspawn = 0_int_p

    end subroutine spawn_null

!--- Attempt spawning based upon random excitation ---

    function nspawn_from_prob(rng, spawn_cutoff, real_factor, probability) result(nspawn)

        ! Generate the number spawned from a probability. If probability is greater than
        ! zero, then number spawned = int(probability) + stochastic{0,1}
        ! where the latter half of the RHS is a stochastic spawning from the remainder
        !
        ! In:
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
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
        integer(int_p), intent(in) :: spawn_cutoff
        integer(int_p), intent(in) :: real_factor
        integer(int_s) :: nspawn
        real(p) :: pspawn

        ! 'Encode' the spawning probability by multiplying by 2^(real_bit_shift).
        ! We then stochastically round this probability either up or down to
        ! the nearest integers. This allows a resolution of 2^(-real_spawning)
        ! when we later divide this factor back out. (See comments for
        ! particle_t%pops).
        pspawn = probability*real_factor

        if (abs(pspawn) < spawn_cutoff) then

            ! If the spawning amplitude is below the minimum spawning event
            ! allowed, stochastically round it either down to zero or up
            ! to the cutoff.
            if (pspawn > get_rand_close_open(rng)*spawn_cutoff) then
                nspawn = spawn_cutoff
            else
                nspawn = 0_int_p
            end if

        else

            ! Need to take into account the possibilty of a spawning attempt
            ! producing multiple offspring...
            ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and
            ! then spawn a particle with probability pspawn-floor(pspawn).
            nspawn = int(pspawn, int_p)
            pspawn = pspawn - nspawn

            if (pspawn > get_rand_close_open(rng)) nspawn = nspawn + 1_int_p

        end if

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
        integer(int_p), intent(in) :: parent_sign
        integer(int_p), intent(inout) :: nspawn

        ! If H_ij is positive, then the spawned walker is of opposite
        ! sign to the parent, otherwise the spawned walkers if of the same
        ! sign as the parent.
        if (hmatel > 0.0_p) then
            nspawn = -sign(nspawn, parent_sign)
        else
            nspawn = sign(nspawn, parent_sign)
        end if

    end subroutine set_child_sign

    function attempt_to_spawn(rng, tau, spawn_cutoff, real_factor, hmatel, pgen, parent_sign) result(nspawn)

        ! In:
        !    tau: timestep being used.
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
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

        integer(int_p) :: nspawn

        real(p), intent(in) :: tau, pgen, hmatel
        integer(int_p), intent(in) :: parent_sign
        type(dSFMT_t), intent(inout) :: rng
        integer(int_p), intent(in) :: spawn_cutoff
        integer(int_p), intent(in) :: real_factor
        real(p) :: pspawn

        ! Calculate probability spawning is successful.
        pspawn = tau*abs(hmatel)/pgen

        ! 'Encode' the spawning probability by multiplying by 2^(real_bit_shift).
        ! We then stochastically round this probability either up or down to
        ! the nearest integers. This allows a resolution of 2^(-real_spawning)
        ! when we later divide this factor back out. (See comments for
        ! particle_t%pops).
        pspawn = pspawn*real_factor

        if (pspawn < spawn_cutoff) then

            ! If the spawning amplitude is below the minimum spawning event
            ! allowed, stochastically round it either down to zero or up
            ! to the cutoff.
            if (pspawn > get_rand_close_open(rng)*spawn_cutoff) then
                nspawn = spawn_cutoff
            else
                nspawn = 0_int_p
            end if

        else

            ! Need to take into account the possibilty of a spawning attempt
            ! producing multiple offspring...
            ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and
            ! then spawn a particle with probability pspawn-floor(pspawn).
            nspawn = int(pspawn, int_p)
            pspawn = pspawn - nspawn

            if (pspawn > get_rand_close_open(rng)) nspawn = nspawn + 1_int_p

        end if

        if (nspawn > 0) then
            ! If H_ij is positive, then the spawned walker is of opposite
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

    subroutine assign_particle_processor(particle_label, nbits, seed, shift, freq, np, particle_proc, slot_pos, proc_map, nslots)

        ! In:
        !    particle_label: bit string which describes the location/basis
        !       function/etc of the particle (ie psip or excip).
        !    nbits: length (in bits) of particle_label.  This allows us to ignore
        !       any additional padding at the end of the bit string for different
        !       sizes of i0 integers.
        !    seed: seed to pass to the hashing function.
        !    shift: value to add to the hash of the label before determining
        !       the processor to which the label is assigned.
        !    freq: frequency over which the result changes exactly once.
        !       See comments below.  Ignored if the shift is 0.  Must be smaller
        !       than 32.
        !    np: number of processors over which the particles are to be
        !       distributed.
        !    proc_map: array which maps determinants to processors.
        !    nslots: number of slots proc_map is divided into.
        ! Out:
        !    particle_proc: processor where determinant resides
        !    slot_pos: position in proc_map for this determinant

        use hashing, only: murmurhash_bit_string

        integer(i0), intent(in) :: particle_label(:)
        integer, intent(in) :: nbits, seed, shift, freq, np
        integer, intent(in) :: proc_map(0:)
        integer, intent(in) :: nslots
        integer, intent(out) :: particle_proc, slot_pos

        integer :: hash, offset, i, tmp1, tmp2
        integer(i0) :: mod_label(size(particle_label))

        ! (Extra credit for parallel calculations)
        ! Hash the label to get a (hopefully uniform) distribution across all
        ! possible particle labels and then modulo it to assign each label in
        ! a (hopefully uniform) fashion.
        hash = murmurhash_bit_string(particle_label, nbits, seed)
        if (shift == 0) then
            ! p = hash(label) % np
            slot_pos = modulo(hash, np*nslots)
            particle_proc = proc_map(slot_pos)
        else
            ! o = [ hash(label) + shift ] >> freq
            ! p = [ hash(label + o) ] % np
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
! [review] - AJWT: I had a recollection that offset only added to the lowest element of the array - or at least that's how I had
! [review] - AJWT: conceived of it.  I don't think it makes much difference to add it to all of the elements if the hash function is any good, so might save some coding and space just to add it to the last element?
            if (i0_length == 32) then
                mod_label = particle_label + offset
            else
                ! Ensure the offset is applied bitwise identically every 32 bits.
                do i = 1, size(particle_label)
                    tmp1 = transfer(particle_label(i), tmp1) + offset
                    tmp2 = transfer(ishft(particle_label(i),-32), tmp1) + offset
                    mod_label(i) = ior(transfer(tmp1, 0_i0), ishft(transfer(tmp2,0_i0),32))
                end do
            end if
            hash = murmurhash_bit_string(mod_label, nbits, seed)
            slot_pos = modulo(hash, np*nslots)
            particle_proc = proc_map(slot_pos)
        end if

    end subroutine assign_particle_processor

    subroutine assign_particle_processor_dmqmc(particle_label, nbits, seed, shift, freq, np, particle_proc, slot_pos, &
                                               proc_map, nslots)

        ! Wrapper around assign_particle_processor to ensure we hash the same
        ! amount of data for tensor labels in DMQMC (which involve two labels)
        ! irrespective of DET_SIZE.

        ! The tensor label is formed by concatenating together the labels for
        ! both determinants.  We therefore need to ensure the same amount of
        ! padding exists between the two labels for DET_SIZE=32 and DET_SIZE=64.
        ! We do this by inserting an additional integer in the DET_SIZE=32 case
        ! if required.

        ! In:
        !    particle_label: bit string which describes the location/basis
        !       function/etc of the particle (ie psip or excip).
        !    nbits: length (in bits) of particle_label.  This allows us to ignore
        !       any additional padding at the end of the bit string for different
        !       sizes of i0 integers.
        !    seed: seed to pass to the hashing function.
        !    shift: value to add to the hash of the label before determining
        !       the processor to which the label is assigned.
        !    freq: frequency over which the result changes exactly once.
        !       See comments below.  Ignored if the shift is 0.  Must be smaller
        !       than 32.
        !    np: number of processors over which the particles are to be
        !       distributed.
        !    proc_map: array which maps determinants to processors.
        !    nslots: number of slots proc_map is divided into.
        ! Out:
        !    particle_proc: processor where determinant resides
        !    slot_pos: position in proc_map for this determinant

        use parallel, only: parent

        integer(i0), intent(in) :: particle_label(:)
        integer, intent(in) :: nbits, seed, shift, freq, np
        integer, intent(in) :: proc_map(0:)
        integer, intent(in) :: nslots
        integer, intent(out) :: particle_proc, slot_pos

        integer(i0) :: particle_label_padded(size(particle_label)+1)
        integer :: label_len

        if (i0_length == 32) then
            label_len = size(particle_label)/2
            if (mod(label_len,2) == 0) then
                call assign_particle_processor(particle_label, nbits, seed, shift, freq, np, &
                                               particle_proc, slot_pos, proc_map, nslots)
            else
                particle_label_padded(:label_len) = particle_label(:label_len)
                particle_label_padded(label_len+1) = 0_i0
                particle_label_padded(label_len+2:) = particle_label(label_len+1:)
                call assign_particle_processor(particle_label_padded, nbits+i0_length, seed, shift, freq, np, &
                                               particle_proc, slot_pos, proc_map, nslots)
            end if
        else
            call assign_particle_processor(particle_label, nbits, seed, shift, freq, np, &
                                           particle_proc, slot_pos, proc_map, nslots)
        end if

    end subroutine assign_particle_processor_dmqmc

    subroutine add_spawned_particle(f_new, nspawn, particle_type, iproc_spawn, spawn)

        ! Add a new particle to a store of spawned particles.

        ! In:
        !    f_new:  determinant on which to spawn.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        !    iproc_spawn: processor to which f_new belongs (see assign_particle_processor).
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use parallel, only: nthreads
        use spawn_data, only: spawn_t
        use omp_lib

        integer(i0), intent(in) :: f_new(:)
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type, iproc_spawn
        type(spawn_t), intent(inout) :: spawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        ! Move to the next position in the spawning array.
        spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

        ! Zero it as not all fields are set.
        spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0_int_s

        ! Set info in spawning array.
        spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = int(f_new, int_s)
        spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = int(nspawn, int_s)

    end subroutine add_spawned_particle

    subroutine add_flagged_spawned_particle(f_new, nspawn, particle_type, flag, iproc_spawn, spawn)

        ! Add a new particle to a store of spawned particles with the flag field
        ! set.

        ! In:
        !    f_new:  determinant on which to spawn.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        !    flag: flag value of the determinant/particle to set in the spawn store.
        !    iproc_spawn: processor to which f_new belongs (see assign_particle_processor).
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use parallel, only: nthreads
        use spawn_data, only: spawn_t
        use omp_lib

        integer(i0), intent(in) :: f_new(:)
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type, flag, iproc_spawn
        type(spawn_t), intent(inout) :: spawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        ! Move to the next position in the spawning array.
        spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

        ! Zero it as not all fields are set.
        spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0_int_s

        ! Set info in spawning array.
        spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = int(f_new, int_s)
        spawn%sdata(spawn%bit_str_len+particle_type,spawn%head(thread_id,iproc_spawn)) = int(nspawn, int_s)
        spawn%sdata(spawn%flag_indx,spawn%head(thread_id,iproc_spawn)) = int(flag, int_s)

    end subroutine add_flagged_spawned_particle

    subroutine add_spawned_particles(f_new, nspawn, iproc_spawn, spawn)

        ! Add a set of particles to a store of spawned particles.

        ! In:
        !    f_new:  determinant on which to spawn.
        !    nspawn: the (signed) number of particles of each particle type to
        !       create on the spawned determinant.
        !    iproc_spawn: processor to which f_new belongs (see assign_particle_processor).
        ! In/Out:
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use parallel, only: nthreads
        use spawn_data, only: spawn_t
        use omp_lib

        integer(i0), intent(in) :: f_new(:)
        integer(int_p), intent(in) :: nspawn(:) ! (spawn%ntypes)
        integer, intent(in) :: iproc_spawn
        type(spawn_t), intent(inout) :: spawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        ! Move to the next position in the spawning array.
        spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

        ! Zero it as not all fields are set.
        spawn%sdata(:,spawn%head(thread_id,iproc_spawn)) = 0_int_s

        ! Set info in spawning array.
        spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,iproc_spawn)) = int(f_new, int_s)
        spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes,spawn%head(thread_id,iproc_spawn)) = int(nspawn, int_s)

    end subroutine add_spawned_particles

    subroutine create_spawned_particle(basis, reference, cdet, connection, nspawn, particle_type, spawn, fexcit)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.  Ignored if fexcit is
        !        given.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        !    fexcit (optional): bit string representation of the determinant spawned onto.
        ! In/Out:
        !    spawn: spawn_t object to which the spawned particle will be added.

        use parallel, only: nprocs

        use basis_types, only: basis_t
        use determinants, only: det_info_t
        use excitations, only: excit_t, create_excited_det
        use spawn_data, only: spawn_t
        use qmc_data, only: reference_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type
        integer(i0), intent(in), target, optional :: fexcit(:)
        type(spawn_t), intent(inout) :: spawn

        integer(i0), target :: f_local(basis%string_len)
        integer(i0), pointer :: f_new(:)
        integer :: iproc_spawn, slot

        if (present(fexcit)) then
            f_new => fexcit
        else
            call create_excited_det(basis, cdet%f, connection, f_local)
            f_new => f_local
        end if

        call assign_particle_processor(f_new, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs, &
                                       iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)

        call add_spawned_particle(f_new, nspawn, particle_type, iproc_spawn, spawn)

    end subroutine create_spawned_particle

    subroutine create_spawned_particle_initiator(basis, reference, cdet, connection, nspawn, particle_type, spawn, fexcit)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.  Ignored if fexcit is
        !        given.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        !    fexcit (optional): bit string representation of the determinant spawned onto.
        ! In/Out:
        !    spawn: spawn_t object to which the spawned particle will be added.

        use parallel, only: nprocs

        use basis_types, only: basis_t
        use determinants, only: det_info_t
        use excitations, only: excit_t, create_excited_det
        use spawn_data, only: spawn_t
        use qmc_data, only: reference_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type
        integer(i0), intent(in), target, optional :: fexcit(:)
        type(spawn_t), intent(inout) :: spawn

        integer(i0), target :: f_local(basis%string_len)
        integer(i0), pointer :: f_new(:)
        integer :: iproc_spawn, slot

        if (present(fexcit)) then
            f_new => fexcit
        else
            call create_excited_det(basis, cdet%f, connection, f_local)
            f_new => f_local
        end if

        call assign_particle_processor(f_new, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs, &
                                       iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)

        call add_flagged_spawned_particle(f_new, nspawn, particle_type, cdet%initiator_flag, iproc_spawn, spawn)

    end subroutine create_spawned_particle_initiator

    subroutine create_spawned_particle_truncated(basis, reference, cdet, connection, nspawn, particle_type, spawn, fexcit)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant that excitation level is calculated from.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.  Ignored if fexcit is
        !        given.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        !    fexcit (optional): bit string representation of the determinant spawned onto.
        ! In/Out:
        !    spawn: spawn_t object to which the spawned particle will be added.

        use parallel, only: nprocs

        use basis_types, only: basis_t
        use determinants, only: det_info_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use spawn_data, only: spawn_t
        use qmc_data, only: reference_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type
        integer(i0), intent(in), target, optional :: fexcit(:)
        type(spawn_t), intent(inout) :: spawn

        integer(i0), target :: f_local(basis%string_len)
        integer(i0), pointer :: f_new(:)
        integer :: iproc_spawn, slot

        if (present(fexcit)) then
            f_new => fexcit
        else
            call create_excited_det(basis, cdet%f, connection, f_local)
            f_new => f_local
        end if

        ! Only accept spawning if it's within the truncation level.
        if (get_excitation_level(reference%hs_f0, f_new) <= reference%ex_level) then

            call assign_particle_processor(f_new, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs, &
                                           iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)

            call add_spawned_particle(f_new, nspawn, particle_type, iproc_spawn, spawn)

        end if

    end subroutine create_spawned_particle_truncated

    subroutine create_spawned_particle_initiator_truncated(basis, reference, cdet, connection, nspawn, particle_type, &
                                                           spawn, fexcit)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.  Ignored if fexcit is
        !        given.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        !    fexcit (optional): bit string representation of the determinant spawned onto.
        ! In/Out:
        !    spawn: spawn_t object to which the spawned particle will be added.

        use parallel, only: nprocs

        use basis_types, only: basis_t
        use determinants, only: det_info_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use spawn_data, only: spawn_t
        use qmc_data, only: reference_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type
        integer(i0), intent(in), target, optional :: fexcit(:)
        type(spawn_t), intent(inout) :: spawn

        integer(i0), target :: f_local(basis%string_len)
        integer(i0), pointer :: f_new(:)
        integer :: iproc_spawn, slot

        if (present(fexcit)) then
            f_new => fexcit
        else
            call create_excited_det(basis, cdet%f, connection, f_local)
            f_new => f_local
        end if

        ! Only accept spawning if it's within the truncation level.
        if (get_excitation_level(reference%hs_f0, f_new) <= reference%ex_level) then

            call assign_particle_processor(f_new, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs, &
                                           iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)

            call add_flagged_spawned_particle(f_new, nspawn, particle_type, cdet%initiator_flag, iproc_spawn, spawn)

        end if

    end subroutine create_spawned_particle_initiator_truncated

    subroutine create_spawned_particle_ras(basis, reference, cdet, connection, nspawn, particle_type, spawn, fexcit)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.  Ignored if fexcit is
        !        given.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        !    fexcit (optional): bit string representation of the determinant spawned onto.
        ! In/Out:
        !    spawn: spawn_t object to which the spawned particle will be added.

        use parallel, only: nprocs

        use basis_types, only: basis_t
        use bit_utils, only: count_set_bits
        use calc, only: ras1, ras3, ras1_min, ras3_max
        use determinants, only: det_info_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level, in_ras
        use spawn_data, only: spawn_t
        use qmc_data, only: reference_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type
        integer(i0), intent(in), target, optional :: fexcit(:)
        type(spawn_t), intent(inout) :: spawn

        integer(i0), target :: f_local(basis%string_len)
        integer(i0), pointer :: f_new(:)
        integer :: iproc_spawn, slot

        if (present(fexcit)) then
            f_new => fexcit
        else
            call create_excited_det(basis, cdet%f, connection, f_local)
            f_new => f_local
        end if

        ! Only accept spawning if it's within the RAS space.
        if (in_ras(ras1, ras3, ras1_min, ras3_max, f_new)) then

            call assign_particle_processor(f_new, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, nprocs, &
                                           iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)

            call add_spawned_particle(f_new, nspawn, particle_type, iproc_spawn, spawn)

        end if

    end subroutine create_spawned_particle_ras

    subroutine create_spawned_particle_initiator_ras(basis, reference, cdet, connection, nspawn, particle_type, spawn, fexcit)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.  Ignored if fexcit is
        !        given.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        !    fexcit (optional): bit string representation of the determinant spawned onto.
        ! In/Out:
        !    spawn: spawn_t object to which the spawned particle will be added.

        use parallel, only: nprocs

        use basis_types, only: basis_t
        use bit_utils, only: count_set_bits
        use calc, only: ras1, ras3, ras1_min, ras3_max
        use determinants, only: det_info_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level, in_ras
        use spawn_data, only: spawn_t
        use qmc_data, only: reference_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: connection
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type
        integer(i0), intent(in), target, optional :: fexcit(:)
        type(spawn_t), intent(inout) :: spawn

        integer(i0), target :: f_local(basis%string_len)
        integer(i0), pointer :: f_new(:)
        integer :: iproc_spawn, slot

        if (present(fexcit)) then
            f_new => fexcit
        else
            call create_excited_det(basis, cdet%f, connection, f_local)
            f_new => f_local
        end if

        ! Only accept spawning if it's within the RAS space.
        if (in_ras(ras1, ras3, ras1_min, ras3_max, f_new)) then

            call assign_particle_processor_dmqmc(f_new, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, &
                                                 nprocs, iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)

            call add_flagged_spawned_particle(f_new, nspawn, particle_type, cdet%initiator_flag, iproc_spawn, spawn)

        end if

    end subroutine create_spawned_particle_initiator_ras

    subroutine create_spawned_particle_density_matrix(basis, reference, f1, f2, connection, nspawn, spawning_end, &
                                                      particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant defining the
        !         accessible region of the Hilbert space.
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
        !    spawn: spawn_t object to which the spawned particle will be added.

        use basis_types, only: basis_t
        use errors, only: stop_all
        use excitations, only: excit_t, create_excited_det
        use parallel, only: nprocs, nthreads
        use qmc_data, only: reference_t
        use spawn_data, only: spawn_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        integer(i0), intent(in) :: f1(basis%string_len), f2(basis%string_len)
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: spawning_end
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn
        type(excit_t), intent(in) :: connection

        integer(i0) :: f_new(basis%string_len)
        integer(i0) :: f_new_tot(basis%tensor_label_len)

        ! DMQMC is not yet OpenMP parallelised.
        integer, parameter :: thread_id = 0

        integer :: iproc_spawn, slot

        ! Create bit string of new determinant. The entire two-ended
        ! bitstring is eventually stored in f_new_tot.
        call create_excited_det(basis, f1, connection, f_new)

        f_new_tot = 0_i0
        if (spawning_end==1) then
            f_new_tot(:basis%string_len) = f_new
            f_new_tot((basis%string_len+1):(basis%tensor_label_len)) = f2
        else
            f_new_tot(:basis%string_len) = f2
            f_new_tot((basis%string_len+1):(basis%tensor_label_len)) = f_new
        end if

        call assign_particle_processor_dmqmc(f_new_tot, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, &
                                             nprocs, iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)

        if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(nthreads-1,iproc_spawn) >= spawn%block_size) &
            call stop_all('create_spawned_particle_density_matrix',&
                           'There is no space left in the spawning array.')

        call add_spawned_particle(f_new_tot, nspawn, particle_type, iproc_spawn, spawn)

    end subroutine create_spawned_particle_density_matrix

    subroutine create_spawned_particle_half_density_matrix(basis, reference, f1, f2, connection, nspawn, spawning_end, &
                                                           particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! If walker tries to spawn in the lower triangle of density matrix
        ! then reflect it to upper triangle by swapping bit strings for the
        ! two ends. The current position in the spawning array is updated.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant defining the
        !         accessible region of the Hilbert space.
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
        !    spawn: spawn_t object to which the spawned particle will be added.

        use bit_utils, only: bit_str_cmp
        use basis_types, only: basis_t
        use errors, only: stop_all
        use excitations, only: excit_t, create_excited_det
        use parallel, only: nprocs, nthreads
        use qmc_data, only: reference_t
        use spawn_data, only: spawn_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        integer(i0), intent(in) :: f1(basis%string_len), f2(basis%string_len)
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: spawning_end
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn
        type(excit_t), intent(in) :: connection

        integer(i0) :: f_new(basis%string_len)
        integer(i0) :: f_new_tot(basis%tensor_label_len)

        ! DMQMC is not yet OpenMP parallelised.
        integer, parameter :: thread_id = 0

        integer :: iproc_spawn, slot

        ! Create bit string of new determinant. The entire two-ended
        ! bitstring is eventually stored in f_new_tot.
        call create_excited_det(basis, f1, connection, f_new)
        f_new_tot = 0_i0

        ! Test to see whether the new determinant resides in the upper
        ! triangle of the density matrix. If so keep bit string ends
        ! as they are. If not then swap bitstring ends so that the
        ! new psips are reflected into the upper triangle of the density
        ! matrix as they try to spawn.
        if (bit_str_cmp(f_new, f2) == -1) then
            f_new_tot(:basis%string_len) = f_new
            f_new_tot((basis%string_len+1):(basis%tensor_label_len)) = f2
        else
            f_new_tot(:basis%string_len) = f2
            f_new_tot((basis%string_len+1):(basis%tensor_label_len)) = f_new
        end if

        call assign_particle_processor_dmqmc(f_new_tot, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, &
                                             nprocs, iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)

        if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(nthreads-1,iproc_spawn) >= spawn%block_size) &
            call stop_all('create_spawned_particle_half_density_matrix',&
                           'There is no space left in the spawning array.')

        call add_spawned_particle(f_new_tot, nspawn, particle_type, iproc_spawn, spawn)

    end subroutine create_spawned_particle_half_density_matrix

    subroutine create_spawned_particle_truncated_half_density_matrix(basis, reference, f1, f2, connection, nspawn, spawning_end, &
                                                                     particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! A spawned walker is only created on (f1', f2) if f1' and f2 do not differ by
        ! more than reference%ex_level basis functions, where f1' is obtained by
        ! applying the connection to f1.

        ! Note: This is the half density matrix version of create_spawned_particle_truncated_density_matrix
        ! It works in an ientical way apart from attempts to spawn on the lower triangle of the
        ! density matrix are reflected so that they are spawned on the upper triangle of the
        ! density matrix.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant defining the
        !         accessible region of the Hilbert space.
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
        !    spawn: spawn_t object to which the spawned particle will be added.

        use bit_utils, only: bit_str_cmp
        use basis_types, only: basis_t
        use errors, only: stop_all
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use parallel, only: nprocs, nthreads
        use qmc_data, only: reference_t
        use spawn_data, only: spawn_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        integer(i0), intent(in) :: f1(basis%string_len), f2(basis%string_len)
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: spawning_end
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn
        type(excit_t), intent(in) :: connection

        integer(i0) :: f_new(basis%string_len)
        integer(i0) :: f_new_tot(basis%tensor_label_len)

        ! DMQMC is not yet OpenMP parallelised.
        integer, parameter :: thread_id = 0

        integer :: iproc_spawn, slot

        ! Create bit string of new determinant. The entire two-ended
        ! bitstring is eventually stored in f_new_tot.
        call create_excited_det(basis, f1, connection, f_new)

        if (get_excitation_level(f2, f_new) <= reference%ex_level) then

            f_new_tot = 0_i0
            ! Test to see whether the new determinant resides in the upper
            ! triangle of the density matrix. If so keep bit string ends
            ! as they are. If not then swap bitstring ends so that the
            ! new psips are reflected into the upper triangle of the density
            ! matrix as they try to spawn.
            if (bit_str_cmp(f_new, f2) == -1) then
                f_new_tot(:basis%string_len) = f_new
                f_new_tot((basis%string_len+1):(basis%tensor_label_len)) = f2
            else
                f_new_tot(:basis%string_len) = f2
                f_new_tot((basis%string_len+1):(basis%tensor_label_len)) = f_new
            end if

            call assign_particle_processor_dmqmc(f_new_tot, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, &
                                                 spawn%move_freq, nprocs, iproc_spawn, slot, spawn%proc_map%map, &
                                                 spawn%proc_map%nslots)

            if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(nthreads-1,iproc_spawn) >= spawn%block_size) &
                call stop_all('create_spawned_particle_truncated_half_density_matrix',&
                               'There is no space left in the spawning array.')

            call add_spawned_particle(f_new_tot, nspawn, particle_type, iproc_spawn, spawn)

        end if

    end subroutine create_spawned_particle_truncated_half_density_matrix

    subroutine create_spawned_particle_truncated_density_matrix(basis, reference, f1, f2, connection, nspawn, spawning_end, &
                                                                particle_type, spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! A spawned walker is only created on (f1', f2) if f1' and f2 do not
        ! differ by more than reference%ex_level basis functions, where f1' is
        ! obtained by applying the connection to f1.

        ! In:
        !    basis: information about the single-particle basis.
        !    reference: current reference determinant defining the
        !         accessible region of the Hilbert space.
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
        !    spawn: spawn_t object to which the spawned particle will be added.

        use basis_types, only: basis_t
        use errors, only: stop_all
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use parallel, only: nprocs, nthreads
        use qmc_data, only: reference_t
        use spawn_data, only: spawn_t

        type(basis_t), intent(in) :: basis
        type(reference_t), intent(in) :: reference
        integer(i0), intent(in) :: f1(basis%string_len), f2(basis%string_len)
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: spawning_end
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn
        type(excit_t), intent(in) :: connection


        integer(i0) :: f_new(basis%string_len)
        integer(i0) :: f_new_tot(basis%tensor_label_len)

        ! DMQMC is not yet OpenMP parallelised.
        integer, parameter :: thread_id = 0

        integer :: iproc_spawn, slot

        ! Create bit string of new determinant. The entire two-ended
        ! bitstring is eventually stored in f_new_tot.
        call create_excited_det(basis, f1, connection, f_new)

        if (get_excitation_level(f2, f_new) <= reference%ex_level) then

            f_new_tot = 0_i0
            if (spawning_end==1) then
                f_new_tot(:basis%string_len) = f_new
                f_new_tot((basis%string_len+1):(basis%tensor_label_len)) = f2
            else
                f_new_tot(:basis%string_len) = f2
                f_new_tot((basis%string_len+1):(basis%tensor_label_len)) = f_new
            end if

            call assign_particle_processor_dmqmc(f_new_tot, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, &
                                                 spawn%move_freq, nprocs, iproc_spawn, slot, spawn%proc_map%map, &
                                                 spawn%proc_map%nslots)

            if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(nthreads-1,iproc_spawn) >= spawn%block_size) &
                call stop_all('create_spawned_particle_truncated_density_matrix', &
                               'There is no space left in the spawning array.')

            call add_spawned_particle(f_new_tot, nspawn, particle_type, iproc_spawn, spawn)

        end if

    end subroutine create_spawned_particle_truncated_density_matrix

    subroutine create_spawned_particle_rdm(rdm, nspawn_in, particle_type, rdm_spawn)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    rdm: the rdm that this psip contributes to.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the index of particle type to be created.
        ! In/Out:
        !    rdm_spawn: rdm_spawn_t object to which the spanwed particle
        !        will be added.

        use bit_utils, only: operator(.bitstrgt.)
        use dmqmc_data, only: rdm_t, rdm_spawn_t
        use errors, only: stop_all
        use parallel, only: iproc, nprocs, nthreads
        use hash_table, only: hash_table_pos_t, lookup_hash_table_entry
        use hash_table, only: assign_hash_table_entry
        use utils, only: int_fmt

        type(rdm_t), intent(in) :: rdm
        integer(int_p), intent(in) :: nspawn_in
        integer, intent(in) :: particle_type
        type(rdm_spawn_t), intent(inout) :: rdm_spawn
        integer(int_p) :: nspawn

        integer(i0) :: f_new_tot(2*rdm%string_len)

        integer :: iproc_spawn, slot
        ! WARNING!  The below algorithm is *not* suitable for conversion to
        ! thread-safety as each thread could be spawning onto the same RDM
        ! element, yet the hash table requires a given element to exist in one
        ! (and only one) location, which is not compatible with how the threaded
        ! spawning arrays work.
        integer, parameter :: thread_id = 0

        type(hash_table_pos_t) :: pos
        logical :: hit
        integer :: err_code

        associate(f1=>rdm%end1, f2=>rdm%end2, rdm_bl=>rdm%string_len, spawn=>rdm_spawn%spawn, &
                  ht=>rdm_spawn%ht, bsl=>rdm_spawn%spawn%bit_str_len)

            nspawn = nspawn_in
            f_new_tot = 0_i0

            ! Symmetry is enforced on the RDM in the following.
            if (f1 .bitstrgt. f2) then
                ! If below the diagonal, swap the bitstrings so that the spawning occurs above it.
                f_new_tot(:rdm_bl) = f2
                f_new_tot(rdm_bl+1:2*rdm_bl) = f1
            else
                f_new_tot(:rdm_bl) = f1
                f_new_tot(rdm_bl+1:2*rdm_bl) = f2
                if (all(f1 == f2)) then
                    ! Because off-diagonal elements have been doubled (elements above the diagonal taking
                    ! contributions from both below and above it), we must double the diagonal elements too.
                    nspawn = nspawn*2
                end if
            end if

            ! Can just hash everything as RDMs don't enter the Markov chain.
            call assign_particle_processor(f_new_tot, 2*i0_length*rdm_bl, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, &
                                          nprocs, iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)

            call lookup_hash_table_entry(ht, f_new_tot, pos, hit)

            if (hit) then
                associate(indx => ht%table(pos%ientry,pos%islot))
                    ! Direct annihilation/cancellation in spawning array.
                    spawn%sdata(bsl+particle_type,indx) = spawn%sdata(bsl+particle_type,indx) + int(nspawn, int_s)
                end associate
            else
                ! Move to the next position in the spawning array.
                call assign_hash_table_entry(ht, pos%islot, pos, err_code)
                if (err_code /= 0) then
                    write(6,'(1X,a9,'//int_fmt(err_code,1)//')') 'err_code:', err_code
                    call stop_all('create_spawned_particle_rdm','Error in assigning hash &
                                  &table entry.')
                end if
                ! Fix hash table to point to the head of the spawn data for this thread/processor.
                spawn%head(thread_id,iproc_spawn) = spawn%head(thread_id,iproc_spawn) + nthreads

                if (spawn%head(thread_id,iproc_spawn) - spawn%head_start(nthreads-1,iproc_spawn) >= spawn%block_size) &
                    call stop_all('create_spawned_particle_rdm','There is no space left in the RDM array.')

                ht%table(pos%ientry,pos%islot) = spawn%head(thread_id,iproc_spawn)
                associate(indx => ht%table(pos%ientry,pos%islot))
                    ! Set info in spawning array.
                    ! Zero it as not all fields are set.
                    spawn%sdata(:,indx) = 0_int_s
                    spawn%sdata(:bsl,indx) = int(f_new_tot, int_s)
                    spawn%sdata(bsl+particle_type,indx) = int(nspawn, int_s)
                end associate
            end if

        end associate

    end subroutine create_spawned_particle_rdm

end module spawning
