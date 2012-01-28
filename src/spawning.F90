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

    subroutine spawn(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant.

        ! This is just a thin wrapper around a system-specific excitation
        ! generator and a utility function.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use excitations, only: excit
        use proc_pointers, only: gen_excit_ptr

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel

        ! 1. Generate random excitation.
        call gen_excit_ptr(cdet, pgen, connection, hmatel)

        ! 2. Attempt spawning.
        nspawn = attempt_to_spawn(hmatel, pgen, parent_sign)

    end subroutine spawn

    subroutine spawn_importance_sampling(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant.

        ! This subroutine applies a transformation to the Hamiltonian to achieve
        ! importance sampling of the stochastic process.

        ! This is just a thin wrapper around a system-specific excitation
        ! generator, trial-function specific transformation routine and
        ! a utility function.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use excitations, only: excit
        use proc_pointers, only: gen_excit_ptr, trial_fn_ptr

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel

        ! 1. Generate random excitation.
        call gen_excit_ptr(cdet, pgen, connection, hmatel)

        ! 2. Transform Hamiltonian matrix element by trial function.
        call trial_fn_ptr(cdet, connection, hmatel)

        ! 3. Attempt spawning.
        nspawn = attempt_to_spawn(hmatel, pgen, parent_sign)

    end subroutine spawn_importance_sampling

    subroutine spawn_lattice_split_gen(cdet, parent_sign, nspawn, connection)

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

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use excitations, only: excit
        use fciqmc_data, only: tau
        use proc_pointers, only: gen_excit_init_ptr, gen_excit_finalise_ptr

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, abs_hmatel, hmatel

        ! 1. Generate enough of a random excitation to determinant the
        ! generation probability and |H_ij|.
        call gen_excit_init_ptr(cdet, pgen, connection, abs_hmatel)

        ! 2. Attempt spawning.
        nspawn = nspawn_from_prob(tau*abs_hmatel/pgen)

        if (nspawn /= 0) then

            ! 3. Complete excitation and find sign of connecting matrix element.
            call gen_excit_finalise_ptr(cdet, connection, hmatel)

            ! 4. Find sign of offspring.
            call set_child_sign(hmatel, parent_sign, nspawn)

        end if

    end subroutine spawn_lattice_split_gen

    subroutine spawn_lattice_split_gen_importance_sampling(cdet, parent_sign, nspawn, connection)

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

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use excitations, only: excit
        use fciqmc_data, only: tau
        use proc_pointers, only: gen_excit_init_ptr, gen_excit_finalise_ptr, gen_excit_ptr, trial_fn_ptr

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, tilde_hmatel, hmatel

        ! 1. Generate enough of a random excitation to determinant the
        ! generation probability and |H_ij|.
        call gen_excit_init_ptr(cdet, pgen, connection, tilde_hmatel)

        ! 2. Transform Hamiltonian matrix element by trial function.
        call trial_fn_ptr(cdet, connection, hmatel)

        ! 3. Attempt spawning.
        nspawn = nspawn_from_prob(tau*abs(tilde_hmatel)/pgen)

        if (nspawn /= 0) then

            ! 4. Complete excitation and find sign of connecting matrix element.
            call gen_excit_finalise_ptr(cdet, connection, hmatel)

            ! 5. Find sign of offspring.
            ! Note that we don't care about the value of H_ij at this step, only
            ! the sign.
            call set_child_sign(tilde_hmatel*hmatel, parent_sign, nspawn)

        end if

    end subroutine spawn_lattice_split_gen_importance_sampling

!--- Attempt spawning based upon random excitation ---

    function nspawn_from_prob(probability) result(number_spawned)

        ! Generate the number spawned from a probability. If probability is greater than
        ! zero, then number spawned = int(probability) + stochastic{0,1}
        ! where the latter half of the RHS is a stochastic spawning from the remainder
        !
        ! In:
        !    probability: the spawning probability
        !
        ! Returns:
        !    number_spawned: the number spawned from this probability

        use dSFMT_interface , only: genrand_real2

        implicit none
        real(p), intent(in) :: probability
        integer             :: number_spawned
        real(p)             :: psuccess, pstochastic

        ! Generate random number
        psuccess = genrand_real2()

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

    function attempt_to_spawn(hmatel, pgen, parent_sign) result(nspawn)

        ! In:
        !    hmatel: Hamiltonian matrix element connecting a determinant and an
        !    excitation of that determinant.
        !    pgen: probability of generating the excitation.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Returns:
        !    number of particles spawned.  0 indicates the spawning attempt was
        !    unsuccessful.

        use dSFMT_interface, only:  genrand_real2
        use fciqmc_data, only: tau

        integer :: nspawn

        real(p), intent(in) :: pgen, hmatel
        integer, intent(in) :: parent_sign

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

        if (pspawn > genrand_real2()) nspawn = nspawn + 1

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

    subroutine create_spawned_particle(cdet, connection, nspawn, particle_type)

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

        use hashing
        use parallel, only: iproc, nprocs

        use basis, only: basis_length
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det
        use fciqmc_data, only: spawned_walkers, spawning_head, spawned_pop

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type

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
        iproc_spawn = modulo(murmurhash_bit_string(f_new, basis_length), nprocs)
#endif

        ! Move to the next position in the spawning array.
        spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawned_walkers(:,spawning_head(iproc_spawn)) = 0
        spawned_walkers(:basis_length,spawning_head(iproc_spawn)) = f_new
        spawned_walkers(particle_type,spawning_head(iproc_spawn)) = nspawn

    end subroutine create_spawned_particle

    subroutine create_spawned_particle_initiator(cdet, connection, nspawn, particle_type)

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

        use hashing
        use parallel, only: iproc, nprocs

        use basis, only: basis_length
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det
        use fciqmc_data, only: spawned_walkers, spawning_head, spawned_pop, spawned_parent

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type

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
        iproc_spawn = modulo(murmurhash_bit_string(f_new, basis_length), nprocs)
#endif

        ! Move to the next position in the spawning array.
        spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawned_walkers(:,spawning_head(iproc_spawn)) = 0
        spawned_walkers(:basis_length,spawning_head(iproc_spawn)) = f_new
        spawned_walkers(particle_type,spawning_head(iproc_spawn)) = nspawn
        ! initiator_flag: flag indicating the staturs of the parent determinant.
        !     initiator_flag = 0 indicates the parent is an initiator.
        !     initiator_flag = 1 indicates the parent is not an initiator.
        spawned_walkers(spawned_parent,spawning_head(iproc_spawn)) = cdet%initiator_flag

    end subroutine create_spawned_particle_initiator

    subroutine create_spawned_particle_truncated(cdet, connection, nspawn, particle_type)

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
        !    particle_type: the type of particle created.  Must correspond to
        !        the desired element in the spawning array (i.e. be spawned_pop
        !        for Hamiltonian particles and spawned_hf_pop for
        !        Hellmann--Feynman particles).

        use hashing
        use parallel, only: iproc, nprocs

        use basis, only: basis_length
        use calc, only: truncation_level
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det, get_excitation_level
        use fciqmc_data, only: spawned_walkers, spawning_head, spawned_pop, f0

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type

        integer(i0) :: f_new(basis_length)
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

        ! Only accept spawning if it's within the truncation level.
        if (get_excitation_level(f0, f_new) <= truncation_level) then

#ifdef PARALLEL
            ! (Extra credit for parallel calculations)
            ! Need to determine which processor the spawned walker should be sent
            ! to.  This communication is done during the annihilation process, after
            ! all spawning and death has occured..
            iproc_spawn = modulo(murmurhash_bit_string(f_new, basis_length), nprocs)
#endif

            ! Move to the next position in the spawning array.
            spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

            ! Set info in spawning array.
            ! Zero it as not all fields are set.
            spawned_walkers(:,spawning_head(iproc_spawn)) = 0
            spawned_walkers(:basis_length,spawning_head(iproc_spawn)) = f_new
            spawned_walkers(particle_type,spawning_head(iproc_spawn)) = nspawn

        end if

    end subroutine create_spawned_particle_truncated

    subroutine create_spawned_particle_initiator_truncated(cdet, connection, nspawn, particle_type)

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
        !    particle_type: the type of particle created.  Must correspond to
        !        the desired element in the spawning array (i.e. be spawned_pop
        !        for Hamiltonian particles and spawned_hf_pop for
        !        Hellmann--Feynman particles).

        use hashing
        use parallel, only: iproc, nprocs

        use basis, only: basis_length
        use calc, only: truncation_level
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det, get_excitation_level
        use fciqmc_data, only: spawned_walkers, spawning_head, spawned_pop, spawned_parent, f0

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type

        integer(i0) :: f_new(basis_length)
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

        ! Only accept spawning if it's within the truncation level.
        if (get_excitation_level(f0, f_new) <= truncation_level) then

#ifdef PARALLEL
            ! (Extra credit for parallel calculations)
            ! Need to determine which processor the spawned walker should be sent
            ! to.  This communication is done during the annihilation process, after
            ! all spawning and death has occured..
            iproc_spawn = modulo(murmurhash_bit_string(f_new, basis_length), nprocs)
#endif

            ! Move to the next position in the spawning array.
            spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

            ! Set info in spawning array.
            ! Zero it as not all fields are set.
            spawned_walkers(:,spawning_head(iproc_spawn)) = 0
            spawned_walkers(:basis_length,spawning_head(iproc_spawn)) = f_new
            spawned_walkers(particle_type,spawning_head(iproc_spawn)) = nspawn
            ! initiator_flag: flag indicating the staturs of the parent determinant.
            !     initiator_flag = 0 indicates the parent is an initiator.
            !     initiator_flag = 1 indicates the parent is not an initiator.
            spawned_walkers(spawned_parent,spawning_head(iproc_spawn)) = cdet%initiator_flag

        end if

    end subroutine create_spawned_particle_initiator_truncated

    subroutine create_spawned_particle_density_matrix(f1, f2, connection, nspawn, spawning_end)

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

        ! Note that particle type, ie specifying hellman-feynman properties, is
        ! not currently considered, as the density matrix formulation (hopefully!)
        ! won't require it.

        use basis, only: basis_length, total_basis_length
        use excitations, only: excit, create_excited_det
        use fciqmc_data, only: spawned_walkers, spawning_head
        use fciqmc_data, only: spawned_parent, spawned_pop
        use hashing
        use parallel, only: iproc, nprocs

        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer, intent(in) :: nspawn
        integer, intent(in) :: spawning_end

        type(excit), intent(in) :: connection
        integer(i0) :: f_new(basis_length)
        integer(i0) :: f_new_tot(total_basis_length)

#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

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

#ifdef PARALLEL
        ! (Extra credit for parallel calculations)
        ! Need to determine which processor the spawned walker should be sent
        ! to.  This communication is done during the annihilation process, after
        ! all spawning and death has occured..
        iproc_spawn = modulo(murmurhash_bit_string(f_new_tot, (total_basis_length)), nprocs)
#endif

        ! Move to the next position in the spawning array.
        spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawned_walkers(:,spawning_head(iproc_spawn)) = 0
        spawned_walkers(:(total_basis_length),spawning_head(iproc_spawn)) = f_new_tot
        spawned_walkers((total_basis_length)+1,spawning_head(iproc_spawn)) = nspawn

    end subroutine create_spawned_particle_density_matrix

    subroutine create_spawned_particle_truncated_density_matrix(f1, f2, connection, nspawn, spawning_end)

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

        ! Note that particle type, ie specifying hellman-feynman properties, is
        ! not currently considered, as the density matrix formulation (hopefully!)
        ! won't require it.

        use basis, only: basis_length, total_basis_length
        use calc, only: truncation_level
        use excitations, only: excit, create_excited_det, get_excitation_level
        use fciqmc_data, only: spawned_walkers, spawning_head
        use fciqmc_data, only: spawned_parent, spawned_pop
        use hashing
        use parallel, only: iproc, nprocs

        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer, intent(in) :: nspawn
        integer, intent(in) :: spawning_end

        type(excit), intent(in) :: connection
        integer(i0) :: f_new(basis_length)
        integer(i0) :: f_new_tot(total_basis_length)

#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

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

#ifdef PARALLEL
            ! (Extra credit for parallel calculations)
            ! Need to determine which processor the spawned walker should be sent
            ! to.  This communication is done during the annihilation process, after
            ! all spawning and death has occured..
            iproc_spawn = modulo(murmurhash_bit_string(f_new_tot, (total_basis_length)), nprocs)
#endif

            ! Move to the next position in the spawning array.
            spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

            ! Set info in spawning array.
            ! Zero it as not all fields are set.
            spawned_walkers(:,spawning_head(iproc_spawn)) = 0
            spawned_walkers(:(total_basis_length),spawning_head(iproc_spawn)) = f_new_tot
            spawned_walkers((total_basis_length)+1,spawning_head(iproc_spawn)) = nspawn

        end if

    end subroutine create_spawned_particle_truncated_density_matrix

end module spawning
