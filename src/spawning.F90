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

use const
implicit none

contains

    subroutine spawn_hub_k(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant for the 
        ! momentum space formulation of the Hubbard model.
        !
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
        use dSFMT_interface, only:  genrand_real2
        use excitations, only: calc_pgen_hub_k, excit
        use hamiltonian, only: slater_condon2_hub_k_excit
        use fciqmc_data, only: tau
        use system, only: hub_k_coulomb

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, psuccess, pspawn, hmatel
        integer :: i, j, a, b, ij_sym

        ! Single excitations are not connected determinants within the 
        ! momentum space formulation of the Hubbard model.

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(cdet%occ_list_alpha, cdet%occ_list_beta, i ,j, ij_sym)

        ! 2. Calculate the generation probability of the excitation.
        ! For two-band systems this depends only upon the orbitals excited from.
        pgen = calc_pgen_hub_k(ij_sym, cdet%f, cdet%unocc_list_alpha, cdet%unocc_list_beta)

        ! The hubbard model in momentum space is a special case. Connected
        ! non-identical determinants have the following properties:
        !     a) They differ by two spin-orbitals.
        !     b) In the (i,j)->(a,b) connecting excitation, the spins of i and
        !     j have to be opposite.  This is because
        !     < ij | ab > = U/N_k \delta_{s_i,s_a} \delta_{s_j,s_b} 
        !     and so < ij || ab > = 0 if s_i = s_a = s_j = s_b.
        !     In fact:
        !     < ij || ab > = 0          if s_i = s_a = s_j = s_b
        !                    U/N        if s_i = s_a & s_j = s_b & s_i /= s_b
        !                   -U/N        if s_i = s_b & s_j = s_a & s_i /= s_a
        ! The FCIQMC method allows us to only generate connected excitations, so
        ! we can actually test whether we accept the excitation before we finish
        ! completing the excitation.

        ! 3. Test that whether the attempted spawning is successful.
        ! As we can only excite from (alpha,beta) or (beta,alpha),
        !   H_ij =  < ij | ab >  *or*
        !        = -< ij | ba >
        ! so
        !   |H_ij| = U/\Omega
        ! for allowed excitations.
        pspawn = tau*abs(hub_k_coulomb)/pgen
        psuccess = genrand_real2()

        ! Need to take into account the possibilty of a spawning attempt
        ! producing multiple offspring...
        ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and 
        ! then spawn a particle with probability pspawn-floor(pspawn).
        nspawn = int(pspawn)
        pspawn = pspawn - nspawn

        if (pspawn > psuccess) nspawn = nspawn + 1

        if (nspawn > 0) then
            ! 4. Well, I suppose we should find out which determinant we're spawning
            ! on...
            call choose_ab_hub_k(cdet%f, cdet%unocc_list_alpha, ij_sym, a, b)

            connection%nexcit = 2
            connection%from_orb(1:2) = (/ i,j /)
            connection%to_orb(1:2) = (/ a,b /)

            ! 5. Is connecting matrix element positive (in which case we spawn with
            ! negative walkers) or negative (in which case we spawn with positive
            ! walkers)?
            call slater_condon2_hub_k_excit(cdet%f, connection, hmatel)

            ! If H_ij is positive, then the spawned walker is of opposite sign
            ! to the parent.
            ! If H_ij is negative, then the spawned walker is of the same sign
            ! as the parent.
            if (hmatel > 0) then
                nspawn = -sign(nspawn, parent_sign)
            else
                nspawn = sign(nspawn, parent_sign)
            end if

        end if
        
    end subroutine spawn_hub_k

    subroutine spawn_hub_k_no_renorm(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant for the 
        ! momentum space formulation of the Hubbard model.
        !
        ! This doesn't use excitation generators which exclude the case where,
        ! having selected 2 occupied orbitals (i and j) and the first virtual
        ! orbital (a), the final orbital is already occupied and so that
        ! excitation is impossible.  Whilst it is somewhat wasteful (generating
        ! excitations which can't be performed), there is a balance between the
        ! cost of generating forbidden excitations and the O(N) cost of
        ! renormalising the generation probabilities.
        !
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
        use dSFMT_interface, only:  genrand_real2
        use excitations, only: excit
        use hamiltonian, only: slater_condon2_hub_k_excit
        use fciqmc_data, only: tau
        use system, only: hub_k_coulomb, nalpha, nbeta, nvirt

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign

        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, psuccess, pspawn, hmatel
        integer :: ij_sym
        logical :: allowed_excitation

        ! Single excitations are not connected determinants within the 
        ! momentum space formulation of the Hubbard model.

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(cdet%occ_list_alpha, cdet%occ_list_beta, connection%to_orb(1), connection%to_orb(2), ij_sym)

        ! 2. Chose a random pair of spin orbitals to excite to.
        call find_ab_hub_k(cdet%f, cdet%unocc_list_alpha, ij_sym, connection%to_orb(1), connection%to_orb(2), allowed_excitation)

        if (allowed_excitation) then

            ! 3. probability of this excitation.
            ! pgen = p(i,j) [ p(a|i,j) p(b|i,j,a) + p(b|i,j) p(a|i,j,b) ]
            ! pgen = 1/(nalpha*nbeta) [ 1/(nbasis-nel) + 1/(basis-nel) ]
            !                    2
            !      =  -------------------------
            !         nalpha*nbeta*(nbasis-nel)
            pgen = 2.0_dp/(nalpha*nbeta*nvirt)

            ! 4. Test that whether the attempted spawning is successful.
            pspawn = tau*abs(hub_k_coulomb)/pgen
            psuccess = genrand_real2()
            nspawn = int(pspawn)
            pspawn = pspawn - nspawn
            if (pspawn > psuccess) nspawn = nspawn + 1

            if (nspawn > 0) then

                connection%nexcit = 2

                ! 5. Is connecting matrix element positive (in which case we spawn with
                ! negative walkers) or negative (in which case we spawn with positive
                ! walkers)?
                call slater_condon2_hub_k_excit(cdet%f, connection, hmatel)

                ! If H_ij is positive, then the spawned walker is of opposite sign
                ! to the parent.
                ! If H_ij is negative, then the spawned walker is of the same sign
                ! as the parent.
                if (hmatel > 0) then
                    nspawn = -sign(nspawn, parent_sign)
                else
                    nspawn = sign(nspawn, parent_sign)
                end if

            end if

        else

            ! Generated a forbidden excitation (b is already occupied)
            nspawn = 0

        end if

    end subroutine spawn_hub_k_no_renorm

    subroutine gen_excit_hub_k(cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate the probability of
        ! selecting that excitation.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use excitations, only: calc_pgen_hub_k, excit
        use fciqmc_data, only: tau
        use system, only: hub_k_coulomb
        use hamiltonian, only: slater_condon2_hub_k_excit

        type(det_info), intent(in) :: cdet
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        integer :: i, j, a, b, ij_sym

        ! See notes in spawn_hub_k for more details regarding this algorithm.
        ! spawn_hub_k contains a somewhat more optimised version for fciqmc and
        ! i-fciqmc, as you can decide whether to spawn (or not) without actually
        ! generating the excitation when working in momentum space.

        ! However, we need pgen and an excitation for use with, e.g. the folded
        ! spectrum algorithm.

        connection%nexcit = 2

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(cdet%occ_list_alpha, cdet%occ_list_beta, connection%from_orb(1), connection%from_orb(2), ij_sym)

        ! 2. Select a random pair of spin orbitals to excite to.
        call choose_ab_hub_k(cdet%f, cdet%unocc_list_alpha, ij_sym, connection%to_orb(1), connection%to_orb(2))

        ! 3. Calculate the generation probability of the excitation.
        ! For two-band systems this depends only upon the orbitals excited from.
        pgen = calc_pgen_hub_k(ij_sym, cdet%f, cdet%unocc_list_alpha, cdet%unocc_list_beta)


        ! 4. find the connecting matrix element.
        call slater_condon2_hub_k_excit(cdet%f, connection, hmatel)

    end subroutine gen_excit_hub_k

    subroutine spawn_hub_real(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant for the 
        ! real space formulation of the Hubbard model.
        !
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
        use dSFMT_interface, only:  genrand_real2
        use excitations, only: excit
        use fciqmc_data, only: tau

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: psuccess, pspawn, pgen, hmatel

        ! 1. Generate random excitation from cdet as well as both the probability 
        ! of spawning there and the connecting matrix element.
        call gen_excit_hub_real(cdet, pgen, connection, hmatel)

        ! 2. Calculate P_gen
        pspawn = tau*abs(hmatel)/pgen

        ! 3. Attempt spawning.
        psuccess = genrand_real2()

        ! Need to take into account the possibilty of a spawning attempt
        ! producing multiple offspring...
        ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and 
        ! then spawn a particle with probability pspawn-floor(pspawn).
        nspawn = int(pspawn)
        pspawn = pspawn - nspawn

        if (pspawn > psuccess) nspawn = nspawn + 1

        if (nspawn > 0) then

            ! 4. If H_ij is positive, then the spawned walker is of opposite
            ! sign to the parent, otherwise the spawned walkers if of the same
            ! sign as the parent.
            if (hmatel > 0.0_p) then
                nspawn = -sign(nspawn, parent_sign)
            else
                nspawn = sign(nspawn, parent_sign)
            end if

        end if

    end subroutine spawn_hub_real

    subroutine spawn_hub_real_no_renorm(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant for the 
        ! real space formulation of the Hubbard model.
        !
        ! This doesn't use excitation generators which exclude the case where,
        ! having selected an occupied orbital (i) there are no possible
        ! connected orbitals which are unoccupied which can be excited into,
        ! and so that excitation is impossible.  Whilst it is somewhat wasteful
        ! (generating excitations which can't be performed), there is a balance
        ! between the cost of generating forbidden excitations and the O(N) cost
        ! of renormalising the generation probabilities.
        !
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
        use dSFMT_interface, only:  genrand_real2
        use excitations, only: excit
        use fciqmc_data, only: tau
        use hamiltonian, only: slater_condon1_hub_real_excit
        use system, only: nel
        use hubbard_real, only: connected_sites
        use basis, only: bit_lookup

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        integer :: i, iel, ipos
        real(p) :: psuccess, pspawn, pgen, hmatel

        ! 1. Generate random excitation from cdet and probability of spawning
        ! there.
        ! Random selection of i.
        i = int(genrand_real2()*nel) + 1
        connection%from_orb(1) = cdet%occ_list(i)
        ! Select a at random from one of the connected orbitals.
        i = int(genrand_real2()*connected_sites(0,connection%from_orb(1)) + 1)
        connection%to_orb(1) = connected_sites(i,connection%from_orb(1))

        ipos = bit_lookup(1, connection%to_orb(1))
        iel = bit_lookup(2, connection%to_orb(1))

        if (btest(cdet%f(iel),ipos)) then

            ! a is occupied; forbidden excitation
            nspawn = 0

        else

            connection%nexcit=1

            ! 2. find the connecting matrix element.
            call slater_condon1_hub_real_excit(cdet%f, connection, hmatel)

            ! 3. Probability of spawning...
            ! For single excitations
            !   pgen = p(i) p(a|i)
            !        = 1/(nel*nconnected_sites)
            pgen = 1.0_dp/(nel*connected_sites(0,i))
            pspawn = tau*abs(hmatel)/pgen

            ! 4. Attempt spawning.
            psuccess = genrand_real2()

            ! Need to take into account the possibilty of a spawning attempt
            ! producing multiple offspring...
            ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and 
            ! then spawn a particle with probability pspawn-floor(pspawn).
            nspawn = int(pspawn)
            pspawn = pspawn - nspawn

            if (pspawn > psuccess) nspawn = nspawn + 1

            if (nspawn > 0) then

                ! 5. If H_ij is positive, then the spawned walker is of opposite
                ! sign to the parent, otherwise the spawned walkers if of the same
                ! sign as the parent.
                if (hmatel > 0.0_p) then
                    nspawn = -sign(nspawn, parent_sign)
                else
                    nspawn = sign(nspawn, parent_sign)
                end if

            end if

        end if

    end subroutine spawn_hub_real_no_renorm

    subroutine gen_excit_hub_real(cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate both the probability 
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D_i^a >, the Hamiltonian matrix element between a 
        !    determinant and a single excitation of it in the real space
        !    formulation of the Hubbard model.

        use determinants, only: det_info
        use excitations, only: calc_pgen_real, excit
        use hamiltonian, only: slater_condon1_hub_real_excit

        type(det_info), intent(in) :: cdet
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        integer :: nvirt_avail

        ! Double excitations are not connected determinants within the 
        ! real space formulation of the Hubbard model.
        connection%nexcit = 1

        ! 1. Chose a random connected excitation.
        call choose_ia_real(cdet%occ_list, cdet%f, connection%from_orb(1), connection%to_orb(1), nvirt_avail)

        ! 2. Find probability of generating this excited determinant.
        pgen = calc_pgen_real(cdet%occ_list, cdet%f, nvirt_avail)

        ! 3. find the connecting matrix element.
        call slater_condon1_hub_real_excit(cdet%f, connection, hmatel)

    end subroutine gen_excit_hub_real
    
    subroutine spawn_heisenberg(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new psip on a connected determinant for the 
        ! Heisenberg model
        !
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
        use dSFMT_interface, only:  genrand_real2
        use excitations, only: calc_pgen_real, excit
        use fciqmc_data, only: tau
        use hamiltonian, only: slater_condon1_hub_real_excit
        use system, only: J_coupling

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, psuccess, pspawn, hmatel
        integer :: i, a, n
        integer :: nvirt_avail

        ! 1. Chose a random connected spin basis vector
        ! (Use real space hubbard model procedure, since condition for
        ! connected 'determinants' is the same)
        call choose_ia_real(cdet%occ_list, cdet%f, i, a, nvirt_avail)

        ! 2. Find probability of generating this excited determinant.
        ! Again, same procedure for Heisenberg as for real space Hubbard
        pgen = calc_pgen_real(cdet%occ_list, cdet%f, nvirt_avail)

        ! 3. Construct the excited determinant and find the connecting matrix
        ! element.
        connection%nexcit = 1
        connection%from_orb(1) = i
        connection%to_orb(1) = a

        ! Non-zero off-diagonal elements are always -2J for Heisenebrg model
        hmatel = -2.0_p*J_coupling

        ! 4. Attempt spawning.
        pspawn = tau*abs(hmatel)/pgen
        psuccess = genrand_real2()

        ! Need to take into account the possibilty of a spawning attempt
        ! producing multiple offspring...
        ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and 
        ! then spawn a particle with probability pspawn-floor(pspawn).
        nspawn = int(pspawn)
        pspawn = pspawn - nspawn

        if (pspawn > psuccess) nspawn = nspawn + 1

        if (nspawn > 0) then

            ! 5. If H_ij is positive, then the spawned walker is of opposite
            ! sign to the parent, otherwise the spawned walkers if of the same
            ! sign as the parent.
            if (hmatel > 0.0_p) then
                nspawn = -sign(nspawn, parent_sign)
            else
                nspawn = sign(nspawn, parent_sign)
            end if

        end if

    end subroutine spawn_heisenberg
    
    subroutine spawn_heisenberg_importance_sampling(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new psip on a connected determinant for the 
        ! Heisenberg model. This subroutine applies a transformed Hamiltonian
        ! to apply importance sampling. It uses the Neel singlet state as a
        ! trial function. In all other ways it is the same as the basic
        ! Heisenberg spawning routine, spawn_heisenberg.
        !
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

        use basis, only: bit_lookup, basis_length
        use determinants, only: det_info, lattice_mask
        use dSFMT_interface, only:  genrand_real2
        use excitations, only: calc_pgen_real, excit
        use fciqmc_data, only: tau
        use fciqmc_data, only: neel_singlet_amp
        use fciqmc_data, only: walker_data, sampling_size
        use hamiltonian, only: slater_condon1_hub_real_excit
        use hamiltonian, only: diagonal_element_heisenberg
        use system, only: J_coupling, guiding_function
        use system, only: neel_singlet_guiding

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, psuccess, pspawn, hmatel, amp_j, amp_i
        integer(i0) :: bitstring_j(basis_length)
        integer :: i, a, n, up_spins_to, up_spins_from
        integer :: bit_position, bit_element
        integer :: nvirt_avail

        ! 1. Chose a random connected spin basis vector
        ! (Use real space hubbard model procedure, since condition for
        ! connected 'determinants' is the same)
        call choose_ia_real(cdet%occ_list, cdet%f, i, a, nvirt_avail)

        ! 2. Find probability of generating this excited determinant.
        ! Again, same procedure for Heisenberg as for real space Hubbard
        pgen = calc_pgen_real(cdet%occ_list, cdet%f, nvirt_avail)

        ! 3. Construct the excited determinant and find the connecting matrix
        ! element.
        connection%nexcit = 1
        connection%from_orb(1) = i
        connection%to_orb(1) = a

        ! Non-zero off-diagonal elements are always -2J for Heisenebrg model
        hmatel = -2.0_p*J_coupling
        
        ! Find the number of up spins on sublattice 1.
        up_spins_from = nint(cdet%data(sampling_size+1))
        ! For the spin up which was flipped to create the connected
        ! basis function, find whether this spin was on sublattice 1 or 2.
        ! If it was on sublattice 1, the basis function we go to has 1 less
        ! up spin on sublattice 1, else it will have one more spin up here.
        bit_position = bit_lookup(1,connection%from_orb(1))
        bit_element = bit_lookup(2,connection%from_orb(1))
        if (btest(lattice_mask(bit_element), bit_position)) then
            up_spins_to = up_spins_from-1
        else
            up_spins_to = up_spins_from+1
        end if
        
        ! When using a trial function |psi_T> = \sum{i} a_i|D_i>, the Hamiltonian
        ! used in importance sampling is 
        ! H_ji^T = a_j * H_ji * (1/a_i)
        ! For a given number of spins up on sublattice 1, n, the corresponding
        ! ampltidue of this basis function in the trial function is stored as
        ! neel_singlet_amp(n), for this particular trial function. Hence we have:
        hmatel = (neel_singlet_amp(up_spins_to)*hmatel)/neel_singlet_amp(up_spins_from)
        
        ! 4. Attempt spawning.
        pspawn = tau*abs(hmatel)/pgen
        psuccess = genrand_real2()

        ! Need to take into account the possibilty of a spawning attempt
        ! producing multiple offspring...
        ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and 
        ! then spawn a particle with probability pspawn-floor(pspawn).
        nspawn = int(pspawn)
        pspawn = pspawn - nspawn

        if (pspawn > psuccess) nspawn = nspawn + 1

        if (nspawn > 0) then

            ! 5. If H_ij is positive, then the spawned walker is of opposite
            ! sign to the parent, otherwise the spawned walkers if of the same
            ! sign as the parent.
            if (hmatel > 0.0_p) then
                nspawn = -sign(nspawn, parent_sign)
            else
                nspawn = sign(nspawn, parent_sign)
            end if

        end if

    end subroutine spawn_heisenberg_importance_sampling

    subroutine gen_excit_heisenberg(cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate both the probability 
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! In:
        !    cdet: info on the current basis function (equivalent to determinant
        !        in electron systems) that we will gen from.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D_i^a >, the Hamiltonian matrix element between a 
        !    determinant and a single excitation of it in the real space
        !    formulation of the Hubbard model.

        use determinants, only: det_info
        use excitations, only: calc_pgen_real, excit
        use system, only: J_coupling

        type(det_info), intent(in) :: cdet
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        integer :: nvirt_avail

        ! Double excitations are not connected.
        connection%nexcit = 1

        ! 1. Chose a random connected excitation.
        ! (Use real space hubbard model procedure, since condition for
        ! connected 'determinants' is the same)
        call choose_ia_real(cdet%occ_list, cdet%f, connection%from_orb(1), connection%to_orb(1), nvirt_avail)

        ! 2. Find probability of generating this excited determinant.
        ! Again, same procedure for Heisenberg as for real space Hubbard
        pgen = calc_pgen_real(cdet%occ_list, cdet%f, nvirt_avail)

        ! 3. find the connecting matrix element.
        ! Non-zero off-diagonal elements are always -2J for Heisenebrg model
        hmatel = -2.0_p*J_coupling

    end subroutine gen_excit_heisenberg

    subroutine choose_ij(occ_list, i ,j, ij_sym, ij_spin)

        ! Randomly choose a pair of spin-orbitals.
        ! See choose_ij_hub_k for a specific procedure for the momentum space
        ! formulation of the hubbard model.
        ! In:
        !    occ_list: integer list of occupied spin-orbitals in a determinant.
        ! Out:
        !    i, j: randomly selected spin-orbitals.
        !    ij_sym: symmetry label of the (i,j) combination.
        !    ij_spin: -1 if (i,j) are both beta, +1 if (i,j) are both alpha and
        !        0 if it's a "mixed" excitation, i.e. alpha, beta or beta,
        !        alpha.

        use basis, only: basis_fns
        use symmetry, only: sym_table
        use system, only: nel
        use dSFMT_interface, only: genrand_real2

        integer, intent(in) :: occ_list(nel)
        integer, intent(out) :: i,j, ij_sym, ij_spin
        integer :: ind, spin_sum
        real(dp) :: r

        ! We use a triangular indexing scheme to compress 2 electron indices
        ! into 1.
        ! For i/=j and (for an arbitrary choice of i>j), a 1D index of 
        ! a strictly lower triangular array is:
        !   p = (i-1)(i-2)/2 + j,
        ! where 1<=j<i and 1<=p<=n(n-1)/2
        ! This maps the indexing scheme as:
        !    .                  .
        !   2,1  .              1  .
        !   3,1 3,2  .      to  2  3  .
        !   4,1 4,2 4,3  .      3  4  5  .
        ! We want to do the reverse process in order to pick 2 electron labels
        ! from one random number.
        ! Consider the case where j=1.  i can trivially be determined from the
        ! quadratic equation:
        !   i = 3/2 + \sqrt(2p-1.75)
        ! As j<i and, for a fixed i, p increases monotonically with j, the
        ! integer part of i given by the above equation can never exceed the
        ! correct value.  Hence i can be found for arbitrary j by taking the
        ! floor of the preceeding equation.  j follows trivially.
        !
        ! See (for lower triangular arrays rather than strictly lower):
        ! Decoding the sequential indices of a (lower) triangular array
        ! SIS-75-1783,  S Rifkin, CERN report (CERN-DD-75-7).

        ! This might seem odd, but it enables us to pick the (i,j) pair to
        ! excite with half the calls to the random number generator, which
        ! represents a substantial saving. :-)

        r = genrand_real2()
        ind = int(r*nel*(nel-1)/2) + 1.

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.
        i = int(1.50_p + sqrt(2*ind-1.750_p))
        j = ind - ((i-1)*(i-2))/2

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.
        i = occ_list(i)
        j = occ_list(j)

        ! Symmetry info is a simple lookup...
        ij_sym = sym_table((i+1)/2,(j+1)/2)

        ! Is mod faster than lookup?  Not sure...
        spin_sum = basis_fns(i)%Ms + basis_fns(j)%Ms
        select case(spin_sum)
        case(2)
            ! alpha, alpha
            ij_spin = 1
        case(0)
            ! alpha, beta 
            ij_spin = 0
        case(-2)
            ! beta, beta 
            ij_spin = -1
        end select

    end subroutine choose_ij

    subroutine choose_ij_hub_k(occ_list_alpha, occ_list_beta, i ,j, ij_sym)

        ! Randomly choose a pair of spin-orbitals.
        !
        ! This is specific to the Hubbard model in momentum space.
        ! Only double excitations which excite from an alpha and a beta
        ! spin orbital are connected, so we return only i,j which meet this
        ! criterion.
        !
        ! In:
        !    occ_list_alpha: Integer list of occupied alpha spin-orbitals.
        !    occ_list_beta: Integer list of occupied beta spin-orbitals.
        ! Out:
        !    i, j: randomly selected spin-orbitals.
        !    ij_sym: symmetry label of the (i,j) combination.

        use symmetry, only: sym_table
        use system, only: nalpha, nbeta
        use dSFMT_interface, only: genrand_real2

        integer, intent(in) :: occ_list_alpha(nalpha), occ_list_beta(nbeta)
        integer, intent(out) :: i,j, ij_sym
        integer :: ind
        real(dp) :: r

        ! We use a similar indexing scheme to choose_ij, except our initial
        ! indices refer to an index in the occupied alpha array and in the
        ! occupied beta array.  This means that rather than be a triangular
        ! indexing scheme, we have a rectangular one:
        !   1  2  3    to    1,1  2,1  1,1
        !   3  5  6          1,2  2,2  3,2
        !   7  8  9          1,3  2,3  3,3
        !  10 11 12          1,4  2,4  3,4
        ! total number of possible combinations is nalpha*nbeta.
        ! The indexing scheme is:
        !  p = (i-1)*n_j + j
        ! Hence to invert this, following a similar method to Rifkin:
        !  i = floor( (p-1)/n_j ) + 1
        !  j = p - (i-1)*n_j

        r = genrand_real2()

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.

        ind = int(r*nalpha*nbeta) + 1
        i = int( (ind-1.0_p)/nbeta ) + 1
        j = ind - (i-1)*nbeta

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.
        i = occ_list_alpha(i)
        j = occ_list_beta(j)

        ! Symmetry info is a simple lookup...
        ij_sym = sym_table((i+1)/2,(j+1)/2)

    end subroutine choose_ij_hub_k

    subroutine choose_ab_hub_k(f, unocc_list_alpha, ij_sym, a, b)

        ! Choose a random pair of (a,b) unoccupied virtual spin-orbitals into
        ! which electrons are excited.
        ! (a,b) are chosen such that the (i,j)->(a,b) excitation is symmetry-
        ! allowed.
        ! In: 
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        !    unocc_alpha: integer list of the unoccupied alpha spin-orbitals.
        !    ij_sym: symmetry spanned by the (i,j) combination of unoccupied
        !        spin-orbitals into which electrons are excited.
        ! Returns:
        !    a,b: virtual spin orbitals involved in the excitation.

        use basis, only: basis_length
        use system, only: nvirt_alpha

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: unocc_list_alpha(nvirt_alpha)
        integer, intent(in) :: ij_sym
        integer, intent(out) :: a, b

        logical :: allowed_excitation

        allowed_excitation = .false.

        do while (.not.allowed_excitation)

            ! Until we find an allowed excitation.

            call find_ab_hub_k(f, unocc_list_alpha, ij_sym, a, b, allowed_excitation)

        end do

    end subroutine choose_ab_hub_k

    subroutine find_ab_hub_k(f, unocc_list_alpha, ij_sym, a, b, allowed_excitation)

        ! Choose a random pair of (a,b) spin-orbitals. 
        ! (a,b) are chosen such that the (i,j)->(a,b) excitation is symmetry-
        ! allowed and a is a virtual spin-orbital.  As (a,b) must be one alpha
        ! orbital and one beta orbital, we can choose a to be an alpha orbital
        ! without loss of generality.
        ! In: 
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        !    unocc_alpha: integer list of the unoccupied alpha spin-orbitals.
        !    ij_sym: symmetry spanned by the (i,j) combination of unoccupied
        !        spin-orbitals into which electrons are excited.
        ! Returns:
        !    a,b: spin orbitals involved in the excitation.
        !    allowed_excitation: is true if the excitation (i,j)->(a,b) is
        !        allowed (i.e. both a and b are indeed spin orbitals).

        use basis, only: basis_length, bit_lookup, nbasis
        use dSFMT_interface, only:  genrand_real2
        use system, only: nvirt_alpha
        use symmetry, only: sym_table, inv_sym

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: unocc_list_alpha(nvirt_alpha)
        integer, intent(in) :: ij_sym
        integer, intent(out) :: a, b
        logical, intent(out) :: allowed_excitation

        integer :: r, b_pos, b_el, ka

        ! The excitation i,j -> a,b is only allowed if k_i + k_j - k_a - k_b = 0
        ! (up to a reciprocal lattice vector).  We store k_i + k_j as ij_sym, so
        ! k_a + k_b must be identical to this.
        ! If we view this in terms of the representations spanned by i,j,a,b
        ! under translational symmetry (which forms an Abelian group) then
        !   \Gamma_i* x \Gamma_j* x \Gamma_a x \Gamma_b = \Gamma_1
        ! is equivalent to the conversation of crystal momentum (where \Gamma_1
        ! is the totally symmetric representation).  
        ! Further, as
        !   \Gamma_i* x \Gamma_i = 1
        ! and direct products in Abelian groups commute, it follows that:
        !   \Gamma_b = \Gamma_i x \Gamma_j x \Gamma_a*
        ! Thus k_b is defined by i,j and a.  As we only consider two-band
        ! systems, b is thus defined by spin conservation.

        ! One electron must be in unocc_list_alpha, so we can use the
        ! random number just to find which unoccupied alpha orbital is in
        ! the excitation.

        r = int(genrand_real2()*(nvirt_alpha)) + 1

        a = unocc_list_alpha(r)
        ! Find corresponding beta orbital which satisfies conservation
        ! of crystal momentum.
        ka = (a+1)/2
        b = 2*sym_table(ij_sym, inv_sym(ka))

        b_pos = bit_lookup(1,b)
        b_el = bit_lookup(2,b)

        ! If b is unoccupied then have found an allowed excitation.
        allowed_excitation = .not.btest(f(b_el), b_pos)

    end subroutine find_ab_hub_k

    subroutine choose_ia_real(occ_list, f, i, a, nvirt_avail)

        ! Find a random connected excitation from a Slater determinant for the
        ! Hubbard model in the real space formulation.
        ! This is also used for the Heisenberg model, which has the same conditions
        ! for two (different) basis functions to be connected. In this context,
        ! the number of electrons is actually the number of spins up (1 in binary).
        ! In: 
        !    f: bit string representation of the Slater determinant.
        !    occ_list: integer list of the occupied spin-orbitals in 
        !        the Slater determinant.
        ! Returns:
        !    i,a: spin orbitals excited from/to respectively.
        !    nvirt_avail: the number of virtual orbitals which can be excited
        !        into from the i-th orbital.

        use basis, only: basis_length, bit_lookup, nbasis
    
        use dSFMT_interface, only:  genrand_real2

        use basis, only: basis_length, basis_lookup
        use bit_utils, only: count_set_bits
        use hubbard_real, only: connected_orbs, connected_sites
        use system, only: nel, ndim

        integer, intent(in) :: occ_list(nel)
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(out) :: i, a, nvirt_avail
        integer(i0) :: virt_avail(basis_length)
        integer :: ivirt, ipos, iel, virt(3*ndim) ! 3*ndim to allow for triangular lattices; minor memory waste for other cases is irrelevant!

        do
            ! Until we find an i orbital which has at least one allowed
            ! excitation.

            ! Random selection of i.
            i = int(genrand_real2()*nel) + 1
            i = occ_list(i)

            ! Does this have at least one allowed excitation?
            ! connected_orbs(:,i) is a bit string with the bits corresponding to
            ! orbials connected to i set.
            ! The complement of the determinant bit string gives the bit string
            ! containing the virtual orbitals and thus taking the and of this
            ! with the relevant connected_orbs element gives the bit string
            ! containing the virtual orbitals which are connected to i.
            ! Neat, huh?
            virt_avail = iand(not(f), connected_orbs(:,i))

            if (any(virt_avail /= 0)) then
                ! Have found an i with at least one available orbital we can
                ! excite into.
                exit
            end if

        end do

        ! Find a.
        nvirt_avail = 0
        ! Now need to find out what orbital this corresponds to...
        do ivirt = 1, connected_sites(0,i)
            ipos = bit_lookup(1, connected_sites(ivirt,i))
            iel = bit_lookup(2, connected_sites(ivirt,i))
            if (btest(virt_avail(iel), ipos)) then
                nvirt_avail = nvirt_avail + 1
                virt(nvirt_avail) = connected_sites(ivirt,i)
            end if
        end do
        a = virt(int(genrand_real2()*nvirt_avail) + 1)

    end subroutine choose_ia_real

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

end module spawning
