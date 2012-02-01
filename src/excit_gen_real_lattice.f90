module excit_gen_real_lattice

! Module for random excitation generators and related routines for lattice model
! systems with real-space orbitals (i.e. Heisenberg model and Hubbard model in
! the local/real-space orbital basis).  These are grouped together due to the
! close relationship between the the systems means that we can re-use many of
! the same routines.

! See top-level comments in spawning about the overall aim of the spawning step.

use const
implicit none

contains

!---- Top level spawning routines ---

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
        use excitations, only: excit
        use fciqmc_data, only: tau
        use spawning, only: attempt_to_spawn

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel

        ! 1. Generate random excitation from cdet as well as both the probability
        ! of spawning there and the connecting matrix element.
        call gen_excit_hub_real(cdet, pgen, connection, hmatel)

        ! 2. Attempt to spawn child.
        nspawn = attempt_to_spawn(hmatel, pgen, parent_sign)

    end subroutine spawn_hub_real

    subroutine spawn_hub_real_no_renorm(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant for the
        ! real space formulation of the Hubbard model.
        !
        ! This uses excitation generators which, having selected the occupied
        ! orbital (i) to excite from, select any orbital (a) connected to
        ! i which conserves spin to excite to.  As a might be occupied, this
        ! is somewhat wasteful (generating excitations which can't be
        ! performed), there is a balance between the cost of generating
        ! forbidden excitations and the O(N) cost of renormalising the
        ! generation probabilities.
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
        use dSFMT_interface, only: genrand_real2
        use excitations, only: excit
        use system, only: nel
        use hamiltonian, only: slater_condon1_hub_real_excit
        use hubbard_real, only: connected_sites
        use basis, only: bit_lookup
        use spawning, only: attempt_to_spawn

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        integer :: i, iel, ipos
        real(p) :: pgen, hmatel

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

            ! 4. Attempt spawning.
            nspawn = attempt_to_spawn(hmatel, pgen, parent_sign)

        end if

    end subroutine spawn_hub_real_no_renorm

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
        use excitations, only: excit
        use spawning, only: attempt_to_spawn

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel

        ! 1. Generate a random excitation
        call gen_excit_heisenberg(cdet, pgen, connection, hmatel)

        ! 2. Attempt spawning.
        nspawn = attempt_to_spawn(hmatel, pgen, parent_sign)

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
        use excitations, only: excit
        use fciqmc_data, only: neel_singlet_amp, sampling_size
        use spawning, only: attempt_to_spawn

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel
        integer :: up_spins_to, up_spins_from
        integer :: bit_position, bit_element

        ! 1. Generate a random excitation.
        call gen_excit_heisenberg(cdet, pgen, connection, hmatel)

        ! 2. When using a trial function |psi_T> = \sum{i} a_i|D_i>, the Hamiltonian
        ! used in importance sampling is H_ji^T = a_j * H_ji * (1/a_i), so we
        ! need to adjust hmatel returned by gen_excit_heisenberg accordingly.

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

        ! For a given number of spins up on sublattice 1, n, the corresponding
        ! ampltidue of this basis function in the trial function is stored as
        ! neel_singlet_amp(n), for this particular trial function. Hence we have:
        hmatel = (neel_singlet_amp(up_spins_to)*hmatel)/neel_singlet_amp(up_spins_from)

        ! 3. Attempt spawning.
        nspawn = attempt_to_spawn(hmatel, pgen, parent_sign)

    end subroutine spawn_heisenberg_importance_sampling

    subroutine spawn_heisenberg_no_renorm(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new psip on a connected determinant for the
        ! Heisenberg model
        !
        ! This uses excitation generators which, having selected the spin-up
        ! site (i), select any site (a) connected to i to exchange spins with i.
        ! As a might be already spin-up, this is somewhat wasteful (generating
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
        use excitations, only: excit
        use system, only: J_coupling
        use spawning, only: attempt_to_spawn

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel

        ! 1. Generate a random excitation
        call gen_excit_heisenberg_no_renorm(cdet, pgen, connection, hmatel)

        ! 2. Attempt spawning if excitation is allowed.
        if (abs(hmatel) > depsilon) then

            nspawn = attempt_to_spawn(hmatel, pgen, parent_sign)

        else

            ! Generated a forbidden excitation (ie selected 2 spin up sites).
            nspawn = 0

        end if

    end subroutine spawn_heisenberg_no_renorm

    subroutine spawn_heisenberg_importance_sampling_no_renorm(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new psip on a connected determinant for the
        ! Heisenberg model. This subroutine applies a transformed Hamiltonian
        ! to apply importance sampling. It uses the Neel singlet state as a
        ! trial function.
        !
        ! This uses excitation generators which, having selected the spin-up
        ! site (i), select any site (a) connected to i to exchange spins with i.
        ! As a might be already spin-up, this is somewhat wasteful (generating
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

        use basis, only: bit_lookup, basis_length
        use determinants, only: det_info, lattice_mask
        use excitations, only: excit
        use fciqmc_data, only: neel_singlet_amp, sampling_size
        use spawning, only: attempt_to_spawn

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, hmatel
        integer :: up_spins_to, up_spins_from
        integer :: bit_position, bit_element

        ! 1. Generate a random excitation.
        call gen_excit_heisenberg_no_renorm(cdet, pgen, connection, hmatel)

        if (abs(hmatel) < depsilon) then

            ! Generated a forbidden excitation (ie selected two spin up sites)
            nspawn = 0

        else

            ! 2. When using a trial function |psi_T> = \sum{i} a_i|D_i>, the Hamiltonian
            ! used in importance sampling is H_ji^T = a_j * H_ji * (1/a_i), so we
            ! need to adjust hmatel returned by gen_excit_heisenberg accordingly.

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

            ! For a given number of spins up on sublattice 1, n, the corresponding
            ! ampltidue of this basis function in the trial function is stored as
            ! neel_singlet_amp(n), for this particular trial function. Hence we have:
            hmatel = (neel_singlet_amp(up_spins_to)*hmatel)/neel_singlet_amp(up_spins_from)

            ! 3. Attempt spawning.
            nspawn = attempt_to_spawn(hmatel, pgen, parent_sign)

        end if

    end subroutine spawn_heisenberg_importance_sampling_no_renorm

!--- Excitation generators: Hubbard model ---

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
        use excitations, only: excit
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

!--- Excitation generators: Heisenberg model ---

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
        use excitations, only: excit
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

    subroutine gen_excit_heisenberg_no_renorm(cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        !
        ! This selects the spin-up site (i) and then selects any site (a)
        ! connected to i to exchange spins with i.  As a might be already
        ! spin-up, this is somewhat wasteful (generating excitations which can't
        ! be performed), there is a balance between the cost of generating
        ! forbidden excitations and the O(N) cost of renormalising the
        ! generation probabilities.
        !
        ! If a forbidden excitation is generated, then hmatel is set to 0 and
        ! pgen to 1.
        !
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

        use basis, only: bit_lookup
        use determinants, only: det_info
        use excitations, only: excit
        use hubbard_real, only: connected_sites
        use system, only: J_coupling, nel

        use dSFMT_interface, only: genrand_real2

        type(det_info), intent(in) :: cdet
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        integer :: i, ipos, iel

        ! 1. Chose a random connected excitation.
        ! Random selection of i (an up spin).
        i = int(genrand_real2()*nel) + 1
        connection%from_orb(1) = cdet%occ_list(i)
        ! Select a at random from one of the connected sites.
        ! nb: a might already be up.
        i = int(genrand_real2()*connected_sites(0,connection%from_orb(1)) + 1)
        connection%to_orb(1) = connected_sites(i,connection%from_orb(1))

        ! Is a already up?  If so, not an allowed excitation: return a null
        ! event.
        ipos = bit_lookup(1, connection%to_orb(1))
        iel = bit_lookup(2, connection%to_orb(1))

        if (btest(cdet%f(iel),ipos)) then

            ! null event means setting hmatel = 0.
            ! also set pgen to avoid potential division by an undefined quantity
            ! causing issues.
            hmatel = 0.0_p
            pgen = 1.0_p

        else

            connection%nexcit = 1

            ! 2. Generation probability
            ! For single excitations
            !   pgen = p(i) p(a|i)
            !        = 1/(nel*nconnected_sites)
            pgen = 1.0_dp/(nel*connected_sites(0,i))

            ! 3. find the connecting matrix element.
            ! Non-zero off-diagonal elements are always -2J for Heisenebrg model
            hmatel = -2.0_p*J_coupling

        end if

    end subroutine gen_excit_heisenberg_no_renorm

!--- Selecting random orbitals ---

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

!--- Excitation generation probabilities ---

    pure function calc_pgen_real(occ_list, f, nvirt_avail) result(pgen)

        ! Calculate the generation probability of a given excitation for the
        ! real space models, the real Hubbard model and the Heisenberg model
        !
        ! Note that all the information required for input should be available
        ! during the FCIQMC algorithm and should not be needed to be calculated.
        !
        ! Note also that electrons and virtual orbitals in the Hubbard model are
        ! equivalent to spins up and spins down respectively in the Heisenberg model.
        !
        ! Further, we assume only allowed excitations are generated.
        !
        ! In:
        !    occ_list: integer list of occupied orbitals in the Slater determinant.
        !    f: bit string representation of the determinant we're exciting
        !        from.
        !    nvirt_avail: the number of available orbitals into which an
        !        electron can be excited, given the choice of the orbital which
        !        is being excited from (i.e. having chosen i, how many
        !        possibilities are there for a, where i is occupied and
        !        a occupied and D_i^a is connected to D.
        !        For Heisenberg model, nvirt_avail is the number spins down
        !        which can be changed to a spin up, given the spin up we have
        !        chosen to try and excite from
        ! Returns:
        !    pgen: the generation probability of the excitation.  See notes in
        !        spawning.

        use basis, only: basis_length
        use system, only: nel
        use hubbard_real, only: connected_orbs

        use errors

        real(p) :: pgen
        integer, intent(in) :: occ_list(nel)
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: nvirt_avail
        integer :: i, no_excit

        ! For single excitations
        !   pgen = p(i) p(a|i) \chi_r
        ! where
        !   p(i) is the probability of choosing the i-th electron to excite
        !   from.
        !   p(i) = 1/nel
        !
        !   p(a|i) is the probability of choosing to excite into orbital a given
        !   that the i-th electron has been chosen.
        !   nvirt_avail is the number of virtual orbitals connected to i and is
        !   calculated when the random excitation is chosen.
        !   p(a|i) = 1/nvirt_avail

        !   \chi_r is a renormalisation to take into account the fact that not
        !   all electrons may be excited from (e.g. no connected orbitals are
        !   vacant).
        !   \chi_r = nel/(nel - no_excit)
        !   where no_excit is the number of occupied orbitals which have no
        !   connected excitations.

        ! \chi_r is a constant for a given determinant, so an optimisation is to
        ! calculate this once per determinant rather than for each walker on the
        ! same determinant.

        no_excit = 0
        do i = 1, nel
            ! See if there are any allowed excitations from this electron
            ! (Or excitations from this spin up for Hesienberg)
            ! (see notes in choose_ia_real for how this works)
            if (all(iand(not(f), connected_orbs(:,occ_list(i))) == 0)) then
                ! none allowed from this orbial
                no_excit = no_excit + 1
            end if
        end do

        pgen = 1.0_p/(nvirt_avail * (nel - no_excit))

    end function calc_pgen_real

end module excit_gen_real_lattice
