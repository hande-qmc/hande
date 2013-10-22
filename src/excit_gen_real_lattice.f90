module excit_gen_real_lattice

! Module for random excitation generators and related routines for sys_global%lattice%lattice model
! systems with real-space orbitals (i.e. Heisenberg model and Hubbard model in
! the local/real-space orbital basis).  These are grouped together due to the
! close relationship between the the systems means that we can re-use many of
! the same routines.

! See top-level comments in spawning about the overall aim of the spawning step.

! TODO: 'split' excitation generators (optimisation already used for Hubbard
! model in momentum space).

use const
implicit none

contains

!--- Excitation generators: Hubbard model ---

    subroutine gen_excit_hub_real(rng, cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D_i^a >, the Hamiltonian matrix element between a
        !    determinant and a single excitation of it in the real space
        !    formulation of the Hubbard model.

        use determinants, only: det_info
        use excitations, only: excit
        use hamiltonian_hub_real, only: slater_condon1_hub_real_excit
        use dSFMT_interface, only: dSFMT_t

        type(det_info), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        integer :: nvirt_avail

        ! Double excitations are not connected determinants within the
        ! real space formulation of the Hubbard model.
        connection%nexcit = 1

        ! 1. Chose a random connected excitation.
        call choose_ia_real(rng, cdet%occ_list, cdet%f, connection%from_orb(1), connection%to_orb(1), nvirt_avail)

        ! 2. Find probability of generating this excited determinant.
        pgen = calc_pgen_real(cdet%occ_list, cdet%f, nvirt_avail)

        ! 3. find the connecting matrix element.
        call slater_condon1_hub_real_excit(cdet%f, connection, hmatel)

    end subroutine gen_excit_hub_real

    subroutine gen_excit_hub_real_no_renorm(rng, cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        !
        ! This selects the occupied orbital (i) and then selects any orbital (a)
        ! connected to i to excite an electron into.  As a might be already
        ! occupied, this is somewhat wasteful (generating excitations which can't
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
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D_i^a >, the Hamiltonian matrix element between a
        !    determinant and a single excitation of it in the real space
        !    formulation of the Hubbard model.

        use determinants, only: det_info
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use excitations, only: excit
        use system
        use hamiltonian_hub_real, only: slater_condon1_hub_real_excit
        use hubbard_real, only: connected_sites
        use basis, only: bit_lookup
        use spawning, only: attempt_to_spawn

        type(det_info), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        type(excit), intent(out) :: connection
        real(p), intent(out) :: pgen, hmatel

        integer :: i, iel, ipos

        ! 1. Generate random excitation from cdet and probability of spawning
        ! there.
        ! Random selection of i.
        i = int(get_rand_close_open(rng)*sys_global%nel) + 1
        connection%from_orb(1) = cdet%occ_list(i)
        ! Select a at random from one of the connected orbitals.
        i = int(get_rand_close_open(rng)*connected_sites(0,connection%from_orb(1)) + 1)
        connection%to_orb(1) = connected_sites(i,connection%from_orb(1))

        ipos = bit_lookup(1, connection%to_orb(1))
        iel = bit_lookup(2, connection%to_orb(1))

        if (btest(cdet%f(iel),ipos)) then

            ! a is occupied; forbidden excitation
            pgen = 1.0_p
            hmatel = 0.0_p

        else

            connection%nexcit=1

            ! 2. find the connecting matrix element.
            call slater_condon1_hub_real_excit(cdet%f, connection, hmatel)

            ! 3. Probability of spawning...
            ! For single excitations
            !   pgen = p(i) p(a|i)
            !        = 1/(sys_global%nel*nconnected_sites)
            pgen = 1.0_dp/(sys_global%nel*connected_sites(0,i))

        end if

    end subroutine gen_excit_hub_real_no_renorm

!--- Excitation generators: Heisenberg model ---

    subroutine gen_excit_heisenberg(rng, cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! In:
        !    cdet: info on the current basis function (equivalent to determinant
        !        in electron systems) that we will gen from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D_i^a >, the Hamiltonian matrix element between a
        !    determinant and a single excitation of it in the real space
        !    formulation of the Hubbard model.

        use determinants, only: det_info
        use excitations, only: excit
        use system
        use dSFMT_interface, only: dSFMT_t

        type(det_info), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        integer :: nvirt_avail

        ! Double excitations are not connected.
        connection%nexcit = 1

        ! 1. Chose a random connected excitation.
        ! (Use real space hubbard model procedure, since condition for
        ! connected 'determinants' is the same)
        call choose_ia_real(rng, cdet%occ_list, cdet%f, connection%from_orb(1), connection%to_orb(1), nvirt_avail)

        ! 2. Find probability of generating this excited determinant.
        ! Again, same procedure for Heisenberg as for real space Hubbard
        pgen = calc_pgen_real(cdet%occ_list, cdet%f, nvirt_avail)

        ! 3. find the connecting matrix element.
        ! Non-zero off-diagonal elements are always -2J for Heisenebrg model
        hmatel = -2.0_p*sys_global%heisenberg%J

    end subroutine gen_excit_heisenberg

    subroutine gen_excit_heisenberg_no_renorm(rng, cdet, pgen, connection, hmatel)

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
        ! In/Out:
        !    rng: random number generator.
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
        use system

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(det_info), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        integer :: i, ipos, iel

        ! 1. Chose a random connected excitation.
        ! Random selection of i (an up spin).
        i = int(get_rand_close_open(rng)*sys_global%nel) + 1
        connection%from_orb(1) = cdet%occ_list(i)
        ! Select a at random from one of the connected sites.
        ! nb: a might already be up.
        i = int(get_rand_close_open(rng)*connected_sites(0,connection%from_orb(1)) + 1)
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
            !        = 1/(sys_global%nel*nconnected_sites)
            pgen = 1.0_dp/(sys_global%nel*connected_sites(0,i))

            ! 3. find the connecting matrix element.
            ! Non-zero off-diagonal elements are always -2J for Heisenebrg model
            hmatel = -2.0_p*sys_global%heisenberg%J

        end if

    end subroutine gen_excit_heisenberg_no_renorm

!--- Selecting random orbitals ---

    subroutine choose_ia_real(rng, occ_list, f, i, a, nvirt_avail)

        ! Find a random connected excitation from a Slater determinant for the
        ! Hubbard model in the real space formulation.
        ! This is also used for the Heisenberg model, which has the same conditions
        ! for two (different) basis functions to be connected. In this context,
        ! the number of electrons is actually the number of spins up (1 in binary).
        ! In:
        !    f: bit string representation of the Slater determinant.
        !    occ_list: integer list of the occupied spin-orbitals in
        !        the Slater determinant.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    i,a: spin orbitals excited from/to respectively.
        !    nvirt_avail: the number of virtual orbitals which can be excited
        !        into from the i-th orbital.

        use basis, only: basis_length, bit_lookup, nbasis

        use dSFMT_interface, only:  dSFMT_t, get_rand_close_open

        use basis, only: basis_length, basis_lookup
        use bit_utils, only: count_set_bits
        use hubbard_real, only: connected_orbs, connected_sites
        use system

        integer, intent(in) :: occ_list(sys_global%nel)
        integer(i0), intent(in) :: f(basis_length)
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: i, a, nvirt_avail
        integer(i0) :: virt_avail(basis_length)
        integer :: ivirt, ipos, iel, virt(3*sys_global%lattice%ndim) ! 3*sys_global%lattice%ndim to allow for triangular lattices; minor memory waste for other cases is irrelevant!

        do
            ! Until we find an i orbital which has at least one allowed
            ! excitation.

            ! Random selection of i.
            i = int(get_rand_close_open(rng)*sys_global%nel) + 1
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
        a = virt(int(get_rand_close_open(rng)*nvirt_avail) + 1)

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
        use system
        use hubbard_real, only: connected_orbs

        use errors

        real(p) :: pgen
        integer, intent(in) :: occ_list(sys_global%nel)
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: nvirt_avail
        integer :: i, no_excit

        ! For single excitations
        !   pgen = p(i) p(a|i) \chi_r
        ! where
        !   p(i) is the probability of choosing the i-th electron to excite
        !   from.
        !   p(i) = 1/sys_global%nel
        !
        !   p(a|i) is the probability of choosing to excite into orbital a given
        !   that the i-th electron has been chosen.
        !   nvirt_avail is the number of virtual orbitals connected to i and is
        !   calculated when the random excitation is chosen.
        !   p(a|i) = 1/nvirt_avail

        !   \chi_r is a renormalisation to take into account the fact that not
        !   all electrons may be excited from (e.g. no connected orbitals are
        !   vacant).
        !   \chi_r = sys_global%nel/(sys_global%nel - no_excit)
        !   where no_excit is the number of occupied orbitals which have no
        !   connected excitations.

        ! \chi_r is a constant for a given determinant, so an optimisation is to
        ! calculate this once per determinant rather than for each walker on the
        ! same determinant.

        no_excit = 0
        do i = 1, sys_global%nel
            ! See if there are any allowed excitations from this electron
            ! (Or excitations from this spin up for Hesienberg)
            ! (see notes in choose_ia_real for how this works)
            if (all(iand(not(f), connected_orbs(:,occ_list(i))) == 0)) then
                ! none allowed from this orbial
                no_excit = no_excit + 1
            end if
        end do

        pgen = 1.0_p/(nvirt_avail * (sys_global%nel - no_excit))

    end function calc_pgen_real

end module excit_gen_real_lattice
