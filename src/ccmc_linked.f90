module ccmc_linked

! Module containing all utility functions only used within linked ccmc.
! For full explanation see top of ccmc.F90.

use const, only: i0, p, dp

implicit none

contains

    pure subroutine linked_excitation(basis, f0, connection, cluster, linked, single_unlinked, excitor)

        ! For Linked Coupled Cluster, the only terms of H that need to
        ! be sampled are those which are connected to (ie have a
        ! creation/annihilation operator in common with) each excitor in the
        ! cluster (by Wick's Theorem). This routine tests which excitors are
        ! connected to the Hamiltonian.

        ! In:
        !    basis: information about the single-particle basis
        !    f0: bit string of the reference
        !    connection: the excitation connecting the current excitor and the
        !       child excitor
        !    cluster: the cluster of excitation operators
        ! Out:
        !    linked: true if the Hamiltonian is connected to the cluster
        !    single_unlinked: true if connection is a single excitation and
        !       cluster contains one excitor that does not share an orbital with
        !       connection. (It is instead connected to the Hamiltonian by the
        !       dummy index in the two-body term \sum_j <aj||ij> a^+_a a^+_j a_j a_i)
        !       Any other excitors in the cluster must share an orbital with connection.
        !    excitor: if single_unlinked is true, this is the bit string of the
        !       excitor not linked to connection (otherwise 0)

        use excitations, only: excit_t, create_excited_det
        use basis_types, only: basis_t
        use ccmc_data, only: cluster_t

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f0(basis%string_len)
        type(excit_t), intent(in) :: connection
        type(cluster_t), intent(in) :: cluster
        logical, intent(out) :: linked, single_unlinked
        integer(i0), intent(out) :: excitor(basis%string_len)

        integer :: i, orb, bit_pos, bit_element, unconnected
        integer(i0) :: excitor_excitation(basis%string_len)
        integer(i0) :: h_excitation(basis%string_len)

        single_unlinked = .false.
        unconnected = 0
        excitor = 0

        ! Get bit string of orbitals in H excitation
        ! (modified from create_excited_det)
        h_excitation = 0
        do i = 1, connection%nexcit
            ! i/j orbital
            orb = connection%from_orb(i)
            bit_pos = basis%bit_lookup(1,orb)
            bit_element = basis%bit_lookup(2,orb)
            h_excitation(bit_element) = ibset(h_excitation(bit_element), bit_pos)
            ! a/b orbital
            orb = connection%to_orb(i)
            bit_pos = basis%bit_lookup(1,orb)
            bit_element = basis%bit_lookup(2,orb)
            h_excitation(bit_element) = ibset(h_excitation(bit_element), bit_pos)
        end do

        do i = 1, cluster%nexcitors
            ! check each cluster operator shares an index with connection
            ! orbitals involved in cluster operator excitation (from reference)
            excitor_excitation = ieor(cluster%excitors(i)%f, f0)
            if (all(iand(h_excitation, excitor_excitation) == 0)) then
                ! no orbitals in common between H and cluster
                unconnected = unconnected + 1
                excitor = cluster%excitors(i)%f
            end if
        end do

        select case(connection%nexcit)
        case(1)
            ! For a single excitation, H contains the sum
            ! \sum_j <ij||aj>a^+j^+ji so can connect to one excitor that has no
            ! orbitals in common with the excitation
            linked = (unconnected <= 1)
            single_unlinked = (unconnected == 1)
        case(2)
            ! Double excitation, H only has the term <ab||ij>a^+b^+ji so must
            ! have an orbital in common with all cluster operators
            linked = (unconnected == 0)
        end select

        if (.not. linked) excitor = 0

    end subroutine linked_excitation

    pure function unlinked_commutator(sys, f0, connection, cluster, cdet, funlinked) result(hmatel)

        ! When H is a single excitation, it can be connected to one of the
        ! excitors in a cluster (a1) by the creation / annihilation operators
        ! corresponding to the dummy index in the two-body term
        !   \sum_j <aj||ij> a^+ j^+ i j.
        ! In that case, both terms of the commutator <D_j|[H,a1]|D_i> (where D_i
        ! is the determinant obtained by applying the other excitors in the
        ! cluster to the reference) are non-zero. <D_j|H a1|D_i> has already been
        ! calculated in the excitation generator; this function evaluates <D_j|a1 H|D_i>.

        ! In:
        !    sys: the system being studied
        !    f0: bit string of the reference
        !    connection: excitation connection between the current excitor
        !        and the child excitor, on which progeny are spawned.
        !    cluster: the cluster of excitors
        !    cdet: the determinant formed by the cluster
        !    funlinked: the excitor in the cluster that is not linked to H
        ! Returns:
        !    <Dj|a H|Di> (to be subtracted from <Dj|H a|Di> calculated when the
        !           excitation was chosen to give the commutator)

        use system, only: sys_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use ccmc_data, only: cluster_t
        use ccmc_utils, only: collapse_cluster, convert_excitor_to_determinant
        use hamiltonian, only: get_hmatel
        use hamiltonian_data

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(sys%basis%string_len)
        type(excit_t), intent(in) :: connection
        type(cluster_t), intent(in) :: cluster
        integer(i0), intent(in) :: cdet(sys%basis%string_len)
        integer(i0), intent(in) :: funlinked(sys%basis%string_len)
        real(p) :: hmatel
        type(hmatel_t) :: dummy_hmatel

        integer(i0) :: deti(sys%basis%string_len), detj(sys%basis%string_len)
        integer(i0) :: temp(sys%basis%string_len)
        integer :: i, found, excitor_level, excitor_sign
        logical :: allowed
        real(p) :: population

        ! Find the determinant obtained by applying all of the cluster operators
        ! linked to H to D0
        population = 0.0_p ! The population doesn't matter as the commutator does not change the amplitude
        found = 0
        if (cluster%nexcitors > 1) then
            do i = 1, cluster%nexcitors
                if (any(cluster%excitors(i)%f /= funlinked)) then
                    ! Linked excitor, needed in cluster
                    found = found + 1
                    if (found == 1) then
                        deti = cluster%excitors(i)%f
                    else
                        call collapse_cluster(sys%basis, f0, cluster%excitors(i)%f, 1.0_p, &
                            deti, population, allowed)
                    end if
                end if
            end do
        else
            ! Only the unlinked excitor
            deti = f0
        end if

        ! Now we want to evaluate <D_i^a|H_i^a|D> ...
        call create_excited_det(sys%basis, deti, connection, detj)
        ! [todo] - general case call is slow.  Improvements: Slater--Condon procedure for
        ! [todo] - the relevant excitation level and system-specific procedures.
        dummy_hmatel = get_hmatel(sys, deti, detj)
        hmatel = dummy_hmatel%r

        ! hmatel will be multiplied by cluster%amplitude and cluster%cluster_to_det_sign which
        ! potentially introduce unwanted sign changes, so we deal with them here
        hmatel = hmatel*cluster%cluster_to_det_sign
        ! Multiplying excitors can give a sign change, which is absorbed into cluster%amplitude
        population = 1.0_p
        temp = deti
        call collapse_cluster(sys%basis, f0, funlinked, 1.0_p, temp, population, allowed)
        hmatel = population*hmatel

        ! Possible sign changes from <D|a|D_0> ...
        if (cluster%nexcitors > 1) then
            excitor_level = get_excitation_level(f0, deti)
            call convert_excitor_to_determinant(deti, excitor_level, excitor_sign, f0)
            if (excitor_sign < 0) hmatel = -hmatel
        end if

        ! ... and <D_k|a_unlinked|D_i^a>
        call create_excited_det(sys%basis, cdet, connection, deti)
        excitor_level = get_excitation_level(deti, detj)
        call convert_excitor_to_determinant(deti, excitor_level, excitor_sign, detj)
        if (excitor_sign < 0) hmatel = -hmatel

    end function unlinked_commutator

    subroutine partition_cluster(rng, sys, f0, cluster, left_cluster, right_cluster, ppart, &
                                 ldet, rdet, allowed, sign_change, part_number)

        ! Divides a cluster into two halves such that any excitors that share an
        ! orbital are in different halves

        ! In:
        !    cluster: the cluster of excitors
        !    sys: the system being studied
        !    f0: bit string of the reference
        !    part_number: (optional) use to enumerate partitions rather than
        !            choosing a random one
        ! In/Out:
        !    rng: random number generator
        ! Out:
        !    left_cluster: the excitors to be applied after the Hamiltonian
        !    right_cluster: the excitors to be applied before the Hamiltonian
        !    ppart: the probability of choosing this partition
        !    ldet: the determinant formed by applying left_cluster to the reference
        !    rdet: the determinant formed by applying right_cluster to the reference
        !    allowed: are left_cluster and right_cluster both valid clusters
        !    sign_change: is there a net sign change on collapsing the two clusters

        use ccmc_data, only: cluster_t
        use ccmc_utils, only: collapse_cluster, convert_excitor_to_determinant
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(sys%basis%string_len)
        type(cluster_t), intent(in) :: cluster
        type(cluster_t), intent(inout) :: left_cluster, right_cluster
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: ppart
        integer(i0), intent(inout) :: ldet(sys%basis%string_len), rdet(sys%basis%string_len)
        logical, intent(out) :: allowed, sign_change
        integer, intent(in), optional :: part_number

        integer :: i, in_left, in_right, side
        real(dp) :: rand
        real(p) :: population

        in_left = 0
        in_right = 0
        ppart = 1.0_p
        allowed = .true.
        population = 1.0_p

        do i = 1, cluster%nexcitors
            if (present(part_number)) then
                ! Use integers 1 to 2^n as bitstrings to consider all
                ! partititions in turn
                if (btest(part_number, i-1)) then
                    side = 1
                else
                    side = -1
                end if
            else
                ! for simplicity of generation probabilities just assign everything
                ! randomly to left or right
                ppart = ppart * 0.5_p
                rand = get_rand_close_open(rng)
                if (rand < 0.5_dp) then
                    side = -1
                else
                    side = 1
                end if
            end if
            if (side == -1) then
                ! add to left_cluster
                in_left = in_left + 1
                left_cluster%excitors(in_left)%f => cluster%excitors(i)%f
                if (in_left == 1) then
                    ldet = cluster%excitors(i)%f
                else
                    call collapse_cluster(sys%basis, f0, cluster%excitors(i)%f, 1.0_p, ldet, population, allowed)
                    if (.not.allowed) exit
                end if
            else
                ! add to right
                in_right = in_right + 1
                right_cluster%excitors(in_right)%f => cluster%excitors(i)%f
                if (in_right == 1) then
                    rdet = cluster%excitors(i)%f
                else
                    call collapse_cluster(sys%basis, f0, cluster%excitors(i)%f, 1.0_p, rdet, population, allowed)
                    if (.not.allowed) exit
                end if
            end if
        end do

        left_cluster%nexcitors = in_left
        right_cluster%nexcitors = in_right

        sign_change = (population < 0)

    end subroutine partition_cluster

    function calc_pgen(sys, excit_gen, excit_gen_data, f, connection, parent_det) result(pgen)

        ! Calculate the probability of an excitation being selected.
        ! Wrapper round system specific functions.

        ! In:
        !    sys: the system being studied
        !    excit_gen: which excitation generator is being used.  Should correspond to a value in the excit_gen_* enum in qmc_data.
        !    qmc_state: input options relating to QMC methods.
        !    f: bit string representation of parent excitor
        !    connection: excitation connection between the current excitor
        !        and the child excitor, on which progeny are spawned.
        !    parent_det: information on the parent excitor
        ! Returns:
        !    The probability of generating the excitation specified by
        !    connection, assuming the excitation is valid

        use bit_utils, only: count_set_bits
        use errors, only: stop_all
        use system
        use excitations, only: excit_t
        use excit_gen_mol, only: calc_pgen_single_mol_no_renorm, calc_pgen_double_mol_no_renorm, &
                                 calc_pgen_single_mol, calc_pgen_double_mol
        use excit_gen_ueg, only: calc_pgen_ueg_no_renorm
        use excit_gen_ringium, only: calc_pgen_ringium
        use read_in_symmetry, only: cross_product_basis_read_in
        use determinants, only: det_info_t
        use qmc_data, only: excit_gen_no_renorm
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: excit_gen
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        integer(i0), intent(in) :: f(sys%basis%string_len)
        type(excit_t), intent(in) :: connection
        type(det_info_t), intent(in) :: parent_det
        real(p) :: pgen

        integer :: spin, ij_sym, max_na, ij_k(3)
        integer(i0) :: poss_a(sys%basis%string_len)

        associate(a=>connection%to_orb(1), b=>connection%to_orb(2), i=>connection%from_orb(1), j=>connection%from_orb(2))
            select case(sys%system)
            case(read_in)
                if (excit_gen == excit_gen_no_renorm) then
                    if (connection%nexcit == 1) then
                        pgen = excit_gen_data%pattempt_single * calc_pgen_single_mol_no_renorm(sys, a)
                    else
                        spin = sys%basis%basis_fns(a)%ms + sys%basis%basis_fns(b)%ms
                        pgen = excit_gen_data%pattempt_double * calc_pgen_double_mol_no_renorm(sys, a, b, spin)
                    end if
                else
                    if (connection%nexcit == 1) then
                        pgen = excit_gen_data%pattempt_single * calc_pgen_single_mol(sys, sys%read_in%pg_sym%gamma_sym, &
                                                                                parent_det%occ_list, parent_det%symunocc, a)
                    else
                        spin = sys%basis%basis_fns(a)%ms + sys%basis%basis_fns(b)%ms
                        ij_sym = sys%read_in%sym_conj_ptr(sys%read_in, &
                                    cross_product_basis_read_in(sys, a, b))
                        pgen = excit_gen_data%pattempt_double * calc_pgen_double_mol(sys, ij_sym, a, b, spin, parent_det%symunocc)
                    end if
                end if
            case(ueg)
                spin = sys%basis%basis_fns(i)%ms + sys%basis%basis_fns(j)%ms
                ij_k = 0
                ij_k(1:sys%lattice%ndim) = sys%basis%basis_fns(i)%l + sys%basis%basis_fns(j)%l
                if (spin == -2) then
                    poss_a = iand(not(f), ishft(excit_gen_data%ueg_ternary_conserve(1:,ij_k(1),ij_k(2),ij_k(3)),1))
                else
                    poss_a = iand(not(f), excit_gen_data%ueg_ternary_conserve(1:,ij_k(1),ij_k(2),ij_k(3)))
                end if
                max_na = sum(count_set_bits(poss_a))
                pgen = calc_pgen_ueg_no_renorm(sys, max_na, spin)
            case(ringium)
                pgen = calc_pgen_ringium(sys)
            case default
                call stop_all('calc_pgen', 'Linked CCMC is not implemented for this system.')
            end select
        end associate

    end function calc_pgen

end module ccmc_linked
