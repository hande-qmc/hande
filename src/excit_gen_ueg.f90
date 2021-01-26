module excit_gen_ueg

! Module for random excitation generators and related routines for the uniform
! electron gas.

! The electron gas is wonderfully simple and only contains double
! excitations (much like, in fact, the Hubbard model).

! NOTE: Besides power-pitzer excitation generators, only unrenormalised
! excitation generators are fully implemented (ie where we explicitly discard
! excitations (i,j)->(a,b) if b is already occupied rather than only generating
! excitations which are allowed). This is because, as there is no integral
! storage bottleneck in the UEG, the number of basis functions tends to be much
! larger than the number of electrons.  Consequently the fraction of such
! forbidden excitations is quite small and so it is not worth the expense of
! renormalising. We still generate excitations (i,j)->(a,b) which are allowed
! both by symmetry (conservation of crystal momentum) and spin.

use const, only: i0, p, dp

implicit none

contains

!--- Excitation generators ---

    subroutine gen_excit_ueg_no_renorm(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! This doesn't exclude the case where, having selected all orbitals
        ! involved in the excitation, the final orbital selected is already
        ! occupied and so cannot be excited into.  Whilst this is somewhat
        ! wasteful (generating excitations which can't be performed), there is
        ! a balance between the cost of generating forbidden excitations and the
        ! O(N) cost of renormalising the generation probabilities.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Data for excitation generator (not used) 
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    rng: random number generator.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !       determinant and a connected determinant in the UEG.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation2
        use excit_gens, only: excit_gen_data_t
        use hamiltonian_ueg, only: slater_condon2_ueg_excit
        use system, only: sys_t
        use hamiltonian_data

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use bit_utils

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(inout) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation

        integer :: ij_k(sys%lattice%ndim), ij_spin, max_na

        ! 1. must have a double excitation.
        connection%nexcit = 2

        ! 2. Select orbitals to excite from and into.
        call choose_ij_k(rng, sys, cdet%occ_list, connection%from_orb(1), connection%from_orb(2), ij_k, ij_spin)
        call find_ab_ueg(rng, sys, cdet%f, ij_k, ij_spin, excit_gen_data%ueg_ternary_conserve, max_na, connection%to_orb(1), &
                         connection%to_orb(2), allowed_excitation)
        
        if (allowed_excitation) then

            ! 3. Return generation probability and conecting matrix element.

            pgen = calc_pgen_ueg_no_renorm(sys, max_na, ij_spin)

            call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)
            hmatel%r = slater_condon2_ueg_excit(sys, connection%from_orb(1), connection%to_orb(1), connection%to_orb(2), &
                connection%perm)


        else

            ! Carelessly selected ij with no available excitations.  Not worth
            ! renormalising.  Discard: return a null excitation.
            hmatel%r = 0.0_p
            pgen = 1.0_p

        end if

    end subroutine gen_excit_ueg_no_renorm

!--- Select random orbitals involved in a valid double excitation ---

    subroutine choose_ij_k(rng, sys, occ_list, i, j, ij_k, ij_spin)

        ! Randomly choose a pair of spin-orbitals for 1-band systems with
        ! Bloch orbitals.
        !
        ! See choose_ij_hub_k for a specific procedure for the momentum space
        ! formulation of the hubbard model.
        !
        ! In:
        !    occ_list: integer list of occupied spin-orbitals in a determinant.
        !        (min length: sys%nel.)
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    i, j: randomly selected spin-orbitals.
        !    ij_spin: -2 if (i,j) are both beta, +2 if (i,j) are both alpha and
        !        0 if it's a "mixed" excitation, i.e. alpha, beta or beta,
        !        alpha.

        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(:)
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: i,j, ij_k(sys%lattice%ndim), ij_spin
        integer :: ind
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

        r = get_rand_close_open(rng)
        ind = int(r*sys%nel*(sys%nel-1)/2) + 1

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.
        ! Note that the indexing scheme for the strictly lower triangular array
        ! assumes j>i.  As occ_list is ordered, this means that we will return
        ! i,j (referring to spin-orbitals) where j>i.  This ordering is
        ! convenient subsequently, e.g. is assumed in the
        ! find_excitation_permutation2 routine.
        j = int(1.50_p + sqrt(2*ind-1.750_p))
        i = ind - ((j-1)*(j-2))/2

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.
        i = occ_list(i)
        j = occ_list(j)

        ij_spin = sys%basis%basis_fns(i)%Ms + sys%basis%basis_fns(j)%Ms
        ij_k = sys%basis%basis_fns(i)%l + sys%basis%basis_fns(j)%l 

    end subroutine choose_ij_k

!--- Select random orbitals in double excitations ---

    subroutine find_ab_ueg(rng, sys, f, ij_k, ij_spin, ternary_conserve, max_na, a, b, allowed_excitation)

        ! Select a random a (and hence b) to create an allowed excitation,
        ! (i,j)->(a,b), given a certain (i,j) pair.  Simply choose a random
        ! spin-orbital (a) and hence find the spin-orbital (b) which conserves
        ! spin and crystal momentum.  Check if the excitation is actually
        ! allowed (ie there exists an unoccupied orbital for which k_i+k_j-k_a
        ! lies within the Fermi sphere and, if so, that the b spin-orbital is
        ! actually unoccupied).

        ! In:
        !    sys: system object being studied.
        !    f: bit string representation of the Slater determinant from which
        !       the excitation is.
        !    ij_k: sum of the wavevectors of the i and j spin-orbitals.
        !    ij_spin: sum of the spins of the i and j spin-orbitals (in units of
        !       1/2).
        !    ternary_conserve: details of allowed excitations to conserve crystal momentum
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    max_na: the maximum number of choices for the a spin-orbital given
        !       the choice of the (i,j) pair of spin-orbitals (ie the number
        !       of unoccupied spin-orbitals in the basis ignoring occupancy such
        !       that there exists at least one other spin-orbital (either
        !       occupied or unoccupied) which allows spin and crystal momentum
        !       to be conserved.
        !    a,b: Indices of the virtual orbitals selected for the excitation.
        !        Meaningless/unset if allowed_excitation is false.
        !    allowed_excitation: if true, then the excitation is allowed.


        use const, only: i0_length
        use bit_utils, only: count_set_bits, decode_bit_string
        use system, only: sys_t
        use ueg_system, only: ueg_basis_index
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use utils, only: tri_ind

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        integer, intent(in) :: ij_k(sys%lattice%ndim), ij_spin
        integer(i0), intent(in), allocatable :: ternary_conserve(:,:,:,:) ! Declared as allocatable to get correct array bounds
        integer, intent(out) :: a, b, max_na
        logical, intent(out) :: allowed_excitation

        integer :: ibp, ibe, n, ind, kb(sys%lattice%ndim), k3(3)
        ! [todo] rewrite using on bit_string_len, as we don't require the additional info_string_len info.
        integer(i0) :: poss_a(sys%basis%tot_string_len)
        integer :: nposs_a(sys%basis%tot_string_len), poss_a_orbs(i0_length)

        ! Let's just check there are possible a,b first!
        ! Adjust for spin.  Allow a to be up (without bias) if possible.
        k3(1:sys%lattice%ndim) = ij_k
        k3(sys%lattice%ndim+1:) = 0
        if (ij_spin == -2) then
            poss_a = iand(not(f), ishft(ternary_conserve(1:,k3(1),k3(2),k3(3)),1))
        else
            poss_a = iand(not(f), ternary_conserve(1:,k3(1),k3(2),k3(3)))
        end if
        nposs_a = count_set_bits(poss_a)
        max_na = sum(nposs_a)

        if (max_na > 0) then

            ! Select a random a.
            a = int(max_na*get_rand_close_open(rng)) + 1

            n = 0
            finda: do ibe = 1, sys%basis%tot_string_len
                if (n + nposs_a(ibe) >= a) then
                    call decode_bit_string(poss_a(ibe), poss_a_orbs)
                    a = sys%basis%basis_lookup(poss_a_orbs(a-n), ibe)
                    exit finda
                else
                    n = n + nposs_a(ibe)
                end if
            end do finda

            ! Ok, b is now completely defined.
            kb = ij_k - sys%basis%basis_fns(a)%l

            ! Implementation detail: b must be down spin if possible as a is up spin if possible.
            select case(ij_spin)
            case(2)
                b = ueg_basis_index(sys%ueg%basis, kb, 1)
            case(0)
                b = ueg_basis_index(sys%ueg%basis, kb, -1)
            case(-2)
                b = ueg_basis_index(sys%ueg%basis, kb, -1)
            end select

            ! Excitation is forbidden if b is already occupied!
            ibp = sys%basis%bit_lookup(1,b)
            ibe = sys%basis%bit_lookup(2,b)
            allowed_excitation = .not.btest(f(ibe), ibp)

            ! It is useful to return a,b ordered (e.g. for the find_excitation_permutation2 routine).
            if (a > b) then
                ind = a
                a = b
                b = ind
            end if

        else

            allowed_excitation = .false.

        end if

    end subroutine find_ab_ueg

!--- Excitation generation probabilities ---

    pure function calc_pgen_ueg_no_renorm(sys, max_na, ij_spin) result(pgen)

        ! NOTE: this assumes that the corresponding excitation (i,j)->(a,b) is
        ! actually allowed.  Cases where this is not true must be explicitly
        ! handled.

        ! In:
        !    sys: system object being studied.
        !    max_na: the maximum number of choices for the a spin-orbital given
        !       the choice of the (i,j) pair of spin-orbitals (ie the number
        !       of unoccupied spin-orbitals in the basis ignoring occupancy such
        !       that there exists at least one other spin-orbital (either
        !       occupied or unoccupied) which allows spin and crystal momentum
        !       to be conserved.
        ! Returns:
        !    pgen: the generation probability of the excitation.

        use system, only: sys_t

        real(p) :: pgen
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: max_na, ij_spin

        ! We explicitly reject excitations i,j->a,b where b is already
        ! occupied.

        ! Having selected i and j, we could either choose a and then b or vice
        ! versa.  Hence:

        ! pgen = p(i,j) [ p(a|i,j) p(b|i,j,a) + p(b|a,j) p(a|i,j,b) ]

        ! p(i,j) = 1/binom(nel,2) = 2/(nel*(nel-1))
        ! As the UEG is a one-band system, the fourth orbital is completely
        ! defined by the selection of the previous three:
        ! p(b|i,j,a) = p(a|i,j,b) = 1.
        ! p(a|i,j) = max_na, the number of choices for the a orbital (see
        ! find_ab_ueg).

        ! **WARNING**
        ! As an implementation detail, we have required a to be up spin if possible.
        ! This introduces an asymmetry into the excitation generator.
        ! If the excitation is spin parallel (up,up->up,up or
        ! down,down->down,down), then we could have chosen the a orbital first
        ! *or* the b orbital first as both will have been in the same list or
        ! available virtual orbitals.
        ! However, if the excitation is spin anti-parallel, then we could not
        ! have chosen b first.
        ! Hence p(b|i,j) = p(a|i,j) \delta(Ms_a, Ms_b).

        pgen = 2.0_p / ( sys%nel*(sys%nel-1) * max_na )
        if (ij_spin /= 0) pgen = pgen * 2

    end function calc_pgen_ueg_no_renorm

    subroutine init_excit_ueg_power_pitzer(sys, pp)

        ! Initialize the pp data for the gen_excit_ueg_power_pitzer excitation generator.
        ! This creates a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        ! Weight the double excitations according the the Power-Pitzer bound
        ! <ij|ab> <= Sqrt(<ia|ai><jb|bj>). 
        ! Since the value of the integral in the UEG is uniquely determined from the i->a 
        ! excitation (momentum conservation), we can use exactly |<ia|ai>| as the upper 
        ! bound for <ij|ab> here. (Note that because momentum is conserved, 
        ! |q| = | k_i - k_a | = | k_j - k_b | and hence for the (non-zero) integral <ij|ab>, 
        ! <ia|ai> = <jb|bj>.)
        ! This version uses the excitation weights from the reference rather than
        ! recalculating for each determinant.

        ! In:
        !    sys: system object being studied.
        ! Out:
        !    pp: excit_gen_power_pitzer_t which will be initialized with 
        !           the alias tables for the excitations.

        use system, only: sys_t
        use sort, only: qsort
        use excit_gens, only: excit_gen_power_pitzer_t, alloc_alias_table_data_t
        use alias, only: generate_alias_tables
        use hamiltonian_ueg, only: create_weighted_excitation_list_ueg 
        type(sys_t), intent(in) :: sys
        type(excit_gen_power_pitzer_t), intent(inout) :: pp

        integer :: i, j, maxv, nbas
        integer :: vlist(sys%basis%nbasis-1)   !List of 'virtuals' for each orbital

        nbas = sys%basis%nbasis    
        maxv = nbas / 2

        call alloc_alias_table_data_t(pp%pp_ia_d, maxv, nbas)

        do i=1, nbas
            do j=1, maxv     !make a temporary array of 'virtuals' for this occ.
               vlist(j) = j*2 - mod(i,2)      !get the virtual of the right spin
            end do
            call create_weighted_excitation_list_ueg(sys, i, vlist, maxv, pp%pp_ia_d%weights(:,i), pp%pp_ia_d%weights_tot(i))
            call generate_alias_tables(maxv, pp%pp_ia_d%weights(:,i), pp%pp_ia_d%weights_tot(i), pp%pp_ia_d%aliasU(:,i), &
                pp%pp_ia_d%aliasK(:,i)) 
        end do

    end subroutine init_excit_ueg_power_pitzer

    subroutine gen_excit_ueg_power_pitzer(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        ! Weight the double excitations according the the Power-Pitzer bound
        ! <ij|ab> <= Sqrt(<ia|ai><jb|bj>)

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Excitation generation data.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    rng: random number generator.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant in molecular systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation2
        use system, only: sys_t
        use hamiltonian_ueg, only: slater_condon2_ueg_excit
        use excit_gens, only: excit_gen_data_t
        use alias, only: select_weighted_value_precalc
        use ueg_system, only: ueg_basis_index
        use hamiltonian_data 

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(inout) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation

        integer :: ij_k(sys%lattice%ndim), ij_spin

      
        integer :: a, b, i, j, a_ind, b_ind, maxv
        integer :: ibp, ibe
        integer :: kb(sys%lattice%ndim)

        associate( pp => excit_gen_data%excit_gen_pp )

            ! 1. must have a double excitation.
            connection%nexcit = 2

            ! 2. Select orbitals to excite from and into.
            call choose_ij_k(rng, sys, cdet%occ_list, i, j,  ij_k, ij_spin)

            ! We now need to select the orbitals to excite into which we do with weighting:
            ! p(ab|ij) = p(a|i) p(b|j) + p(a|j) p(b|i)

            ! We actually choose a|i then b|j, but since we could have also generated the excitation b from i and a from j, we need to include that prob too.
            
            maxv = sys%basis%nbasis / 2

            ! Just use electron i
            a_ind = select_weighted_value_precalc(rng, maxv, pp%pp_ia_d%aliasU(:,i), pp%pp_ia_d%aliasK(:,i))
            ! Use the alias method to select i with the appropriate probability
            ! Map those >=i to the one after ( we're not allowed to select i)
            ! convert from spatial orbital back to spin orbital

            a = 2*a_ind - mod(i,2)

            ! Excitation is forbidden if a is already occupied!
            ibp = sys%basis%bit_lookup(1,a)
            ibe = sys%basis%bit_lookup(2,a)
            allowed_excitation = .not.btest(cdet%f(ibe), ibp)

            if (allowed_excitation) then
                !b is now fully determined (but we have to work out which one it is)
                ! Ok, b is now completely defined.
                kb = ij_k - sys%basis%basis_fns(a)%l

                select case(ij_spin)
                case(2)
                    b = ueg_basis_index(sys%ueg%basis, kb, 1)
                case(0)
                    b = ueg_basis_index(sys%ueg%basis, kb, -sys%basis%basis_fns(a)%Ms)
                case(-2)
                    b = ueg_basis_index(sys%ueg%basis, kb, -1)
                end select

                if (b<=0) then
                    allowed_excitation = .false.
                else
                !convert to an ind
                    b_ind = (b+1)/2

                    ! Excitation is forbidden if b is already occupied!
                    ibp = sys%basis%bit_lookup(1,b)
                    ibe = sys%basis%bit_lookup(2,b)
                    allowed_excitation = .not.btest(cdet%f(ibe), ibp)
                end if
            end if
            
            if (allowed_excitation) then

                ! 3b. Probability of generating this excitation.

                ! Calculate p(ab|ij) = p(a|i) p(b|ija) + p(b|i)p(a|ijb)
              
                if (ij_spin==0) then 
                    ! Not possible to have chosen the reversed excitation.
                    pgen=pp%pp_ia_d%weights(a_ind,i) / pp%pp_ia_d%weights_tot(i)   ! p(j|b)=1
                else
                    ! i and j have same spin, so could have been selected in the other order.
                    pgen= ( pp%pp_ia_d%weights(a_ind, i) + pp%pp_ia_d%weights(b_ind, i) )  / (pp%pp_ia_d%weights_tot(i))
                end if

                pgen = pgen * 2.0_p/(sys%nel*(sys%nel-1)) ! pgen(ij)

                allowed_excitation = (a/=b)

            end if
            if (allowed_excitation) then
                connection%from_orb(1) = i
                connection%from_orb(2) = j
                if (a<b) then
                    connection%to_orb(1) = a
                    connection%to_orb(2) = b 
                else
                    connection%to_orb(2) = a
                    connection%to_orb(1) = b 
                end if

                ! 4b. Parity of permutation required to line up determinants.
                ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

                ! 5b. Find the connecting matrix element.
                hmatel%r = slater_condon2_ueg_excit(sys, connection%from_orb(1), connection%to_orb(1), connection%to_orb(2), &
                    connection%perm)
            else
                ! Carelessly selected ij with no possible excitations.  Such
                ! events are not worth the cost of renormalising the generation
                ! probabilities.
                ! Return a null excitation.
                hmatel%r = 0.0_p
                pgen = 1.0_p
            end if
        end associate

    end subroutine gen_excit_ueg_power_pitzer

end module excit_gen_ueg
