module spawning

! Module for procedures involved in the spawning step of the FCIQMC algorithm.

use const
implicit none

contains

    subroutine spawn_hub_k(spos, cdet)

        ! Attempt to spawn a new particle on a connected determinant.
        ! In:
        !    spos: current position within the spawned walkers array.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.

        use determinants, only: det_info
        use dSFMT_interface, only:  genrand_real2
        use excitations, only: calc_pgen_hub_k
        use fciqmc_data, only: tau
        use system, only: hub_k_coulomb
        use symmetry, only: inv_sym

        integer, intent(in) :: spos
        type(det_info), intent(in) :: cdet
        real(dp) :: pgen, psuccess, pspawn
        integer :: i, j, a, b, ij_sym, ab_sym, nparticles

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(cdet%occ_list_alpha, cdet%occ_list_beta, i ,j, ij_sym)

        ! 2. Calculate the generation probability of the excitation.
        ! For two-band systems this depends only upon the orbitals excited from.
        pgen = calc_pgen_hub_k(ij_sym, cdet%f, cdet%unocc_list_alpha, cdet%unocc_list_beta)

        ! The excitation i,j -> a,b is only allowed if k_i + k_j - k_a - k_b = 0
        ! (up to a reciprocal lattice vector).  We store k_i + k_j as ij_sym, so
        ! k_a + k_b must be the inverse of this.
        ab_sym = inv_sym(ij_sym)

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
        ! We wish to spawn with probability
        !   tau |H_ij|,
        ! where tau is the timestep of the simulation.
        ! The probability of doing something is the probability of selecting to
        ! attempt to do it multiplied by the probability of actually doing it,
        ! hence:
        !   p_spawn = p_select tau*|H_ij|/p_gen
        ! and p_select = p_gen for normalised probabilities.
        ! p_select is included intrinsically in the algorithm, so we compare
        ! to the probability tau*|H_ij|/p_gen.
        pspawn = tau*hub_k_coulomb/pgen
        psuccess = genrand_real2()

        ! Need to take into account the possibilty of a spawning attempt
        ! producing multiple offspring...
        ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and 
        ! then spawn a particle with probability pspawn-floor(pspawn).
        nparticles = int(pspawn)
        pspawn = pspawn - nparticles

        if (pspawn > psuccess) nparticles = nparticles + 1

        if (nparticles > 0) then
            ! 4. Well, I suppose we should find out which determinant we're spawning
            ! on...
            call choose_ab_hub_k(cdet%f, cdet%unocc_list_alpha, cdet%unocc_list_beta, ab_sym, a, b)

            ! 5. Is connecting matrix element positive (in which case we spawn with
            ! negative walkers) or negative (in which case we spawn with positive
            ! walkers)?
            if (mod(i-a,2) == 0) then
                ! i,a are both alpha spin-orbitals or both beta spin-orbitals.
                ! Thus the connecting matrix element comes from the Coulomb
                ! integral and is positive.
                nparticles = -nparticles
            end if

            ! 6. Set info in spawning array.

        end if
        
    end subroutine spawn_hub_k

    subroutine choose_ij(occ_list, i ,j, ij_sym, ij_spin)

        ! Randomly choose a pair of spin-orbitals.
        ! See choose_ij_hub_k for a specific procedure for the momentum space
        ! formulation of the hubbard model.
        ! In:
        !    occ_list: Integer list of occupied spin-orbitals.
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
        i = int(1.50_dp + sqrt(2*ind-1.750_dp))
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
        ! occupied beta array.  This means that:
        !   * i and j can be identical.
        !   * we need to handle the case where nalpha/=nbeta (in choose_ij we
        !     just consider a single list of nel spin-orbitals, whereas here
        !     we need to consider a pick once from the list of of nalpha
        !     spin-orbitals and once form the list of nbeta spin-orbitals.
        !   * the number of possible unique i,j pairs is given by
        !       n1*(n1+1)/2 + n1*(n2-n1)
        !     where n1 = min(nalpha,nbeta) and n2 = max(nalpha,nbeta).
        !     The first term accounts for the "triangular" part of the array,
        !     The second part accounts for the case where we can choose
        !     a j which is greater than the maximum i.

        r = genrand_real2()

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.
        if (nalpha == nbeta) then
            ! Identical case to that given by Rifkin apart from the fact that
            ! our indices go from 1 rather than 0.
            ind = int(r*nalpha*(nalpha+1)/2) + 1
            call decode_tri_ind(ind, i, j)
        else if (nalpha > nbeta) then
            ind = int(r*( nbeta*(nbeta+1)/2 + nbeta*(nalpha-nbeta)))
            if (ind <= nbeta*(nbeta+1)/2) then
                ! "triangular" part of array, where i,j can be generated either
                ! way round.
                call decode_tri_ind(ind, i, j)
            else
                ! Rectangular part.
                call decode_rect_ind(ind, nbeta, i, j)
            end if
        else
            ! Same as nalpha > nbeta case but with labels the other way round.
            ind = int(r*( nalpha*(nalpha+1)/2 + nalpha*(nbeta-nalpha)))
            if (ind <= nalpha*(nalpha+1)/2) then
                ! "triangular" part of array, where i,j can be generated either
                ! way round.
                call decode_tri_ind(ind, j, i)
            else
                ! Rectangular part.
                call decode_rect_ind(ind, nalpha, j, i)
            end if
        end if

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.
        i = occ_list_alpha(i)
        j = occ_list_beta(j)

        ! Symmetry info is a simple lookup...
        ij_sym = sym_table((i+1)/2,(j+1)/2)

        contains

            ! The following are routines for decoding the combined index of an array
            ! like:
            !   1  .  .          1,1
            !   2  3  .          1,2  2,2
            !   3  4  5    to    1,3  2,3  3,3
            !   6  7  8          1,4  2,4  3,4
            !   9 10 11          1,5  2,5  3,5
            ! ie we can't distinguish between i and j apart from the ranges:
            !   1 <= i <= 5
            !   1 <= j <= min(i,3)
            ! the triangular part is decoded in decode_tri_ind, the rectangular part
            ! in decode_rect_ind.

            subroutine decode_tri_ind(ind, i, j)

                ! Decode a combined index.
                ! The combined index, p, is given by:
                !   p = i(i-1)/2 + j
                ! Following the procedure by Rifkin, it is easy to show:
                !   i = int(0.5 + sqrt( 2p - 7/4 )
                ! and j hence follows.
                !
                ! In:
                !    ind: combined index based upon i,j.
                ! Out:
                !    i: i is in the range 1 <= i.
                !    j: j is in the range 1 <= j <= i.

                integer, intent(in) :: ind
                integer, intent(out) :: i, j

                i = int(0.50_dp + sqrt(2*ind-1.750_dp))
                j = ind - (i*(i-1))/2

            end subroutine decode_tri_ind

            subroutine decode_rect_ind(ind, n, i, j)

                ! Decode a combined index.
                ! The combined index, p, is given by:
                !   p = p = n(n+1)/2 + (i-n-1)*n + j
                ! Following the procedure by Rifkin, it is easy to show:
                !   i = int[ 1/n (p - n(n+1)/2 - 1 )] + n + 1
                ! and j hence follows.
                !
                ! In:
                !    ind: combined index based upon i,j.
                !    n: max value of j.
                ! Out:
                !    i: i is in the range n+1 <= i.
                !    j: j is in the range 1 <= j <= n.

                integer, intent(in) :: ind, n
                integer, intent(out) :: i, j

                i = int( (ind - (n*(n+1))/2 - 1.0_dp)/n ) + n + 1
                j = ind - (i*(i-1))/2

            end subroutine decode_rect_ind

    end subroutine choose_ij_hub_k

    subroutine choose_ab_hub_k(f, unocc_list_alpha, unocc_list_beta, ab_sym, a, b)

        ! Choose a random pair of (a,b) unoccupied virtual spin-orbitals into
        ! which electrons are excited.
        ! (a,b) are chosen such that the (i,j)->(a,b) excitation is symmetry-
        ! allowed.
        ! In: 
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        !    unocc_alpha, unocc_beta: integer list of the unoccupied alpha and
        !        beta (respectively) spin-orbitals.
        !    ab_sym: symmetry spanned by the (a,b) combination of unoccupied
        !        spin-orbitals into which electrons are excited.
        ! Returns:
        !    a,b: virtual spin orbitals involved in the excitation.

        use basis, only: basis_length, bit_lookup, nbasis
        use dSFMT_interface, only:  genrand_real2
        use system, only: nvirt_alpha, nvirt_beta
        use symmetry, only: sym_table

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: unocc_list_alpha(nvirt_alpha), unocc_list_beta(nvirt_beta)
        integer, intent(in) :: ab_sym
        integer, intent(out) :: a, b

        integer :: r, b_pos, b_el, ka

        do

            ! Until we find an allowed excitation.

            r = int(genrand_real2()*(nvirt_alpha+nvirt_beta)) + 1

            if (r <= nvirt_alpha) then

                a = unocc_list_alpha(r)
                ! Find corresponding beta orbital which satisfies conservation
                ! of crystal momentum.
                ka = (a+1)/2
                b = 2*sym_table(ab_sym, ka)

            else

                a = unocc_list_beta(r)
                ! Find corresponding alpha orbital which satisfies conservation
                ! of crystal momentum.
                ka = a/2
                b = 2*sym_table(ab_sym, ka) - 1

            end if

            b_pos = bit_lookup(1,b)
            b_el = bit_lookup(2,b)

            ! If b is unoccupied then have found the excitation.
            if (.not.btest(f(b_el), b_pos)) exit

        end do

    end subroutine choose_ab_hub_k

end module spawning
