module ccmc_utils

! Module containing all utility functions only used within ccmc.
! For full explanation see top of ccmc.F90.

use const, only: i0, p, int_p, dp, int_64

implicit none

interface update_noncumulative_dist_int
    module procedure update_noncumulative_dist_int_32
    module procedure update_noncumulative_dist_int_64
end interface

contains

    subroutine find_D0(psip_list, f0, D0_pos)

        ! Find the reference determinant in the list of walkers

        ! In:
        !    psip_list: particle_t object containing current excip distribution on
        !       this processor.
        !    f0: bit string representing the reference.
        ! In/Out:
        !    D0_pos: on input, the position of the reference in
        !       particle_t%states in the previous iteration (or -1 if it was
        !       not on this processor).  On output, the current position.

        use bit_utils, only: bit_str_cmp
        use search, only: binary_search
        use qmc_data, only: particle_t
        use errors, only: stop_all

        type(particle_t), intent(in) :: psip_list
        integer(i0), intent(in) :: f0(:)
        integer, intent(inout) :: D0_pos

        integer :: D0_pos_old
        logical :: hit

        if (D0_pos == -1) then
            ! D0 was just moved to this processor.  No idea where it might be...
            call binary_search(psip_list%states, f0, 1, psip_list%nstates, hit, D0_pos)
        else
            D0_pos_old = D0_pos
            select case(bit_str_cmp(f0, psip_list%states(:,D0_pos)))
            case(0)
                ! D0 hasn't moved.
                hit = .true.
            case(1)
                ! D0 < psip_list%states(:,D0_pos) -- it has moved to earlier in
                ! the list and the old D0_pos is an upper bound.
                call binary_search(psip_list%states, f0, 1, D0_pos_old, hit, D0_pos)
            case(-1)
                ! D0 > psip_list%states(:,D0_pos) -- it has moved to later in
                ! the list and the old D0_pos is a lower bound.
                call binary_search(psip_list%states, f0, D0_pos_old, psip_list%nstates, hit, D0_pos)
            end select
        end if
        if (.not.hit) call stop_all('find_D0', 'Cannot find reference!')

    end subroutine find_D0

    pure subroutine collapse_cluster(basis, f0, excitor, excitor_population, cluster_excitor, cluster_population, allowed)

        ! Collapse two excitors.  The result is returned in-place.

        ! In:
        !    basis: information about the single-particle basis.
        !    f0: bit string representation of the reference determinant.
        !    excitor: bit string of the Slater determinant formed by applying
        !        the excitor, e1, to the reference determinant.
        !    excitor_population: number of excips on the excitor e1.
        ! In/Out:
        !    cluster_excitor: bit string of the Slater determinant formed by applying
        !        the excitor, e2, to the reference determinant.
        !    cluster_population: number of excips on the 'cluster' excitor, e2.
        ! Out:
        !    allowed: true if excitor e1 can be applied to excitor e2 (i.e. e1
        !       and e2 do not involve exciting from/to the any identical
        !       spin-orbitals).

        ! On input, cluster excitor refers to an existing excitor, e2.  On output,
        ! cluster excitor refers to the excitor formed from applying the excitor
        ! e1 to the cluster e2.
        ! ***WARNING***: if allowed is false then cluster_excitor is *not* updated.

        use basis_types, only: basis_t, reset_extra_info_bit_string

        use bit_utils, only: count_set_bits
        use const, only: i0_end

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f0(basis%tot_string_len)
        integer(i0), intent(in) :: excitor(basis%tot_string_len)
        complex(p), intent(in) :: excitor_population
        integer(i0), intent(inout) :: cluster_excitor(basis%tot_string_len)
        complex(p), intent(inout) :: cluster_population
        logical,  intent(out) :: allowed

        integer(i0) :: excitor_loc(basis%tot_string_len)

        integer :: ibasis, ibit
        integer(i0) :: excitor_excitation(basis%tot_string_len)
        integer(i0) :: excitor_annihilation(basis%tot_string_len)
        integer(i0) :: excitor_creation(basis%tot_string_len)
        integer(i0) :: cluster_excitation(basis%tot_string_len)
        integer(i0) :: cluster_annihilation(basis%tot_string_len)
        integer(i0) :: cluster_creation(basis%tot_string_len)
        integer(i0) :: permute_operators(basis%tot_string_len)

        excitor_loc = excitor
        call reset_extra_info_bit_string(basis, excitor_loc)

        ! Apply excitor to the cluster of excitors.

        ! orbitals involved in excitation from reference
        excitor_excitation = ieor(f0, excitor_loc)
        cluster_excitation = ieor(f0, cluster_excitor)
        ! annihilation operators (relative to the reference)
        excitor_annihilation = iand(excitor_excitation, f0)
        cluster_annihilation = iand(cluster_excitation, f0)
        ! creation operators (relative to the reference)
        excitor_creation = iand(excitor_excitation, excitor_loc)
        cluster_creation = iand(cluster_excitation, cluster_excitor)

        ! First, let's find out if the excitor is valid...
        if (any(iand(excitor_creation,cluster_creation) /= 0) &
                .or. any(iand(excitor_annihilation,cluster_annihilation) /= 0)) then
            ! excitor attempts to excite from an orbital already excited from by
            ! the cluster or into an orbital already excited into by the
            ! cluster.
            ! => not valid
            allowed = .false.

            ! We still use the cluster in linked ccmc so need its amplitude
            cluster_population = cluster_population*excitor_population
        else
            ! Applying the excitor to the existing cluster of excitors results
            ! in a valid cluster.
            allowed = .true.

            ! Combine amplitudes.
            ! Might need a sign change as well...see below!
            cluster_population = cluster_population*excitor_population

            ! Now apply the excitor to the cluster (which is, in its own right,
            ! an excitor).
            ! Consider a cluster, e.g. t_i^a = a^+_a a_i (which corresponds to i->a).
            ! We wish to collapse two excitors, e.g. t_i^a t_j^b, to a single
            ! excitor, e.g. t_{ij}^{ab} = a^+_a a^+_b a_j a_i (where i<j and
            ! a<b).  However t_i^a t_j^b = a^+_a a_i a^+_b a_j.  We thus need to
            ! permute the creation and annihilation operators.  Each permutation
            ! incurs a sign change.

            do ibasis = 1, basis%bit_string_len
                do ibit = 0, i0_end
                    if (btest(excitor_excitation(ibasis),ibit)) then
                        if (btest(f0(ibasis),ibit)) then
                            ! Exciting from this orbital.
                            cluster_excitor(ibasis) = ibclr(cluster_excitor(ibasis),ibit)
                            ! We need to swap it with every annihilation
                            ! operator and every creation operator referring to
                            ! an orbital with a higher index already in the
                            ! cluster.
                            ! Note that an orbital cannot be in the list of
                            ! annihilation operators and the list of creation
                            ! operators.
                            ! First annihilation operators:
                            permute_operators = iand(basis%excit_mask(:,basis%basis_lookup(ibit,ibasis)),cluster_annihilation)
                            ! Now add the creation operators:
                            permute_operators = ior(permute_operators,cluster_creation)
                        else
                            ! Exciting into this orbital.
                            cluster_excitor(ibasis) = ibset(cluster_excitor(ibasis),ibit)
                            ! Need to swap it with every creation operator with
                            ! a lower index already in the cluster.
                            permute_operators = iand(not(basis%excit_mask(:,basis%basis_lookup(ibit,ibasis))),cluster_creation)
                            permute_operators(ibasis) = ibclr(permute_operators(ibasis),ibit)
                        end if
                        if (mod(sum(count_set_bits(permute_operators)),2) == 1) &
                            cluster_population = -cluster_population
                    end if
                end do
            end do

        end if

    end subroutine collapse_cluster

    pure subroutine convert_excitor_to_determinant(excitor, excitor_level, excitor_sign, f)

        ! We usually consider an excitor as a bit string representation of the
        ! determinant formed by applying the excitor (a group of annihilation
        ! and creation operators) to the reference determinant; indeed the
        ! determinant form is required when constructing Hamiltonian matrix
        ! elements.  However, the resulting determinant might actually contain
        ! a sign change, which needs to be absorbed into the (signed) population
        ! of excips on the excitor.

        ! This results from the fact that a determinant, |D>, is defined as:
        !   |D> = a^+_i a^+_j ... a^+_k |0>,
        ! where |0> is the vacuum, a^+_i creates an electron in the i-th
        ! orbital, i<j<...<k and |0> is the vacuum.  An excitor is defined as
        !   t_{ij...k}^{ab...c} = a^+_a a^+_b ... a^+_c a_k ... a_j a_i
        ! where i<j<...<k and a<b<...<c.  (This definition is somewhat
        ! arbitrary; the key thing is to be consistent.)  Hence applying an
        ! excitor to the reference might result in a change of sign, i.e.
        ! t_{ij}^{ab} |D_0> = -|D_{ij}^{ab}>.  As a more concrete example,
        ! consider a set of spin states (as the extension to fermions is
        ! irrelevant to the argument), with the reference:
        !   |D_0> = | 1 2 3 >
        ! and the excitor
        !   t_{13}^{58}
        ! Thus, using |0> to denote the vacuum:
        !   t_{13}^{58} | 1 2 3 > = + a^+_5 a^+_8 a_3 a_1 a^+_1 a^+_2 a^+_3 |0>
        !                         = + a^+_5 a^+_8 a_3 a^+_2 a^+_3 |0>
        !                         = - a^+_5 a^+_8 a_3 a^+_3 a^+_2 |0>
        !                         = - a^+_5 a^+_8 a^+_2 |0>
        !                         = + a^+_5 a^+_2 a^+_8 |0>
        !                         = - a^+_2 a^+_5 a^+_8 |0>
        !                         = - | 2 5 8 >
        ! Similarly
        !   t_{12}^{58} | 1 2 3 > = + a^+_5 a^+_8 a_2 a_1 a^+_1 a^+_2 a^+_3 |0>
        !                         = + a^+_5 a^+_8 a_2 a^+_2 a^+_3 |0>
        !                         = + a^+_5 a^+_8 a^+_3 |0>
        !                         = - a^+_5 a^+_3 a^+_8 |0>
        !                         = + a^+_3 a^+_5 a^+_8 |0>
        !                         = + | 3 5 8 >

        ! This potential sign change must be taken into account; we do so by
        ! absorbing the sign into the signed population of excips on the
        ! excitor.

        ! Essentially taken from Alex Thom's original implementation.

        ! In:
        !    excitor: bit string of the Slater determinant formed by applying
        !        the excitor to the reference determinant.
        !    excitor_level: excitation level, relative to the determinant f,
        !        of the excitor.  Equal to the number of
        !        annihilation (or indeed creation) operators in the excitor.
        !    f: bit string of determinant to which excitor is
        !       applied to generate a new determinant.
        ! Out:
        !    excitor_sign: sign due to applying the excitor to the
        !       determinant f to form a Slater determinant, i.e. < D_i | a_i D_f >,
        !       which is +1 or -1, where D_i is the determinant formed from
        !       applying the cluster of excitors, a_i, to D_f

        use const, only: i0_end

        integer(i0), intent(in) :: excitor(:)
        integer, intent(in) :: excitor_level
        integer, intent(inout) :: excitor_sign
        integer(i0), intent(in) :: f(:)

        integer(i0) :: excitation(size(excitor))
        integer :: ibasis, ibit, ncreation, nannihilation

        ! Bits involved in the excitation from the reference determinant.
        excitation = ieor(f, excitor)

        nannihilation = excitor_level
        ncreation = excitor_level

        excitor_sign = 1

        ! Obtain sign change by (implicitly) constructing the determinant formed
        ! by applying the excitor to the reference determinant.
        do ibasis = 1, size(excitor)
            do ibit = 0, i0_end
                if (btest(f(ibasis),ibit)) then
                    ! Occupied orbital in reference.
                    if (btest(excitation(ibasis),ibit)) then
                        ! Orbital excited from...annihilate electron.
                        ! This amounts to one fewer operator in the cluster through
                        ! which the other creation operators in the determinant must
                        ! permute.
                        nannihilation = nannihilation - 1
                    else
                        ! Orbital is occupied in the reference and once the
                        ! excitor has been applied.
                        ! Permute the corresponding creation operator through
                        ! the remaining creation and annihilation operators of
                        ! the excitor (which operate on orbitals with a higher
                        ! index than the current orbital).
                        ! If the permutation is odd, then we incur a sign
                        ! change.
                        if (mod(nannihilation+ncreation,2) == 1) &
                            excitor_sign = -excitor_sign
                    end if
                else if (btest(excitation(ibasis),ibit)) then
                    ! Orbital excited to...create electron.
                    ! This amounts to one fewer operator in the cluster through
                    ! which the creation operators in the determinant must
                    ! permute.
                    ! Note due to definition of the excitor, it is guaranteed
                    ! that this is created in the correct place, ie there are no
                    ! other operators in the excitor it needs to be interchanged
                    ! with.
                    ncreation = ncreation - 1
                end if
            end do
        end do

    end subroutine convert_excitor_to_determinant

    subroutine zero_estimators_t(estimators)

        use qmc_data, only: estimators_t

        type(estimators_t), intent(inout) :: estimators

        estimators%D0_population = 0.0_p
        estimators%proj_energy = 0.0_p
        estimators%D0_population_comp = cmplx(0.0, 0.0, p)
        estimators%proj_energy_comp = cmplx(0.0, 0.0, p)

    end subroutine zero_estimators_t

    subroutine cumulative_population(pops, ex_lvls, nactive, D0_proc, D0_pos, real_factor, calc_dist, complx, &
                                    cumulative_pops, tot_pop, ex_lvl_dist)

        ! Calculate the cumulative population, i.e. the number of psips/excips
        ! residing on a determinant/an excitor and all determinants/excitors which
        ! occur before it in the determinant/excitor list.

        ! This is primarily so in CCMC we can select clusters of excitors with each
        ! excip being equally likely to be part of a cluster.  (If we just select
        ! each occupied excitor with equal probability, then we get wildy
        ! fluctuating selection probabilities and hence large population blooms.)
        ! As 'excips' on the reference cannot be part of a cluster, then the
        ! population on the reference is treated as 0 if required.

        ! In:
        !    pops: list of populations on each determinant/excitor.  Must have
        !       minimum length of nactive.
        !    nactive: number of occupied determinants/excitors (ie pops(:,1:nactive)
        !       contains the population(s) on each currently "active"
        !       determinant/excitor.
        !    D0_proc: processor on which the reference resides.
        !    D0_pos: position in the pops list of the reference.  Only relevant if
        !       1<=D0_pos<=nactive and the processor holds the reference.
        !    real_factor: the encoding factor by which the stored populations are multiplied
        !       to enable non-integer populations.
        !    calc_dist: whether to update excitation level distributions. NB. this is only
        !       possible if the walker list is adjusted to be sorted by excitation level.
        ! Out:
        !    cumulative_pops: running total of excitor population, i.e.
        !        cumulative_pops(i) = sum(abs(pops(1:i))), excluding the
        !        population on the reference if appropriate.
        !    tot_pop: total population (possibly excluding the population on the
        !       reference).
        !    ex_lvl_dist: derived types containing distributions of states
        !       and populations between different excitation levels.

        ! NOTE: currently only the populations in the first psip/excip space are
        ! considered.  This should be changed if we do multiple simulations at
        ! once/Hellmann-Feynman sampling/etc.

        use parallel, only: iproc
        use ccmc_data, only: ex_lvl_dist_t

        integer(int_p), intent(in) :: pops(:,:), real_factor
        integer, intent(in) :: nactive, D0_proc, D0_pos
        real(p), allocatable, intent(inout) :: cumulative_pops(:)
        real(p), intent(out) :: tot_pop
        type(ex_lvl_dist_t), intent(inout) :: ex_lvl_dist
        logical, intent(in) :: complx, calc_dist
        integer(i0), intent(in) :: ex_lvls(:)

        integer :: i
        integer(i0) :: j, ex_lvl


        ! First need to set values to account for reference correctly in ex_lvl_dist.
        if (calc_dist) then
            if (D0_proc==iproc) then
                ex_lvl_dist%cumulative_nstates_ex_lvl(0:ex_lvls(1)) = 1
            else
                ex_lvl_dist%cumulative_nstates_ex_lvl(0:ex_lvls(1)) = 0
            end if
            ex_lvl_dist%cumulative_pop_ex_lvl(0) = 0
            ex_lvl = 0_i0
        end if

        ! Need to combine spaces if doing complex; we choose combining in quadrature.
        cumulative_pops(1) = get_pop_contrib(pops(:,1), real_factor, complx)
        if (D0_proc == iproc) then
            ! Let's be a bit faster: unroll loops and skip over the reference
            ! between the loops.
            do i = 2, d0_pos-1
                cumulative_pops(i) = cumulative_pops(i-1) + &
                                        get_pop_contrib(pops(:,i), real_factor, complx)
                if (calc_dist) then
                    if (ex_lvls(i-1) < ex_lvls(i)) then
                        ! i-1 is last entry of it's excitation level.
                        ! Need to update all intervening excitation levels
                         do j = ex_lvls(i-1), ex_lvls(i) - 1
                            ex_lvl_dist%cumulative_nstates_ex_lvl(j) = i-1
                            ex_lvl_dist%cumulative_pop_ex_lvl(j) = cumulative_pops(i-1)
                        end do
                        ex_lvl = ex_lvls(i)
                    end if
                end if
            end do
            ! Set cumulative on the reference to be the running total merely so we
            ! can continue accessing the running total from the i-1 element in the
            ! loop over excitors in slots above the reference.
            if (d0_pos == 1) then
                cumulative_pops(d0_pos) = 0
            end if
            if (d0_pos > 1) cumulative_pops(d0_pos) = cumulative_pops(d0_pos-1)
            do i = d0_pos+1, nactive
                cumulative_pops(i) = cumulative_pops(i-1) + &
                                        get_pop_contrib(pops(:,i), real_factor, complx)
                if (calc_dist) then
                    if (ex_lvls(i-1) < ex_lvls(i)) then
                        ! i-1 is last entry of it's excitation level.
                        ! Need to update all intervening excitation levels
                         do j = ex_lvls(i-1), ex_lvls(i) - 1
                            ex_lvl_dist%cumulative_nstates_ex_lvl(j) = i-1
                            ex_lvl_dist%cumulative_pop_ex_lvl(j) = cumulative_pops(i-1)
                        end do
                        ex_lvl = ex_lvls(i)
                    end if
                end if
            end do
        else
            ! V simple on other processors: no reference to get in the way!
            do i = 2, nactive
                cumulative_pops(i) = cumulative_pops(i-1) + &
                                        get_pop_contrib(pops(:,i), real_factor, complx)
                if (calc_dist) then
                    if (ex_lvls(i-1) < ex_lvls(i)) then
                        ! i-1 is last entry of it's excitation level.
                        ! Need to update all intervening excitation levels
                         do j = ex_lvls(i-1), ex_lvls(i) - 1
                            ex_lvl_dist%cumulative_nstates_ex_lvl(j) = i-1
                            ex_lvl_dist%cumulative_pop_ex_lvl(j) = cumulative_pops(i-1)
                        end do
                        ex_lvl = ex_lvls(i)
                    end if
                end if
            end do
        end if
        if (nactive > 0) then
            tot_pop = cumulative_pops(nactive)
        else
            tot_pop = 0.0_p
        end if

        if (calc_dist) then
            ex_lvl_dist%cumulative_nstates_ex_lvl(ex_lvl:) = nactive
            ex_lvl_dist%cumulative_pop_ex_lvl(ex_lvl:) = tot_pop
        end if

    end subroutine cumulative_population

    pure function get_pop_contrib(pops, real_factor, complx) result(contrib)

        ! Get contribution from a given position in a population list,
        ! given real or complex.

        ! In:
        !    pops: population on determinant/excitor position.
        !    real_factor: the encoding factor by which the stored populations
        !       are multiplied to enable non-integer populations.
        !    complx: logical, true if complex populations, false otherwise.
        ! Out:
        !    contrib: contribution to cumulative population accumulation from
        !       given position in pop list.

        integer(int_p), intent(in) :: pops(:), real_factor
        logical, intent(in) :: complx

        real(p) :: contrib

        if (complx) then
            contrib = abs(cmplx(pops(1), pops(2), p))/real(real_factor, p)
        else
            contrib = abs(real(pops(1), p))/real(real_factor, p)
        end if

    end function get_pop_contrib

    subroutine init_contrib(sys, cluster_size, linked, contrib_info)

        ! Allocates cdet and cluster, and their components.

        ! In:
        !    sys: system being studied
        !    cluster_size: the maximum number of excitors allowed in a cluster
        !    linked: whether we are performing propogation through the linked
        !       coupled cluster equations.
        ! Out:
        !    contrib_info: Array of wfn_contrib_t variables, one for each thread,
        !       with appropriate cluster and det components allocated.

        use parallel, only: nthreads
        use determinants, only: alloc_det_info_t
        use ccmc_data, only: wfn_contrib_t
        use system, only: sys_t
        use checking, only: check_allocate

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: cluster_size
        logical, intent(in) :: linked
        type(wfn_contrib_t), allocatable, intent(out) :: contrib_info(:)

        integer :: i, ierr
        integer :: cluster_size_loc

        cluster_size_loc = cluster_size
        if (linked) cluster_size_loc = 4

        ! Allocate arrays
        allocate(contrib_info(0:nthreads-1), stat=ierr)
        call check_allocate('contrib_info', nthreads, ierr)

        do i = 0, nthreads-1
            ! Allocate det_info_t and cluster_t components
            call alloc_det_info_t(sys, contrib_info(i)%cdet)
            allocate(contrib_info(i)%cluster%excitors(cluster_size_loc), stat=ierr)
            call check_allocate('contrib_info%cluster%excitors', cluster_size_loc, ierr)
            ! Only require right/left components if performing linked propogation.
            if (linked) then
                call alloc_det_info_t(sys, contrib_info(i)%ldet)
                call alloc_det_info_t(sys, contrib_info(i)%rdet)

                allocate(contrib_info(i)%left_cluster%excitors(cluster_size_loc), stat=ierr)
                allocate(contrib_info(i)%right_cluster%excitors(cluster_size_loc), stat=ierr)
                call check_allocate('contrib_info%left_cluster%excitors', cluster_size_loc, ierr)
                call check_allocate('contrib_info%right_cluster%excitors', cluster_size_loc, ierr)
            end if
        end do

    end subroutine init_contrib

    subroutine dealloc_contrib(contrib, linked)

        ! Deallocates cdet and cluster, and their components.

        ! In:
        !    linked: whether we are performing propogation through the linked
        !       coupled cluster equations.
        ! In/Out:
        !    contrib_info: Array of wfn_contrib_t variables, one for each thread,
        !       with all components deallocated.

        use ccmc_data, only: wfn_contrib_t
        use checking, only: check_deallocate
        use determinants, only: dealloc_det_info_t

        type(wfn_contrib_t), allocatable, intent(inout) :: contrib(:)
        logical, intent(in) :: linked
        integer :: i, ierr

        do i = lbound(contrib,dim=1), ubound(contrib, dim=1)
            call dealloc_det_info_t(contrib(i)%cdet)
            deallocate(contrib(i)%cluster%excitors, stat=ierr)
            call check_deallocate('contrib%cluster%excitors', ierr)
            if (linked) then
                call dealloc_det_info_t(contrib(i)%ldet)
                call dealloc_det_info_t(contrib(i)%rdet)
                deallocate(contrib(i)%left_cluster%excitors, stat=ierr)
                call check_deallocate('contrib%left_cluster%excitors', ierr)
                deallocate(contrib(i)%right_cluster%excitors, stat=ierr)
                call check_deallocate('contrib%right_cluster%excitors', ierr)
            end if
        end do

        deallocate(contrib, stat=ierr)
        call check_deallocate('contrib', ierr)

    end subroutine dealloc_contrib

    subroutine init_ex_lvl_dist_t(max_ex_lvl, ex_lvl_dist)

        use ccmc_data, only: ex_lvl_dist_t
        use checking, only: check_allocate

        integer, intent(in) :: max_ex_lvl
        type(ex_lvl_dist_t), intent(inout) :: ex_lvl_dist

        integer :: ierr

        associate(eld => ex_lvl_dist)
            allocate(eld%nstates_ex_lvl(0:max_ex_lvl), stat=ierr)
            call check_allocate('eld%nstates_ex_lvl', max_ex_lvl+1, ierr)

            allocate(eld%cumulative_nstates_ex_lvl(0:max_ex_lvl), stat=ierr)
            call check_allocate('eld%cumulative_nstates_ex_lvl', max_ex_lvl+1, ierr)

            allocate(eld%pop_ex_lvl(0:max_ex_lvl), stat=ierr)
            call check_allocate('eld%pop_ex_lvl', max_ex_lvl+1, ierr)

            allocate(eld%cumulative_pop_ex_lvl(0:max_ex_lvl), stat=ierr)
            call check_allocate('eld%cumulative_pop_ex_lvl', max_ex_lvl+1, ierr)
        end associate

    end subroutine init_ex_lvl_dist_t

    subroutine end_ex_lvl_dist_t(ex_lvl_dist)

        use ccmc_data, only: ex_lvl_dist_t
        use checking, only: check_deallocate

        type(ex_lvl_dist_t), intent(inout) :: ex_lvl_dist

        integer :: ierr

        associate(eld => ex_lvl_dist)
            deallocate(eld%nstates_ex_lvl, stat=ierr)
            call check_deallocate('eld%nstates_ex_lvl', ierr)

            deallocate(eld%cumulative_nstates_ex_lvl, stat=ierr)
            call check_deallocate('eld%cumulative_nstates_ex_lvl', ierr)

            deallocate(eld%pop_ex_lvl, stat=ierr)
            call check_deallocate('eld%pop_ex_lvl', ierr)

            deallocate(eld%cumulative_pop_ex_lvl, stat=ierr)
            call check_deallocate('eld%cumulative_pop_ex_lvl', ierr)
        end associate

    end subroutine end_ex_lvl_dist_t

    subroutine add_ex_level_bit_string_calc(basis, f0, f)

        ! Sets bits within bit string to give excitation level at end of bit strings.
        ! This routine sets ex level from provided reference.

        use excitations, only: get_excitation_level
        use basis_types, only: basis_t

        type(basis_t), intent(in) :: basis
        integer(i0), intent(inout) :: f(:)
        integer(i0), intent(in) :: f0(:)

        integer(i0) :: ex_lvl

        if (basis%info_string_len/=0) then
            ex_lvl = int(get_excitation_level(f, f0), kind=i0)

            f(basis%bit_string_len+1) = ex_lvl
        end if

    end subroutine add_ex_level_bit_string_calc

    subroutine add_ex_level_bit_string_provided(basis, ex_lvl, f)

        ! Sets bits within bit string to give excitation level at end of bit strings.
        ! This routine uses a provided excitation level.

        use basis_types, only: basis_t

        type(basis_t), intent(in) :: basis
        integer, intent(in) :: ex_lvl
        integer(i0) :: ex_lvl_loc

        integer(i0), intent(inout) :: f(:)

        if (basis%info_string_len/=0) then
            ex_lvl_loc = int(ex_lvl, kind=i0)

            f(basis%bit_string_len+1) = ex_lvl
        end if

    end subroutine add_ex_level_bit_string_provided

    subroutine update_ex_lvl_dist(ex_lvl_dist)

        ! Takes ex_lvl_dist derived type with updated cumulative distributions and
        ! sets non-cumulative distributions correspondingly.

        use ccmc_data, only: ex_lvl_dist_t

        type(ex_lvl_dist_t), intent(inout) :: ex_lvl_dist

        call update_noncumulative_dist_int(ex_lvl_dist%cumulative_nstates_ex_lvl, ex_lvl_dist%nstates_ex_lvl)
        call update_noncumulative_dist_realp(ex_lvl_dist%cumulative_pop_ex_lvl, ex_lvl_dist%pop_ex_lvl)

    end subroutine update_ex_lvl_dist

! ---- Helper functions to re-construct cumulative and non-cumulative distributions from each other ----

    subroutine update_cumulative_dist_real(dist, cumulative_dist, normalised)

        ! Updates cumulative_dist to be consistent with provided dist.
        ! If normalised is set the distribution is set to have total amplitude
        ! of 1.0_dp by adjusting the final non-zero entry of dist.

        use const, only: depsilon, dp

        real(dp), intent(inout), allocatable :: dist(:)
        real(dp), intent(inout), allocatable :: cumulative_dist(:)
        logical, intent(in) :: normalised
        integer :: i
        real(dp) :: diff

        if (normalised .and. sum(dist) > depsilon) dist = dist / sum(dist)

        cumulative_dist(lbound(cumulative_dist,dim=1)) = dist(lbound(dist,dim=1))
        do i = lbound(dist, dim=1) + 1, ubound(dist, dim=1)
            cumulative_dist(i) = cumulative_dist(i-1) + dist(i)
        end do

        if (normalised) then
            diff = 1.0_dp - cumulative_dist(ubound(dist,dim=1))
            do i = ubound(dist, dim=1), lbound(dist,dim=1), -1
                if (dist(i) > depsilon) then
                    dist(i) = dist(i) + diff
                    cumulative_dist(i:) = 1.0_dp
                    exit
                end if
            end do
        end if

    end subroutine update_cumulative_dist_real

    subroutine update_noncumulative_dist_int_32(cumulative_dist, dist)

        ! Uses the provided cumulative distribution to update dist in a
        ! consistent manner.

        use const, only: int_32

        integer(int_32), intent(inout), allocatable :: dist(:)
        integer(int_32), intent(in), allocatable :: cumulative_dist(:)
        integer :: i

        dist(lbound(cumulative_dist,dim=1)) = cumulative_dist(lbound(dist,dim=1))
        do i = lbound(cumulative_dist, dim=1) + 1, ubound(cumulative_dist, dim=1)
            dist(i) = cumulative_dist(i) - cumulative_dist(i-1)
        end do

    end subroutine update_noncumulative_dist_int_32

    subroutine update_noncumulative_dist_int_64(cumulative_dist, dist)

        ! Uses the provided cumulative distribution to update dist in a
        ! consistent manner.

        use const, only: int_64

        integer(int_64), intent(inout), allocatable :: dist(:)
        integer(int_64), intent(inout), allocatable :: cumulative_dist(:)
        integer :: i

        dist(lbound(cumulative_dist,dim=1)) = cumulative_dist(lbound(dist,dim=1))
        do i = lbound(cumulative_dist, dim=1) + 1, ubound(cumulative_dist, dim=1)
            dist(i) = cumulative_dist(i) - cumulative_dist(i-1)
        end do

    end subroutine update_noncumulative_dist_int_64

    subroutine update_noncumulative_dist_realp(cumulative_dist, dist)

        ! Uses the provided cumulative distribution to update dist in a
        ! consistent manner.

        use const, only: int_64

        real(p), intent(inout), allocatable :: dist(:)
        real(p), intent(inout), allocatable :: cumulative_dist(:)
        integer :: i

        dist(lbound(cumulative_dist,dim=1)) = cumulative_dist(lbound(dist,dim=1))
        do i = lbound(cumulative_dist, dim=1) + 1, ubound(cumulative_dist, dim=1)
            dist(i) = cumulative_dist(i) - cumulative_dist(i-1)
        end do

    end subroutine update_noncumulative_dist_realp

    subroutine regenerate_ex_levels_psip_list(basis, qs)

        ! Regenerates excitation level information stored at start of bit string
        ! within states in psip list. For use when restarting from a restart file
        ! not containing this information.
        ! Also sorts the list, as ordering will change.
        ! Hashing only uses nbasis bits, so should be unaffected by the additional
        ! information at the start of the bit string.
        ! [todo] figure out a way to double check this is the case.

        ! In:
        !   basis: information on single-particle basis in use.
        ! In/Out:
        !   qmc_state: information on current state of calculation. We update and
        !       reorder the bit strings within the psip list, using the reference
        !       determinant bit string stored within qs%ref%f0.

        use basis_types, only: basis_t
        use qmc_data, only: qmc_state_t
        use sort, only: qsort

        type(basis_t), intent(in) :: basis
        type(qmc_state_t), intent(inout) :: qs

        integer :: istate

        do istate = 1, qs%psip_list%nstates
            call add_ex_level_bit_string_calc(basis, qs%ref%f0, qs%psip_list%states(:,istate))
        end do

        associate(pl=>qs%psip_list)
             call qsort(pl%nstates, pl%states, pl%pops, pl%dat)
        end associate


    end subroutine regenerate_ex_levels_psip_list

end module ccmc_utils
