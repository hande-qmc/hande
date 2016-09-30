module ccmc_utils

! Module containing all utility functions only used within ccmc.
! For full explanation see top of ccmc.F90.

use const, only: i0, p

implicit none

contains

    subroutine init_cluster(sys, cluster_size, cdet, cluster)

        ! Allocates cdet and cluster, and their components.

        ! In:
        !    sys: system being studied
        !    cluster_size: the maximum number of excitors allowed in a cluster
        ! Out:
        !    cdet: Array of det_info_t variables, one for each thread, with
        !       components allocated
        !    cluster: Array of cluster_t variables, one for each thread, with
        !       components allocated

        use parallel, only: nthreads
        use determinants, only: det_info_t, alloc_det_info_t
        use ccmc_data, only: cluster_t
        use system, only: sys_t
        use checking, only: check_allocate

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: cluster_size
        type(det_info_t), allocatable, intent(out) :: cdet(:)
        type(cluster_t), allocatable, intent(out) :: cluster(:)

        integer :: i, ierr

        ! Allocate arrays
        allocate(cdet(0:nthreads-1), stat=ierr)
        call check_allocate('cdet', nthreads, ierr)
        allocate(cluster(0:nthreads-1), stat=ierr)
        call check_allocate('cluster', nthreads, ierr)

        do i = 0, nthreads-1
            ! Allocate det_info_t and cluster_t components
            call alloc_det_info_t(sys, cdet(i))
            allocate(cluster(i)%excitors(cluster_size), stat=ierr)
            call check_allocate('cluster%excitors', cluster_size, ierr)
        end do

    end subroutine init_cluster

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

        use basis_types, only: basis_t

        use bit_utils, only: count_set_bits
        use const, only: i0_end

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f0(basis%string_len)
        integer(i0), intent(in) :: excitor(basis%string_len)
        real(p), intent(in) :: excitor_population
        integer(i0), intent(inout) :: cluster_excitor(basis%string_len)
        real(p), intent(inout) :: cluster_population
        logical,  intent(out) :: allowed

        integer :: ibasis, ibit
        integer(i0) :: excitor_excitation(basis%string_len)
        integer(i0) :: excitor_annihilation(basis%string_len)
        integer(i0) :: excitor_creation(basis%string_len)
        integer(i0) :: cluster_excitation(basis%string_len)
        integer(i0) :: cluster_annihilation(basis%string_len)
        integer(i0) :: cluster_creation(basis%string_len)
        integer(i0) :: permute_operators(basis%string_len)

        ! Apply excitor to the cluster of excitors.

        ! orbitals involved in excitation from reference
        excitor_excitation = ieor(f0, excitor)
        cluster_excitation = ieor(f0, cluster_excitor)
        ! annihilation operators (relative to the reference)
        excitor_annihilation = iand(excitor_excitation, f0)
        cluster_annihilation = iand(cluster_excitation, f0)
        ! creation operators (relative to the reference)
        excitor_creation = iand(excitor_excitation, excitor)
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

            do ibasis = 1, basis%string_len
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

end module ccmc_utils
