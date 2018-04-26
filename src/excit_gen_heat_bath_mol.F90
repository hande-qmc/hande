module excit_gen_heat_bath_mol

! A module containing excitations generators which weight excitations according to the heat bath algorithm
! by Holmes et al. (Holmes, A. A.; Changlani, H. J.; Umrigar, C. J. J. Chem. Theory Comput. 2016, 12, 1561–1571).

use const, only: i0, p, depsilon

implicit none

contains

    subroutine init_excit_mol_heat_bath(sys, hb, original)

        ! Generate excitation tables for all spinorbitals for the gen_excit_mol_heat_bath
        ! excitation generator.

        ! In:
        !    sys: system object being studied.
        !    original: boolean, true if original heat bath, false if single excitations are sampled uniformly/exactly.
        ! In/Out:
        !    hb: an empty excit_gen_heat_bath_t object which gets filled with
        !           the alias tables required to generate excitations.

        use checking, only: check_allocate
        use system, only: sys_t
        use sort, only: qsort
        use proc_pointers, only: create_weighted_excitation_list_ptr, slater_condon2_excit_ptr
        use proc_pointers, only: abs_hmatel_ptr
        use excit_gens, only: excit_gen_heat_bath_t, alloc_alias_table_data_t
        use alias, only: generate_alias_tables
        use read_in_symmetry, only: cross_product_basis_read_in
        use hamiltonian_data, only: hmatel_t
        use errors, only: stop_all, warning
#ifdef PARALLEL
        use parallel

        integer :: displs_nbasis(0:nprocs-1)
        integer :: sizes_nbasis(0:nprocs-1)
        integer :: ierr
#endif
        type(sys_t), intent(in) :: sys
        type(excit_gen_heat_bath_t), intent(inout) :: hb
        logical, intent(in) :: original
        type(hmatel_t) :: hmatel

        integer :: i, j, a, b, ij_sym, isymb, ims, isyma
        integer :: i_tmp, j_tmp, a_tmp, b_tmp
        integer :: ierr_alloc
        real(p) :: i_weight, ij_weight, ija_weight, ijab_weight, ija_weights_tot, i_weight_extra
        integer :: iproc_nbasis_start, iproc_nbasis_end
        
        integer :: j_nonzero(sys%basis%nbasis,sys%basis%nbasis) 
        ! Temp storage
        allocate(hb%i_weights(sys%basis%nbasis), stat=ierr_alloc) ! This will hold S_p in the Holmes JCTC paper.
        call check_allocate('hb%i_weights', sys%basis%nbasis, ierr_alloc)
        allocate(hb%ij_weights(sys%basis%nbasis,sys%basis%nbasis), stat=ierr_alloc) ! This stores D_pq.
        call check_allocate('hb%ij_weights', (sys%basis%nbasis*sys%basis%nbasis), ierr_alloc)
        call alloc_alias_table_data_t(hb%hb_ija, sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis)
        call alloc_alias_table_data_t(hb%hb_ijab, sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis,sys%basis%nbasis)

        ! [todo] - consider setting hmatel%r and hmatel%c to zero to avoid undefined behaviour.
        ! [todo] - although it is probably fine, abs_hmatel_ptr ignores the irrelevant bit.

#ifdef PARALLEL
        ! Initialise do-loop range for each processor, [iproc_nbasis_start,iproc_nbasis_end].
        ! [todo] - get_proc_loop_range can also assign in serial mode. Disadvantage: would need to define displs_nbasis
        ! [todo] - and sizes_nbasis for serial mode.
        call get_proc_loop_range(sys%basis%nbasis, iproc_nbasis_start, iproc_nbasis_end, displs_nbasis, sizes_nbasis)
#else
        iproc_nbasis_start = 1
        iproc_nbasis_end = sys%basis%nbasis
#endif
        ! The following big for-loop/do-loop structure pre-computes weights to select i (hb%i_weights), to select j given i
        ! (hb%ij_weights), to select a given ij (hb%hb_ija%weights) and to select b given ija (hb%hb_ijab%weights).
        ! The weights are (as described in Holmes et al.):
        !   hb%i_weights(:) = \frac{\sum_jab h_ijab}{\sum_ijab h_ijab}
        !   hb%ij_weights(:,i) = \frac{\sum_ab h_ijab}{\sum_jab h_ijab}
        !   hb%hb_ija%weights(:,j,i) = \frac{\sum_b h_ijab}{\sum_ab h_ijab}
        !   hb%hb_ijab%weights(:,a,j,i) = \frac{h_ijab}{\sum_b h_ijab}
        ! where the orbitals ijab are all different and chosen from the set of all spinorbitals.
        !   h_ijab = <D_0| H |D_{ij}^{ab}>, i.e. the Hamiltonian matrix element representing the double excitation
        ! ij -> ab. Note that for double excitations, the value of h_ijab does not depend on wether we excite from
        ! the reference or not. That is useful as at the time of initialisation we do not know the determinant
        ! we will look at at a certain spawn attempt.

        ! We also check for bias in the loop if we are doing the "original" heat bath algorithm (original == .true.)
        ! and stop the calculation if we have to. There is a bias if for certain determinants to spawn from, a valid
        ! single excitation i -> a cannot be chosen as there is no occupied spinorbital (different to i) that could
        ! play the role of j when selecting ija.
        ! Reminder: We only decide whether we should do a single excitation i->a or a double excitation ij->ab
        ! after having selected ija. A single excitation i->a can therefore only be selected if there is a j
        ! we can select so that we can select ija.

        j_nonzero = 0

        do i = iproc_nbasis_start, iproc_nbasis_end
            i_weight = 0.0_p
            do j = 1, sys%basis%nbasis
                ij_weight = 0.0_p
                hb%hb_ija%weights_tot(j,i) = 0.0_p
                if (i /= j) then
                    ! [todo] - is sorting really necessary?
                    if (j < i) then
                        i_tmp = j
                        j_tmp = i
                    else
                        i_tmp = i
                        j_tmp = j
                    end if
                    ! The symmetry of b, isymb, is given by
                    ! (sym_i* x sym_j* x sym_a)* = sym_b
                    ! (at least for Abelian point groups)
                    ! ij_sym: symmetry conjugate of the irreducible representation spanned by the codensity
                    !        \phi_i*\phi_j. (We assume that ij is going to be in the bra of the excitation.)
                    ! [todo] - Check whether order of i and j matters here.
                    ij_sym = sys%read_in%sym_conj_ptr(sys%read_in, cross_product_basis_read_in(sys, i_tmp, j_tmp))
                    ija_weights_tot = 0.0_p
                    i_weight_extra = 0.0_p
            
                    !$omp parallel do default(none) &
                    !$omp shared(sys,i,j,i_tmp,j_tmp,hb,ij_sym,j_nonzero,slater_condon2_excit_ptr,abs_hmatel_ptr) &
                    !$omp private(a,ija_weight,isymb,ijab_weight,b,a_tmp,b_tmp,hmatel) &
                    !$omp reduction(+:i_weight_extra,ij_weight,ija_weights_tot)
                    do a = 1, sys%basis%nbasis
                        ija_weight = 0.0_p
                        hb%hb_ijab%weights_tot(a,j,i) = 0.0_p
                        if ((a /= i) .and. (a /= j)) then
                            isymb = sys%read_in%sym_conj_ptr(sys%read_in, &
                                        sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, sys%basis%basis_fns(a)%sym))
                            do b = 1, sys%basis%nbasis
                                ijab_weight = 0.0_p
                                ! Check spin conservation and symmetry conservation.
                                if ((((sys%basis%basis_fns(i_tmp)%Ms == sys%basis%basis_fns(a)%Ms) .and. &
                                    (sys%basis%basis_fns(j_tmp)%Ms == sys%basis%basis_fns(b)%Ms)) .or. &
                                    ((sys%basis%basis_fns(i_tmp)%Ms == sys%basis%basis_fns(b)%Ms) .and. &
                                    (sys%basis%basis_fns(j_tmp)%Ms == sys%basis%basis_fns(a)%Ms))) .and. &
                                    (sys%basis%basis_fns(b)%sym == isymb) .and. ((b /= a) .and. (b /= i) .and. &
                                    (b /= j))) then
                                    if (b < a) then
                                        a_tmp = b
                                        b_tmp = a
                                    else
                                        a_tmp = a
                                        b_tmp = b
                                    end if
                                    hmatel = slater_condon2_excit_ptr(sys, i_tmp, j_tmp, a_tmp, b_tmp, .false.)
                                    i_weight_extra = i_weight_extra + abs_hmatel_ptr(hmatel)
                                    ij_weight = ij_weight + abs_hmatel_ptr(hmatel)
                                    ija_weight = ija_weight + abs_hmatel_ptr(hmatel)
                                    ijab_weight = ijab_weight + abs_hmatel_ptr(hmatel)
                                end if
                                hb%hb_ijab%weights(b,a,j,i) = ijab_weight
                                hb%hb_ijab%weights_tot(a,j,i) = hb%hb_ijab%weights_tot(a,j,i) + ijab_weight
                            end do
                        else
                            do b = 1, sys%basis%nbasis
                                hb%hb_ijab%weights(b,a,j,i) = 0.0_p
                            end do
                        end if
                        hb%hb_ija%weights(a,j,i) = ija_weight
                        ija_weights_tot = ija_weights_tot + ija_weight
                        ! [todo] - depsilon or 0.0_p?
                        if (ija_weight > depsilon) then
                            j_nonzero(a,i) = j_nonzero(a,i) + 1
                        end if
                    end do
                    !$omp end parallel do

                    i_weight = i_weight + i_weight_extra
                    hb%hb_ija%weights_tot(j,i) = hb%hb_ija%weights_tot(j,i) + ija_weights_tot
                else
                    !$omp parallel do default(none) &
                    !$omp shared(hb,i,j,sys) &
                    !$omp private(a,b)
                    do a = 1, sys%basis%nbasis
                        hb%hb_ija%weights(a,j,i) = 0.0_p
                        hb%hb_ijab%weights_tot(a,j,i) = 0.0_p
                        do b = 1, sys%basis%nbasis
                            hb%hb_ijab%weights(b,a,j,i) = 0.0_p
                        end do
                    end do
                    !$omp end parallel do
                end if
                hb%ij_weights(j,i) = ij_weight
            end do
            hb%i_weights(i) = i_weight

            ! Test for bias and stop calculation if a bias is found.
            ! Test that all single excitation i -> a that are allowed have some non zero weight ija for some j.
            do a = 1, sys%basis%nbasis
                if ((j_nonzero(a,i) < (sys%basis%nbasis - sys%nel)) .and. (original) .and. (i /= a)) then
                    ! no non zero weight ija for some j and (i /= a). Now check that i -> a is valid.
                    ims = sys%basis%basis_fns(i)%ms
                    isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(i)%sym, &
                                                    sys%read_in%pg_sym%gamma_sym)
                    if ((sys%basis%basis_fns(a)%sym == isyma) .and. (sys%basis%basis_fns(a)%ms == ims)) then
                        call stop_all('init_excit_mol_heat_bath','Not all possible single excitations can be accounted for. &
                            &Use another excitation generator.')
                    end if
                end if
            end do
        end do

#ifdef PARALLEL
        ! note how FORTRAN stores arrays: array(2,1) comes before array(1,2) in memory.
        associate(nb=>sys%basis%nbasis)
            call mpi_allgatherv(MPI_IN_PLACE, sizes_nbasis(iproc), &
                mpi_preal, hb%i_weights, sizes_nbasis, displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, nb*sizes_nbasis(iproc), &
                mpi_preal, hb%ij_weights, nb*sizes_nbasis, nb*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, nb*sizes_nbasis(iproc), &
                mpi_preal, hb%hb_ija%weights_tot, nb*sizes_nbasis, nb*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, nb*nb*sizes_nbasis(iproc), &
                mpi_preal, hb%hb_ija%weights, nb*nb*sizes_nbasis, nb*nb*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, nb*nb*sizes_nbasis(iproc), &
                mpi_preal, hb%hb_ijab%weights_tot, nb*nb*sizes_nbasis, nb*nb*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, nb*nb*nb*sizes_nbasis(iproc), &
                mpi_preal, hb%hb_ijab%weights, nb*nb*nb*sizes_nbasis, nb*nb*nb*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
        end associate
#endif

        ! Generate alias tables for hb%hb_ija%weights and hb%hb_ijab%weights. hb%i_weights and hb%ij_weights do not
        ! get pre-computed alias tables, they will be computed on-the-fly during the run, just including occupied
        ! spinorbitals.
        do i = iproc_nbasis_start, iproc_nbasis_end
            do j = 1, sys%basis%nbasis
                if (abs(hb%hb_ija%weights_tot(j,i)) > 0.0_p) then
                    call generate_alias_tables(sys%basis%nbasis, hb%hb_ija%weights(:,j,i), hb%hb_ija%weights_tot(j,i), &
                                                hb%hb_ija%aliasU(:,j,i), hb%hb_ija%aliasK(:,j,i))
                    do a = 1, sys%basis%nbasis
                        if (abs(hb%hb_ijab%weights_tot(a,j,i)) > 0.0_p) then
                            call generate_alias_tables(sys%basis%nbasis, hb%hb_ijab%weights(:,a,j,i), &
                                hb%hb_ijab%weights_tot(a,j,i), hb%hb_ijab%aliasU(:,a,j,i), hb%hb_ijab%aliasK(:,a,j,i))
                        end if
                    end do
                end if
            end do
        end do

#ifdef PARALLEL
        ! note how FORTRAN stores arrays: array(2,1) comes before array(1,2) in memory.
        associate(nb=>sys%basis%nbasis)
            call mpi_allgatherv(MPI_IN_PLACE, nb*nb*sizes_nbasis(iproc), &
                mpi_preal, hb%hb_ija%aliasU, nb*nb*sizes_nbasis, nb*nb*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, nb*nb*sizes_nbasis(iproc), &
                MPI_INTEGER, hb%hb_ija%aliasK, nb*nb*sizes_nbasis, nb*nb*displs_nbasis, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, nb*nb*nb*sizes_nbasis(iproc), &
                mpi_preal, hb%hb_ijab%aliasU, nb*nb*nb*sizes_nbasis, nb*nb*nb*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, nb*nb*nb*sizes_nbasis(iproc), &
                MPI_INTEGER, hb%hb_ijab%aliasK, nb*nb*nb*sizes_nbasis, nb*nb*nb*displs_nbasis, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        end associate
#endif

    end subroutine init_excit_mol_heat_bath

    subroutine gen_excit_mol_heat_bath(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! This is the main heat bath algorithm described in Holmes et al. (2016)

        ! A short description (note that excit_gen_data%hb => hb):
        ! --------------------------------------------------------
        ! 1a: We first select i using the pre-calculated weights hb%i_weights(:), first restricting
        !     this weight array to occupied orbitals i. We calculate alias tables on-the-fly.
        ! 1b: Then, we find j using the pre-calculated weights hb%ij_weights(:,i), again restricting
        !     the weight array to occupied orbitals j/=i and calculate alias tables on-the-fly.
        ! 2:  We select a using pre-computed alias tables hb%hb_ija%alias{K,U}.
        ! 3:  We decide whether we should do a single excitation i->a or double excitation ij->ab,
        !     with probability 0.5 if (|h_ia| > hb%hb_ijab%weights_tot(a,j,i)) where h_ia is the
        !     Hamiltonian matrix element of the single excitation i->a from the current determinant.
        !     Else, we do a single excitation with probability
        !     |h_ia|/(hb%hb_ijab%weights_tot(b,i,j) + |h_ia|) and a double excitation otherwise.
        !     This step is slightly different to Holmes et al. who might choose to do a single
        !     and a double excitation at the same time under certain conditions.
        ! [todo] - double check what to do if |h_ia| == hb%hb_ijab%weights_tot(a,j,i)
        ! 4:  If we do a double excitation, we now select b using pre-computed alias tables
        !     hb%hb_ijab%alias{K,U}.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Excitation generation data.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant in molecular systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2, get_excitation_locations
        use proc_pointers, only: slater_condon1_excit_ptr, slater_condon2_excit_ptr, abs_hmatel_ptr
        use proc_pointers, only: create_weighted_excitation_list_ptr
        use checking, only: check_allocate, check_deallocate
        use system, only: sys_t
        use excit_gen_mol, only: gen_single_excit_mol_no_renorm
        use excit_gens, only: excit_gen_heat_bath_t, excit_gen_data_t
        use excit_gen_utils, only: select_ij_heat_bath
        use alias, only: select_weighted_value_precalc, select_weighted_value
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use hamiltonian_data, only: hmatel_t
        use read_in_symmetry, only: cross_product_basis_read_in
        use search, only: binary_search
        use sort, only: qsort

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation
        real(p) :: ij_weights_occ(sys%nel)
        real(p) :: ji_weights_occ(sys%nel)

        integer :: pos_q, i_ind, j_ind, i, j, a, b, ims, isyma
        
        real(p) :: ij_weights_occ_tot, ji_weights_occ_tot, hmatel_mod_ia, x, psingle
        real(p) :: pgen_ija, pgen_ijb, pgen_jia, pgen_jib, psingle_q
        logical :: double
        type(hmatel_t) :: hmatel_single

        associate( hb => excit_gen_data%excit_gen_hb )
            ! 1: Select orbitals i and j to excite from.
            
            call select_ij_heat_bath(rng, sys%nel, hb%ij_weights, cdet, i, j, i_ind, j_ind, ij_weights_occ, ij_weights_occ_tot, &
                ji_weights_occ, ji_weights_occ_tot, allowed_excitation)
            
            ! Is there a possible a?
            if ((allowed_excitation) .and. (abs(hb%hb_ija%weights_tot(j, i)) > 0.0_p)) then
                allowed_excitation = .true.
            else
                allowed_excitation = .false.
            end if
            
            ! Important: We do not order i and j here as the original order matters before we decide whether to do a single
            ! excitation i -> a or a double excitation ij -> ab.

            ! If allowed_excitation == .true. at this stage, we will have checked that
            ! (abs(hb%hb_ija%weights_tot(j, i)) > 0.0_p) in the subroutine select_ij_heat_bath, i.e. there exists an a we
            ! can find.
            if (allowed_excitation) then
                ! 2: Find a.
                a = select_weighted_value_precalc(rng, sys%basis%nbasis, hb%hb_ija%aliasU(:, j, i), hb%hb_ija%aliasK(:, j, i))
                ! Check a is not occupied.
                if (.not.btest(cdet%f(sys%basis%bit_lookup(2,a)), sys%basis%bit_lookup(1,a))) then
                    ! 3: Determine whether to do a single or a double excitation.
                    ims = sys%basis%basis_fns(i)%ms
                    isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(i)%sym, &
                                                    sys%read_in%pg_sym%gamma_sym)    
                    if ((sys%basis%basis_fns(a)%sym == isyma) .and. (sys%basis%basis_fns(a)%ms == ims)) then
                        connection%from_orb(1) = i
                        connection%to_orb(1) = a
                        connection%nexcit = 1
                        call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)
                        
                        hmatel_single = slater_condon1_excit_ptr(sys, cdet%occ_list, i, a, connection%perm)
                        hmatel_mod_ia = abs_hmatel_ptr(hmatel_single)
                        x = get_rand_close_open(rng)                        
                        if (hmatel_mod_ia < hb%hb_ijab%weights_tot(a,j,i)) then
                            psingle = hmatel_mod_ia/(hb%hb_ijab%weights_tot(a,j,i) + hmatel_mod_ia)
                        else
                            psingle = 0.5_p
                        end if
                        if (x < psingle) then
                            double = .false.
                        else
                            double = .true.
                        end if
                    else
                        ! i -> a not possible without j and b.
                        double = .true.
                        psingle = 0.0_p
                    end if
                else
                    allowed_excitation = .false.
                end if
            else
                allowed_excitation = .false.
            end if

            if (allowed_excitation) then
                ! 4: We have decided to do a double excitation. Find b.
                if (double) then
                    b = select_weighted_value_precalc(rng, sys%basis%nbasis, hb%hb_ijab%aliasU(:, a, j, i), &
                                                        hb%hb_ijab%aliasK(:, a, j, i))
                    ! Is b occupied?
                    if (.not.btest(cdet%f(sys%basis%bit_lookup(2,b)), sys%basis%bit_lookup(1,b))) then
                        allowed_excitation = .true.
                        
                        ! Calculate all contributions to pgen (need to consider different orders of having selected i and j
                        ! and different combinations of a and b. Sub-pgens are pgen_ija, pgen_ijb, pgen_jia, pgen_jib.
                        pgen_ija = ((cdet%i_d_weights_occ(i_ind)/cdet%i_d_weights_occ_tot) * & 
                                   (ij_weights_occ(j_ind)/ij_weights_occ_tot)) * &
                            (hb%hb_ija%weights(a,j,i)/hb%hb_ija%weights_tot(j,i)) * (1.0_p - psingle) * &
                            (hb%hb_ijab%weights(b,a,j,i)/hb%hb_ijab%weights_tot(a,j,i))

                        ims = sys%basis%basis_fns(i)%ms
                        isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(i)%sym, &
                                                                            sys%read_in%pg_sym%gamma_sym)

                        connection%from_orb(1) = i
                        connection%to_orb(1) = b
                        connection%nexcit = 1

                        if ((sys%basis%basis_fns(b)%sym == isyma) .and. (sys%basis%basis_fns(b)%ms == ims)) then
                            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)
                            hmatel_single = slater_condon1_excit_ptr(sys, cdet%occ_list, i, b, connection%perm)
                            hmatel_mod_ia = abs_hmatel_ptr(hmatel_single)
                            if (hmatel_mod_ia < hb%hb_ijab%weights_tot(b,j,i)) then
                                psingle = hmatel_mod_ia/(hb%hb_ijab%weights_tot(b,j,i) + hmatel_mod_ia)
                            else
                                psingle = 0.5_p
                            end if
                        else
                            psingle = 0.0_p
                        end if

                        pgen_ijb = ((cdet%i_d_weights_occ(i_ind)/cdet%i_d_weights_occ_tot) * &
                             (ij_weights_occ(j_ind)/ij_weights_occ_tot)) * &
                            (hb%hb_ija%weights(b,j,i)/hb%hb_ija%weights_tot(j,i)) * (1.0_p - psingle) * &
                            (hb%hb_ijab%weights(a,b,j,i)/hb%hb_ijab%weights_tot(b,j,i))

                        ims = sys%basis%basis_fns(j)%ms
                        isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(j)%sym, &
                                                    sys%read_in%pg_sym%gamma_sym)

                        connection%from_orb(1) = j
                        connection%to_orb(1) = a
                        connection%nexcit = 1

                        if ((sys%basis%basis_fns(a)%sym == isyma) .and. (sys%basis%basis_fns(a)%ms == ims)) then
                            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)
                            hmatel_single = slater_condon1_excit_ptr(sys, cdet%occ_list, j, a, connection%perm)
                            hmatel_mod_ia = abs_hmatel_ptr(hmatel_single)
                            if (hmatel_mod_ia < hb%hb_ijab%weights_tot(a,i,j)) then
                                psingle = hmatel_mod_ia/(hb%hb_ijab%weights_tot(a,i,j) + hmatel_mod_ia)
                            else
                                psingle = 0.5_p
                            end if
                        else
                            psingle = 0.0_p
                        end if

                        pgen_jia = ((cdet%i_d_weights_occ(j_ind)/cdet%i_d_weights_occ_tot) * &
                            (ji_weights_occ(i_ind)/ji_weights_occ_tot)) * &
                            (hb%hb_ija%weights(a,i,j)/hb%hb_ija%weights_tot(i,j)) * (1.0_p - psingle) * &
                            (hb%hb_ijab%weights(b,a,i,j)/hb%hb_ijab%weights_tot(a,i,j))

                        connection%from_orb(1) = j
                        connection%to_orb(1) = b
                        connection%nexcit = 1
                        
                        if ((sys%basis%basis_fns(b)%sym == isyma) .and. (sys%basis%basis_fns(b)%ms == ims)) then
                            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)
                            hmatel_single = slater_condon1_excit_ptr(sys, cdet%occ_list, j, b, connection%perm)
                            hmatel_mod_ia = abs_hmatel_ptr(hmatel_single)
                            if (hmatel_mod_ia < hb%hb_ijab%weights_tot(b,i,j)) then
                                psingle = hmatel_mod_ia/(hb%hb_ijab%weights_tot(b,i,j) + hmatel_mod_ia)
                            else
                                psingle = 0.5_p
                            end if
                        else
                            psingle = 0.0_p
                        end if

                        pgen_jib = ((cdet%i_d_weights_occ(j_ind)/cdet%i_d_weights_occ_tot) * &
                            (ji_weights_occ(i_ind)/ji_weights_occ_tot)) * &
                            (hb%hb_ija%weights(b,i,j)/hb%hb_ija%weights_tot(i,j)) * (1.0_p - psingle) * &
                            (hb%hb_ijab%weights(a,b,i,j)/hb%hb_ijab%weights_tot(b,i,j))

                        pgen = pgen_ija + pgen_ijb + pgen_jia + pgen_jib

                        if (i < j) then
                            connection%from_orb(1) = i
                            connection%from_orb(2) = j
                        else
                            connection%from_orb(1) = j
                            connection%from_orb(2) = i
                        end if
                        if (a < b) then
                            connection%to_orb(1) = a
                            connection%to_orb(2) = b
                        else
                            connection%to_orb(2) = a
                            connection%to_orb(1) = b
                        end if
                        connection%nexcit = 2

                        call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)
                        hmatel = slater_condon2_excit_ptr(sys, connection%from_orb(1), connection%from_orb(2), &
                                                    connection%to_orb(1), connection%to_orb(2), connection%perm)
                    else
                        allowed_excitation = .false.
                        ! We have not found a valid excitation.
                        ! Return a null excitation.
                        hmatel%c = cmplx(0.0_p, 0.0_p, p)
                        hmatel%r = 0.0_p
                        pgen = 1.0_p
                    end if
                else 
                    ! We have decided to do a single excitation. Need to calculate pgen. Given that j is now a dummy variable,
                    ! we need to consider that we could have achieved the same single excitation i->a with a different (occupied)
                    ! j/=i.
                    connection%from_orb(1) = i
                    connection%to_orb(1) = a
                    connection%nexcit = 1
                    hmatel = slater_condon1_excit_ptr(sys, cdet%occ_list, i, a, connection%perm)
                    pgen = 0.0_p
                    do pos_q = 1, sys%nel
                        if ((i /= cdet%occ_list(pos_q)) .and. (a /= cdet%occ_list(pos_q))) then
                            if (hmatel_mod_ia < hb%hb_ijab%weights_tot(a,cdet%occ_list(pos_q),i)) then
                                psingle_q = hmatel_mod_ia/(hb%hb_ijab%weights_tot(a,cdet%occ_list(pos_q),i) + hmatel_mod_ia)
                            else
                                psingle_q = 0.5_p
                            end if
                            pgen = pgen + (psingle_q * (ij_weights_occ(pos_q)/ij_weights_occ_tot) * &
                                (hb%hb_ija%weights(a, cdet%occ_list(pos_q), i)/hb%hb_ija%weights_tot(cdet%occ_list(pos_q), i)))
                        end if
                    end do
                    pgen = pgen * (cdet%i_d_weights_occ(i_ind)/cdet%i_d_weights_occ_tot)
                end if
            else
                ! We have not found a valid excitation.
                ! Return a null excitation.
                hmatel%c = cmplx(0.0_p, 0.0_p, p)
                hmatel%r = 0.0_p
                pgen = 1.0_p
            end if

        end associate

    end subroutine gen_excit_mol_heat_bath

    subroutine gen_excit_mol_heat_bath_uniform(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! This is the heat bath algorithm that selects single excitations uniformly or exactly.
        ! This is based on Holmes et al. (2016) (uniform) and a recommendation by Pablo Rios Lopez and Ali Alavi (exact).
        ! How we handle double excitations is the same as in the "original" version above.
        
        ! A short description (note that excit_gen_data%hb => hb):
        ! --------------------------------------------------------
        ! 0:  We decide whether to do a single or a double excitation. Either we use the renorm
        !     excitation generator to find i->a if doing singles (if we selected the uniform option)
        !     or we calculate the weights exactly (single option). For more information on the latter,
        !     see subroutine below. Points 1-3 are for double excitations only.
        ! 1a: If we decide to do a double excitation, we first select i using the pre-calculated weights
        !     hb%i_weights(:), first restricting
        !     this weight array to occupied orbitals i. We calculate alias tables on-the-fly.
        ! 1b: Then, we find j using the pre-calculated weights hb%ij_weights(:,i), again restricting
        !     the weight array to occupied orbitals j/=i and calculate alias tables on-the-fly.
        ! 2:  We select a using pre-computed alias tables hb%hb_ija%alias{K,U}.
        ! 3:  We now select b using pre-computed alias tables hb%hb_ijab%alias{K,U}.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Excitation generation data.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant in molecular systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2, get_excitation_locations
        use proc_pointers, only: slater_condon2_excit_ptr
        use proc_pointers, only: create_weighted_excitation_list_ptr
        use checking, only: check_allocate, check_deallocate
        use system, only: sys_t
        use excit_gen_mol, only: gen_single_excit_mol
        use excit_gens, only: excit_gen_heat_bath_t, excit_gen_data_t
        use excit_gen_utils, only: select_ij_heat_bath
        use alias, only: select_weighted_value_precalc, select_weighted_value
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use hamiltonian_data, only: hmatel_t
        use qmc_data, only: excit_gen_heat_bath_uniform
        use read_in_symmetry, only: cross_product_basis_read_in
        use search, only: binary_search
        use sort, only: qsort

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation
        real(p) :: ij_weights_occ(sys%nel)
        real(p) :: ji_weights_occ(sys%nel)

        integer :: i_ind, j_ind, i, j, a, b, j_tmp
        
        real(p) :: ij_weights_occ_tot, ji_weights_occ_tot

        ! 0: Select single or double.
        if (get_rand_close_open(rng) < excit_gen_data%pattempt_single) then  
            ! We have a single
            if (excit_gen_data%excit_gen == excit_gen_heat_bath_uniform) then
                ! The user has chosen the uniform option. Call renorm single excitation generator.
                call gen_single_excit_mol(rng, sys, excit_gen_data%pattempt_single, cdet, pgen, connection, hmatel, &
                                                allowed_excitation)
            else
                ! The user has chosen the single option (exact weights). Call subroutine below.
                call gen_single_excit_heat_bath_exact(rng, sys, excit_gen_data%pattempt_single, cdet, pgen, connection, hmatel,&
                                                    allowed_excitation)
            end if
        else
            ! We have a double
            associate( hb => excit_gen_data%excit_gen_hb )
                ! 1: Select orbitals i and j to excite from.
                
                call select_ij_heat_bath(rng, sys%nel, hb%ij_weights, cdet, i, j, i_ind, j_ind, ij_weights_occ, &
                    ij_weights_occ_tot, ji_weights_occ, ji_weights_occ_tot, allowed_excitation)
            
                ! Is there a possible a?
                if ((allowed_excitation) .and. (abs(hb%hb_ija%weights_tot(j, i)) > 0.0_p)) then
                    allowed_excitation = .true.
                else
                    allowed_excitation = .false.
                end if
                
                pgen = ((cdet%i_d_weights_occ(i_ind)/cdet%i_d_weights_occ_tot) * (ij_weights_occ(j_ind)/ij_weights_occ_tot)) + &
                        ((cdet%i_d_weights_occ(j_ind)/cdet%i_d_weights_occ_tot) * (ji_weights_occ(i_ind)/ji_weights_occ_tot))

                ! [todo] - this is technically not necessary at this stage, the weights are symmetric in i and j (and in fact
                ! [todo] - when precalculating the weights, we sort i and j such that i < j).
                ! sort i and j.
                if (j < i) then
                    j_tmp = j
                    j = i
                    i = j_tmp
                end if
                
                ! If allowed_excitation == .true. at this stage, we will have checked that
                ! (abs(hb%hb_ija%weights_tot(j, i)) > 0.0_p) in the subroutine select_ij_heat_bath, i.e. there exists an a we
                ! can find.
                if (allowed_excitation) then
                    ! 2: Find a.
                    a = select_weighted_value_precalc(rng, sys%basis%nbasis, hb%hb_ija%aliasU(:, j, i), &
                                                    hb%hb_ija%aliasK(:, j, i))
                    ! Check that a possible b exists and that a is not occupied.
                    if ((abs(hb%hb_ijab%weights_tot(a, j, i)) > 0.0_p) .and. (.not.btest(cdet%f(sys%basis%bit_lookup(2,a)), &
                            sys%basis%bit_lookup(1,a)))) then
                        ! 3: Select b.
                        b = select_weighted_value_precalc(rng, sys%basis%nbasis, hb%hb_ijab%aliasU(:, a, j, i), &
                                                            hb%hb_ijab%aliasK(:, a, j, i))
                        if (.not.btest(cdet%f(sys%basis%bit_lookup(2,b)), sys%basis%bit_lookup(1,b))) then
                            allowed_excitation = .true.               
                        else
                            allowed_excitation = .false.
                        end if
                    else
                        allowed_excitation = .false.
                    end if
                end if

                if (allowed_excitation) then
                    ! i and j were sorted above.
                    connection%from_orb(1) = i
                    connection%from_orb(2) = j
                    if (a < b) then
                        connection%to_orb(1) = a
                        connection%to_orb(2) = b
                    else
                        connection%to_orb(2) = a
                        connection%to_orb(1) = b
                    end if
                    connection%nexcit = 2
                    ! 4b. Parity of permutation required to line up determinants.
                    ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                    call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

                    ! 5b. Find the connecting matrix element.
                    hmatel = slater_condon2_excit_ptr(sys, connection%from_orb(1), connection%from_orb(2), &
                                                      connection%to_orb(1), connection%to_orb(2), connection%perm)
                    pgen = (1.0_p - excit_gen_data%pattempt_single) * pgen * &
                        (((hb%hb_ija%weights(a, j, i)/hb%hb_ija%weights_tot(j, i)) * &
                        (hb%hb_ijab%weights(b, a, j, i)/hb%hb_ijab%weights_tot(a, j, i))) + &
                        ((hb%hb_ija%weights(b, j, i)/hb%hb_ija%weights_tot(j, i)) * &
                        (hb%hb_ijab%weights(a, b, j, i)/hb%hb_ijab%weights_tot(b, j, i))))
                else
                    ! We have not found a valid excitation ij -> ab.  Such
                    ! events are not worth the cost of renormalising the generation
                    ! probabilities.
                    ! Return a null excitation.
                    hmatel%c = cmplx(0.0_p, 0.0_p, p)
                    hmatel%r = 0.0_p
                    pgen = 1.0_p
                end if
            end associate
        end if

    end subroutine gen_excit_mol_heat_bath_uniform

    subroutine gen_single_excit_heat_bath_exact(rng, sys, pattempt_single, cdet, pgen, connection, hmatel, &
            allowed_excitation)

        ! Create a random single excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! This uses weights for i and a calculated as exact as possibly. Called if user selected
        ! heat bath single option and we are doing a single excitation.
        ! For i, the weight is \sum_a H_ia, for a given i (a|i) it is H_ia.

        ! In:
        !    sys: system object being studied.
        !    pattempt_single: probability of having chosen to attempt a single and not a double
        !                    excitation.
        !    cdet: determinant to attempt spawning from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !       determinant and a connected determinant in molecular systems.
        !    allowed_excitation: false if a valid symmetry allowed excitation was not generated

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation1
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t
        use alias, only: select_weighted_value
        use hamiltonian_data, only: hmatel_t
        use proc_pointers, only: slater_condon1_excit_ptr

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: pattempt_single
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation

        integer :: i_ind, a_ind, i, a

        allowed_excitation = .true.
        connection%nexcit = 1
        
        if (cdet%i_s_weights_occ_tot < depsilon) then
            ! no allowed single excitations.
            allowed_excitation = .false.
            hmatel%c = cmplx(0.0_p, 0.0_p, p)
            hmatel%r = 0.0_p
            pgen = 1.0_p ! Avoid any dangerous division by pgen by returning a sane (but cheap) value.
        else
            ! Select i.
            i_ind = select_weighted_value(rng, sys%nel, cdet%i_s_weights_occ, cdet%i_s_weights_occ_tot)
            i = cdet%occ_list(i_ind)

            ! Select a.
            a_ind = select_weighted_value(rng, sys%nvirt, cdet%ia_s_weights_occ(:, i_ind), cdet%i_s_weights_occ(i_ind))
            a = cdet%unocc_list(a_ind)

            connection%from_orb(1) = i
            connection%to_orb(1) = a

            ! Calculate pgen and hmatel.
            pgen = pattempt_single * (cdet%i_s_weights_occ(i_ind)/cdet%i_s_weights_occ_tot) * &
                                    (cdet%ia_s_weights_occ(a_ind, i_ind)/cdet%i_s_weights_occ(i_ind))

            ! Parity of permutation required to line up determinants.
            call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)

            ! Find the connecting matrix element.
            hmatel = slater_condon1_excit_ptr(sys, cdet%occ_list, connection%from_orb(1), connection%to_orb(1), connection%perm)
        end if

    end subroutine gen_single_excit_heat_bath_exact
    
end module excit_gen_heat_bath_mol
