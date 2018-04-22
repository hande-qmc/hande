module excit_gen_power_pitzer_mol

! A module containing excitations generators for molecules which weight excitations according to the exchange matrix elements.

use const, only: i0, p, depsilon

implicit none

! [review] - JSS: some code-style uniformity please.  e.g. end do/end if instead of end do and end if

contains

    subroutine init_excit_mol_power_pitzer_occ_ref(sys, ref, pp)

        ! Generate excitation tables from the reference for the
        ! gen_excit_mol_power_pitzer_occ_ref excitation generator.
        ! This creates a random excitation from a det and calculates both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        ! Weight the double excitations according the the Power-Pitzer bound
        ! <ij|ab> <= Sqrt(<ia|ai><jb|bj>), see J.D. Power, R.M. Pitzer, Chem. Phys. Lett.,
        ! 478-483 (1974).
        ! This is an O(M) version which pretends the determinant excited from is the reference,
        ! then modifies the selected orbitals to be those of the determinant given.
        ! Each occupied and virtual not present in the det we're actually given will be
        ! mapped to the one of the equivalent free numerical index in the reference.

        ! In:
        !    sys: system object being studied.
        !    ref: the reference from which we are exciting.
        ! In/Out:
        !    pp: an empty excit_gen_power_pitzer_t object which gets filled with
        !           the alias tables required to generate excitations.

        use checking, only: check_allocate
        use system, only: sys_t
        use qmc_data, only: reference_t
        use sort, only: qsort
        use proc_pointers, only: create_weighted_excitation_list_ptr
        use excit_gens, only: excit_gen_power_pitzer_t, alloc_alias_table_data_t
        use alias, only: generate_alias_tables
        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: ref
        type(excit_gen_power_pitzer_t), intent(inout) :: pp

        integer :: i, j, ind_a, ind_b, maxv, nv, bsym, ierr_alloc
        
        ! Temp storage
        maxv = max(sys%nvirt_alpha,sys%nvirt_beta)
        call alloc_alias_table_data_t(pp%pp_ia_d, maxv, sys%nel)
        call alloc_alias_table_data_t(pp%pp_jb_d, maxval(sys%read_in%pg_sym%nbasis_sym_spin), &
            (/sys%sym0_tot,sys%sym_max_tot/), sys%nel)
        allocate(pp%occ_list(sys%nel + 1), stat=ierr_alloc)  ! The +1 is a pad to allow loops to look better
        call check_allocate('pp%occ_list', (sys%nel + 1), ierr_alloc)
        allocate(pp%virt_list_alpha(sys%nvirt_alpha), stat=ierr_alloc)
        call check_allocate('pp%virt_list_alpha', sys%nvirt_alpha, ierr_alloc)
        allocate(pp%virt_list_beta(sys%nvirt_beta), stat=ierr_alloc)
        call check_allocate('pp%virt_list_beta', sys%nvirt_beta, ierr_alloc)
        
        pp%occ_list(:sys%nel) = ref%occ_list0(:sys%nel)
        pp%occ_list(sys%nel + 1) = sys%basis%nbasis*2  ! A pad
        ! [todo] - Consider testing j <= sys%nel below instead of having this pad.
        
        ! Now sort this, just in case we have an old restart file and the reference was not sorted then.
        call qsort(pp%occ_list,sys%nel)

        ! Make the unocc list.
        j = 1     ! The next occ to look at
        ind_a = 0 ! The present position in the virt_list we're making
        ind_b = 0 ! The present position in the virt_list we're making
        
        do i = 1, sys%basis%nbasis
            if (i==pp%occ_list(j)) then ! Our basis fn is in the ref
                ! Due to the +1 pad in occ_list, there is not danger of going past array boundaries here.
                j = j + 1
            else ! Need to store it as a virt
                if (sys%basis%basis_fns(i)%Ms == -1) then ! beta
                    ind_b = ind_b + 1
                    pp%virt_list_beta(ind_b) = i
                else
                    ind_a = ind_a + 1
                    pp%virt_list_alpha(ind_a) = i
                end if
            end if
        end do

        ! pp%virt_list_*(:) now contains the lists of virtual alpha and beta with respect to ref.

        ! Now generate the occ->virtual weighting lists and alias tables.
        do i = 1, sys%nel
            j = pp%occ_list(i)  ! The elec we're looking at
            if (sys%basis%basis_fns(j)%Ms == -1) then ! beta
                nv = sys%nvirt_beta
                if (nv > 0) then
                    call create_weighted_excitation_list_ptr(sys, j, 0, pp%virt_list_beta, nv, pp%pp_ia_d%weights(:,i), &
                                                            pp%pp_ia_d%weights_tot(i))
                    call check_min_weight_ratio(pp%pp_ia_d%weights(:,i), pp%pp_ia_d%weights_tot(i), nv, &
                                               pp%power_pitzer_min_weight)
                    call generate_alias_tables(nv, pp%pp_ia_d%weights(:,i), pp%pp_ia_d%weights_tot(i), pp%pp_ia_d%aliasU(:,i), &
                                               pp%pp_ia_d%aliasK(:,i))
                end if
                do bsym = sys%sym0_tot, sys%sym_max_tot
                    if (sys%read_in%pg_sym%nbasis_sym_spin(1,bsym) > 0) then
                        call create_weighted_excitation_list_ptr(sys, j, 0, sys%read_in%pg_sym%sym_spin_basis_fns(:,1,bsym), &
                            sys%read_in%pg_sym%nbasis_sym_spin(1,bsym), pp%pp_jb_d%weights(:,bsym,i), &
                            pp%pp_jb_d%weights_tot(bsym,i))
                        call check_min_weight_ratio(pp%pp_jb_d%weights(:,bsym,i), pp%pp_jb_d%weights_tot(bsym,i), &
                                                   sys%read_in%pg_sym%nbasis_sym_spin(1,bsym), pp%power_pitzer_min_weight)
                        call generate_alias_tables(sys%read_in%pg_sym%nbasis_sym_spin(1,bsym), pp%pp_jb_d%weights(:,bsym,i), &
                            pp%pp_jb_d%weights_tot(bsym,i), pp%pp_jb_d%aliasU(:,bsym,i), pp%pp_jb_d%aliasK(:,bsym,i))
                    end if
                end do
            else ! alpha
                nv = sys%nvirt_alpha
                if (nv > 0) then
                    call create_weighted_excitation_list_ptr(sys, j, 0, pp%virt_list_alpha, nv, pp%pp_ia_d%weights(:,i), &
                                                             pp%pp_ia_d%weights_tot(i))
                    call check_min_weight_ratio(pp%pp_ia_d%weights(:,i), pp%pp_ia_d%weights_tot(i), nv, &
                        pp%power_pitzer_min_weight)
                    call generate_alias_tables(nv, pp%pp_ia_d%weights(:,i), pp%pp_ia_d%weights_tot(i), pp%pp_ia_d%aliasU(:,i), &
                                               pp%pp_ia_d%aliasK(:,i))
                end if
                do bsym = sys%sym0_tot, sys%sym_max_tot
                    if (sys%read_in%pg_sym%nbasis_sym_spin(2,bsym) > 0) then
                        call create_weighted_excitation_list_ptr(sys, j, 0, sys%read_in%pg_sym%sym_spin_basis_fns(:,2,bsym), &
                            sys%read_in%pg_sym%nbasis_sym_spin(2,bsym), pp%pp_jb_d%weights(:,bsym,i), &
                            pp%pp_jb_d%weights_tot(bsym,i))
                        call check_min_weight_ratio(pp%pp_jb_d%weights(:,bsym,i), pp%pp_jb_d%weights_tot(bsym,i), &
                                                   sys%read_in%pg_sym%nbasis_sym_spin(2,bsym), pp%power_pitzer_min_weight)
                        call generate_alias_tables(sys%read_in%pg_sym%nbasis_sym_spin(2,bsym), pp%pp_jb_d%weights(:,bsym,i), &
                            pp%pp_jb_d%weights_tot(bsym,i), pp%pp_jb_d%aliasU(:,bsym,i), pp%pp_jb_d%aliasK(:,bsym,i))
                    end if
                end do
            end if
        end do
    end subroutine init_excit_mol_power_pitzer_occ_ref

    subroutine check_min_weight_ratio(weights, weights_tot, weights_len, min_ratio)
        
        ! Restrict the minimum ratio of weights(i)/weights_tot to be min_ratio/number of orbitals with nonzero
        ! weight.

        ! In:
        !   weights_len: number of elements in weights list
        !   min_ratio: minimum value of weights(i)/weights_tot
        ! In/Out:
        !   weights: list of weights
        !   weights_tot: sum of weights
    
        use const, only: depsilon

        integer, intent(in) :: weights_len
        real(p), intent(in) :: min_ratio
        real(p), intent(inout) :: weights(:), weights_tot

        integer :: i, j, k, min_weights_counter, non_zero_weights_len
        real(p) :: weights_keep, min_weight, min_weight_tmp

        min_weight_tmp = 0.0_p
        non_zero_weights_len = 0

        if ((weights_tot > 0.0_p) .and. (min_ratio > 0.0_p)) then
            
            ! Find the number of non zero weights.
            do i = 1, weights_len
                if (weights(i) > 0.0_p) then
                    non_zero_weights_len = non_zero_weights_len + 1
                end if
            end do

            min_weight = (min_ratio/real(non_zero_weights_len)) * weights_tot

            ! Find the correct min_weight such that min_ratio*weights_tot/non_zero_weights_len is min_weight at the end.            
            do while (abs(min_weight_tmp - min_weight) > depsilon)
                min_weight_tmp = min_weight
                weights_keep = 0.0_p
                min_weights_counter = 0
                
                do k = 1, weights_len
                    if ((weights(k) > 0.0_p) .and. (weights(k) < min_weight)) then
                        min_weights_counter = min_weights_counter + 1
                    else
                        weights_keep = weights_keep + weights(k)
                    end if            
                end do
                if (min_weights_counter == non_zero_weights_len) then
                    exit
                end if
                ! The weights_tot of the next (future) iteration is calculated as
                ! weights_tot (next iteration) = weights_keep + 
                ! min_ratio*min_weights_counter*weights_tot (next iteration) /non_zero_weights_len,
                ! which is then solved to be 
                ! weights_tot (next iteration) = weights_keep / (1-min_ratio*min_weights_counter/non_zero_weights_len).
                ! Substitute into min_weight = weights_tot * min_ratio/non_zero_weights_len
                ! and get the expression below.
                min_weight = (min_ratio/real(non_zero_weights_len)) * &
                    (weights_keep/(1-(min_ratio*real(min_weights_counter)/real(non_zero_weights_len))))
            end do

            ! Set the new weights if required and update weights_tot.
            weights_tot = 0.0_p
            do j = 1, weights_len
                if ((weights(j) > 0.0_p) .and. (weights(j) < min_weight)) then
                    weights(j) = min_weight
                end if
                weights_tot = weights_tot + weights(j)
            end do
                
        end if

    end subroutine check_min_weight_ratio
 
    subroutine init_excit_mol_power_pitzer_orderN(sys, ref, pp)

        ! Generate excitation tables for all spinorbitals for the gen_excit_mol_power_pitzer_orderN
        ! excitation generator.

        ! In:
        !    sys: system object being studied.
        !    ref: the reference from which we are exciting.
        ! In/Out:
        !    pp: an empty excit_gen_power_pitzer_t object which gets filled with
        !           the alias tables required to generate excitations.

        use checking, only: check_allocate
        use system, only: sys_t
        use qmc_data, only: reference_t
        use sort, only: qsort
        use proc_pointers, only: create_weighted_excitation_list_ptr, slater_condon2_excit_ptr
        use proc_pointers, only: abs_hmatel_ptr, single_excitation_weight_ptr
        use excit_gens, only: excit_gen_power_pitzer_t, alloc_alias_table_data_t
        use excit_gen_utils, only: init_double_weights_ab
        use alias, only: generate_alias_tables
        use read_in_symmetry, only: cross_product_basis_read_in
        use hamiltonian_data, only: hmatel_t
#ifdef PARALLEL
        use parallel

        integer :: displs_nel(0:nprocs-1), displs_nbasis(0:nprocs-1)
        integer :: sizes_nel(0:nprocs-1), sizes_nbasis(0:nprocs-1)
        integer :: ierr, sr
#endif
        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: ref
        type(excit_gen_power_pitzer_t), intent(inout) :: pp

        integer :: i, j, a, ind_a, ind_b, bsym, isyma, ims, imsa
        integer :: nall
        integer :: iproc_nel_start, iproc_nel_end, iproc_nbasis_start, iproc_nbasis_end
        integer :: ierr_alloc
        real(p) :: i_weight, ij_weight, ia_s_weights_tot, ij_d_weights_tot
        
        ! Store weights and alias tables.
        ! pp%ppn_i_d%weights(:) selects i from orbitals occupied in the reference in a double excitation.
        call alloc_alias_table_data_t(pp%ppn_i_d, sys%nel)
        ! pp%ppn_i_s%weights(:) selects i from orbitals occupied in the reference in a single excitation.
        call alloc_alias_table_data_t(pp%ppn_i_s, sys%nel)
        ! pp%ppn_ij_d%weights(:,i) selects j from orbitals occupied in the reference in a double excitation.
        call alloc_alias_table_data_t(pp%ppn_ij_d, sys%nel, sys%basis%nbasis)
        ! pp%ppn_ia_d%weights(:,i) select a from all orbitals in a double excitation.
        ! Note that the alias table values for this are allocated further below.
        ! pp%ppn_ia_s%weights(:,i) selects a from all orbitals in a single excitation.
        call alloc_alias_table_data_t(pp%ppn_ia_s, maxval(sys%read_in%pg_sym%nbasis_sym_spin), sys%basis%nbasis)
        ! pp%ppn_jb_d%weights(:,symb,j) selects b from orbitals with j's spin and symmetry symb in a double excitation.
        call alloc_alias_table_data_t(pp%ppn_jb_d, maxval(sys%read_in%pg_sym%nbasis_sym_spin), (/sys%sym0_tot,sys%sym_max_tot/),&
            sys%basis%nbasis)
        allocate(pp%occ_list(sys%nel + 1), stat=ierr_alloc)  ! The +1 is a pad to allow loops to look better
        call check_allocate('pp%occ_list', (sys%nel + 1), ierr_alloc)
        allocate(pp%all_list_alpha(sys%basis%nbasis), stat=ierr_alloc)
        call check_allocate('pp%all_list_alpha', sys%basis%nbasis, ierr_alloc)
        allocate(pp%all_list_beta(sys%basis%nbasis), stat=ierr_alloc)
        call check_allocate('pp%all_list_beta', sys%basis%nbasis, ierr_alloc)

        pp%occ_list(:sys%nel) = ref%occ_list0(:sys%nel)
        pp%occ_list(sys%nel + 1) = sys%basis%nbasis*2  ! A pad
        ! [todo] - Consider testing j <= sys%nel below instead of having this pad.

        ! Now sort this, just in case we have an old restart file and the reference was not sorted then.
        call qsort(pp%occ_list,sys%nel)

#ifdef PARALLEL
        ! Initialise do-loop range for each processor, e.g. [iproc_nel_start,iproc_nel_end], in the case for
        ! a do-loop over sys%nel.
        call get_proc_loop_range(sys%nel, iproc_nel_start, iproc_nel_end, displs_nel, sizes_nel)
        call get_proc_loop_range(sys%basis%nbasis, iproc_nbasis_start, iproc_nbasis_end, displs_nbasis, sizes_nbasis)
#else
        iproc_nel_start = 1
        iproc_nbasis_start = 1
        iproc_nel_end = sys%nel
        iproc_nbasis_end = sys%basis%nbasis
#endif

        ! Fill the list with alpha and the list with betas.
        ind_a = 0
        ind_b = 0

        do i = 1, sys%basis%nbasis
            if (sys%basis%basis_fns(i)%Ms == -1) then ! beta
                ind_b = ind_b + 1
                pp%all_list_beta(ind_b) = i
            else
                ind_a = ind_a + 1
                pp%all_list_alpha(ind_a) = i
            end if
        end do

        pp%n_all_alpha = ind_a
        pp%n_all_beta = ind_b

        call alloc_alias_table_data_t(pp%ppn_ia_d, max(pp%n_all_alpha,pp%n_all_beta), sys%basis%nbasis)

        ! Now set up weight tables.
        ! Single Excitations.
        do i = iproc_nel_start, iproc_nel_end
            i_weight = 0.0_p
            ims = sys%basis%basis_fns(pp%occ_list(i))%ms
            isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(pp%occ_list(i))%sym, &
                                                    sys%read_in%pg_sym%gamma_sym)
            do a = 1, sys%basis%nbasis
            ! Check that a is not i and has the required spin and symmetry.
                if ((a /= pp%occ_list(i)) .and. (sys%basis%basis_fns(a)%sym == isyma) .and. &
                    (sys%basis%basis_fns(a)%ms == ims)) then
                    i_weight = i_weight + single_excitation_weight_ptr(sys, ref, pp%occ_list(i), a)
                end if
            end do
            if (i_weight < depsilon) then
                ! because we map i later, even if i->a is not allowed/ has been wrongly assigned a zero weight,
                ! it needs a weight. It will be assigned a minimum weight by the routine check_in_weight_ratio
                ! called after this loop. The finite weight here should be very small (and will be raised in
                ! the min. weight function) but still detectable (> depsilon).
                ! [todo] - factor in front of depsilon is a bit arbitary. 10depsilon is assumed to be less than
                ! [todo] - the min. weight.
                i_weight = 10.0_p*depsilon
            end if
            pp%ppn_i_s%weights(i) = i_weight
        end do

#ifdef PARALLEL
        call mpi_allgatherv(MPI_IN_PLACE, sizes_nel(iproc), &
            mpi_preal, pp%ppn_i_s%weights, sizes_nel, displs_nel, mpi_preal, MPI_COMM_WORLD, ierr)
#endif
        pp%ppn_i_s%weights_tot = sum(pp%ppn_i_s%weights)

        ! The i that we find with i_weight is i_ref which is mapped to i_cdet later. If i_weight is very close to 0
        ! (i.e. if there is not a with a valid excitaiton i->a), we need to make that weight finite as the mapped
        ! i_cdet might have allowed exciations. This is to prevent a bias where a valid i cannot be selected.
        call check_min_weight_ratio(pp%ppn_i_s%weights(:), pp%ppn_i_s%weights_tot, sys%nel, pp%power_pitzer_min_weight)
        
        ! Generate alias tables.
        call generate_alias_tables(sys%nel, pp%ppn_i_s%weights(:), pp%ppn_i_s%weights_tot, pp%ppn_i_s%aliasU(:), &
                                pp%ppn_i_s%aliasK(:))
        
        do i = iproc_nbasis_start, iproc_nbasis_end
            ia_s_weights_tot = 0.0_p
            ims = sys%basis%basis_fns(i)%ms
            ! Convert ims which is in {-1, +1} notation to imsa which is {1, 2}.
            imsa = (3 + ims)/2
            isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(i)%sym, sys%read_in%pg_sym%gamma_sym)
        
            !$omp parallel do default(none) &
            !$omp shared(sys,i,imsa,isyma,pp,ref,single_excitation_weight_ptr) &
            !$omp private(a) reduction(+:ia_s_weights_tot)
            do a = 1, sys%read_in%pg_sym%nbasis_sym_spin(imsa,isyma)
                pp%ppn_ia_s%weights(a, i) = 0.0_p
                if (sys%read_in%pg_sym%sym_spin_basis_fns(a,imsa, isyma) /= i) then
                    pp%ppn_ia_s%weights(a, i) = single_excitation_weight_ptr(sys, ref, i, &
                        sys%read_in%pg_sym%sym_spin_basis_fns(a,imsa,isyma))
                    if (pp%ppn_ia_s%weights(a, i) < depsilon) then
                        ! i-> a spin and symmetry allowed but weight assigned zero by our assumptions made when calculating
                        ! the weight. fix this by giving it a very small, but finite weight. The check_min_weight_ratio
                        ! routine called after this loop will convert nonzero finite weights that are less than the minimum
                        ! weight to a minimum weight. The min. weight depends on the other weights in the array, which is why
                        ! we cannot just set it here. It is important that this weight is non zero, otherwise it will not
                        ! be touched. We want i->a which is not symmetry and spin allowed to have a zero weight.
                        ! [todo] - factor in front of depsilon is a bit arbitary. 10depsilon is assumed to be less than
                        ! [todo] - the min. weight.
                        pp%ppn_ia_s%weights(a, i) = 10.0_p*depsilon
                    end if
                end if
                ia_s_weights_tot = ia_s_weights_tot + pp%ppn_ia_s%weights(a, i)
            end do
            !$omp end parallel do
            
            ! Use of tmp variable ia_s_weights_tot to make openmp happy.
            pp%ppn_ia_s%weights_tot(i) = ia_s_weights_tot
        
            ! Make sure that very small non zero weights get raised to the minimum weight.
            call check_min_weight_ratio(pp%ppn_ia_s%weights(:,i), pp%ppn_ia_s%weights_tot(i), &
                sys%read_in%pg_sym%nbasis_sym_spin(imsa,isyma), pp%power_pitzer_min_weight)
            call generate_alias_tables(sys%read_in%pg_sym%nbasis_sym_spin(imsa,isyma), pp%ppn_ia_s%weights(:,i), &
                                    pp%ppn_ia_s%weights_tot(i), pp%ppn_ia_s%aliasU(:,i), pp%ppn_ia_s%aliasK(:,i))
        end do
#ifdef PARALLEL
        ! note how FORTRAN stores arrays: array(2,1) comes before array(1,2) in memory.
        associate(mv=>maxval(sys%read_in%pg_sym%nbasis_sym_spin))
            call mpi_allgatherv(MPI_IN_PLACE, mv*sizes_nbasis(iproc), &
                mpi_preal, pp%ppn_ia_s%weights, mv*sizes_nbasis, mv*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, sizes_nbasis(iproc), &
                mpi_preal, pp%ppn_ia_s%weights_tot, sizes_nbasis, displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, mv*sizes_nbasis(iproc), &
                mpi_preal, pp%ppn_ia_s%aliasU, mv*sizes_nbasis, mv*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            ! [todo] - this is not the safest thing in the universe: Once someone changes whether aliasK is int_32 or int_64
            ! [todo] - in excit_gens.f90, this needs to be changed too.
            ! [todo] -  One could do a test on int_bas to be sure.
            call mpi_allgatherv(MPI_IN_PLACE, mv*sizes_nbasis(iproc), &
                MPI_INTEGER, pp%ppn_ia_s%aliasK, mv*sizes_nbasis, mv*displs_nbasis, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        end associate
#endif
        ! Double Excitations.
        ! Generate the i weighting lists and alias tables for possible (spin conserved - symmetry cannot be checked due to latter
        ! mapping which only conserves spin) double excitations.
        ! i and j are drawn from spinorbitals occupied in the reference and a and b can be all orbitals. i!=j and a!=b.
        ! i and j can equal a and b because we later map i. Single excitations are not considered. [todo]
        ! subroutine init_double_weights_ab is OpenMP parallelised.
        do i = iproc_nel_start, iproc_nel_end
            i_weight = 0.0_p
            
            do j = 1, sys%nel
                if (i /= j) then
                    call init_double_weights_ab(sys, pp%occ_list(i), pp%occ_list(j), i_weight)
                end if
            end do
            
            if (i_weight < depsilon) then
                ! because we map i later, even if i->a is not allowed/ has been wrongly assigned a zero weight,
                ! it needs a weight. It will be assigned a minimum weight by the routine check_in_weight_ratio
                ! called after this loop. The finite weight here should be very small (and will be raised in
                ! the min. weight function) but still detectable (> depsilon).
                ! [todo] - factor in front of depsilon is a bit arbitary. 10depsilon is assumed to be less than
                ! [todo] - the min. weight.
                i_weight = 10.0_p*depsilon
            end if
            pp%ppn_i_d%weights(i) = i_weight
        end do


#ifdef PARALLEL
        call mpi_allgatherv(MPI_IN_PLACE, sizes_nel(iproc), &
            mpi_preal, pp%ppn_i_d%weights, sizes_nel, displs_nel, mpi_preal, MPI_COMM_WORLD, ierr)
#endif
        
        pp%ppn_i_d%weights_tot = sum(pp%ppn_i_d%weights)

        ! The i that we find with i_weight is i_ref which is mapped to i_cdet later. If i_weight is very close to 0,
        ! we need to make that weight finite as the mapped i_cdet might have allowed excitations. This is to prevent
        ! a bias where a valid i cannot be selected.
        call check_min_weight_ratio(pp%ppn_i_d%weights(:), pp%ppn_i_d%weights_tot, sys%nel, pp%power_pitzer_min_weight)
        ! Generate alias tables for i.
        call generate_alias_tables(sys%nel, pp%ppn_i_d%weights(:), pp%ppn_i_d%weights_tot, pp%ppn_i_d%aliasU(:), &
                                pp%ppn_i_d%aliasK(:))

        ! Generate the j given i weighting lists and alias tables. a and b cannot equal i (they are drawn from the same set).
        ! subroutine init_double_weights_ab is OpenMP parallelised.
        do i = iproc_nbasis_start, iproc_nbasis_end
            ij_d_weights_tot = 0.0_p
            do j = 1, sys%nel
                ij_weight = 0.0_p
                if (pp%occ_list(j) /= i) then
                    call init_double_weights_ab(sys, i, pp%occ_list(j), ij_weight)
                end if
                if (ij_weight < depsilon) then
                    ! because we map i later, even if ij->ab is never allowed for any ab/ has been wrongly assigned
                    ! a zero weight, it needs a weight. It will be assigned a minimum weight by the routine
                    ! check_in_weight_ratio called after this loop. The finite weight here should be very small
                    ! (and will be raised in the min. weight function) but still detectable (> depsilon).
                    ! if i==j and that is why ij_weight is zero, set it to a nonzero finite weight since i=i_cdet but
                    ! j=j_ref here and in general i_cdet/=j_cdet if i==j.
                    ! [todo] - factor in front of depsilon is a bit arbitary. 10depsilon is assumed to be less than
                    ! [todo] - the min. weight.
                    ij_weight = 10.0_p*depsilon
                end if
                pp%ppn_ij_d%weights(j,i) = ij_weight
                ij_d_weights_tot = ij_d_weights_tot + ij_weight
            end do
            
            ! Use of tmp variable to keep OpenMP happy.
            pp%ppn_ij_d%weights_tot(i) = ij_d_weights_tot

            ! The j that we find with ij_weight is j_ref which is mapped to j_cdet later. If ij_weight is very close to 0,
            ! we need to make that weight finite as the mapped j_cdet might have allowed excitations. This is to prevent
            ! a bias where a valid j cannot be selected.
            call check_min_weight_ratio(pp%ppn_ij_d%weights(:,i), pp%ppn_ij_d%weights_tot(i), sys%nel, &
                                    pp%power_pitzer_min_weight)
            ! Generate alias tables.
            call generate_alias_tables(sys%nel, pp%ppn_ij_d%weights(:,i), pp%ppn_ij_d%weights_tot(i), pp%ppn_ij_d%aliasU(:,i), &
                                    pp%ppn_ij_d%aliasK(:,i))
        end do

#ifdef PARALLEL
        ! note how FORTRAN stores arrays: array(2,1) comes before array(1,2) in memory.
        associate(ne=>sys%nel)
            call mpi_allgatherv(MPI_IN_PLACE, ne*sizes_nbasis(iproc), &
                mpi_preal, pp%ppn_ij_d%weights, ne*sizes_nbasis, ne*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, sizes_nbasis(iproc), &
                mpi_preal, pp%ppn_ij_d%weights_tot, sizes_nbasis, displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, ne*sizes_nbasis(iproc), &
                mpi_preal, pp%ppn_ij_d%aliasU, ne*sizes_nbasis, ne*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, ne*sizes_nbasis(iproc), &
                MPI_INTEGER, pp%ppn_ij_d%aliasK, ne*sizes_nbasis, ne*displs_nbasis, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        end associate
#endif

        ! Now generate the occ->virtual weighting lists and alias tables.
        ! This breaks in the case where there is only one spinorbital in the system with a certain spin. This is unlikely
        ! and will be ignored. [todo] - check for case where sys%nbeta and sys%nalpha are 1?
        ! [todo] - do we need min weights here? probably not because we do not map and if <ai|ia> = 0 that is probably deserved?
        do i = iproc_nbasis_start, iproc_nbasis_end
            if (sys%basis%basis_fns(i)%Ms == -1) then ! beta
                nall = pp%n_all_beta
                if (nall > 0) then
                    call create_weighted_excitation_list_ptr(sys, i, i, pp%all_list_beta, nall, pp%ppn_ia_d%weights(:,i), &
                                                            pp%ppn_ia_d%weights_tot(i))
                    call generate_alias_tables(nall, pp%ppn_ia_d%weights(:,i), pp%ppn_ia_d%weights_tot(i), &
                                            pp%ppn_ia_d%aliasU(:,i), pp%ppn_ia_d%aliasK(:,i))
                end if
                do bsym = sys%sym0_tot, sys%sym_max_tot
                    if (sys%read_in%pg_sym%nbasis_sym_spin(1,bsym) > 0) then
                        call create_weighted_excitation_list_ptr(sys, i, i, sys%read_in%pg_sym%sym_spin_basis_fns(:,1,bsym), &
                            sys%read_in%pg_sym%nbasis_sym_spin(1,bsym), pp%ppn_jb_d%weights(:,bsym,i), &
                            pp%ppn_jb_d%weights_tot(bsym,i))
                        call generate_alias_tables(sys%read_in%pg_sym%nbasis_sym_spin(1,bsym), pp%ppn_jb_d%weights(:,bsym,i), &
                            pp%ppn_jb_d%weights_tot(bsym,i), pp%ppn_jb_d%aliasU(:,bsym,i), pp%ppn_jb_d%aliasK(:,bsym,i))
                    end if
                end do
            else ! alpha
                nall = pp%n_all_alpha
                if (nall > 0) then
                    call create_weighted_excitation_list_ptr(sys, i, i, pp%all_list_alpha, nall, pp%ppn_ia_d%weights(:,i), &
                                                             pp%ppn_ia_d%weights_tot(i))

                    call generate_alias_tables(nall, pp%ppn_ia_d%weights(:,i), pp%ppn_ia_d%weights_tot(i), &
                                            pp%ppn_ia_d%aliasU(:,i), pp%ppn_ia_d%aliasK(:,i))
                end if
                do bsym = sys%sym0_tot, sys%sym_max_tot
                    if (sys%read_in%pg_sym%nbasis_sym_spin(2,bsym) > 0) then
                        call create_weighted_excitation_list_ptr(sys, i, i, sys%read_in%pg_sym%sym_spin_basis_fns(:,2,bsym), &
                            sys%read_in%pg_sym%nbasis_sym_spin(2,bsym), pp%ppn_jb_d%weights(:,bsym,i), &
                            pp%ppn_jb_d%weights_tot(bsym,i))
                        call generate_alias_tables(sys%read_in%pg_sym%nbasis_sym_spin(2,bsym), pp%ppn_jb_d%weights(:,bsym,i), &
                            pp%ppn_jb_d%weights_tot(bsym,i), pp%ppn_jb_d%aliasU(:,bsym,i), pp%ppn_jb_d%aliasK(:,bsym,i))
                    end if
                end do
            end if
        end do
#ifdef PARALLEL
        sr = sys%sym_max_tot - sys%sym0_tot + 1
        associate(mv=>maxval(sys%read_in%pg_sym%nbasis_sym_spin), na=>max(pp%n_all_alpha,pp%n_all_beta))
            call mpi_allgatherv(MPI_IN_PLACE, na*sizes_nbasis(iproc), &
                mpi_preal, pp%ppn_ia_d%weights, na*sizes_nbasis, na*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, sizes_nbasis(iproc), &
                mpi_preal, pp%ppn_ia_d%weights_tot, sizes_nbasis, displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, na*sizes_nbasis(iproc), &
                mpi_preal, pp%ppn_ia_d%aliasU, na*sizes_nbasis, na*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, na*sizes_nbasis(iproc), &
                MPI_INTEGER, pp%ppn_ia_d%aliasK, na*sizes_nbasis, na*displs_nbasis, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        
            call mpi_allgatherv(MPI_IN_PLACE, &
                sr*mv*sizes_nbasis(iproc), mpi_preal, pp%ppn_jb_d%weights(:,sys%sym0_tot:sys%sym_max_tot,:), &
                sr*mv*sizes_nbasis, sr*mv*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, &
                sr*sizes_nbasis(iproc), mpi_preal, pp%ppn_jb_d%weights_tot(sys%sym0_tot:sys%sym_max_tot,:), sr*sizes_nbasis, &
                sr*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, &
                sr*mv*sizes_nbasis(iproc), mpi_preal, pp%ppn_jb_d%aliasU(:,sys%sym0_tot:sys%sym_max_tot,:), sr*mv*sizes_nbasis, &
                sr*mv*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, &
                sr*mv*sizes_nbasis(iproc), MPI_INTEGER, pp%ppn_jb_d%aliasK(:,sys%sym0_tot:sys%sym_max_tot,:), &
                sr*mv*sizes_nbasis, sr*mv*displs_nbasis, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        end associate
#endif

    end subroutine init_excit_mol_power_pitzer_orderN

    subroutine init_excit_mol_power_pitzer_orderM_ij(sys, ref, pp)

        ! Generate weights for i and j for the gen_excit_mol_power_pitzer_orderM_ij
        ! excitation generator based on heat bath approach (Holmes et al.).

        ! In:
        !    sys: system object being studied.
        !    ref: the reference from which we are exciting.
        ! In/Out:
        !    pp: an empty excit_gen_power_pitzer_t object which gets filled with
        !           the alias tables required to generate excitations.

        use checking, only: check_allocate
        use system, only: sys_t
        use qmc_data, only: reference_t
        use sort, only: qsort
        use proc_pointers, only: create_weighted_excitation_list_ptr, slater_condon2_excit_ptr
        use proc_pointers, only: abs_hmatel_ptr, single_excitation_weight_ptr
        use excit_gens, only: excit_gen_power_pitzer_t, alloc_alias_table_data_t
        use excit_gen_utils, only: init_double_weights_ab
        use alias, only: generate_alias_tables
        use read_in_symmetry, only: cross_product_basis_read_in
        use hamiltonian_data, only: hmatel_t
#ifdef PARALLEL
        use parallel

        integer :: displs_nbasis(0:nprocs-1)
        integer :: sizes_nbasis(0:nprocs-1)
        integer :: ierr
#endif
        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: ref
        type(excit_gen_power_pitzer_t), intent(inout) :: pp

        integer :: i, j
        integer :: iproc_nbasis_start, iproc_nbasis_end
        integer :: ierr_alloc

        real(p) :: i_weight_extra

        ! Store weights and alias tables.
        allocate(pp%ppm_i_d_weights(sys%basis%nbasis), stat=ierr_alloc)
        call check_allocate('pp%ppm_i_d_weights', sys%basis%nbasis, ierr_alloc)
        allocate(pp%ppm_ij_d_weights(sys%basis%nbasis,sys%basis%nbasis), stat=ierr_alloc)
        call check_allocate('pp%ppm_ij_d_weights', (sys%basis%nbasis*sys%basis%nbasis), ierr_alloc)
        
#ifdef PARALLEL
        ! Initialise do-loop range for each processor, [iproc_nbasis_start,iproc_nbasis_end].
        call get_proc_loop_range(sys%basis%nbasis, iproc_nbasis_start, iproc_nbasis_end, displs_nbasis, sizes_nbasis)
#else
        iproc_nbasis_start = 1
        iproc_nbasis_end = sys%basis%nbasis
#endif

        ! subroutine init_double_weights_ab is OpenMP parallelised.
        do i = iproc_nbasis_start, iproc_nbasis_end
            i_weight_extra = 0.0_p
            do j = 1, sys%basis%nbasis
                pp%ppm_ij_d_weights(j,i) = 0.0_p
                if (i /= j) then
                    call init_double_weights_ab(sys, i, j, pp%ppm_ij_d_weights(j,i))
                end if
                i_weight_extra = i_weight_extra + pp%ppm_ij_d_weights(j,i)
            end do
            pp%ppm_i_d_weights(i) = i_weight_extra
        end do

#ifdef PARALLEL
        ! note how FORTRAN stores arrays: array(2,1) comes before array(1,2) in memory.
        associate(nb=>sys%basis%nbasis)
            call mpi_allgatherv(MPI_IN_PLACE, sizes_nbasis(iproc), &
                mpi_preal, pp%ppm_i_d_weights, sizes_nbasis, displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
            call mpi_allgatherv(MPI_IN_PLACE, nb*sizes_nbasis(iproc), &
                mpi_preal, pp%ppm_ij_d_weights, nb*sizes_nbasis, nb*displs_nbasis, mpi_preal, MPI_COMM_WORLD, ierr)
        end associate
#endif

    end subroutine init_excit_mol_power_pitzer_orderM_ij

    subroutine gen_excit_mol_power_pitzer_occ_ref(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        ! Weight the double excitations according the the Power-Pitzer bound
        ! <ij|ab> <= Sqrt(<ia|ai><jb|bj>), see J.D. Power, R.M. Pitzer, Chem. Phys. Lett.,
        ! 478-483 (1974).
        ! This is an O(M/64) version which pretends the determinant excited from is the reference,
        ! then modifies the selected orbitals to be those of the determinant given.
        ! Each occupied and virtual not present in the det we're actually given will be
        ! mapped to the one of the equivalent free numerical index in the reference.

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Excitation generation data.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
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
        use excitations, only: find_excitation_permutation2,get_excitation_locations
        use proc_pointers, only: slater_condon2_excit_ptr
        use system, only: sys_t
        use excit_gen_mol, only: gen_single_excit_mol_no_renorm, choose_ij_mol
        use excit_gens, only: excit_gen_power_pitzer_t, excit_gen_data_t
        use alias, only: select_weighted_value_precalc
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use hamiltonian_data, only: hmatel_t
        use read_in_symmetry, only: cross_product_basis_read_in
        use search, only: binary_search

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation


        integer ::  ij_spin, ij_sym, i_ref_j_ref_sym, imsb, isyma, isymb
        
        ! We distinguish here between a_ref and a_cdet. a_ref is a in the world where the reference is fully occupied
        ! and we excite from the reference. We then map a_ref onto a_cdet which is a in the world where cdet is fully
        ! occupied (which is the world we are actually in). This mapping is a one-to-one mapping.

        ! a_ind_ref is the index of a in the world where the reference is fully occupied, etc.

        ! a_ind_rev_cdet is the index of a_cdet if a had been chosen after b (relevant when both have the same spin
        ! and the reverse selection has to be considered).
        
        integer :: i_ind_ref, j_ind_ref, a_ind_ref, b_ind_cdet
        integer :: a_ind_rev_cdet, b_ind_rev_ref
        integer :: i_ref, j_ref, a_ref, b_ref, i_cdet, j_cdet, a_cdet, b_cdet

        integer :: nex
        integer :: cdet_store(sys%nel)
        integer :: ref_store(sys%nel)
        integer :: ii, jj, t
        real(p) :: pgen_ij

        logical :: found, a_found

        ! 1. Select single or double.
        if (get_rand_close_open(rng) < excit_gen_data%pattempt_single) then

            call gen_single_excit_mol_no_renorm(rng, sys, excit_gen_data%pattempt_single, cdet, pgen, connection, hmatel, &
                                            allowed_excitation)

        else
            ! We have a double
            associate( pp => excit_gen_data%excit_gen_pp )

                ! 2b. Select orbitals to excite from
                
                call choose_ij_mol(rng, sys, pp%occ_list, i_ind_ref, j_ind_ref, i_ref, j_ref, i_ref_j_ref_sym, ij_spin, pgen_ij)

                ! At this point we pretend we're the reference, and fix up mapping ref's orbitals to cdet's orbitals later.
                i_ref = pp%occ_list(i_ind_ref)
                j_ref = pp%occ_list(j_ind_ref)

                ! We now need to select the orbitals to excite into which we do with weighting:
                ! p(ab|ij) = p(a|i) p(b|j) + p(a|j) p(b|i)
                ! We actually choose a|i then b|j, but since we could have also generated the excitation b from i and a from j,
                ! we need to include that prob too.

                ! Given i_ref, use the alias method to select a_ref with appropriate probability from the set of orbitals
                ! of the same spin as i_ref that are unoccupied if all electrons are in the reference.
                if (sys%basis%basis_fns(i_ref)%Ms < 0) then
                    if (sys%nvirt_beta > 0) then
                        a_ind_ref = select_weighted_value_precalc(rng, sys%nvirt_beta, pp%pp_ia_d%aliasU(:,i_ind_ref), &
                                                               pp%pp_ia_d%aliasK(:,i_ind_ref))
                        a_ref = pp%virt_list_beta(a_ind_ref)
                        a_found = .true.
                    else
                        a_found = .false.
                    end if
                else
                    if (sys%nvirt_alpha > 0) then
                        a_ind_ref = select_weighted_value_precalc(rng, sys%nvirt_alpha, pp%pp_ia_d%aliasU(:,i_ind_ref), &
                                                               pp%pp_ia_d%aliasK(:,i_ind_ref))
                        a_ref = pp%virt_list_alpha(a_ind_ref)
                        a_found = .true.
                    else
                        a_found = .false.
                    end if
                end if

                if (a_found) then
                    ! To conserve total spin, b and j will have the same spin, as a and i have the same spin.
                    ! To find what symmetry b should have, we first have to map i,j and a to what they correspond to
                    ! had we considered cdet and not the reference.
                
                    ! We need a list of the substitutions in cdet vs the ref.  i.e. for each orbital in the ref-lined-up
                    ! cdet which is not in ref, we need the location.
                    ! This is currently done with an O(N) step, but might be sped up at least.

                    call get_excitation_locations(pp%occ_list, cdet%occ_list, ref_store, cdet_store, sys%nel, nex)
                    ! These orbitals might not be aligned in the most efficient way:
                    !  They may not match in spin, so first deal with this

                    ! ref store (e.g.) contains the indices within pp%occ_list of the orbitals
                    ! which have been excited from.
                    ! [todo] - Consider having four lists separated by spin instead of two, i.e. ref_store_alpha and
                    ! [todo] - ref_store_beta instead of ref_store and the same for cdet_store.
                    ! [todo] - Consider calculating these lists once per determinants/excitor instead of once per excitation
                    ! [todo] - generator call.
                    do ii=1, nex
                        associate(bfns=>sys%basis%basis_fns)
                            if (bfns(pp%occ_list(ref_store(ii)))%Ms /= bfns(cdet%occ_list(cdet_store(ii)))%Ms) then
                                jj = ii + 1
                                do while (bfns(pp%occ_list(ref_store(ii)))%Ms /= bfns(cdet%occ_list(cdet_store(jj)))%Ms)
                                    jj = jj + 1
                                end do
                                ! det's jj now points to an orb of the same spin as ref's ii, so swap cdet_store's ii and jj.
                                t = cdet_store(ii)
                                cdet_store(ii) = cdet_store(jj)
                                cdet_store(jj) = t
                            end if
                        end associate
                    end do
                
                    ! Now see if i_ref and j_ref are in ref_store or a_ref in cdet_store and map appropriately
                    i_cdet = i_ref
                    j_cdet = j_ref
                    a_cdet = a_ref
                    ! ref_store  contains the indices within pp%occ_list of the orbitals which are excited out of ref into cdet
                    ! cdet_store  contains the indices within cdet%occ_list of the orbitals which are in cdet (excited out
                    ! of ref).  i_ind_ref and j_ind_ref are the indices of the orbitals in pp%occ_list which we're exciting from.
                    do ii=1, nex
                        if (ref_store(ii) == i_ind_ref) then  
                            ! i_ref isn't actually in cdet, so we assign i_cdet to the orb that is
                            i_cdet = cdet%occ_list(cdet_store(ii))
                        else if (ref_store(ii) == j_ind_ref) then 
                            ! j_ref isn't actually in cdet, so we assign j_cdet to the orb that is
                            j_cdet = cdet%occ_list(cdet_store(ii))
                        end if
                        if (cdet%occ_list(cdet_store(ii)) == a_ref) then
                            ! a_ref is occupied in cdet, assign a_cdet to the orb that is not
                            a_cdet = pp%occ_list(ref_store(ii)) 
                        end if
                    end do

                    ! The symmetry of b (=b_cdet), isymb, is given by
                    ! (sym_i_cdet* x sym_j_cdet* x sym_a_cdet)* = sym_b_cdet
                    ! (at least for Abelian point groups)
                    ! ij_sym: symmetry conjugate of the irreducible representation spanned by the codensity
                    !        \phi_i_cdet*\phi_j_cdet. (We assume that ij is going to be in the bra of the excitation.)
                    ij_sym = sys%read_in%sym_conj_ptr(sys%read_in, cross_product_basis_read_in(sys, i_cdet, j_cdet))

                    isymb = sys%read_in%sym_conj_ptr(sys%read_in, &
                                sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, sys%basis%basis_fns(a_cdet)%sym))
                    ! Ms_i + Ms_j = Ms_a + Ms_b => Ms_b = Ms_i + Ms_j - Ms_a.
                    ! Given that Ms_a = Ms_i, Ms_b = Ms_j.
                    ! Ms_k is +1 if up and -1 if down but imsb is +2 if up and +1 if down,
                    ! therefore a conversion is necessary.
                    imsb = (sys%basis%basis_fns(j_ref)%Ms + 3)/2
                end if

                ! Check whether an orbital (occupied or not) with required spin and symmetry exists.
                if (a_found .and. (sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb) > 0)) then
                    ! Using alias tables based on the reference, find b_cdet out of the set of (all) orbitals that have
                    ! the required spin and symmetry. Note that these orbitals might be occupied in the reference and/or
                    ! cdet (although we only care about whether they are occupied in cdet which we deal with later).

                    b_ind_cdet = select_weighted_value_precalc(rng, sys%read_in%pg_sym%nbasis_sym_spin(imsb, isymb), &
                                        pp%pp_jb_d%aliasU(:, isymb, j_ind_ref), pp%pp_jb_d%aliasK(:, isymb, j_ind_ref))
                    b_cdet = sys%read_in%pg_sym%sym_spin_basis_fns(b_ind_cdet, imsb, isymb)

                    ! Check that a_cdet /= b_cdet and that b_cdet is not occupied in cdet:
                    if (a_cdet /= b_cdet .and. .not.btest(cdet%f(sys%basis%bit_lookup(2,b_cdet)), &
                        sys%basis%bit_lookup(1,b_cdet))) then
                        
                        ! 3b. Probability of generating this excitation.

                        ! Calculate p(ab|ij) = p(a|i) p(j|b) + p(b|i)p(a|j)
                        if (ij_spin == 0) then
                            ! Not possible to have chosen the reversed excitation.
                            pgen = pp%pp_ia_d%weights(a_ind_ref, i_ind_ref) / pp%pp_ia_d%weights_tot(i_ind_ref) &
                                    * pp%pp_jb_d%weights(b_ind_cdet, isymb, j_ind_ref) / pp%pp_jb_d%weights_tot(isymb, j_ind_ref)
                        else
                            ! i and j have same spin, so could have been selected in the other order.
                            ! Need to find b_ref, the orbital b would have been in the world where we focus on the reference.
                            b_ref = b_cdet

                            do ii=1, nex
                                if (pp%occ_list(ref_store(ii)) == b_cdet) then
                                    ! b_cdet is occupied in ref, assign b_ref to the orb that is not
                                    b_ref = cdet%occ_list(cdet_store(ii))
                                    exit
                                end if
                            end do

                            if (imsb == 1) then
                                ! find index b as if we had it selected first and as a from list of unoccupied virtual orbitals.
                                call binary_search(pp%virt_list_beta, b_ref, 1, sys%nvirt_beta, found, b_ind_rev_ref)
                            else
                                call binary_search(pp%virt_list_alpha, b_ref, 1, sys%nvirt_alpha, found, b_ind_rev_ref)
                            end if
                            isyma = sys%read_in%sym_conj_ptr(sys%read_in, &
                                        sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, isymb))
                            ! imsa = imsb
                            call binary_search(sys%read_in%pg_sym%sym_spin_basis_fns(:,imsb,isyma), a_cdet, 1, &
                                    sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma), found, a_ind_rev_cdet)

                            pgen = pp%pp_ia_d%weights(a_ind_ref, i_ind_ref) / pp%pp_ia_d%weights_tot(i_ind_ref) &
                                * pp%pp_jb_d%weights(b_ind_cdet, isymb, j_ind_ref) / pp%pp_jb_d%weights_tot(isymb, j_ind_ref) &
                                +  pp%pp_ia_d%weights(b_ind_rev_ref, i_ind_ref) / pp%pp_ia_d%weights_tot(i_ind_ref) &
                                * pp%pp_jb_d%weights(a_ind_rev_cdet, isyma, j_ind_ref) / pp%pp_jb_d%weights_tot(isyma, j_ind_ref)
                        end if

                        pgen = excit_gen_data%pattempt_double * pgen * pgen_ij ! pgen(ab)
                        connection%nexcit = 2
                        allowed_excitation = .true.
                    else
                        allowed_excitation = .false.
                    end if
                else
                    allowed_excitation = .false.
                end if

                if (allowed_excitation) then

                    if (i_cdet<j_cdet) then
                        connection%from_orb(1) = i_cdet
                        connection%from_orb(2) = j_cdet
                    else
                        connection%from_orb(2) = i_cdet
                        connection%from_orb(1) = j_cdet
                    end if
                    if (a_cdet<b_cdet) then
                        connection%to_orb(1) = a_cdet
                        connection%to_orb(2) = b_cdet
                    else
                        connection%to_orb(2) = a_cdet
                        connection%to_orb(1) = b_cdet
                    end if

                    ! 4b. Parity of permutation required to line up determinants.
                    ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                    call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

                    ! 5b. Find the connecting matrix element.
                    hmatel = slater_condon2_excit_ptr(sys, connection%from_orb(1), connection%from_orb(2), &
                                                      connection%to_orb(1), connection%to_orb(2), connection%perm)
                else
                    ! Carelessly selected ij with no possible excitations.  Such
                    ! events are not worth the cost of renormalising the generation
                    ! probabilities.
                    ! Return a null excitation.
                    hmatel%c = cmplx(0.0_p, 0.0_p, p)
                    hmatel%r = 0.0_p
                    pgen = 1.0_p
                end if
            end associate
        end if

    end subroutine gen_excit_mol_power_pitzer_occ_ref

    subroutine gen_excit_mol_power_pitzer_orderN(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        ! Weight the double excitations according the the Power-Pitzer bound
        ! <ij|ab> <= Sqrt(<ia|ai><jb|bj>), see J.D. Power, R.M. Pitzer, Chem. Phys. Lett.,
        ! 478-483 (1974).

        ! A short overview over the method:
        ! ---------------------------------
        ! We have a reference and a simulation frame of reference. In the reference frame, the
        ! reference determinant is occupied and in the simulation frame, cdet is occupied (which
        ! is actually correct). i_ref is i in the reference frame which can be mapped to i_cdet in
        ! simulation frame.
        ! 1. Decide whether to do a single or a double excitation using p_single.

        ! Single excitation:
        ! 2. Select i_ref from the set of occupied orbitals in the reference using pre-computed
        !    alias tables pp%ppn_i_s%alias{K,U}(:).
        ! 3. Done in decoder: Map i_ref to i_cdet. This is a one-to-one mapping. All spinorbitals that are occupied
        !    in one of the two frames but not the other are separated by spin and then ordered.
        !    If i_ref is occupied in both frames, i_ref = i_cdet. If not, a spinorbital of the same
        !    spin is found that is at the same position in the lists of differing orbitals.
        ! 4. Find a_cdet from all orbitals of required spin and symmetry using pre-calculated alias
        !    information pp%ppn_ia_s%alias{K,U}(:, i_cdet).
        ! 5. Calculate pgen.

        ! Double excitation:
        ! 2. Select i_ref from occupied orbitals in the reference using pre-calculated alias tables
        !    pp%ppn_i_d%alias{K,U}.
        ! 3. Done in decoder: Map i_ref to i_cdet. This is a one-to-one mapping. All spinorbitals that are occupied
        !    in one of the two frames but not the other are separated by spin and then ordered.
        !    If i_ref is occupied in both frames, i_ref = i_cdet. If not, a spinorbital of the same
        !    spin is found that is at the same position in the lists of differing orbitals.
        ! 4. Find j_ref from occupied orbitals in the reference using pre-calculated alias tables
        !    pp%ppn_ij_d%alias{K,U}(:,i_cdet).
        ! 5. Done in decoder: Map j_ref to j_cdet.
        ! 6. Select a_cdet from all orbitals of the same spin as i using pre-calculated alias tables
        !    pp%ppn_ia_d%alias{K,U}(:,i_cdet). If a_cdet is occupied, excitation is forbidden.
        ! 7. Select b_cdet from all orbitals of required spin and symmetry isymb using pre-calculated
        !    alias tables pp%ppn_jb_d%alias{K,U}(:, isymb, j_cdet). If b is occupied, excitation is
        !    forbidden.
        ! 8. Calculate pgen, considering that if ij are of the same spin. b could have been selected
        !    from alias tables pp%ppn_ia_d%alias{K,U}(:,i_cdet) and a from
        !    pp%ppn_jb_d%alias{K,U}(:, isymb', j_cdet).

        ! In:
        !    sys: system object being studied.
        !    excit_gen_data: Excitation generation data.
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
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
        use excitations, only: find_excitation_permutation1, find_excitation_permutation2
        use proc_pointers, only: slater_condon1_excit_ptr, slater_condon2_excit_ptr
        use system, only: sys_t
        use excit_gens, only: excit_gen_power_pitzer_t, excit_gen_data_t
        use alias, only: select_weighted_value_precalc
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use hamiltonian_data, only: hmatel_t
        use read_in_symmetry, only: cross_product_basis_read_in
        use search, only: binary_search

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation


        integer ::  ij_spin, ij_sym, imsa, imsb, isyma, isymb
        
        ! We distinguish here between a_ref and a_cdet. a_ref is a in the world where the reference is fully occupied
        ! and we excite from the reference. We then map a_ref onto a_cdet which is a in the world where cdet is fully
        ! occupied (which is the world we are actually in). This mapping is a one-to-one mapping.

        ! a_ind_ref is the index of a in the world where the reference is fully occupied, etc.

        ! a_ind_rev_cdet is the index of a_cdet if a had been chosen after b (relevant when both have the same spin
        ! and the reverse selection has to be considered).
        
        integer :: i_ind_ref, j_ind_ref, a_ind_cdet, b_ind_cdet
        integer :: a_ind_rev_cdet, b_ind_rev_cdet
        integer :: i_cdet, j_cdet, a_cdet, b_cdet, j_cdet_tmp

        logical :: found

        ! 1. Select single or double.
        if (get_rand_close_open(rng) < excit_gen_data%pattempt_single) then  
            ! We have a single
            associate( pp => excit_gen_data%excit_gen_pp )
                ! 2. Select i_ref.
                
                i_ind_ref = select_weighted_value_precalc(rng, sys%nel, pp%ppn_i_s%aliasU(:), pp%ppn_i_s%aliasK(:))
                i_cdet = cdet%ref_cdet_occ_list(i_ind_ref)
                
                ! 4. Select a_cdet.
                ! Convert ims which is in {-1, +1} notation to imsa which is {1, 2}.
                imsa = (3 + sys%basis%basis_fns(i_cdet)%ms)/2
                isyma = sys%read_in%cross_product_sym_ptr(sys%read_in, sys%basis%basis_fns(i_cdet)%sym, &
                                                        sys%read_in%pg_sym%gamma_sym)
                if (sys%read_in%pg_sym%nbasis_sym_spin(imsa,isyma) > 0) then
                    a_ind_cdet = select_weighted_value_precalc(rng, sys%read_in%pg_sym%nbasis_sym_spin(imsa,isyma), &
                                                    pp%ppn_ia_s%aliasU(:, i_cdet), pp%ppn_ia_s%aliasK(:, i_cdet))
                    a_cdet = sys%read_in%pg_sym%sym_spin_basis_fns(a_ind_cdet, imsa, isyma)
                    if (.not. btest(cdet%f(sys%basis%bit_lookup(2,a_cdet)), sys%basis%bit_lookup(1,a_cdet))) then
                        allowed_excitation = .true.
                    else
                        allowed_excitation = .false.
                    end if
                else
                    allowed_excitation = .false.
                end if

                if (allowed_excitation) then
                    ! 5. Calculate pgen.
                    pgen = (pp%ppn_i_s%weights(i_ind_ref)/pp%ppn_i_s%weights_tot) * &
                            (pp%ppn_ia_s%weights(a_ind_cdet, i_cdet)/pp%ppn_ia_s%weights_tot(i_cdet))
                    pgen = excit_gen_data%pattempt_single * pgen ! pgen(ab)
                    connection%nexcit = 1
                    connection%to_orb(1) = a_cdet
                    connection%from_orb(1) = i_cdet
                    call find_excitation_permutation1(sys%basis%excit_mask, cdet%f, connection)
                    hmatel = slater_condon1_excit_ptr(sys, cdet%occ_list, connection%from_orb(1), &
                                          connection%to_orb(1), connection%perm)
                else
                    ! We have not found a valid excitation i -> a.  Such
                    ! events are not worth the cost of renormalising the generation
                    ! probabilities.
                    ! Return a null excitation.
                    hmatel%c = cmplx(0.0_p, 0.0_p, p)
                    hmatel%r = 0.0_p
                    pgen = 1.0_p
                end if

            end associate

        else
            ! We have a double
            associate( pp => excit_gen_data%excit_gen_pp )
                ! 2. Select i_ref.
                
                i_ind_ref = select_weighted_value_precalc(rng, sys%nel, pp%ppn_i_d%aliasU(:), pp%ppn_i_d%aliasK(:))
                i_cdet = cdet%ref_cdet_occ_list(i_ind_ref)
                
                ! 4. Find j_ref.
                ! j is chosen from set of occupied spinorbitals in reference but with weights calculated for i_cdet.
                j_ind_ref = select_weighted_value_precalc(rng, sys%nel, pp%ppn_ij_d%aliasU(:,i_cdet), &
                                                    pp%ppn_ij_d%aliasK(:,i_cdet))
                j_cdet = cdet%ref_cdet_occ_list(j_ind_ref)

                if (j_cdet /= i_cdet) then
                    allowed_excitation = .true. ! for now - that can change.
                    pgen = (pp%ppn_i_d%weights(i_ind_ref) / pp%ppn_i_d%weights_tot) * &
                            (pp%ppn_ij_d%weights(j_ind_ref,i_cdet) / pp%ppn_ij_d%weights_tot(i_cdet))
                    ij_spin = sys%basis%basis_fns(i_cdet)%Ms + sys%basis%basis_fns(j_cdet)%Ms
                    ! Could have selected i and j the other way around. 
                    pgen = pgen + ((pp%ppn_i_d%weights(j_ind_ref) / pp%ppn_i_d%weights_tot) * &
                            (pp%ppn_ij_d%weights(i_ind_ref,j_cdet) / pp%ppn_ij_d%weights_tot(j_cdet)))
                    ! Order i and j such that i<j.
                    if (j_cdet < i_cdet) then
                        j_cdet_tmp = j_cdet
                        j_cdet = i_cdet
                        i_cdet = j_cdet_tmp
                    end if
                else
                    allowed_excitation = .false.
                end if
                if (allowed_excitation) then
                    ! We now need to select the orbitals to excite into which we do with weighting:
                    ! p(ab|ij) = p(a|i) p(b|j) + p(a|j) p(b|i).
                    ! We actually choose a|i then b|j, but since we could have also generated the excitation b from i and a from
                    ! j, we need to include that too if they have the same spin.

                    ! 6. Find a_cdet.
                    ! Given i_cdet, use the alias method to select a_cdet with appropriate probability from the set of orbitals
                    ! of the same spin as i_cdet.
                    if (sys%basis%basis_fns(i_cdet)%Ms < 0) then
                        a_ind_cdet = select_weighted_value_precalc(rng, pp%n_all_beta, pp%ppn_ia_d%aliasU(:,i_cdet), &
                                                               pp%ppn_ia_d%aliasK(:,i_cdet))
                        a_cdet = pp%all_list_beta(a_ind_cdet)
                    else
                        a_ind_cdet = select_weighted_value_precalc(rng, pp%n_all_alpha, pp%ppn_ia_d%aliasU(:,i_cdet), &
                                                               pp%ppn_ia_d%aliasK(:,i_cdet))
                        a_cdet = pp%all_list_alpha(a_ind_cdet)
                    end if
                    ! Make sure a_cdet is a virtual orbital (not occupied in c_det).
                    if (btest(cdet%f(sys%basis%bit_lookup(2,a_cdet)), sys%basis%bit_lookup(1,a_cdet))) then
                        allowed_excitation = .false.
                    end if
                end if
                if (allowed_excitation) then
                    ! 7. Select b_cdet.
                    ! To conserve total spin, b and j will have the same spin, as a and i have the same spin.

                    ! The symmetry of b (=b_cdet), isymb, is given by
                    ! (sym_i_cdet* x sym_j_cdet* x sym_a_cdet)* = sym_b_cdet
                    ! (at least for Abelian point groups)
                    ! ij_sym: symmetry conjugate of the irreducible representation spanned by the codensity
                    !        \phi_i_cdet*\phi_j_cdet. (We assume that ij is going to be in the bra of the excitation.)
                    ij_sym = sys%read_in%sym_conj_ptr(sys%read_in, cross_product_basis_read_in(sys, i_cdet, j_cdet))

                    isymb = sys%read_in%sym_conj_ptr(sys%read_in, &
                                sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, sys%basis%basis_fns(a_cdet)%sym))
                    ! Ms_i + Ms_j = Ms_a + Ms_b => Ms_b = Ms_i + Ms_j - Ms_a.
                    ! Given that Ms_a = Ms_i, Ms_b = Ms_j.
                    ! Ms_k is +1 if up and -1 if down but imsb is +2 if up and +1 if down,
                    ! therefore a conversion is necessary.
                    imsb = (sys%basis%basis_fns(j_cdet)%Ms + 3)/2

                    if (sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb) == 0) then
                        allowed_excitation = .false.
                    end if
                end if
                ! Check whether an orbital (occupied or not) with required spin and symmetry exists.
                if (allowed_excitation) then
                    ! Find b_cdet out of the set of (all) orbitals that have the required spin and symmetry for given j_cdet. 
                    ! Note that these orbitals might be occupied.
                    b_ind_cdet = select_weighted_value_precalc(rng, sys%read_in%pg_sym%nbasis_sym_spin(imsb, isymb), &
                                        pp%ppn_jb_d%aliasU(:, isymb, j_cdet), pp%ppn_jb_d%aliasK(:, isymb, j_cdet))
                    b_cdet = sys%read_in%pg_sym%sym_spin_basis_fns(b_ind_cdet, imsb, isymb)

                    ! Check that a_cdet /= b_cdet and that b_cdet is not occupied in cdet:
                    if ((a_cdet /= b_cdet) .and. .not. btest(cdet%f(sys%basis%bit_lookup(2,b_cdet)), &
                        sys%basis%bit_lookup(1,b_cdet))) then
                        
                        ! 8. Probability of generating this excitation.

                        ! Calculate p(ab|ij) = p(a|i) p(j|b) + p(b|i)p(a|j)
                        if (ij_spin == 0) then
                            ! Not possible to have chosen the reversed excitation.
                            pgen = pgen * (pp%ppn_ia_d%weights(a_ind_cdet, i_cdet) / pp%ppn_ia_d%weights_tot(i_cdet) &
                                    * pp%ppn_jb_d%weights(b_ind_cdet, isymb, j_cdet) / pp%ppn_jb_d%weights_tot(isymb, j_cdet))
                        else
                            ! i and j have same spin, so could have been selected in the other order.

                            if (imsb == 1) then
                                ! find index b as if we had it selected first and as a from list of beta orbitals.
                                call binary_search(pp%all_list_beta, b_cdet, 1, pp%n_all_beta, found, b_ind_rev_cdet)
                            else
                                call binary_search(pp%all_list_alpha, b_cdet, 1, pp%n_all_alpha, found, b_ind_rev_cdet)
                            end if
                            isyma = sys%read_in%sym_conj_ptr(sys%read_in, &
                                        sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, isymb))
                            ! imsa = imsb
                            call binary_search(sys%read_in%pg_sym%sym_spin_basis_fns(:,imsb,isyma), a_cdet, 1, &
                                    sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma), found, a_ind_rev_cdet)

                            pgen = pgen * (pp%ppn_ia_d%weights(a_ind_cdet, i_cdet) / pp%ppn_ia_d%weights_tot(i_cdet) &
                                * pp%ppn_jb_d%weights(b_ind_cdet, isymb, j_cdet) / pp%ppn_jb_d%weights_tot(isymb, j_cdet) &
                                +  pp%ppn_ia_d%weights(b_ind_rev_cdet, i_cdet) / pp%ppn_ia_d%weights_tot(i_cdet) &
                                * pp%ppn_jb_d%weights(a_ind_rev_cdet, isyma, j_cdet) / pp%ppn_jb_d%weights_tot(isyma, j_cdet))
                        end if

                        pgen = excit_gen_data%pattempt_double * pgen ! pgen(ab)
                        connection%nexcit = 2
                    else
                        allowed_excitation = .false.
                    end if
                end if


                if (allowed_excitation) then
                    ! i and j are ordered.
                    connection%from_orb(1) = i_cdet
                    connection%from_orb(2) = j_cdet
                    if (a_cdet<b_cdet) then
                        connection%to_orb(1) = a_cdet
                        connection%to_orb(2) = b_cdet
                    else
                        connection%to_orb(2) = a_cdet
                        connection%to_orb(1) = b_cdet
                    end if
                    ! 4b. Parity of permutation required to line up determinants.
                    ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                    call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

                    ! 5b. Find the connecting matrix element.
                    hmatel = slater_condon2_excit_ptr(sys, connection%from_orb(1), connection%from_orb(2), &
                                                      connection%to_orb(1), connection%to_orb(2), connection%perm)
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

    end subroutine gen_excit_mol_power_pitzer_orderN

    subroutine gen_excit_mol_power_pitzer_occ(rng, sys, excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.
        ! Weight the double excitations according the the Power-Pitzer bound
        ! <ij|ab> <= Sqrt(<ia|ai><jb|bj>), see J.D. Power, R.M. Pitzer, Chem. Phys. Lett.,
        ! 478-483 (1974).
        ! This requires a lookup of O(M) two-electron integrals in its setup.
        ! This calculates weights on-the-fly. Single excitations are currently treated uniformly and ij are also
        ! selected uniformly currently or according to heat bath.

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

        use checking, only: check_deallocate
        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use excitations, only: find_excitation_permutation2
        use proc_pointers, only: slater_condon2_excit_ptr, create_weighted_excitation_list_ptr
        use system, only: sys_t
        use excit_gen_mol, only: gen_single_excit_mol, choose_ij_mol
        use hamiltonian_data, only: hmatel_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use search, only: binary_search
        use checking, only: check_allocate, check_deallocate
        use excit_gens, only: excit_gen_power_pitzer_t, excit_gen_data_t
        use excit_gen_utils, only: select_ij_heat_bath
        use alias, only: select_weighted_value
        use read_in_symmetry, only: cross_product_basis_read_in
        use qmc_data, only: excit_gen_power_pitzer_occ_ij

        type(sys_t), intent(in) :: sys
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen
        type(hmatel_t), intent(out) :: hmatel
        type(excit_t), intent(out) :: connection
        logical, intent(out) :: allowed_excitation

        integer ::  ij_spin, ij_sym, ierr
        logical :: found, a_found
        real(p), allocatable :: ia_weights(:), ja_weights(:), jb_weights(:)
        real(p) :: ia_weights_tot, ja_weights_tot, jb_weights_tot
        real(p) :: ij_weights_occ(sys%nel), ji_weights_occ(sys%nel)
        real(p) :: ij_weights_occ_tot, ji_weights_occ_tot
        real(p) :: pgen_ij
        integer :: a, b, i, j, j_tmp, a_ind, b_ind, a_ind_rev, b_ind_rev, i_ind, j_ind, isymb, imsb, isyma
        integer :: ierr_dealloc

        ! 1. Select single or double.
        if (get_rand_close_open(rng) < excit_gen_data%pattempt_single) then
            call gen_single_excit_mol(rng, sys, excit_gen_data%pattempt_single, cdet, pgen, connection, hmatel, &
                                        allowed_excitation)
        else
            ! We have a double

            ! 2. Select orbitals to excite from
            if (excit_gen_data%excit_gen == excit_gen_power_pitzer_occ_ij) then
                ! Select ij using heat bath excit. gen. techniques.

                call select_ij_heat_bath(rng, sys%nel, excit_gen_data%excit_gen_pp%ppm_ij_d_weights, cdet, i, j, i_ind, j_ind, &
                    ij_weights_occ, ij_weights_occ_tot, ji_weights_occ, ji_weights_occ_tot, allowed_excitation)

                ij_spin = sys%basis%basis_fns(i)%Ms + sys%basis%basis_fns(j)%Ms
                ! ij_sym: symmetry conjugate of the irreducible representation spanned by the codensity
                !        \phi_i_cdet*\phi_j_cdet. (We assume that ij is going to be in the bra of the excitation.)
                ij_sym = sys%read_in%sym_conj_ptr(sys%read_in, cross_product_basis_read_in(sys, i, j))
                
                pgen_ij = ((cdet%i_d_weights_occ(i_ind)/cdet%i_d_weights_occ_tot) * (ij_weights_occ(j_ind)/ij_weights_occ_tot)) + &
                    ((cdet%i_d_weights_occ(j_ind)/cdet%i_d_weights_occ_tot) * (ji_weights_occ(i_ind)/ji_weights_occ_tot))

                ! Sort i and j such that j>i.
                if (j < i) then
                    j_tmp = j
                    j = i
                    i = j_tmp
                end if

            else
                ! Select ij uniformly.
                call choose_ij_mol(rng, sys, cdet%occ_list, i_ind, j_ind, i, j, ij_sym, ij_spin, pgen_ij)
                allowed_excitation = .true.
            end if

            ! Now we've chosen i and j. They are ordered, j>i.
 
            ! We now need to select the orbitals to excite into which we do with weighting:
            ! p(ab|ij) = p(a|i) p(b|j) + p(a|j) p(b|i)
            
            ! We actually choose a|i then b|j, but since we could have also generated the excitation b from i and a from j,
            ! we need to include that prob too.

            ! 3. Find a.
            ! Given i, construct the weights of all possible a
            if (allowed_excitation) then
                if (sys%basis%basis_fns(i)%Ms < 0) then
                    ! [todo] - Consider doing a binary/linear search instead of using the alias method.
                    if (sys%nvirt_beta > 0) then
                        allocate(ia_weights(1:sys%nvirt_beta), stat=ierr)
                        call check_allocate('ia_weights', sys%nvirt_beta, ierr)
                        call create_weighted_excitation_list_ptr(sys, i, 0, cdet%unocc_list_beta, sys%nvirt_beta, ia_weights, &
                                                             ia_weights_tot)
                        if (ia_weights_tot > 0.0_p) then
                            ! Use the alias method to select i with the appropriate probability
                            a_ind = select_weighted_value(rng, sys%nvirt_beta, ia_weights, ia_weights_tot)
                            a = cdet%unocc_list_beta(a_ind)
                            a_found = .true.
                        else
                            a_found = .false.
                        end if
                    else
                        a_found = .false.
                    end if
                else
                    if (sys%nvirt_alpha > 0) then
                        allocate(ia_weights(1:sys%nvirt_alpha), stat=ierr)
                        call check_allocate('ia_weights', sys%nvirt_alpha, ierr)

                        call create_weighted_excitation_list_ptr(sys, i, 0, cdet%unocc_list_alpha, sys%nvirt_alpha, ia_weights, &
                                                             ia_weights_tot)
                        if (ia_weights_tot > 0.0_p) then
                            ! Use the alias method to select i with the appropriate probability
                            a_ind = select_weighted_value(rng, sys%nvirt_alpha, ia_weights, ia_weights_tot)
                            a = cdet%unocc_list_alpha(a_ind)
                            a_found = .true.
                        else
                            a_found = .false.
                        end if
                    else
                        a_found = .false.
                    end if
                end if
            end if

            if ((a_found) .and. (allowed_excitation)) then
                ! Given i,j,a construct the weights of all possible b
                ! This requires that total symmetry and spin are conserved.
                ! The symmetry of b (isymb) is given by
                ! (sym_i* x sym_j* x sym_a)* = sym_b
                ! (at least for Abelian point groups)
                isymb = sys%read_in%sym_conj_ptr(sys%read_in, &
                        sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, sys%basis%basis_fns(a)%sym))
                ! Ms_i + Ms_j = Ms_a + Ms_b => Ms_b = Ms_i + Ms_j - Ms_a
                ! Ms_k is +1 if up and -1 if down but imsb is +2 if up and +1 if down,
                ! therefore a conversion is necessary.
                imsb = (ij_spin - sys%basis%basis_fns(a)%Ms + 3)/2
            
                if (sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb) > 0) then
                    allocate(jb_weights(1:sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb)), stat=ierr)
                    call check_allocate('jb_weights', sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb), ierr)
                    call create_weighted_excitation_list_ptr(sys, j, a, sys%read_in%pg_sym%sym_spin_basis_fns(:,imsb,isymb), &
                                                             sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb), jb_weights, &
                                                             jb_weights_tot)
                else
                    jb_weights_tot = 0.0_p
                end if
            end if
            ! Test whether at least one possible b given i,j,a exists.
            ! Note that we did not need a btest for orbital a because we only considered
            ! virtual orbitals there.

            if (a_found .and. (jb_weights_tot > 0.0_p) .and. (allowed_excitation)) then
                ! 4. Use the alias method to select b with the appropriate probability
                b_ind = select_weighted_value(rng, sys%read_in%pg_sym%nbasis_sym_spin(imsb,isymb), jb_weights, jb_weights_tot)
                b = sys%read_in%pg_sym%sym_spin_basis_fns(b_ind,imsb,isymb)

                if (.not.btest(cdet%f(sys%basis%bit_lookup(2,b)), sys%basis%bit_lookup(1,b))) then
         
                    ! 5. Probability of generating this excitation.

                    ! Calculate p(ab|ij) = p(a|i) p(j|b) + p(b|i)p(a|j)
                    if (ij_spin==0) then
                        ! not possible to have chosen the reversed excitation
                        pgen=ia_weights(a_ind)/ia_weights_tot*jb_weights(b_ind)/jb_weights_tot
                    else
                        ! i and j have same spin, so could have been selected in the other order.
                        if (imsb == 1) then
                            ! find index b as if we had it selected first and as a from list of unoccupied virtual orbitals.
                            ! This uses the fact that the position of b in cdet%unocc_list_{beta,alpha} is the same
                            ! position index as in b's index in ia_weights is.
                            call binary_search(cdet%unocc_list_beta, b, 1, sys%nvirt_beta, found, b_ind_rev)
                        else
                            call binary_search(cdet%unocc_list_alpha, b, 1, sys%nvirt_alpha, found, b_ind_rev)
                        end if

                        isyma = sys%read_in%sym_conj_ptr(sys%read_in, &
                                sys%read_in%cross_product_sym_ptr(sys%read_in, ij_sym, isymb))
                        ! imsa = imsb
                        allocate(ja_weights(1:sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma)), stat=ierr)
                        call check_allocate('ja_weights', sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma), ierr)
                        call create_weighted_excitation_list_ptr(sys, j, b, sys%read_in%pg_sym%sym_spin_basis_fns(:,imsb,isyma),&
                                            sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma), ja_weights, ja_weights_tot)
                        call binary_search(sys%read_in%pg_sym%sym_spin_basis_fns(:,imsb,isyma), a, 1, &
                                    sys%read_in%pg_sym%nbasis_sym_spin(imsb,isyma), found, a_ind_rev)
                        pgen=       (ia_weights(a_ind)*jb_weights(b_ind)) / (ia_weights_tot*jb_weights_tot) + &
                            (ia_weights(b_ind_rev)*ja_weights(a_ind_rev) ) / (ia_weights_tot*ja_weights_tot)
                    end if
                    pgen = excit_gen_data%pattempt_double * pgen * pgen_ij ! pgen(ab)
                    connection%nexcit = 2

                    ! if create_weighted_excitation_list checks this then this check is redundant
                    ! allowed_excitation = a /= b

                    if (i<j) then
                        connection%from_orb(1) = i
                        connection%from_orb(2) = j
                    else
                        connection%from_orb(2) = i
                        connection%from_orb(1) = j
                    end if
                    if (a<b) then
                        connection%to_orb(1) = a
                        connection%to_orb(2) = b
                    else
                        connection%to_orb(2) = a
                        connection%to_orb(1) = b
                    end if
                    allowed_excitation = .true.
                else
                    allowed_excitation = .false.
                end if
            else
                allowed_excitation = .false.
            end if
           
            if (allowed_excitation) then

                ! Parity of permutation required to line up determinants.
                ! NOTE: connection%from_orb and connection%to_orb *must* be ordered.
                call find_excitation_permutation2(sys%basis%excit_mask, cdet%f, connection)

                ! Find the connecting matrix element.
                hmatel = slater_condon2_excit_ptr(sys, connection%from_orb(1), connection%from_orb(2), &
                                            connection%to_orb(1), connection%to_orb(2), connection%perm)
            else
                ! Carelessly selected ij with no possible excitations.  Such
                ! events are not worth the cost of renormalising the generation
                ! probabilities.
                ! Return a null excitation.
                hmatel%c = cmplx(0.0_p, 0.0_p, p)
                hmatel%r = 0.0_p
                pgen = 1.0_p
            end if

            ! deallocate weight arrays if allocated
            if (allocated(ia_weights)) then
                deallocate(ia_weights, stat=ierr_dealloc)
                call check_deallocate('ia_weights', ierr_dealloc)
            end if
            if (allocated(jb_weights)) then
                deallocate(jb_weights, stat=ierr_dealloc)
                call check_deallocate('jb_weights', ierr_dealloc)
            end if
            if (allocated(ja_weights)) then
                deallocate(ja_weights, stat=ierr_dealloc)
                call check_deallocate('ja_weights', ierr_dealloc)
            end if
        end if

    end subroutine gen_excit_mol_power_pitzer_occ

end module excit_gen_power_pitzer_mol
