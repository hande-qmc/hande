module replica_rdm

! Procedures relating to the accumulation of the 2-RDM
! from replica sampling with FCIQMC and CCMC.

use const

implicit none

contains

    pure subroutine update_rdm(sys, f1, f2, occ_list_1, pop1, pop2, prob, rdm)

        ! Add contribution from a pair of determinants to the 2-RDM

        ! In:
        !   sys: system being studied
        !   f1, f2: bit strings of the two determinants
        !   occ_list_1: a list of orbitals occupied in determinant f1
        !   pop1, pop2: populations of the determinants in *different* replicas
        !   prob: Probability of this spawning event, to weight stochastic contributions
        ! In/Out:
        !   rdm: The 2-RDM, with indices in physiacl notation

        use system, only: sys_t
        use excitations, only: excit_t, get_excitation
        use determinant_data, only: det_info_t
        use utils, only: tri_ind_distinct_reorder

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(:), f2(:)
        integer, intent(in) :: occ_list_1(:)
        real(p), intent(in) :: pop1, pop2, prob
        real(p), intent(inout) :: rdm(:,:)

        type(excit_t) :: excit
        real(p) :: matel
        integer :: i, j

        excit = get_excitation(sys%nel, sys%basis, f1, f2)
        ! Contribution to matrix element
        matel = pop1*pop2/prob
        if (excit%perm) matel = -matel

        select case(excit%nexcit)
        case(0)
            ! Diagonal elements.
            do i = 1, sys%nel
                do j = i+1, sys%nel
                    associate(ind => tri_ind_distinct_reorder(occ_list_1(i), occ_list_1(j)))
                        rdm(ind,ind) = rdm(ind,ind) + matel
                    end associate
                end do
            end do
        case(1)
            ! Single excitation contributes to one term for each orbital in common
            do i = 1, sys%nel
                associate(p => excit%to_orb(1), q => excit%from_orb(1), orb => occ_list_1(i))
                    if (orb == q) cycle
                    if (orb<p .and. orb<q) then
                        rdm(tri_ind_distinct_reorder(orb,p),tri_ind_distinct_reorder(orb,q)) = &
                            rdm(tri_ind_distinct_reorder(orb,p),tri_ind_distinct_reorder(orb,q)) + matel
                    else if (orb>p .and. orb>q) then
                        rdm(tri_ind_distinct_reorder(orb,p),tri_ind_distinct_reorder(q,orb)) = &
                            rdm(tri_ind_distinct_reorder(orb,p),tri_ind_distinct_reorder(q,orb)) + matel
                    else
                        rdm(tri_ind_distinct_reorder(p,orb),tri_ind_distinct_reorder(orb,q)) = &
                            rdm(tri_ind_distinct_reorder(p,orb),tri_ind_distinct_reorder(orb,q)) - matel
                    end if
                end associate
            end do
        case(2)
            ! Double excitation contributes to one term
            associate(p => excit%to_orb(1), q => excit%to_orb(2), r => excit%from_orb(1), s => excit%from_orb(2))
                rdm(tri_ind_distinct_reorder(p,q),tri_ind_distinct_reorder(r,s)) = &
                    rdm(tri_ind_distinct_reorder(p,q),tri_ind_distinct_reorder(r,s)) + matel
            end associate
        end select

    end subroutine update_rdm

    subroutine update_rdm_from_spawns(sys, psip_list, rdm_spawn, rdm)

        ! Use spawning events from rdm_spawn to sample rdm
        ! This is similar to annihilate_main_list

        ! In:
        !   sys: system being studied
        !   psip_list: populations of determinants
        ! In.Out:
        !   rdm_spawn: spawn_t containing events used to sample RDM
        !   rdm: sampled RDM

        use system, only: sys_t
        use qmc_data, only: particle_t
        use spawn_data, only: spawn_t, annihilate_wrapper_spawn_t
        use search, only: binary_search
        use determinants, only: decode_det

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(in) :: psip_list
        type(spawn_t), intent(inout) :: rdm_spawn
        real(p), intent(inout) :: rdm(:,:)

        logical :: hit
        integer, parameter :: thread_id = 0
        integer :: i, pos, istart, iend
        integer(i0) :: f_parent(sys%basis%bit_string_len), f_child(sys%basis%bit_string_len)
        integer :: occ_list(sys%nel)
        integer(int_p) :: nspawned(psip_list%nspaces)

        istart = 1
        iend = psip_list%nstates

        call annihilate_wrapper_spawn_t(rdm_spawn, .false., .false.)

        associate(bsl=>sys%basis%bit_string_len)
            do i = 1, rdm_spawn%head(thread_id, 0)
                f_parent = int(rdm_spawn%sdata(:bsl,i), i0)
                f_child = int(rdm_spawn%sdata(bsl+1:2*bsl,i), i0)
                call binary_search(psip_list%states, f_child, istart, iend, hit, pos, .false.)
                ! Can ignore spawning if child not present in psip_list (population 0)
                if (hit) then
                    nspawned = int(rdm_spawn%sdata(bsl+1:bsl+psip_list%nspaces,i), int_p)
                    call decode_det(sys%basis, f_child, occ_list)
                    ! [todo] - use both sets of spawnings to update the rdm.  It should be a trivial change...
                    ! [review] - JSS: do i = 1,2; pops(i,pos) and nspawned(3-i)?
                    call update_rdm(sys, f_child, f_parent, occ_list, &
                                    real(psip_list%pops(1,pos),p)/psip_list%pop_real_factor, &
                                    real(nspawned(2),p)/psip_list%pop_real_factor, 1.0_p, rdm)
                end if
                ! Next spawn can't be to an earlier determinant (but can be to the same one from a
                ! different parent)
                istart = pos
            end do
        end associate

        ! Diagonal elements
        do i = 1, psip_list%nstates
            call decode_det(sys%basis, psip_list%states(:,i), occ_list)
            call update_rdm(sys, psip_list%states(:,i), psip_list%states(:,i), occ_list, &
                            real(psip_list%pops(1,i),p)/psip_list%pop_real_factor, &
                            real(psip_list%pops(2,i),p)/psip_list%pop_real_factor, 1.0_p, rdm)
        end do


        rdm_spawn%head = rdm_spawn%head_start

    end subroutine update_rdm_from_spawns

    subroutine normalise_rdm(nel, rdm)

        ! Enforce the trace of the 2-RDM as N(N-1)/2

        ! In:
        !   nel: number of electrons
        ! In/Out:
        !   rdm: unnormailsed 2-RDM

        integer, intent(in) :: nel
        real(p), intent(inout) :: rdm(:,:)

        integer :: i
        real(p) :: trace

        trace = 0.0_p

        do i = 1, size(rdm,dim=1)
            trace = trace + rdm(i,i)
        end do

        rdm = rdm * nel * (nel-1) * 0.5 / trace

    end subroutine normalise_rdm

    subroutine check_hermiticity(rdm, comment_unit)

        ! Calculate the deviation from Hermiticity of the RDM and make it Hermitian

        ! In/Out:
        !   rdm: the 2-RDM.
        ! In:
        !   comment_unit: io unit to write any comments to.

        real(p), intent(inout) :: rdm(:,:)
        integer, intent(in) :: comment_unit

        integer :: i,j
        real(p) :: max_abs_error, mean_abs_error, err

        max_abs_error = 0.0_p
        mean_abs_error = 0.0_p

        do i = 1, size(rdm,dim=1)
            do j = i+1, size(rdm,dim=1)
                err = abs(rdm(i,j)-rdm(j,i))
                mean_abs_error = mean_abs_error + err
                max_abs_error = max(max_abs_error,err)
                rdm(i,j) = 0.5*(rdm(i,j)+rdm(j,i))
                rdm(j,i) = rdm(i,j)
            end do
        end do

        mean_abs_error = mean_abs_error/(size(rdm,dim=1)*(size(rdm,dim=1)-1)/2)
        write (comment_unit,'(1X,"#",1X,"Maximum deviation from Hermiticity is:",1X,es17.10)') max_abs_error
        write (comment_unit,'(1X,"#",1X,"Average deviation from Hermiticity is:",1X,es17.10)') mean_abs_error

    end subroutine check_hermiticity

    pure subroutine calc_rdm_energy(sys, reference, rdm, energy, trace)

        ! Calculate the energy from the rdm

        ! In:
        !   sys: system being studied
        !   rdm: the 2 particle reduced density matrix

        use proc_pointers, only: get_one_e_int_ptr, get_two_e_int_ptr
        use system, only: sys_t, read_in
        use qmc_data, only: reference_t
        use utils, only: orbs_from_index

        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: reference
        real(p), intent(in) :: rdm(:,:)
        real(p), intent(out) :: energy, trace

        integer :: i, j, nrows, a, b, c, d

        nrows = size(rdm, dim=1)

        ! Calculate the trace to correct the energy as rdm is not necessarily normalised
        trace = 0.0_p
        do i = 1, nrows
            trace = trace + rdm(i,i)
        end do

        energy = 0.0_p

        do i = 1, nrows
            ! Get orbital indices from combined index
            call orbs_from_index(i, a, b)
            do j = 1, nrows
                if (abs(rdm(i,j)) < depsilon) cycle
                call orbs_from_index(j, c, d)
                energy = energy + rdm(i,j)*get_two_e_int_ptr(sys,a,b,c,d)
                ! Contributions from 1 RDM
                if (a == c) energy = energy + rdm(i,j)*get_one_e_int_ptr(sys,b,d)/(sys%nel-1)
                if (b == d) energy = energy + rdm(i,j)*get_one_e_int_ptr(sys,a,c)/(sys%nel-1)
                if (a == d) energy = energy - rdm(i,j)*get_one_e_int_ptr(sys,b,c)/(sys%nel-1)
                if (b == c) energy = energy - rdm(i,j)*get_one_e_int_ptr(sys,a,d)/(sys%nel-1)
            end do
        end do

        energy = energy*sys%nel*(sys%nel-1)*0.5

        if (sys%system == read_in) then
            energy = energy + (sys%read_in%Ecore - reference%H00)*trace
        else
            energy = energy - reference%H00*trace
        end if

    end subroutine calc_rdm_energy

    subroutine write_final_rdm(rdm, nel, nbasis, filename, comment_unit)

        ! Write (normalised, Hermitian part of) RDM to a file

        ! In:
        !   nel: number of electrons
        !   nbasis: Number of basis fns
        !   filename: file to write to
        !   comment_unit: io unit to write any additional comments to.
        ! In/Out:
        !   rdm: Sampled two particle reduced density matrix. On exit it has been
        !        normalised and made Hermitian and the root processor holds the sum
        !        over all processors

        use parallel
        use utils, only: tri_ind_distinct_reorder

        real(p), intent(inout) :: rdm(:,:)
        integer, intent(in) :: nel, nbasis, comment_unit
        character(*), intent(in) :: filename

        integer :: i, j, k, l, fileunit
        

#ifdef PARALLEL
        integer :: ierr

        real(p), allocatable :: rdm_total(:,:)

        if (parent) allocate(rdm_total(size(rdm, dim=1),size(rdm,dim=2)))
        call mpi_reduce(rdm, rdm_total, size(rdm), MPI_preal, MPI_SUM, root, MPI_COMM_WORLD, ierr)
        if (parent) then
            rdm = rdm_total
            deallocate(rdm_total)
        end if
#endif

        if (parent) then
            call normalise_rdm(nel, rdm)
            call check_hermiticity(rdm, comment_unit)

            open(file=filename, newunit=fileunit)
            do i = 1, nbasis
                do j = i+1, nbasis
                    do k = 1, nbasis
                        do l = k+1, nbasis
                            if (abs(rdm(tri_ind_distinct_reorder(i,j),tri_ind_distinct_reorder(k,l))) > depsilon) then
                                write (fileunit, '(4i4,2x,es13.6)') i, j, k, l, &
                                    rdm(tri_ind_distinct_reorder(i,j),tri_ind_distinct_reorder(k,l))
                            end if
                        end do
                    end do
                end do
            end do
            close(fileunit)
        end if

    end subroutine write_final_rdm

end module replica_rdm
