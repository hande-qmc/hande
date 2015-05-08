module replica_rdm

! Procedures relating to the accumulation of the 1- and 2-RDM
! from replica sampling with FCIQMC.

use const

implicit none

contains

    pure function rdm_ind(p, q)

        ! We require p<q and r<s for \Gamma_pqrs so can use a rank two array
        ! for the 2-RDM to save a factor of 4 in storage.

        ! In:
        !   p, q: orbital indices
        ! Returns:
        !   Composite index for the RDM

        integer, intent(in) :: p, q
        integer :: rdm_ind

        ! The same indexing is used as in choose_ij_k

        rdm_ind = (q-1)*(q-2)/2 + p

    end function rdm_ind

    pure subroutine orbs_from_index(ind, p, q)

        ! Does the reverse of rdm_ind - gets a pair of orbitals from a single index

        ! In:
        !    ind: index to the RDM
        ! Out:
        !    p, q: orbitals, with p<q

        integer, intent(in) :: ind
        integer, intent(out) :: p, q

        q = int(1.50 + sqrt(2*ind-1.750))
        p = ind - ((q-1)*(q-2))/2

    end subroutine orbs_from_index

    pure subroutine update_rdm(sys, f1, f2, pop1, pop2, prob, rdm)

        ! Add contribution from a pair of determinants to the 2-RDM

        ! In:
        !   sys: system being studied
        !   f1, f2: bit strings of the two determinants
        !   pop1, pop2: populations of the determinants in *different* replicas
        !   prob: Probability of this spawning event, to weight stochastic contributions
        ! In/Out:
        !   rdm: The 2-RDM, with indices in physiacl notation

        use system, only: sys_t
        use excitations, only: excit_t, get_excitation
        use determinants, only: decode_det

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(:), f2(:)
        real(p), intent(in) :: pop1, pop2, prob
        real(p), intent(inout) :: rdm(:,:)

        type(excit_t) :: excit
        real(p) :: matel
        integer :: occ_list(sys%nel)
        integer :: i, j

        excit = get_excitation(sys%nel, sys%basis, f1, f2)
        ! Contribution to matrix element
        matel = pop1*pop2/prob
        if (excit%perm) matel = -matel

        select case(excit%nexcit)
        case(0)
            ! Diagonal elements.
            call decode_det(sys%basis, f1, occ_list)
            do i = 1, sys%nel
                do j = i+1, sys%nel
                    associate(ind => rdm_ind(occ_list(i),occ_list(j)))
                        rdm(ind,ind) = rdm(ind,ind) + matel
                    end associate
                end do
            end do
        case(1)
            ! Single excitation contributes to one term for each orbital in common
            call decode_det(sys%basis, iand(f1,f2), occ_list)
            do i = 1, sys%nel-1
                associate(p => excit%to_orb(1), q => excit%from_orb(1), orb => occ_list(i))
                    if (orb<p .and. orb<q) then
                        rdm(rdm_ind(orb,p),rdm_ind(orb,q)) = rdm(rdm_ind(orb,p),rdm_ind(orb,q)) + matel
                    else if (orb<p .and. orb>q) then
                        rdm(rdm_ind(orb,p),rdm_ind(q,orb)) = rdm(rdm_ind(orb,p),rdm_ind(q,orb)) - matel
                    else if (orb>p .and. orb<q) then
                        rdm(rdm_ind(p,orb),rdm_ind(orb,q)) = rdm(rdm_ind(p,orb),rdm_ind(orb,q)) - matel
                    else
                        rdm(rdm_ind(p,orb),rdm_ind(q,orb)) = rdm(rdm_ind(p,orb),rdm_ind(q,orb)) + matel
                    end if
                end associate
            end do
        case(2)
            ! Double excitation contributes to one term
            associate(p => excit%to_orb(1), q => excit%to_orb(2), r => excit%from_orb(1), s => excit%from_orb(2))
                rdm(rdm_ind(p,q),rdm_ind(r,s)) = rdm(rdm_ind(p,q),rdm_ind(r,s)) + matel
            end associate
        end select

    end subroutine update_rdm

    subroutine normalise_rdm(nel, rdm)

        ! Enforce the trace of the 2-RDM as N(N-1)/2

        ! In:
        !   nel: number of electrons
        ! In/Out:
        !   rdm: unnormailsed 2-RDM

        integer, intent(in) :: nel
        real(p), intent(inout) :: rdm(:,:)

        integer :: nbasis, i, j
        real(p) :: trace

        trace = 0.0_p

        do i = 1, size(rdm,dim=1)
            trace = trace + rdm(i,i)
        end do

        rdm = rdm * nel * (nel-1) * 0.5 / trace

    end subroutine normalise_rdm

    subroutine check_hermiticity(rdm)

        ! Calculate the deviation from Hermiticity of the RDM and make it Hermitian

        ! In/Out:
        !   rdm: the 2-RDM.

        real(p), intent(inout) :: rdm(:,:)

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
        write (6,'(1X,"#",1X,"Maximum deviation from Hermiticity is:",1X,es17.10)') max_abs_error
        write (6,'(1X,"#",1X,"Average deviation from Hermiticity is:",1X,es17.10)') mean_abs_error

    end subroutine check_hermiticity

    pure subroutine calc_rdm_energy(sys, reference, rdm, energy, trace)

        ! Calculate the energy from the rdm

        ! In:
        !   sys: system being studied
        !   rdm: the 2 particle reduced density matrix

        use proc_pointers, only: get_one_e_int_ptr, get_two_e_int_ptr
        use system, only: sys_t, read_in
        use qmc_data, only: reference_t

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

    subroutine write_final_rdm(rdm, nel, nbasis, filename)

        ! Write (normalised, Hermitian part of) RDM to a file

        ! In:
        !   nel: number of electrons
        !   nbasis: Number of basis fns
        !   filename: file to write to
        ! In/Out:
        !   rdm: Sampled two particle reduced density matrix. On exit it has been
        !        normalised and made Hermitian

        use parallel

        use parallel

        real(p), intent(inout) :: rdm(:,:)
        integer, intent(in) :: nel, nbasis
        character(*), intent(in) :: filename

        integer :: i, j, k, l, fileunit, ierr

        real(p), allocatable :: rdm_total(:,:)

#ifdef PARALLEL
        if (parent) allocate(rdm_total(size(rdm, dim=1),size(rdm,dim=2)))
        call mpi_reduce(rdm, rdm_total, size(rdm), MPI_preal, MPI_SUM, root, MPI_COMM_WORLD, ierr)
        if (parent) then
            rdm = rdm_total
            deallocate(rdm_total)
        end if
#endif

        if (parent) then
            call normalise_rdm(nel, rdm)
            call check_hermiticity(rdm)

            open(file=filename, newunit=fileunit)
            do i = 1, nbasis
                do j = i+1, nbasis
                    do k = 1, nbasis
                        do l = k+1, nbasis
                            if (abs(rdm(rdm_ind(i,j),rdm_ind(k,l))) > depsilon) &
                                write (fileunit, '(4i4,2x,es13.6)') i, j, k, l, rdm(rdm_ind(i,j),rdm_ind(k,l))
                        end do
                    end do
                end do
            end do
            close(fileunit)
        end if

    end subroutine write_final_rdm

end module replica_rdm
