module diagonalisation

! Module for coordinating complete and Lanczos diagonalisations.

use const
use parallel, only: iproc, nprocs, parent
use calc

implicit none

contains

    subroutine diagonalise()

        ! Construct and diagonalise the Hamiltonian matrix using Lanczos and/or
        ! exact diagonalisation.

        use checking, only: check_allocate, check_deallocate
        use basis, only: nbasis
        use errors
        use system
        use determinant_enumeration, only: enumerate_determinants, ndets, sym_space_size
        use determinants, only: tot_ndets, set_spin_polarisation, spin_orb_list
        use fciqmc_data, only: occ_list0
        use dmqmc_procedures, only: setup_rdm_arrays
        use lanczos
        use fciqmc_data, only: doing_exact_rdm_eigv, reduced_density_matrix
        use full_diagonalisation
        use hamiltonian, only: get_hmatel_dets
        use reference_determinant
        use symmetry, only: symmetry_orb_list

        use errors, only: stop_all
        use utils, only: int_fmt
        use ranking, only: insertion_rank_rp

        type soln
            integer :: ms
            real(p) :: energy
        end type soln

        integer :: ms, ms_min, ms_max, ierr, nlanczos, nexact, nfound, i, ind, state, isym, isym_min, isym_max

        type(soln), allocatable :: lanczos_solns(:), exact_solns(:)

        ! For communication with Lanczos and exact diagonalisers.
        real(dp), allocatable :: lanczos_eigv(:)
        real(p), allocatable :: exact_eigv(:)

        ! For sorting the solutions by energy rather than by energy within each
        ! spin block.
        integer, allocatable :: eigv_rank(:)

        logical :: spin_flip

        ! This stores the reduced density matrix eigenvalues, when being calculated.
        real(p), allocatable :: rdm_eigenvalues(:)
        integer :: rdm_size

        character(50) :: fmt1

        if (parent) write (6,'(1X,a15,/,1X,15("="),/)') 'Diagonalisation'

        ! The Hamiltonian can be written in block diagonal form using the spin
        ! quantum number.  < D_1 | H | D_2 > /= 0 only if D_1 and D_2 have the
        ! same total spin.

        ! For a given number of electrons, n, the total spin can range from -n
        ! to +n in steps of 2.
        if (ms_in == huge(1)) then
            if (sys_global%read_in%uhf) then
                ms_min = -min(sys_global%nel, nbasis-sys_global%nel)
                ms_max = min(sys_global%nel, nbasis-sys_global%nel)
            else
                ! The -ms blocks are degenerate with the +ms blocks so only need to
                ! solve for ms >= 0.
                ms_min = mod(sys_global%nel,2)
                ms_max = min(sys_global%nel, nbasis-sys_global%nel)
            end if
        else
            ! ms was set in input
            ms_min = ms_in
            ms_max = ms_in
        end if

        if (sym_in == huge(1)) then
            isym_min = sys_global%sym0
            isym_max = sys_global%sym_max
        else
            ! sym was set in input
            isym_min = sym_in
            isym_max = sym_in
        end if

        spin_flip = .false.
        if (allocated(occ_list0)) then
            ! If a reference determinant was supplied, then we need to only
            ! consider that spin and symmetry.
            if (isym_min /= isym_max) then
                isym_min = symmetry_orb_list(occ_list0)
                isym_max = isym_min
            else
                spin_flip = .true.
            end if
            if (ms_min /= ms_max) then
                ms_min = spin_orb_list(occ_list0)
                ms_max = ms_min
            else
                spin_flip = .true.
            end if
        endif

        if (truncate_space) then
            if (.not. allocated(occ_list0)) then
                if (isym_min /= isym_max .or. ms_min /= ms_max) then
                    call stop_all('diagonalise', 'Symmetry and spin must be &
                        &specified or a reference determinant supplied for &
                        &truncated CI calculations.')
                end if
                ! Create reference det based upon symmetry labels.
                call set_spin_polarisation(ms_min)
                allocate(occ_list0(sys_global%nel), stat=ierr)
                call check_allocate('occ_list0', sys_global%nel, ierr)
                call set_reference_det(occ_list0, .true., isym_min)
            end if
        end if

        if (doing_calc(lanczos_diag)) then
            allocate(lanczos_solns(nlanczos_eigv*(sys_global%nel/2+1)*nbasis/2), stat=ierr)
            call check_allocate('lanczos_solns',nlanczos_eigv*(sys_global%nel/2+1)*nbasis/2,ierr)
        end if
        if (doing_calc(exact_diag)) then
            allocate(exact_solns(tot_ndets), stat=ierr)
            call check_allocate('exact_solns',tot_ndets,ierr)
        end if

        ! Number of lanczos and exact solutions found.
        nlanczos = 0
        nexact = 0

        do ms = ms_min, ms_max, 2

            if (parent) write (6,'(1X,a35,'//int_fmt(ms)//',/)') 'Considering determinants with spin:', ms

            ! Find and set information about the space.
            call set_spin_polarisation(ms)
            if (allocated(occ_list0)) then
                call enumerate_determinants(.true., spin_flip, occ_list0=occ_list0)
            else
                call enumerate_determinants(.true., spin_flip)
            end if

            ! Diagonalise each symmetry block in turn.
            do isym = isym_min, isym_max

                if (sym_space_size(isym)==0) then
                    if (parent) then
                        associate(nsym=>sys_global%nsym)
                            fmt1 = '(1X,a25,'//int_fmt(isym,1)//',1X,a17,'//int_fmt(nsym,1)//',1X,a6,'//int_fmt(nsym,1)//')'
                            write (6,fmt1) 'No determinants with spin',ms,'in symmetry block',isym,'out of',nsym
                            write (6,'(/,1X,15("-"),/)')
                        end associate
                    end if
                    cycle
                end if

                if (parent) then
                    fmt1 = '(1X,a28,'//int_fmt(isym,1)//',1X,a6,'//int_fmt(sys_global%nsym,1)//')'
                    write (6,fmt1) 'Diagonalising symmetry block',isym,'out of',sys_global%nsym
                end if

                ! Find all determinants with this spin.
                if (allocated(occ_list0)) then
                    call enumerate_determinants(.false., spin_flip, isym, occ_list0)
                else
                    call enumerate_determinants(.false., spin_flip, isym)
                end if

                if (ndets == 1) then

                    if (parent) then
                        write (6,'(/,1X,a35,/)') 'Performing trivial diagonalisation.'
                    end if

                    ! The trivial case seems to trip up TRLan and scalapack in
                    ! parallel.
                    if (doing_calc(lanczos_diag)) then
                        lanczos_solns(nlanczos+1)%energy = get_hmatel_dets(1,1)
                        lanczos_solns(nlanczos+1)%ms = ms
                        nlanczos = nlanczos + 1
                    end if
                    if (doing_calc(exact_diag)) then
                        exact_solns(nexact+1)%energy = get_hmatel_dets(1,1)
                        exact_solns(nexact+1)%ms = ms
                        nexact = nexact + 1
                    end if

                else

                    if (nprocs == 1 .and. (doing_calc(exact_diag) .or. .not.direct_lanczos)) call generate_hamil(distribute_off)

                    ! Lanczos.
                    if (doing_calc(lanczos_diag)) then
                        allocate(lanczos_eigv(ndets), stat=ierr)
                        call check_allocate('lanczos_eigv',ndets,ierr)
                        ! Construct the Hamiltonian matrix distributed over the processors
                        ! if running in parallel.
                        if (nprocs > 1 .and. .not.direct_lanczos) call generate_hamil(distribute_cols)
                        call lanczos_diagonalisation(nfound, lanczos_eigv)
                        lanczos_solns(nlanczos+1:nlanczos+nfound)%energy = lanczos_eigv(:nfound)
                        lanczos_solns(nlanczos+1:nlanczos+nfound)%ms = ms
                        nlanczos = nlanczos + nfound
                        deallocate(lanczos_eigv)
                        call check_deallocate('lanczos_eigv',ierr)
                    end if

                    ! Exact diagonalisation.
                    ! Warning: this destroys the Hamiltonian matrix...
                    if (doing_calc(exact_diag)) then
                        allocate(exact_eigv(ndets), stat=ierr)
                        call check_allocate('exact_eigv',ndets,ierr)
                        ! Construct the Hamiltonian matrix distributed over the processors
                        ! if running in parallel.
                        if (nprocs > 1) call generate_hamil(distribute_blocks)
                        call exact_diagonalisation(exact_eigv)
                        exact_solns(nexact+1:nexact+ndets)%energy = exact_eigv
                        exact_solns(nexact+1:nexact+ndets)%ms = ms
                        nexact = nexact + ndets
                        deallocate(exact_eigv)
                        call check_deallocate('exact_eigv',ierr)
                    end if

                end if

                ! Reduced density matrix eigenvalues calculation.
                if (doing_exact_rdm_eigv) then
                    if (nprocs > 1) then
                        call warning('diagonalise','RDM eigenvalue calculation is only implemented in serial. &
                                      &Skipping this calculation.', 3)
                    else
                        if (doing_calc(lanczos_diag)) then
                            call warning('diagonalise','RDM eigenvalues cannot be calculated with the Lanczos &
                                                        &routine. Skipping this calculation.', 3)
                        else if (doing_calc(exact_diag)) then
                            write(6,'(1x,a46)') "Performing reduced density matrix calculation."
                            call setup_rdm_arrays()
                            rdm_size = size(reduced_density_matrix, 1)

                            allocate(rdm_eigenvalues(rdm_size), stat=ierr)
                            call check_allocate('rdm_eigenvalues',rdm_size,ierr)

                            call get_rdm_eigenvalues(rdm_eigenvalues)
                        end if
                    end if
                end if

                if (parent) write (6,'(1X,15("-"),/)')

            end do

        end do

        ! Output results.
        if (doing_calc(lanczos_diag) .and. parent) then
            write (6,'(1X,a31,/,1X,31("-"),/)') 'Lanczos diagonalisation results'
            allocate(eigv_rank(nlanczos), stat=ierr)
            call check_allocate('eigv_rank',nlanczos,ierr)
            call insertion_rank_rp(lanczos_solns(:nlanczos)%energy, eigv_rank, tolerance=depsilon)
            write (6,'(1X,a8,3X,a4,3X,a12)') 'State','Spin','Total energy'
            state = 0
            do i = 1, nlanczos
                ind = eigv_rank(i)
                state = state + 1
                write (6,'(1X,i6,2X,i6,1X,f18.12)') state, lanczos_solns(ind)%ms, lanczos_solns(ind)%energy
                if (.not.sys_global%read_in%uhf .and. lanczos_solns(ind)%ms /= 0) then
                    ! Also a degenerate -ms solution.
                    state = state + 1
                    write (6,'(1X,i6,2X,i6,1X,f18.12)') state, -lanczos_solns(ind)%ms, lanczos_solns(ind)%energy
                end if
            end do
            write (6,'(/,1X,a21,f18.12,/)') 'Lanczos ground state:', lanczos_solns(eigv_rank(1))%energy
            deallocate(eigv_rank, stat=ierr)
            call check_deallocate('eigv_rank',ierr)
        end if

        if (doing_calc(exact_diag) .and. parent) then
            write (6,'(1X,a29,/,1X,29("-"),/)') 'Exact diagonalisation results'
            allocate(eigv_rank(nexact), stat=ierr)
            call check_allocate('eigv_rank',nexact,ierr)
            call insertion_rank_rp(exact_solns(:nexact)%energy, eigv_rank, tolerance=depsilon)
            write (6,'(1X,a8,3X,a4,3X,a12)') 'State','Spin','Total energy'
            state = 0
            do i = 1, nexact
                ind = eigv_rank(i)
                state = state + 1
                write (6,'(1X,i6,2X,i6,1X,f18.12)') state, exact_solns(ind)%ms, exact_solns(ind)%energy
                if (.not.sys_global%read_in%uhf .and. exact_solns(ind)%ms /= 0) then
                    ! Also a degenerate -ms solution.
                    state = state + 1
                    write (6,'(1X,i6,2X,i6,1X,f18.12)') state, -exact_solns(ind)%ms, exact_solns(ind)%energy
                end if
            end do
            write (6,'(/,1X,a19,f18.12,/)') 'Exact ground state:', exact_solns(eigv_rank(1))%energy
            deallocate(eigv_rank, stat=ierr)
            call check_deallocate('eigv_rank',ierr)
            if (doing_exact_rdm_eigv) then
                write(6, '(1x,a35,/,1X,35("-"))') 'Reduced density matrix eigenvalues:'
                do i = 1, size(reduced_density_matrix, 1)
                    write(6, '(1x,i6,1x,es15.8)') i, rdm_eigenvalues(i)
                end do
                write(6,'()')
                deallocate(reduced_density_matrix, stat=ierr)
                call check_deallocate('reduced_density_matrix',ierr)
                deallocate(rdm_eigenvalues, stat=ierr)
                call check_deallocate('rdm_eigenvalues',ierr)
            end if
        end if

        if (doing_calc(lanczos_diag)) then
            deallocate(lanczos_solns, stat=ierr)
            call check_deallocate('lanczos_solns',ierr)
        end if
        if (doing_calc(exact_diag)) then
            deallocate(exact_solns, stat=ierr)
            call check_deallocate('exact_solns',ierr)
        end if

    end subroutine diagonalise

    subroutine generate_hamil(distribute_mode)

        ! Generate a symmetry block of the Hamiltonian matrix, H = < D_i | H | D_j >.
        ! The list of determinants, {D_i}, is grouped by symmetry and contains
        ! only determinants of a specified spin.
        ! Only generate the upper diagonal for use with (sca)lapack and Lanczos routines.
        ! In:
        !    distribute_mode (optional): flag for determining how the
        !        Hamiltonian matrix is distributed among processors.  It is
        !        irrelevant if only one processor is used: the distribution schemes
        !        all reduce to the storing the entire matrix on the single
        !        processor.  Can take the values given by the distribute_off,
        !        distribute_blocks and distribute_cols parameters.  See above
        !        for descriptions of the different behaviours.

        use checking, only: check_allocate, check_deallocate
        use csr, only: init_csrpsy, end_csrpsy
        use utils, only: get_free_unit
        use errors
        use parallel

        use hamiltonian, only: get_hmatel_dets
        use hubbard_real
        use determinant_enumeration, only: ndets

        integer, intent(in), optional :: distribute_mode
        integer :: ierr, iunit, n1, n2, ind_offset
        integer :: i, j, ii, jj, ilocal, iglobal, jlocal, jglobal, nnz, imode
        real(p) :: hmatel

        ! Store how many determinants have already been considered.
        ! This allows for the correct indexing scheme to be used in
        ! output of the Hamiltonian matrix.
        integer, save :: ndets_prev = 0

        ! scratch array for printing distributed matrix
        real(p), allocatable :: work_print(:)

        if (allocated(hamil)) then
            deallocate(hamil, stat=ierr)
            call check_deallocate('hamil',ierr)
        end if
        if (allocated(hamil_csr%mat)) call end_csrpsy(hamil_csr)

        if (present(distribute_mode)) then
            distribute = distribute_mode
        else
            distribute = distribute_off
        end if

        ! Find dimensions of local array.
        select case(distribute)
        case(distribute_off)
            ! Useful to have a dummy proc_blacs_info to refer to.  Everything
            ! apart from the matrix descriptors are valid in the dummy instance.
            proc_blacs_info = get_blacs_info(ndets, (/1, nprocs/))
            n1 = ndets
            n2 = ndets
        case(distribute_blocks)
            ! Use as square a processor grid as possible.
            proc_blacs_info = get_blacs_info(ndets)
            n1 = proc_blacs_info%nrows
            n2 = proc_blacs_info%ncols
        case(distribute_cols)
            ! TRLan assumes that the only the rows of the matrix
            ! are distributed.
            ! Furthermore it seems TRLan assumes all processors store at least
            ! some of the matrix.
            if (nprocs*block_size > ndets) then
                if (parent) then
                    call warning('generate_hamil','Reducing block size so that all processors contain at least a single row.', 3)
                    write (6,'(1X,a69)') 'Consider running on fewer processors or reducing block size in input.'
                    write (6,'(1X,a19,i4)') 'Old block size was:',block_size
                end if
                block_size = ndets/nprocs
                if (parent) write (6,'(1X,a18,i4,/)') 'New block size is:',block_size
            end if
            proc_blacs_info = get_blacs_info(ndets, (/1, nprocs/))
            n1 = proc_blacs_info%nrows
            n2 = proc_blacs_info%ncols
        case default
            call stop_all('generate_hamil','Unknown distribution scheme.')
        end select

        if (use_sparse_hamil) then
            if (distribute /= distribute_off) then
                write (6,'(1X,a58,/,1X,a26,/,1X,a83)') &
                    'Sparse distributed matrices are not currently implemented.', &
                    'Not using sparse matrices.', &
                    'If this is disagreeable to you, please contribute patches resolving this situation.'
            end if
            if (doing_calc(exact_diag)) then
                write (6,'(1X,a47,1X,a26,/,1X,a83)') &
                    'Sparse matrices are not compatible with LAPACK.', &
                    'Not using sparse matrices.', &
                    'If this is disagreeable to you, please contribute patches resolving this situation.'
            end if
        end if

        if (distribute /= distribute_off .or. doing_calc(exact_diag) .or. .not. use_sparse_hamil) then
            allocate(hamil(n1,n2), stat=ierr)
            call check_allocate('hamil',n1*n2,ierr)
        end if

        ! index offset for the symmetry block compared to the index of the
        ! determinant list.
        !ind_offset = ???
        ind_offset = 0

        ! Form the Hamiltonian matrix < D_i | H | D_j >.
        select case(distribute)
        case(distribute_off)
            if (use_sparse_hamil .and. .not. doing_calc(exact_diag)) then
                ! First, find out how many non-zero elements there are.
                ! We'll be naive and not just test for non-zero by symmetry, but
                ! actually calculate all matrix elements for now.
                ! Then, store the Hamiltonian.
                do imode = 1, 2
                    nnz = 0
                    !$omp parallel
                    do i = 1, ndets
                        ! OpenMP chunk size determined completely empirically
                        ! from a single test.  Please feel free to improve...
                        !$omp do private(j, hmatel) schedule(dynamic, 200)
                        do j = i, ndets
                            hmatel = get_hmatel_dets(i+ind_offset, j+ind_offset)
                            if (abs(hmatel) > depsilon) then
                                !$omp critical
                                nnz = nnz + 1
                                if (imode == 2) then
                                    hamil_csr%mat(nnz) = hmatel
                                    hamil_csr%col_ind(nnz) = j
                                    if (hamil_csr%row_ptr(i) == 0) hamil_csr%row_ptr(i) = nnz
                                end if
                                !$omp end critical
                            end if
                        end do
                        !$omp end do
                    end do
                    !$omp end parallel
                    if (imode == 1) then
                        write (6,'(1X,a50,i8/)') 'Number of non-zero elements in Hamiltonian matrix:', nnz
                        call init_csrpsy(hamil_csr, ndets, nnz)
                        hamil_csr%row_ptr(1:ndets) = 0
                    else
                        ! Any element not set in row_ptr means that the
                        ! corresponding row has no non-zero elements.
                        ! Set it to be identical to the next row (this avoids
                        ! looping over the zero-row).
                        do i = ndets, 1, -1
                            if (hamil_csr%row_ptr(i) == 0) hamil_csr%row_ptr(i) = hamil_csr%row_ptr(i+1)
                        end do
                    end if
                end do
            else
                !$omp parallel
                do i = 1, ndets
                    !$omp do private(j) schedule(dynamic, 200)
                    do j = i, ndets
                        hamil(i,j) = get_hmatel_dets(i+ind_offset,j+ind_offset)
                    end do
                    !$omp end do
                end do
                !$omp end parallel
            end if
        case(distribute_blocks, distribute_cols)
            ! blacs divides the matrix up into sub-matrices of size block_size x block_size.
            ! The blocks are distributed in a cyclic fashion amongst the
            ! processors.
            ! Each processor stores a total of nrows of the matrix.
            ! i gives the index of the first row in the current block.
            ! ii gives the index of the current row within the current block.
            ! The local i index is thus i-1+ii.  This is used to refer to the
            ! matrix element as stored (continuously) on the processor.
            ! The global i index is given by the sum of the rows held on
            ! preceeding proecssors for the previous loops over processors (as
            ! done in the loop over i), the rows held on other processors during
            ! the corrent loop over processors and the position within the
            ! current block.
            ! Similarly for the other indices.
            do i = 1, proc_blacs_info%nrows, block_size
                do ii = 1, min(block_size, proc_blacs_info%nrows - i + 1)
                    ilocal = i - 1 + ii
                    iglobal = (i-1)*nproc_rows + proc_blacs_info%procx*block_size + ii + ind_offset
                    do j = 1, proc_blacs_info%ncols, block_size
                        do jj = 1, min(block_size, proc_blacs_info%ncols - j + 1)
                            jlocal = j - 1 + jj
                            jglobal = (j-1)*nproc_cols + proc_blacs_info%procy*block_size + jj + ind_offset
                            hamil(ilocal, jlocal) = get_hmatel_dets(iglobal, jglobal)
                        end do
                    end do
                end do
            end do
        end select

        if (write_hamiltonian) then
            iunit = get_free_unit()
            open(iunit, file=hamiltonian_file, status='unknown')
            if (nprocs > 1) then
                ! Note that this uses a different format to the serial case...
                allocate(work_print(block_size**2), stat=ierr)
                call check_allocate('work_print', block_size**2, ierr)
#ifdef SINGLE_PRECISION
                call pslaprnt(ndets, ndets, hamil, 1, 1, proc_blacs_info%desc_m, 0, 0, 'hamil', iunit, work_print)
#else
                call pdlaprnt(ndets, ndets, hamil, 1, 1, proc_blacs_info%desc_m, 0, 0, 'hamil', iunit, work_print)
#endif
                deallocate(work_print)
                call check_deallocate('work_print', ierr)
            else
                if (allocated(hamil)) then
                    do i=1, ndets
                        write (iunit,*) i,i,hamil(i,i)
                        do j=i+1, ndets
                            if (abs(hamil(i,j)) > depsilon) write (iunit,*) i+ndets_prev,j+ndets_prev,hamil(i,j)
                        end do
                    end do
                else if (allocated(hamil_csr%mat)) then
                    j = 1
                    do i = 1, size(hamil_csr%mat)
                        if (abs(hamil_csr%mat(i)) > depsilon) then
                            if (hamil_csr%row_ptr(j+1) <= i) j = j+1
                            write (iunit,*) j+ndets_prev, hamil_csr%col_ind(i)+ndets_prev, hamil_csr%mat(i)
                        end if
                    end do
                end if
                ndets_prev = ndets
            end if
            close(iunit, status='keep')
        end if

    end subroutine generate_hamil

    subroutine get_rdm_eigenvalues(rdm_eigenvalues)

        use basis, only: basis_length
        use checking, only: check_allocate, check_deallocate
        use determinant_enumeration, only: ndets, dets_list
        use dmqmc_procedures, only: decode_dm_bitstring, rdms
        use fciqmc_data, only: reduced_density_matrix

        real(p), intent(out) :: rdm_eigenvalues(size(reduced_density_matrix,1))
        integer(i0) :: f1(basis_length), f2(basis_length)
        integer(i0) :: f3(2*basis_length)
        integer :: i, j, rdm_size, info, ierr, lwork
        integer(i0) :: end1, end2
        real(p), allocatable :: work(:)
        real(p) :: rdm_element

        write(6,'(1x,a36)') "Setting up reduced density matrix..."

        ! Loop over all elements of the density matrix and add all contributing elements to the RDM.
        do i = 1, ndets
            do j = 1, ndets
                f1 = iand(rdms(1)%B_masks(:,1),dets_list(:,i))
                f2 = iand(rdms(1)%B_masks(:,1),dets_list(:,j))
                ! If the two bitstrings are the same after bits corresponding to subsystem B have
                ! been unset, then these two bitstrings contribute to the RDM.
                if (sum(abs(f1-f2)) == 0) then
                    ! In f3, concatenate the two bitstrings.
                    f3(1:basis_length) = dets_list(:,i)
                    f3(basis_length+1:basis_length*2) = dets_list(:,j)

                    ! Get the position in the RDM of this density matrix element.
                    call decode_dm_bitstring(f3,1,1)
                    rdms(1)%end1 = rdms(1)%end1 + 1
                    rdms(1)%end2 = rdms(1)%end2 + 1

                    ! The ground state wave function is stored in hamil(:,1).
                    rdm_element = hamil(i,1)*hamil(j,1)
                    ! Finally add in the contribution from this density matrix element.
                    reduced_density_matrix(rdms(1)%end1(1),rdms(1)%end2(1)) = &
                        reduced_density_matrix(rdms(1)%end1(1),rdms(1)%end2(1)) + rdm_element
                end if
            end do
        end do

        ! Now the RDM is completley calculated.

        ! Calculate the eigenvalues:
        write(6,'(1x,a39,/)') "Diagonalising reduced density matrix..."

        rdm_size = size(reduced_density_matrix, 1)
        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        call check_allocate('work',1,ierr)
#ifdef SINGLE_PRECISION
        call ssyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, rdm_eigenvalues, work, -1, info)
#else
        call dsyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, rdm_eigenvalues, work, -1, info)
#endif
        lwork = nint(work(1))
        deallocate(work)
        call check_deallocate('work',ierr)

        ! Perform the diagonalisation.
        allocate(work(lwork), stat=ierr)
        call check_allocate('work',lwork,ierr)

#ifdef SINGLE_PRECISION
        call ssyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, rdm_eigenvalues, work, lwork, info)
#else
        call dsyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, rdm_eigenvalues, work, lwork, info)
#endif

    end subroutine get_rdm_eigenvalues

    subroutine end_hamil()

        ! Clean up hamiltonian module.

        use checking, only: check_deallocate

        integer :: ierr

        if (allocated(hamil)) then
            deallocate(hamil, stat=ierr)
            call check_deallocate('hamil',ierr)
        end if

    end subroutine end_hamil

end module diagonalisation
