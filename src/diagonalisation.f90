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

        use basis, only: nbasis
        use system, only: nel
        use determinants, only: enumerate_determinants, find_sym_space_size, set_spin_polarisation
        use determinants, only: tot_ndets, ndets, sym_space_size
        use lanczos
        use full_diagonalisation
        use hamiltonian, only: get_hmatel
        use symmetry, only: nsym

        use utils, only: int_fmt
        use m_mrgref, only: mrgref

        type soln
            integer :: ms
            real(p) :: energy
        end type soln

        integer :: ms, ms_min, ms_max, ierr, nlanczos, nexact, nfound, i, ind, state, isym, sym_min, sym_max

        type(soln), allocatable :: lanczos_solns(:), exact_solns(:)

        ! For communication with Lanczos and exact diagonalisers.
        real(dp), allocatable :: lanczos_eigv(:)
        real(p), allocatable :: exact_eigv(:)

        ! For sorting the solutions by energy rather than by energy within each
        ! spin block.
        integer, allocatable :: ranking(:)

        character(50) :: fmt1

        if (parent) write (6,'(1X,a15,/,1X,15("="),/)') 'Diagonalisation'

        ! The Hamiltonian can be written in block diagonal form using the spin
        ! quantum number.  < D_1 | H | D_2 > /= 0 only if D_1 and D_2 have the
        ! same total spin.

        ! For a given number of electrons, n, the total spin can range from -n
        ! to +n in steps of 2.
        ! The -ms blocks are degenerate with the +ms blocks so only need to
        ! solve for ms >= 0.
        
        if (ms_in == huge(1)) then
            ms_min = mod(nel,2)
            ms_max = min(nel, nbasis-nel)
        else
            ! ms was set in input
            ms_min = ms_in
            ms_max = ms_in
        end if
        
        if (sym_in == huge(1)) then
            sym_min = 1
            sym_max = nsym
        else
            ! sym was set in input
            sym_min = sym_in
            sym_max = sym_in
        end if

        if (doing_calc(lanczos_diag)) allocate(lanczos_solns(nlanczos_eigv*(nel/2+1)*nbasis/2), stat=ierr)
        if (doing_calc(exact_diag)) allocate(exact_solns(tot_ndets), stat=ierr)

        ! Number of lanczos and exact solutions found.
        nlanczos = 0
        nexact = 0

        do ms = ms_min, ms_max, 2

            if (parent) write (6,'(1X,a35,'//int_fmt(ms)//',/)') 'Considering determinants with spin:', ms

            ! Find and set information about the space.
            call set_spin_polarisation(ms)
            call find_sym_space_size()

            ! Diagonalise each symmetry block in turn.
            do isym = sym_min, sym_max

                if (sym_space_size(isym)==0) then
                    if (parent) then
                        fmt1 = '(1X,a25,'//int_fmt(isym,1)//',1X,a17,'//int_fmt(nsym,1)//',1X,a6,'//int_fmt(nsym,1)//')'
                        write (6,fmt1) 'No determinants with spin',ms,'in symmetry block',isym,'out of',nsym
                        write (6,'(/,1X,15("-"),/)')
                    end if
                    cycle
                end if

                if (parent) then
                    fmt1 = '(1X,a28,'//int_fmt(isym,1)//',1X,a6,'//int_fmt(nsym,1)//')'
                    write (6,fmt1) 'Diagonalising symmetry block',isym,'out of',nsym
                end if

                ! Find all determinants with this spin.
                call enumerate_determinants(isym)

                if (ndets == 1) then

                    if (parent) then
                        write (6,'(/,1X,a35,/)') 'Performing trivial diagonalisation.'
                    end if

                    ! The trivial case seems to trip up TRLan and scalapack in
                    ! parallel.
                    if (doing_calc(lanczos_diag)) then
                        lanczos_solns(nlanczos+1)%energy = get_hmatel(1,1)
                        lanczos_solns(nlanczos+1)%ms = ms 
                        nlanczos = nlanczos + 1
                    end if
                    if (doing_calc(exact_diag)) then
                        exact_solns(nexact+1)%energy = get_hmatel(1,1)
                        exact_solns(nexact+1)%ms = ms 
                        nexact = nexact + 1
                    end if

                else

                    if (nprocs == 1) call generate_hamil(distribute_off)

                    ! Lanczos.
                    if (doing_calc(lanczos_diag)) then
                        allocate(lanczos_eigv(ndets), stat=ierr)
                        ! Construct the Hamiltonian matrix distributed over the processors
                        ! if running in parallel.
                        if (nprocs > 1 .and. .not.direct_lanczos) call generate_hamil(distribute_cols)
                        call lanczos_diagonalisation(nfound, lanczos_eigv)
                        lanczos_solns(nlanczos+1:nlanczos+nfound)%energy = lanczos_eigv(:nfound)
                        lanczos_solns(nlanczos+1:nlanczos+nfound)%ms = ms 
                        nlanczos = nlanczos + nfound
                        deallocate(lanczos_eigv)
                    end if

                    ! Exact diagonalisation.
                    ! Warning: this destroys the Hamiltonian matrix...
                    if (doing_calc(exact_diag)) then
                        allocate(exact_eigv(ndets), stat=ierr)
                        ! Construct the Hamiltonian matrix distributed over the processors
                        ! if running in parallel.
                        if (nprocs > 1) call generate_hamil(distribute_blocks)
                        call exact_diagonalisation(exact_eigv)
                        exact_solns(nexact+1:nexact+ndets)%energy = exact_eigv
                        exact_solns(nexact+1:nexact+ndets)%ms = ms 
                        nexact = nexact + ndets
                        deallocate(exact_eigv)
                    end if

                end if

                if (parent) write (6,'(1X,15("-"),/)')

            end do

        end do

        ! Output results.
        if (doing_calc(lanczos_diag) .and. parent) then
            write (6,'(1X,a31,/,1X,31("-"),/)') 'Lanczos diagonalisation results'
            allocate(ranking(nlanczos), stat=ierr)
            call mrgref(lanczos_solns(:nlanczos)%energy, ranking)
            write (6,'(1X,a8,3X,a4,3X,a12)') 'State','Spin','Total energy'
            state = 0
            do i = 1, nlanczos
                ind = ranking(i)
                state = state + 1
                write (6,'(1X,i6,2X,i6,1X,f18.12)') state, lanczos_solns(ind)%ms, lanczos_solns(ind)%energy
                if (lanczos_solns(ind)%ms /= 0) then
                    ! Also a degenerate -ms solution.
                    state = state + 1
                    write (6,'(1X,i6,2X,i6,1X,f18.12)') state, -lanczos_solns(ind)%ms, lanczos_solns(ind)%energy
                end if
            end do
            write (6,'(/,1X,a21,f18.12,/)') 'Lanczos ground state:', lanczos_solns(ranking(1))%energy
            deallocate(ranking, stat=ierr)
        end if

        if (doing_calc(exact_diag) .and. parent) then
            write (6,'(1X,a29,/,1X,29("-"),/)') 'Exact diagonalisation results'
            allocate(ranking(nexact), stat=ierr)
            call mrgref(exact_solns(:nexact)%energy, ranking)
            write (6,'(1X,a8,3X,a4,3X,a12)') 'State','Spin','Total energy'
            state = 0
            do i = 1, nexact
                ind = ranking(i)
                state = state + 1
                write (6,'(1X,i6,2X,i6,1X,f18.12)') state, exact_solns(ind)%ms, exact_solns(ind)%energy
                if (exact_solns(ind)%ms /= 0) then
                    ! Also a degenerate -ms solution.
                    state = state + 1
                    write (6,'(1X,i6,2X,i6,1X,f18.12)') state, -exact_solns(ind)%ms, exact_solns(ind)%energy
                end if
            end do
            write (6,'(/,1X,a19,f18.12,/)') 'Exact ground state:', exact_solns(ranking(1))%energy
            deallocate(ranking, stat=ierr)
        end if

        if (doing_calc(lanczos_diag)) deallocate(lanczos_solns, stat=ierr)
        if (doing_calc(exact_diag)) deallocate(exact_solns, stat=ierr)

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

        use utils, only: get_free_unit
        use errors
        use parallel

        use hamiltonian, only: get_hmatel
        use hubbard_real
        use determinants, only: ndets

        integer, intent(in), optional :: distribute_mode
        integer :: ierr, iunit, n1, n2, ind_offset
        integer :: i, j, ii, jj, ilocal, iglobal, jlocal, jglobal

        ! Store how many determinants have already been considered.
        ! This allows for the correct indexing scheme to be used in
        ! output of the Hamiltonian matrix.
        integer, save :: ndets_prev = 0

        if (allocated(hamil)) then
            deallocate(hamil, stat=ierr)
        end if

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

        allocate(hamil(n1,n2), stat=ierr)

        ! index offset for the symmetry block compared to the index of the
        ! determinant list.
        !ind_offset = ???
        ind_offset = 0

        ! Form the Hamiltonian matrix < D_i | H | D_j >.
        select case(distribute)
        case(distribute_off)
            forall (i=1:ndets) 
                forall (j=i:ndets) hamil(i,j) = get_hmatel(i+ind_offset,j+ind_offset)
            end forall
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
                            hamil(ilocal, jlocal) = get_hmatel(iglobal, jglobal)
                        end do
                    end do
                end do
            end do
        end select

        if (write_hamiltonian) then
            if (nprocs > 1) then
                if (parent) call warning('generate_hamil','Output of hamiltonian not implemented in parallel.',2)
            else
                iunit = get_free_unit()
                open(iunit, file=hamiltonian_file, status='unknown')
                do i=1, ndets
                    write (iunit,*) i,i,hamil(i,i)
                    do j=i+1, ndets
                        if (abs(hamil(i,j)) > depsilon) write (iunit,*) i+ndets_prev,j+ndets_prev,hamil(i,j)
                    end do
                end do
                close(iunit, status='keep')
                ndets_prev = ndets
            end if
        end if

    end subroutine generate_hamil

    subroutine end_hamil()

        ! Clean up hamiltonian module.

        integer :: ierr

        if (allocated(hamil)) deallocate(hamil, stat=ierr)

    end subroutine end_hamil

end module diagonalisation
