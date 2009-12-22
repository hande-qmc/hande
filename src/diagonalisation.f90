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
        use determinants, only: enumerate_determinants, tot_ndets, ndets
        use lanczos
        use full_diagonalisation

        use utils, only: int_fmt
        use m_mrgref, only: mrgref

        type soln
            integer :: ms
            real(dp) :: energy
        end type soln

        integer :: ms, ms_min, ierr, nlanczos, nexact, nfound, i, ind, state, isym

        type(soln), allocatable :: lanczos_solns(:), exact_solns(:)

        ! For communication with Lanczos and exact diagonalisers.
        real(dp), allocatable :: eigv(:)

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
        
        ms_min = mod(nel,2)

        if (t_lanczos) allocate(lanczos_solns(nlanczos_eigv*(nel/2+1)*nbasis/2), stat=ierr)
        if (t_exact) allocate(exact_solns(tot_ndets), stat=ierr)

        ! Number of lanczos and exact solutions found.
        nlanczos = 0
        nexact = 0

        do ms = ms_min, nel, 2

            if (parent) write (6,'(1X,a58,'//int_fmt(ms)//')') 'Diagonalising Hamiltonian block of determinants with spin:', ms

            ! Find all determinants with this spin.
            call enumerate_determinants(ms)

            ! Symmetry blocks within the spin block.
            call find_sym_blocks()

            ! Diagonalise each symmetry block in turn.
            do isym = 1, nsym_blocks

                if (parent) then
                    fmt1 = '(1X,a28,'//int_fmt(isym,1)//',1X,a6,'//int_fmt(nsym_blocks,1)//')'
                    write (6,fmt1) 'Diagonalising symmetry block',isym,'out of',nsym_blocks
                end if

                if (isym == nsym_blocks) then
                    nhamil = ndets - sym_blocks(isym) + 1
                else
                    nhamil = sym_blocks(isym+1) - sym_blocks(isym)
                end if

                allocate(eigv(nhamil), stat=ierr)

                if (nprocs == 1) call generate_hamil(isym, distribute_off)

                ! Lanczos.
                if (t_lanczos) then
                    ! Construct the Hamiltonian matrix distributed over the processors
                    ! if running in parallel.
                    if (nprocs > 1 .and. .not.direct_lanczos) call generate_hamil(isym, distribute_cols)
                    call lanczos_diagonalisation(nfound, eigv)
                    lanczos_solns(nlanczos+1:nlanczos+nfound)%energy = eigv(:nfound)
                    lanczos_solns(nlanczos+1:nlanczos+nfound)%ms = ms 
                    nlanczos = nlanczos + nfound
                end if

                ! Exact diagonalisation.
                ! Warning: this destroys the Hamiltonian matrix...
                if (t_exact) then
                    ! Construct the Hamiltonian matrix distributed over the processors
                    ! if running in parallel.
                    if (nprocs > 1) call generate_hamil(isym, distribute_blocks)
                    call exact_diagonalisation(eigv)
                    exact_solns(nexact+1:nexact+nhamil)%energy = eigv
                    exact_solns(nexact+1:nexact+nhamil)%ms = ms 
                    nexact = nexact + nhamil
                end if

                deallocate(eigv)

                if (parent) write (6,'(1X,15("-"),/)')

            end do

        end do

        ! Output results.
        if (t_lanczos .and. parent) then
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

        if (t_exact .and. parent) then
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

        if (t_lanczos) deallocate(lanczos_solns, stat=ierr)
        if (t_exact) deallocate(exact_solns, stat=ierr)

    end subroutine diagonalise

    subroutine find_sym_blocks()

        ! Identify the symmetry blocks in the determinant list.
        ! This sets up sym_blocks, an array which contains the index of the
        ! first determinant in each symmetry block.
        ! Currently only relevant to momentum space.

        use basis, only: nbasis
        use determinants, only: ndets, dets
        use kpoints, only: is_reciprocal_lattice_vector 
        use system, only: system_type, hub_real

        integer :: i, iblk, ierr

        ! This is the maximum size of sym_blocks.  Allocating it to be this size
        ! means we can use it throughout an entire calculation.  The memory
        ! wastage is insignificant.
        if (.not.allocated(sym_blocks)) allocate(sym_blocks(nbasis/2), stat=ierr)

        if (system_type == hub_real) then
            nsym_blocks = 1
            sym_blocks(1) = 1
        else
            ! Hamiltonian matrix elements, < D_i | H | D_j >, are only non-zero
            ! if k_i and k_j differ by (at most) a reciprocal lattice vector of
            ! the primitive cell.
            ! The list of determinants are helpfully sorted in this order, so we
            ! just need to find which determinants are the first in each
            ! symmetry block.
            iblk = 1
            sym_blocks(iblk) = 1
            do i = 2, ndets
                if (is_reciprocal_lattice_vector(dets(i)%k - dets(i-1)%k)) then
                    cycle
                else
                    ! New block!
                    iblk = iblk + 1
                    sym_blocks(iblk) = i
                end if
            end do
            nsym_blocks = iblk
        end if

    end subroutine find_sym_blocks

    subroutine generate_hamil(isym, distribute_mode)

        ! Generate a symmetry block of the Hamiltonian matrix, H = < D_i | H | D_j >.
        ! The list of determinants, {D_i}, is grouped by symmetry and contains
        ! only determinants of a specified spin.
        ! Only generate the upper diagonal for use with (sca)lapack and Lanczos routines.
        ! In:
        !    isym: index of the symmetry block of the Hamiltonian matrix to
        !    generate.
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

        use determinants, only: dets
        use hamiltonian, only: get_hmatel
        use hubbard_real

        integer, intent(in) :: isym
        integer, intent(in), optional :: distribute_mode
        integer :: ierr, iunit, n1, n2, ind_offset
        integer :: i, j, ii, jj, ilocal, iglobal, jlocal, jglobal

        ! Store how many determinants have already been considered.
        ! This allows for the correct indexing scheme to be used in
        ! output of the Hamiltonian matrix.
        integer, save :: nhamil_prev = 0

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
            proc_blacs_info = get_blacs_info(nhamil, (/1, nprocs/))
            n1 = nhamil
            n2 = nhamil
        case(distribute_blocks)
            ! Use as square a processor grid as possible.
            proc_blacs_info = get_blacs_info(nhamil)
            n1 = proc_blacs_info%nrows
            n2 = proc_blacs_info%ncols
        case(distribute_cols)
            ! TRLan assumes that the only the rows of the matrix
            ! are distributed.
            ! Furthermore it seems TRLan assumes all processors store at least
            ! some of the matrix.
            if ((nprocs-1)*block_size > nhamil) then
                if (parent) then
                    call warning('generate_hamil','Reducing block size so that all processors contain at least a single row.', 3)
                    write (6,'(1X,a69)') 'Consider running on fewer processors or reducing block size in input.'
                    write (6,'(1X,a19,i4)') 'Old block size was:',block_size
                end if
                block_size = nhamil/nprocs
                if (parent) write (6,'(1X,a18,i4,/)') 'New block size is:',block_size
            end if
            proc_blacs_info = get_blacs_info(nhamil, (/1, nprocs/))
            n1 = proc_blacs_info%nrows
            n2 = proc_blacs_info%ncols
        case default
            call stop_all('generate_hamil','Unknown distribution scheme.')
        end select

        allocate(hamil(n1,n2), stat=ierr)

        ! index offset for the symmetry block compared to the index of the
        ! determinant list.
        ind_offset = sym_blocks(isym) - 1

        ! Form the Hamiltonian matrix < D_i | H | D_j >.
        select case(distribute)
        case(distribute_off)
            forall (i=1:nhamil) 
                forall (j=i:nhamil) hamil(i,j) = get_hmatel(i+ind_offset,j+ind_offset)
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
                do i=1, nhamil
                    write (iunit,*) i,i,hamil(i,i)
                    do j=i+1, nhamil
                        if (abs(hamil(i,j)) > depsilon) write (iunit,*) i+nhamil_prev,j+nhamil_prev,hamil(i,j)
                    end do
                end do
                close(iunit, status='keep')
                nhamil_prev = nhamil
            end if
        end if

    end subroutine generate_hamil

    subroutine end_hamil()

        ! Clean up hamiltonian module.

        integer :: ierr

        if (allocated(hamil)) deallocate(hamil, stat=ierr)
        if (allocated(sym_blocks)) deallocate(sym_blocks, stat=ierr)

    end subroutine end_hamil

end module diagonalisation
