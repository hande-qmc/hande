module hilbert_space

! Estimating the size of the Hilbert space.

implicit none

integer :: nhilbert_cycles

contains

    subroutine estimate_hilbert_space()

        ! Based on Appendix A in GHB's thesis.

        ! Find the size of the Hilbert space which is of the same symmetry as
        ! the reference determinant.

        ! See find_sym_space_size for a dumb but exact enumeration of the size
        ! of the space (which is needed for FCI calculations).

        use basis, only: basis_length, bit_lookup, write_basis_fn, basis_fns
        use calc, only: ms_in
        use const, only: dp
        use determinants, only: decode_det, set_spin_polarisation
        use dSFMT_interface, only: genrand_real2
        use fciqmc_data, only: occ_list0, set_reference_det
        use symmetry, only: sym_table
        use system
        use parallel
        use utils, only: binom_r

        integer :: i, iel, icycle, naccept
        integer :: a, a_el, a_pos, b, b_el, b_pos
        integer :: ref_sym, det_sym
        integer(i0) :: f(basis_length)
        integer :: occ_list(nel)
        real(dp) :: space_size
#ifdef PARALLEL
        integer :: ierr
        real(dp) :: proc_space_size(nprocs), sd_space_size
#endif

        if (parent) write (6,'(1X,a13,/,1X,13("-"),/)') 'Hilbert space'

        if (system_type == hub_real) then

            ! Symmetry not currently implemented for the real space Hubbard
            ! code.
            if (parent) write (6,'(1X,a,g8.4,/)') 'Size of space is', binom_r(nsites, nalpha)*binom_r(nsites, nbeta)

        else if (system_type == hub_k) then

            ! Perform a Monte Carlo sampling of the space.

            call set_spin_polarisation(ms_in)
            call set_reference_det()

            ! Symmetry of the reference determinant.
            ref_sym = 1
            do iel = 1, nel
                ref_sym = sym_table((occ_list0(iel)+1)/2, ref_sym)
            end do

            if (parent) then
                write (6,'(1X,a34)',advance='no') 'Symmetry of reference determinant:'
                call write_basis_fn(basis_fns(2*ref_sym), new_line=.true., print_full=.false.)
            end if

            naccept = 0

            do icycle = 1, nhilbert_cycles
                ! Generate a random determinant.
                ! Alpha electrons.
                f = 0
                do i = 1, nalpha
                    do
                        a = 2*nint(genrand_real2()*(nsites-1))+1
                        a_pos = bit_lookup(1,a)
                        a_el = bit_lookup(2,a)
                        if (.not.btest(f(a_el), a_pos)) then
                            ! found unoccupied alpha orbital.
                            f(a_el) = ibset(f(a_el), a_pos)
                            exit
                        end if
                    end do
                end do
                ! Beta electrons.
                do i = 1, nbeta
                    do
                        b = 2*nint(genrand_real2()*(nsites-1))+2
                        b_pos = bit_lookup(1,b)
                        b_el = bit_lookup(2,b)
                        if (.not.btest(f(b_el), b_pos)) then
                            ! found unoccupied beta orbital.
                            f(b_el) = ibset(f(b_el), b_pos)
                            exit
                        end if
                    end do
                end do
                ! Find the symmetry of the determinant.
                call decode_det(f, occ_list)
                det_sym = 1
                do iel = 1, nel
                    det_sym = sym_table((occ_list(iel)+1)/2, det_sym)
                end do
                ! Is this the same symmetry as the reference determinant?
                if (det_sym == ref_sym) naccept = naccept + 1
            end do

            ! Size of the Hilbert space in the desired symmetry block is given
            ! by
            !   C(nalpha_orbitals, nalpha_electrons)*C(nbeta_orbitals, nbeta_electrons)*naccept/nattempts
            space_size = (binom_r(nsites, nalpha) * binom_r(nsites, nbeta) * naccept) / nhilbert_cycles

#ifdef PARALLEL
            ! If we did this on multiple processors then we can get an estimate
            ! of the error as well as a better mean...
            call mpi_gather(space_size, 1, mpi_real8, proc_space_size, 1, mpi_real8, root, mpi_comm_world, ierr)
            space_size = sum(proc_space_size)/nprocs
            sd_space_size = sqrt(sum((proc_space_size-space_size)**2))/(nprocs-1)
            if (parent) then
                write (6,'(1X,a41,1X,es10.4)',advance='no') 'Monte-Carlo estimate of size of space is:', space_size
                if (nprocs > 1) then
                    write (6,'(1X,a3,1X,es10.4)') '+/-', sd_space_size
                else
                    write (6,'()')
                end if
            end if
#else
            if (parent) write (6,'(1X,a41,1X,es10.4)') 'Monte-Carlo estimate of size of space is:', space_size
#endif

            if (parent) write (6,'()')

        end if

    end subroutine estimate_hilbert_space

end module hilbert_space
