module hilbert_space

! Estimating the size of the Hilbert space.

implicit none

integer :: nhilbert_cycles

contains

    subroutine estimate_hilbert_space()

        ! Based on Appendix A in George Booth's thesis.

        ! Find the size of the Hilbert space which is of the same symmetry as
        ! the reference determinant.

        ! See find_sym_space_size for a dumb but exact enumeration of the size
        ! of the space (which is needed for FCI calculations).

        use basis, only: basis_length, bit_lookup, write_basis_fn, basis_fns, nbasis
        use calc, only: ms_in, truncate_space, truncation_level
        use const, only: dp
        use determinants, only: set_spin_polarisation, encode_det
        use excitations, only: get_excitation_level
        use dSFMT_interface, only: genrand_real2
        use fciqmc_data, only: occ_list0
        use reference_determinant, only: set_reference_det
        use symmetry, only: symmetry_orb_list
        use system
        use parallel
        use utils, only: binom_r

        integer :: iel, icycle, naccept
        integer :: a, a_el, a_pos, b, b_el, b_pos
        integer :: ref_sym, det_sym
        integer(i0) :: f(basis_length), f0(basis_length)
        integer :: occ_list(nel)
        real(dp) :: space_size
#ifdef PARALLEL
        integer :: ierr
        real(dp) :: proc_space_size(nprocs), sd_space_size
#endif

        if (parent) write (6,'(1X,a13,/,1X,13("-"),/)') 'Hilbert space'

        call set_spin_polarisation(ms_in)

        select case(system_type)

        case(heisenberg)

            ! Symmetry not currently implemented for the Heisenberg code.
            ! There is one spin per site, so it's just a case of how many ways
            ! there are to arrange the nalpha spins across the nsites (or
            ! equivalently the nbeta spins across the nsites).
            ! See comments in system for how nel and nvirt are used in the
            ! Heisenberg model.
            if (truncate_space) then
                space_size = binom_r(nsites-(nel-truncation_level),truncation_level)
            else
                space_size = binom_r(nsites, nel)
            end if
            if (parent) write (6,'(1X,a,g12.4,/)') 'Size of space is', space_size

        case default

            if (system_type == hub_real .and. .not.truncate_space) then
                ! Symmetry not currently implemented for the real space Hubbard
                ! code.
                ! Just a case of how we arrange the alpha and beta electrons across
                ! the alpha orbitals and beta orbitals.  As we're dealing with the
                ! simplest possible Hubbard model, the number of orbitals of each
                ! spin is equal to the number of sites.
                if (parent) write (6,'(1X,a,g12.4,/)') 'Size of space is', binom_r(nsites, nalpha)*binom_r(nsites, nbeta)
            else

                ! Perform a Monte Carlo sampling of the space.

                call set_reference_det(occ_list0, .false.)
                call encode_det(occ_list0, f0)

                ! Symmetry of the reference determinant.
                ref_sym = symmetry_orb_list(occ_list0)

                if (parent) then
                    write (6,'(1X,a34)',advance='no') 'Symmetry of reference determinant:'
                    if (momentum_space) then
                        call write_basis_fn(basis_fns(2*ref_sym), new_line=.true., print_full=.false.)
                    else
                        write (6,'(1X,i2)') ref_sym
                    end if
                end if

                naccept = 0

                do icycle = 1, nhilbert_cycles
                    ! Generate a random determinant.
                    ! Alpha electrons.
                    f = 0
                    iel = 0
                    do
                        ! generate random number 1,3,5,...
                        a = 2*int(genrand_real2()*(nbasis/2))+1
                        a_pos = bit_lookup(1,a)
                        a_el = bit_lookup(2,a)
                        if (.not.btest(f(a_el), a_pos)) then
                            ! found unoccupied alpha orbital.
                            f(a_el) = ibset(f(a_el), a_pos)
                            iel = iel + 1
                            occ_list(iel) = a
                            if (iel == nalpha) exit
                        end if
                    end do
                    ! Beta electrons.
                    do
                        ! generate random number 2,4,6,...
                        b = 2*int(genrand_real2()*(nbasis/2))+2
                        b_pos = bit_lookup(1,b)
                        b_el = bit_lookup(2,b)
                        if (.not.btest(f(b_el), b_pos)) then
                            ! found unoccupied beta orbital.
                            f(b_el) = ibset(f(b_el), b_pos)
                            iel = iel + 1
                            occ_list(iel) = b
                            if (iel == nel) exit
                        end if
                    end do
                    ! Find the symmetry of the determinant.
                    det_sym = symmetry_orb_list(occ_list)
                    ! Is this the same symmetry as the reference determinant?
                    if (det_sym == ref_sym) then
                        if (truncate_space) then
                            if (get_excitation_level(f,f0) <= truncation_level) &
                                                                naccept = naccept + 1
                        else
                            naccept = naccept + 1
                        end if
                    end if
                end do

                ! Size of the Hilbert space in the desired symmetry block is given
                ! by
                !   C(nalpha_orbitals, nalpha_electrons)*C(nbeta_orbitals, nbeta_electrons)*naccept/nattempts
                space_size = (binom_r(nbasis/2, nalpha) * binom_r(nbasis/2, nbeta) * naccept) / nhilbert_cycles

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

        end select

    end subroutine estimate_hilbert_space

end module hilbert_space
