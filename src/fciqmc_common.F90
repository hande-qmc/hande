module fciqmc_common

! Module containing routines common to different fciqmc algorithms.

use fciqmc_data

implicit none

contains

! --- Utility routines ---

    subroutine select_ref_det()

        ! Change the reference determinant to be the determinant with the
        ! greatest population if it exceeds some threshold relative to the
        ! current reference determinant.

        ! TODO: make compatible with HFS.

        use basis, only: basis_length
        use calc, only: doing_calc, hfs_fciqmc_calc
        use determinants, only: decode_det, write_det

        use parallel
        use errors, only: stop_all

        integer, parameter :: particle_type = 1
        integer :: i, fmax(basis_length), max_pop
#ifdef PARALLEL
        integer :: in_data(2), out_data(2), ierr
#endif
        real(p) :: H00_max, H00_old
        logical :: updated

        H00_old = H00

        updated = .false.
        ! Find determinant with largest population.
        max_pop = 0
        do i = 1, tot_walkers
            if (abs(walker_population(particle_type,i)) > abs(max_pop)) then
                max_pop = walker_population(particle_type,i)
                fmax = walker_dets(:,i)
                H00_max = walker_data(particle_type, i)
            end if
        end do

        ! Only change reference determinant if the population is larger than the
        ! reference determinant by a given factor to avoid switching
        ! continuously between degenerate determinants.

        ! Note we don't broadcast the population of the new reference det as
        ! that is reset at the start of the next report loop anyway (and this
        ! routine should only be called at the end of the report loop).

#ifdef PARALLEL

        if (abs(max_pop) > ref_det_factor*abs(D0_population)) then
            in_data = (/ max_pop, iproc /)
        else if (iproc == D0_proc) then
            ! Ensure that D0_proc has the correct (average) population.
            in_data = (/ nint(D0_population), iproc /)
        else
            ! No det with sufficient population to become reference det on this
            ! processor.
            in_data = (/ 0, iproc /)
        end if

        call mpi_allreduce(in_data, out_data, 1, MPI_2INTEGER, MPI_MAXLOC, MPI_COMM_WORLD, ierr)

        if (out_data(1) /= nint(D0_population) .and. all(fmax /= f0)) then
            write (6,*) 'MPI_MAXLOC', out_data
            max_pop = out_data(1)
            updated = .true.
            D0_proc = out_data(2)
            f0 = fmax
            H00 = H00_max
            ! Broadcast updated data
            call mpi_bcast(f0, basis_length, mpi_det_integer, D0_proc, MPI_COMM_WORLD, ierr)
            call mpi_bcast(H00, 1, mpi_preal, D0_proc, MPI_COMM_WORLD, ierr)
        end if

#else

        if (abs(max_pop) > ref_det_factor*abs(D0_population) .and. all(fmax /= f0)) then
            updated = .true.
            f0 = fmax
            H00 = H00_max
        end if

#endif

        if (updated) then
            ! Update occ_list.
            call decode_det(f0, occ_list0)
            ! walker_data(1,i) holds <D_i|H|D_i> - H00_old.  Update.
            ! H00 is currently <D_0|H|D_0> - H00_old.
            ! Want walker_data(1,i) to be <D_i|H|D_i> - <D_0|H|D_0>
            ! We'll fix H00 later and avoid an extra tot_walkers*additions.
            do i = 1, tot_walkers
                walker_data(1,i) = walker_data(1,i) - H00
            end do
            ! The fold line in the folded spectrum approach is set (during
            ! initialisation) relative to the reference.
            fold_line = fold_line - H00
            ! Now set H00 = <D_0|H|D_0> so that future references to it are
            ! correct.
            H00 = H00 + H00_old
            if (doing_calc(hfs_fciqmc_calc)) call stop_all('select_ref_det', 'Not implemented for HFS.')
            if (parent) then
                write (6,'(1X,"#",1X,62("-"))')
                write (6,'(1X,"#",1X,"Changed reference det to:",1X)',advance='no')
                call write_det(f0, new_line=.true.)
                write (6,'(1X,"#",1X,"Population on old reference det (averaged over report loop):",f10.2)') D0_population
                write (6,'(1X,"#",1X,"Population on new reference det:",27X,i8)') max_pop
                write (6,'(1X,"#",1X,"E0 = <D0|H|D0> = ",f20.12)') H00
                write (6,'(1X,"#",1X,"Care should be taken with accumulating statistics before this point.")')
                write (6,'(1X,"#",1X,62("-"))')
            end if
        end if

    end subroutine select_ref_det

    subroutine find_single_double_prob(occ_list, psingle, pdouble)

        ! Calculate the probabilities of selecting a single or double excitation
        ! from a given determinant.  We assume all possible excitations (i.e.
        ! those with Hamiltonian matrix elements which are not zero by symmetry)
        ! are equally likely, so this amounts to finding the number of possible
        ! (symmetry-allowed) single and double excitations.
        !
        ! In:
        !    occ_list: integer list of occupied spin-orbitals in a determinant, D.
        ! Out:
        !    psingle: probability of attempting to spawn on a determinant
        !             connected to D by a single excitation.
        !    pdouble: probability of attempting to spawn on a determinant
        !             connected to D by a double excitation.

        use basis, only: basis_fns
        use system, only: nel, system_type, hub_k, hub_real, heisenberg, read_in, sym0, sym_max
        use point_group_symmetry, only: cross_product_pg_basis, cross_product_pg_sym, nbasis_sym_spin

        integer, intent(in) :: occ_list(nel)
        real(p), intent(out) :: psingle, pdouble

        integer :: i, j, virt_syms(2, sym0:sym_max), nsingles, ndoubles, isyma, isymb, ims1, ims2

        select case(system_type)
        case(hub_k)
            ! Only double excitations
            psingle = 0.0_p
            pdouble = 1.0_p
        case(hub_real,heisenberg)
            ! Only single excitations
            psingle = 1.0_p
            pdouble = 0.0_p
        case(read_in)

            ! Count number of basis functions in each symmetry.
            virt_syms = nbasis_sym_spin
            do i = 1, nel
                ! Convert -1->1 and 1->2 for spin index in arrays.
                ims1 = (basis_fns(occ_list(i))%ms+3)/2
                virt_syms(ims1,basis_fns(occ_list(i))%sym) = virt_syms(ims1,basis_fns(occ_list(i))%sym) - 1
            end do

            ! Count number of possible single excitations from the supplied
            ! determinant.
            ! Symmetry and spin must be conserved.
            nsingles = 0
            do i = 1, nel
                ! Convert -1->1 and 1->2 for spin index in arrays.
                ims1 = (basis_fns(occ_list(i))%ms+3)/2
                ! Can't excite into already occupied orbitals.
                nsingles = nsingles + virt_syms(ims1,basis_fns(occ_list(i))%sym)
            end do

            ! Count number of possible double excitations from the supplied
            ! determinant.
            ndoubles = 0
            do i = 1, nel
                ! Convert -1->1 and 1->2 for spin index in arrays.
                ims1 = (basis_fns(occ_list(i))%ms+3)/2
                do j = i+1, nel
                    ! Convert -1->1 and 1->2 for spin index in arrays.
                    ims2 = (basis_fns(occ_list(j))%ms+3)/2
                    do isyma = sym0, sym_max
                        ! Symmetry of the final orbital is determined (for Abelian
                        ! symmetries) from the symmetry of the first three.
                        isymb = cross_product_pg_sym(isyma, cross_product_pg_basis(occ_list(i),occ_list(j)))
                        if (isyma == isymb) then
                            if (ims1 == ims2) then
                                ! Cannot excit 2 electrons into the same spin-orbital.
                                ! Need to avoid double counting.
                                !  => number of unique pairs is identical to
                                !  number of elements in the strictly lower
                                !  triangle of a square matrix.
                                ndoubles = ndoubles + (virt_syms(ims1,isyma)*(virt_syms(ims2,isymb)-1))/2
                            else
                                ndoubles = ndoubles + virt_syms(ims1,isyma)*virt_syms(ims2,isymb)
                            end if
                        else if (isyma < isymb) then
                            ! isyma < isymb to avoid double counting.
                            ndoubles = ndoubles + virt_syms(ims1,isyma)*virt_syms(ims2,isymb)
                            ! can also have the opposite spin structure of
                            ! occupied orbitals have different spins.
                            if (ims1 /= ims2) ndoubles = ndoubles + virt_syms(ims2,isyma)*virt_syms(ims1,isymb)
                        end if
                    end do
                end do
            end do

            psingle = real(nsingles,p)/(nsingles+ndoubles)
            pdouble = real(ndoubles,p)/(nsingles+ndoubles)

        end select

    end subroutine find_single_double_prob

    subroutine load_balancing_report()

        ! Print out a load-balancing report when run in parallel showing how
        ! determinants and walkers/particles are distributed over the processors.

#ifdef PARALLEL
        use annihilation, only: annihilation_comms_time
        use parallel

        integer(lint) :: load_data(nprocs)
        integer :: ierr
        real(dp) :: comms(nprocs)

        if (nprocs > 1) then
            if (parent) then
                write (6,'(1X,a14,/,1X,14("^"),/)') 'Load balancing'
                write (6,'(1X,a77,/)') "The final distribution of walkers and determinants across the processors was:"
            endif
            call mpi_gather(nparticles, 1, mpi_integer8, load_data, 1, mpi_integer8, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                write (6,'(1X,a34,6X,i8)') 'Min # of particles on a processor:', minval(load_data)
                write (6,'(1X,a34,6X,i8)') 'Max # of particles on a processor:', maxval(load_data)
                write (6,'(1X,a35,5X,f11.2)') 'Mean # of particles on a processor:', real(sum(load_data), p)/nprocs
            end if
            call mpi_gather(tot_walkers, 1, mpi_integer, load_data, 1, mpi_integer8, 0, MPI_COMM_WORLD, ierr)
            call mpi_gather(annihilation_comms_time, 1, mpi_real8, comms, 1, mpi_real8, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                write (6,'(1X,a37,3X,i8)') 'Min # of determinants on a processor:', minval(load_data)
                write (6,'(1X,a37,3X,i8)') 'Max # of determinants on a processor:', maxval(load_data)
                write (6,'(1X,a38,2X,f11.2)') 'Mean # of determinants on a processor:', real(sum(load_data), p)/nprocs
                write (6,'()')
                write (6,'(1X,a38,5X,f8.2,a1)') 'Min time take by walker communication:', minval(comms),'s'
                write (6,'(1X,a38,5X,f8.2,a1)') 'Max time take by walker communication:', maxval(comms),'s'
                write (6,'(1X,a39,4X,f8.2,a1)') 'Mean time take by walker communication:', real(sum(comms), p)/nprocs,'s'
                write (6,'()')
            end if
        end if
#endif

    end subroutine load_balancing_report

! --- Output routines ---

    subroutine initial_fciqmc_status()

        ! Calculate the projected energy based upon the initial walker
        ! distribution (either via a restart or as set during initialisation)
        ! and print out.

        use parallel
        use proc_pointers, only: update_proj_energy_ptr
        integer :: idet
        integer(lint) :: ntot_particles
#ifdef PARALLEL
        integer :: ierr
        real(p) :: proj_energy_sum
#endif

        ! Calculate the projected energy based upon the initial walker
        ! distribution.  proj_energy and D0_population are both accumulated in
        ! update_proj_energy.
        proj_energy = 0.0_p
        D0_population = 0
        do idet = 1, tot_walkers
            call update_proj_energy_ptr(idet)
        end do

#ifdef PARALLEL
        call mpi_allreduce(proj_energy, proj_energy_sum, 1, mpi_preal, MPI_SUM, MPI_COMM_WORLD, ierr)
        proj_energy = proj_energy_sum
        call mpi_allreduce(nparticles, ntot_particles, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        ntot_particles = nparticles(1)
#endif

        proj_energy = proj_energy/D0_population

        if (parent) then
            ! See also the format used in write_fciqmc_report if this is changed.
            ! We prepend a # to make it easy to skip this point when do data
            ! analysis.
            write (6,'(1X,"#",3X,i8,2X,2(es17.10,2X),es17.10,4X,i11,6X,a3,3X,a3)') &
                    mc_cycles_done, shift, proj_energy, D0_population, ntot_particles,'n/a','n/a'
        end if

    end subroutine initial_fciqmc_status

end module fciqmc_common
