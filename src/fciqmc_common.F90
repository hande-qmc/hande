module fciqmc_common

! Module containing routines common to different fciqmc algorithms.

use fciqmc_data
implicit none

contains

    subroutine init_fciqmc()

        ! Initialisation for fciqmc calculations.
        ! Setup the spin polarisation for the system, initialise the RNG,
        ! allocate the required memory for the list of walkers and set the
        ! initial walker.

        use errors, only: stop_all
        use hashing, only: murmurhash_bit_string
        use parallel, only: iproc, nprocs, parent
        use utils, only: int_fmt

        use basis, only: basis_length
        use calc, only: sym_in, ms_in, initiator_fciqmc, hfs_fciqmc_calc, doing_calc
        use determinants, only: encode_det, set_spin_polarisation, write_det
        use hamiltonian, only: get_hmatel_real, slater_condon0_hub_real, slater_condon0_hub_k
        use fciqmc_restart, only: read_restart
        use system, only: nel, system_type, hub_real, hub_k

        integer :: ierr
        integer :: i, iproc_ref
        integer :: step

        if (parent) write (6,'(1X,a6,/,1X,6("-"),/)') 'FCIQMC'

        ! Array sizes depending upon FCIQMC algorithm.
        sampling_size = 1
        spawned_size = basis_length + 1
        spawned_pop = spawned_size
        if (doing_calc(hfs_fciqmc_calc)) then
            spawned_size = spawned_size + 1
            spawned_hf_pop = spawned_size
            sampling_size = sampling_size + 1
        else
            spawned_hf_pop = spawned_size
        end if
        if (doing_calc(initiator_fciqmc)) then
            spawned_size = spawned_size + 1
            spawned_parent = spawned_size
        end if

        ! Allocate main walker lists.
        allocate(nparticles(sampling_size), stat=ierr)
        allocate(walker_dets(basis_length,walker_length), stat=ierr)
        allocate(walker_population(sampling_size,walker_length), stat=ierr)
        allocate(walker_energies(sampling_size,walker_length), stat=ierr)

        ! Allocate spawned walker lists.
        if (mod(spawned_walker_length, nprocs) /= 0) then
            if (parent) write (6,'(1X,a68)') 'spawned_walker_length is not a multiple of the number of processors.'
            spawned_walker_length = ceiling(real(spawned_walker_length)/nprocs)*nprocs
            if (parent) write (6,'(1X,a35,'//int_fmt(spawned_walker_length,1)//',a1,/)') &
                                        'Increasing spawned_walker_length to',spawned_walker_length,'.'
        end if
        allocate(spawned_walkers1(spawned_size,spawned_walker_length), stat=ierr)
        spawned_walkers => spawned_walkers1
        ! Allocate scratch space for doing communication.
        allocate(spawned_walkers2(spawned_size,spawned_walker_length), stat=ierr)
        spawned_walkers_recvd => spawned_walkers2

        ! Set spawning_head to be the same size as spawning_block_start.
        allocate(spawning_head(0:max(1,nprocs-1)), stat=ierr)

        ! Find the start position within the spawned walker lists for each
        ! processor.
        ! spawning_block_start(1) should contain the number of elements allocated
        ! for each processor so we allow it to be accessible even if the number
        ! of processors is 1.
        allocate(spawning_block_start(0:max(1,nprocs-1)), stat=ierr)
        step = spawned_walker_length/nprocs
        forall (i=0:nprocs-1) spawning_block_start(i) = i*step

        ! Set spin variables.
        call set_spin_polarisation(ms_in)

        ! Set initial walker population.
        ! occ_list could be set and allocated in the input.
        allocate(f0(basis_length), stat=ierr)
        if (restart) then
            if (.not.allocated(occ_list0)) allocate(occ_list0(nel), stat=ierr)
            call read_restart()
        else
            tot_walkers = 1
            walker_population(1,tot_walkers) = D0_population

            ! Reference det
            ! Set the reference determinant to be the spin-orbitals with the lowest
            ! kinetic energy which satisfy the spin polarisation.
            ! Note: this is for testing only!  The symmetry input is currently
            ! ignored.
            call set_reference_det()

            call encode_det(occ_list0, walker_dets(:,tot_walkers))

            ! Reference det
            f0 = walker_dets(:,tot_walkers)
            ! Energy of reference determinant.
            select case(system_type)
            case(hub_k)
                H00 = slater_condon0_hub_k(f0)
            case(hub_real)
                H00 = slater_condon0_hub_real(f0)
            end select
            ! By definition:
            walker_energies(1,tot_walkers) = 0.0_p

            ! Finally, we need to check if the reference determinant actually
            ! belongs on this processor.
            ! If it doesn't, set the walkers array to be empty.
            if (nprocs > 1) then
                iproc_ref = modulo(murmurhash_bit_string(f0, basis_length), nprocs)
                if (iproc_ref /= iproc) tot_walkers = 0
                D0_proc = iproc_ref
            else
                D0_proc = iproc
            end if
        end if

        ! Total number of particles on processor.
        ! Probably should be handled more simply by setting it to be either 0 or
        ! D0_population or obtaining it from the restart file, as appropriate.
        nparticles = sum(abs(walker_population(:,:tot_walkers)))

        if (parent) then
            write (6,'(1X,a29,1X)',advance='no') 'Reference determinant, |D0> ='
            call write_det(f0, new_line=.true.)
            write (6,'(1X,a16,f20.12)') 'E0 = <D0|H|D0> =',H00
            write (6,'(1X,a44,'//int_fmt(D0_population,1)//',/)') &
                              'Initial population on reference determinant:',D0_population
            write (6,'(1X,a68,/)') 'Note that FCIQMC calculates the correlation energy relative to |D0>.'
            if (doing_calc(initiator_fciqmc)) then
                write (6,'(1X,a24)') 'Initiator method in use.'
                write (6,'(1X,a36,1X,"(",'//int_fmt(CAS(1),0)//',",",'//int_fmt(CAS(2),0)//'")")')  &
                    'CAS space of initiator determinants:',CAS
                write (6,'(1X,a66,'//int_fmt(initiator_population,1)//',/)') &
                    'Population for a determinant outside CAS space to be an initiator:', initiator_population
            end if
            write (6,'(1X,a49,/)') 'Information printed out every FCIQMC report loop:'
            write (6,'(1X,a66)') 'Instant shift: the shift calculated at the end of the report loop.'
            write (6,'(1X,a88)') 'Average shift: the running average of the shift from when the shift was allowed to vary.'
            write (6,'(1X,a98)') 'Proj. Energy: projected energy averaged over the report loop. &
                                 &Calculated at the end of each cycle.'
            write (6,'(1X,a53)') 'Av. Proj. E: running average of the projected energy.'
            write (6,'(1X,a54)') '# D0: current population at the reference determinant.'
            write (6,'(1X,a49)') '# particles: current total population of walkers.'
            write (6,'(1X,a56,/)') 'R_spawn: average rate of spawning across all processors.'
        end if
        
    end subroutine init_fciqmc

    subroutine initial_fciqmc_status(update_proj_energy)

        ! Calculate the projected energy based upon the initial walker
        ! distribution (either via a restart or as set during initialisation)
        ! and print out.

        ! In:
        !    update_proj_energy: relevant subroutine to update the projected
        !        energy.  See the energy_evaluation module.

        use parallel

        interface
            subroutine update_proj_energy(idet, inst_proj_energy)
                use const, only: p
                implicit none
                integer, intent(in) :: idet
                real(p), intent(inout) :: inst_proj_energy
            end subroutine update_proj_energy
        end interface

        integer :: idet, ierr
        integer :: ntot_particles
        real(p) :: inst_proj_energy

        ! Calculate the projected energy based upon the initial walker
        ! distribution.
        inst_proj_energy = 0.0_p
        do idet = 1, tot_walkers 
            call update_proj_energy(idet, inst_proj_energy)
        end do 

#ifdef PARALLEL
        call mpi_allreduce(inst_proj_energy, proj_energy, 1, mpi_preal, MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(nparticles, ntot_particles, sampling_size, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        proj_energy = inst_proj_energy
        ntot_particles = nparticles(1)
#endif 
        
        proj_energy = proj_energy/D0_population

        if (parent) then
            ! See also the format used in write_fciqmc_report if this is changed.
            ! We prepend a # to make it easy to skip this point when do data
            ! analysis.
            write (6,'(1X,"#",3X,i8,2X,4(f14.10,2X),i11,2X,i11,6X,a3,3X,a3)') &
                    mc_cycles_done, shift, 0.0_p, proj_energy, 0.0_p, D0_population, ntot_particles,'n/a','n/a'
        end if

    end subroutine initial_fciqmc_status

    subroutine load_balancing_report()

        ! Print out a load-balancing report when run in parallel showing how
        ! determinants and walkers/particles are distributed over the processors.

#ifdef PARALLEL
        use parallel

        integer :: load_data(nprocs), ierr

        if (nprocs > 1) then
            if (parent) then
                write (6,'(1X,a14,/,1X,14("^"),/)') 'Load balancing'
                write (6,'(1X,a77,/)') "The final distribution of walkers and determinants across the processors was:"
            endif
            call mpi_gather(nparticles, 1, mpi_integer, load_data, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                write (6,'(1X,a34,6X,i8)') 'Min # of particles on a processor:', minval(load_data)
                write (6,'(1X,a34,6X,i8)') 'Max # of particles on a processor:', maxval(load_data)
                write (6,'(1X,a35,5X,f11.2)') 'Mean # of particles on a processor:', real(sum(load_data), p)/nprocs
            end if
            call mpi_gather(tot_walkers, 1, mpi_integer, load_data, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
            if (parent) then
                write (6,'(1X,a37,3X,i8)') 'Min # of determinants on a processor:', minval(load_data)
                write (6,'(1X,a37,3X,i8)') 'Max # of determinants on a processor:', maxval(load_data)
                write (6,'(1X,a38,2X,f11.2)') 'Mean # of determinants on a processor:', real(sum(load_data), p)/nprocs
                write (6,'()')
            end if
        end if
#endif

    end subroutine load_balancing_report

end module fciqmc_common
