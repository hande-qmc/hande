module dmqmc_estimators

use const

implicit none

contains

   subroutine update_dmqmc_estimators(ireport, ntot_particles_old)

        ! Update the shift and average the shift and estimators

        ! Should be called every report loop in an DMQMC calculation.

        ! In:
        !    ireport: index of the current report loop.
        ! Inout:
        !    ntot_particles_old: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop.

        use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use energy_evaluation, only: update_shift
        use fciqmc_data, only: nparticles, sampling_size, target_particles, rspawn
        use fciqmc_data, only: shift, av_shift, vary_shift, start_vary_shift
        use fciqmc_data, only: nreport, ncycles, trace, thermal_energy, thermal_staggered_mag
        use fciqmc_data, only: total_trace, total_estimator_numerators, number_dmqmc_estimators
        use parallel

        integer, intent(in) :: ireport
        integer, intent(inout) :: ntot_particles_old(sampling_size)
        integer :: counter

#ifdef PARALLEL
        real(dp) :: ir(sampling_size+1+((number_dmqmc_estimators+1)*ncycles))
        real(dp) :: ir_sum(sampling_size+1+((number_dmqmc_estimators+1)*ncycles))
        integer :: ntot_particles(sampling_size), ierr

        ! Need to sum the number of particles and other quantites over all processors.
        ir(1:sampling_size) = nparticles
        ir(sampling_size+1) = rspawn
        ir(sampling_size+2:sampling_size+1+ncycles) = trace
        counter = 0
        if (doing_dmqmc_calc(dmqmc_energy)) then
            counter = counter + 1
            ir( sampling_size+2+(counter*ncycles):sampling_size+1+((counter+1)*ncycles) )=thermal_energy
        end if
        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
            counter = counter + 1
            ir( sampling_size+2+(counter*ncycles):sampling_size+1+((counter+1)*ncycles) )=thermal_staggered_mag
        end if
        call mpi_allreduce(ir, ir_sum, size(ir), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        ntot_particles = nint(ir_sum(1:sampling_size))
        rspawn = ir_sum(sampling_size+1)
        if (parent) then
            total_trace = ir_sum(sampling_size+2:sampling_size+1+ncycles)
            do counter = 1, number_dmqmc_estimators
               total_estimator_numerators(:,counter) = &
                  + ir_sum( sampling_size+2+(counter*ncycles):sampling_size+1+((counter+1)*ncycles) )
            end do
        end if

        if (vary_shift) then
            call update_shift(ntot_particles_old(1), ntot_particles(1), ncycles)
        end if
        ntot_particles_old = ntot_particles
        if (ntot_particles(1) > target_particles .and. .not.vary_shift) then
            vary_shift = .true.
            start_vary_shift = ireport
        end if

#else
        total_trace = trace
        counter = 0
        if (doing_dmqmc_calc(dmqmc_energy)) then
            counter = counter+1
            total_estimator_numerators(:,counter) = thermal_energy
        end if
        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
            counter = counter+1
            total_estimator_numerators(:,counter) = thermal_staggered_mag
        end if

        if (vary_shift) call update_shift(ntot_particles_old(1), nparticles(1), ncycles)
        ntot_particles_old = nparticles
        if (nparticles(1) > target_particles .and. .not.vary_shift) then
            vary_shift = .true.
            start_vary_shift = ireport
        end if        
#endif

            rspawn = rspawn/(ncycles*nprocs)

   end subroutine update_dmqmc_estimators

   subroutine call_dmqmc_estimators(idet, beta_index)
       
       use basis, only: basis_length, total_basis_length
       use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
       use excitations, only: get_excitation, excit
       use fciqmc_data, only: walker_dets, walker_population, trace
       use proc_pointers, only: update_dmqmc_energy_ptr, update_dmqmc_stag_mag_ptr

       integer, intent(in) :: idet, beta_index
       type(excit) :: excitation     
   
       excitation = get_excitation(walker_dets(:basis_length,idet), &
                        walker_dets((1+basis_length):total_basis_length,idet))

       if (excitation%nexcit == 0) trace(beta_index)=trace(beta_index) + walker_population(1,idet)
       if (doing_dmqmc_calc(dmqmc_energy)) call update_dmqmc_energy_ptr(idet, beta_index, excitation)
       if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) &
                                     call update_dmqmc_stag_mag_ptr(idet, beta_index, excitation)

   end subroutine call_dmqmc_estimators

   subroutine dmqmc_energy_heisenberg(idet, beta_index, excitation)

       use basis, only: basis_length, total_basis_length
       use basis, only: bit_lookup
       use excitations, only: excit
       use fciqmc_data, only: walker_dets, walker_population
       use fciqmc_data, only: walker_energies, thermal_energy, H00
       use hubbard_real, only: connected_orbs
       use system, only: J_coupling

       integer, intent(in) :: idet, beta_index
       type(excit), intent(in) :: excitation
       integer :: bit_element, bit_position

       if (excitation%nexcit == 0) then
           thermal_energy(beta_index)=thermal_energy(beta_index) + &
                             (walker_energies(1,idet)+H00)*walker_population(1,idet)
       else if (excitation%nexcit == 1) then
           bit_position = bit_lookup(1,excitation%from_orb(1))
           bit_element = bit_lookup(2,excitation%from_orb(1))
           if (btest(connected_orbs(bit_element, excitation%to_orb(1)), bit_position)) &
                 thermal_energy(beta_index)=thermal_energy(beta_index) - &
                             (2.0*J_coupling*walker_population(1,idet))
       end if

   end subroutine dmqmc_energy_heisenberg

   subroutine dmqmc_stag_mag_heisenberg(idet, beta_index, excitation)

       use basis, only: basis_length, total_basis_length
       use bit_utils, only: count_set_bits
       use determinants, only: lattice_mask
       use excitations, only: excit
       use fciqmc_data, only: walker_dets, walker_population, thermal_staggered_mag
       use system, only: nel

       integer, intent(in) :: idet, beta_index
       type(excit), intent(in) :: excitation
       integer :: bit_element, bit_position, stag_mag_el, sublattice1_up_spins
       integer(i0) :: f_mask(basis_length)

       if (excitation%nexcit == 0) then
           f_mask = iand(walker_dets(:basis_length,idet), lattice_mask)
           sublattice1_up_spins = sum(count_set_bits(f_mask))
           stag_mag_el = (4*sublattice1_up_spins - 2*nel)
           thermal_staggered_mag(beta_index)=thermal_staggered_mag(beta_index) + &
                             stag_mag_el*walker_population(1,idet)

       end if

   end subroutine dmqmc_stag_mag_heisenberg

end module dmqmc_estimators
