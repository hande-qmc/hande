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

        use energy_evaluation, only: update_shift
        use fciqmc_data, only: nparticles, sampling_size, target_particles, rspawn
        use fciqmc_data, only: shift, av_shift, vary_shift, start_vary_shift
        use fciqmc_data, only: nreport, ncycles, trace, thermal_energy
        use fciqmc_data, only: total_trace, total_thermal_energy
        use parallel

        integer, intent(in) :: ireport
        integer, intent(inout) :: ntot_particles_old(sampling_size)

#ifdef PARALLEL
        real(dp) :: ir(sampling_size+1+(2*ncycles)), ir_sum(sampling_size+1+(2*ncycles))
        integer :: ntot_particles(sampling_size), ierr

            ! Need to sum the number of particles and other quantites over all processors.
            ir(1:sampling_size) = nparticles
            ir(sampling_size+1) = rspawn
            ir(sampling_size+2:sampling_size+1+ncycles) = trace
            ir(sampling_size+2+ncycles:sampling_size+1+(2*ncycles)) = thermal_energy
            call mpi_allreduce(ir, ir_sum, size(ir), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
            ntot_particles = nint(ir_sum(1:sampling_size))
            rspawn = ir_sum(sampling_size+1)
            if (parent) then
                total_trace = total_trace + ir_sum(sampling_size+2:sampling_size+1+ncycles)
                total_thermal_energy = total_thermal_energy + ir_sum(sampling_size+2+ncycles:&
                                                                     sampling_size+1+(2*ncycles))
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
       use basis, only: bit_lookup
       use excitations, only: get_excitation, excit
       use fciqmc_data, only: walker_dets, walker_population, trace
       use fciqmc_data, only: walker_energies, thermal_energy, H00
       use hubbard_real, only: connected_orbs
       use system, only: J_coupling

       integer, intent(in) :: idet, beta_index
       type(excit) :: excitation
       integer :: bit_element, bit_position

       excitation = get_excitation(walker_dets(:basis_length,idet), &
                        walker_dets((1+basis_length):total_basis_length,idet))

       if (excitation%nexcit == 0) then
           thermal_energy(beta_index)=thermal_energy(beta_index) + &
                             (walker_energies(1,idet)+H00)*walker_population(1,idet)
           trace(beta_index)=trace(beta_index) + walker_population(1,idet)
       else if (excitation%nexcit == 1) then
           bit_position = bit_lookup(1,excitation%from_orb(1))
           bit_element = bit_lookup(2,excitation%from_orb(1))
           if (btest(connected_orbs(bit_element, excitation%to_orb(1)), bit_position)) &
                 thermal_energy(beta_index)=thermal_energy(beta_index) - &
                             (2.0*J_coupling*walker_population(1,idet))
       end if


   end subroutine call_dmqmc_estimators

end module dmqmc_estimators
