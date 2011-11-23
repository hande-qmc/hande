module dmqmc_procedures

use const
implicit none

contains

    subroutine random_distribution_heisenberg()

        ! For the Heisenberg model only. Distribute the initial number of psips
        ! along the main diagonal. Each diagonal element should be chosen
        ! with the same probability.

        ! Currently this creates psips with Ms = ms_in only.

        ! If we have number of sites = nsites,
        ! and total spin value = ms_in,
        ! then number of up spins is equal to up_spins = (ms_in + nsites)/2.

        ! Start from state with all spins down, then choose the above number of
        ! spins to flip up with equal probability.

        use basis, only: nbasis, basis_length, bit_lookup
        use calc, only: ms_in
        use dSFMT_interface, only:  genrand_real2
        use fciqmc_data, only: dmqmc_npsips
        use parallel
        use system, only: nsites

        integer :: i, up_spins, rand_basis, bits_set
        integer :: bit_element, bit_position, npsips
        integer(i0) :: f(basis_length)
        real :: rand_num

        up_spins = (ms_in+nsites)/2
        npsips = int(dmqmc_npsips/nprocs)
        
        do i = 1, npsips
          
            ! Start with all spins down.
            f = 0
            bits_set = 0

            loop1: do
                ! Choose a random spin to flip.
                rand_num = genrand_real2()
                rand_basis = ceiling(rand_num*nbasis)
                ! Find the corresponding positions for this spin.
                bit_position = bit_lookup(1,rand_basis)
                bit_element = bit_lookup(2,rand_basis)
                if (btest(f(bit_element),bit_position)) then
                    ! If already flipped up, flip back down.
                    f(bit_element) = ibclr(f(bit_element),bit_position)
                    bits_set = bits_set - 1
                else
                    ! If not flipped up, flip the spin up.
                    f(bit_element) = ibset(f(bit_element),bit_position)
                    bits_set = bits_set + 1
                end if
                ! If half the spins are now flipped up, we have our basis
                ! function fully created, so exit the loop.
                if (bits_set==up_spins) exit loop1
            end do loop1

        ! Now call a routine to add the corresponding diagonal element to
        ! the spawned walkers list.
        call create_diagonal_particle(f)

        end do

    end subroutine random_distribution_heisenberg

    subroutine create_diagonal_particle(f_new)

        ! Create a psip on a diagonal element of the density
        ! matrix by adding it to the spawned walkers list. This
        ! list can then be sorted correctly by the direct_annihilation
        ! routine

        ! In:
        !    f_new: Bit string representation of index of the diagonal
        !           element upon which a new psip shall be placed.

        use hashing
        use basis, only: basis_length
        use fciqmc_data, only: spawned_walkers, spawning_head, spawned_pop
        use parallel

        integer(i0), intent(in) :: f_new(basis_length)
        integer(i0) :: f_new_diagonal(basis_length*2)
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

        ! Create the bitstring of a psip on a diagonal element.
        f_new_diagonal = 0
        f_new_diagonal(:basis_length) = f_new
        f_new_diagonal((basis_length+1):(2*basis_length)) = f_new

#ifdef PARALLEL
        ! Need to determine which processor the spawned walker should be sent to.
        iproc_spawn = modulo(murmurhash_bit_string(f_new_diagonal, &
                                (2*basis_length)), nprocs)
#endif

        ! Move to the next position in the spawning array.
        spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawned_walkers(:,spawning_head(iproc_spawn)) = 0
        ! indices 1 to total_basis_length store the bitstring.
        spawned_walkers(:(2*basis_length),spawning_head(iproc_spawn)) = f_new_diagonal
        ! The final index stores the number of psips created, always 1 in this situation.
        spawned_walkers((2*basis_length)+1,spawning_head(iproc_spawn)) = 1

    end subroutine create_diagonal_particle

    subroutine initial_dmqmc_status()

        ! Calculate the projected energy based upon the initial walker
        ! distribution (either via a restart or as set during initialisation)
        ! and print out.

        use parallel
        use proc_pointers, only: update_dmqmc_energy_ptr
        use fciqmc_data, only: trace, thermal_energy, nparticles, mc_cycles_done, shift
        use fciqmc_data, only: tot_walkers
        integer :: idet
        integer :: ntot_particles
        integer(i0) :: combined_trace
        real(p) :: combined_energy
#ifdef PARALLEL
        integer :: ierr
#endif

        ! Calculate the projected energy based upon the initial walker
        ! distribution.  proj_energy and D0_population are both accumulated in
        ! update_proj_energy.
        trace = 0
        thermal_energy = 0
        do idet = 1, tot_walkers 
            call update_dmqmc_energy_ptr(idet,1)
        end do 

#ifdef PARALLEL
        call mpi_allreduce(trace(1), combined_trace, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(thermal_energy(1), combined_energy, 1, mpi_preal, MPI_SUM, MPI_COMM_WORLD, ierr)
        call mpi_allreduce(nparticles, ntot_particles, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        combined_trace = trace(1)
        combined_energy = thermal_energy(1)
        ntot_particles = nparticles(1)
#endif 
        
        if (parent) then
            ! See also the format used in write_fciqmc_report if this is changed.
            ! We prepend a # to make it easy to skip this point when do data
            ! analysis. 
            write (6,'(1X,"#",3X,i8,2(2X,es17.10),4X,i8,4x,es17.10,4X,i11,6X,a3,3X,a3)') &
                                             mc_cycles_done, shift,&
                                             0.0_p, combined_trace, combined_energy,&
                                             ntot_particles, 'n/a', 'n/a'
        end if

    end subroutine initial_dmqmc_status

end module dmqmc_procedures
