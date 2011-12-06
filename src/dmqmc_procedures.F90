module dmqmc_procedures

use const
implicit none

contains

    subroutine random_distribution_entire_space()

        ! Distribute the initial number of psips along the main diagonal.
        ! Each diagonal element should be chosen with the same probability.
        ! Note that this does not only create psips within a particular
        ! subspace, but can create psips anywhere in the entire Hilbert space.

        ! To pick each configuration with equal probability, we first choose
        ! a number of bits to be set, i, out of a total of N = nbasis bits. Then
        ! we set choose i bits to set, each with equal probability. So the
        ! probability that a particular configuration is chosen is
 
        ! prob(config) = prob(config with i bits set) * prob(config | config with i bits set)
        !              = fraction of configs with i bits set * (1/number of configs with i bits set)
        !              = (number of configs with i bits set/total number of configs) 
        !                            * (1/number of configs with i bits set)
        ! => prob(config) =  1/total number of configs

        ! So each configuartion will be chosen with equal probability, as desired.

        use basis, only: nbasis, basis_length, bit_lookup, basis_lookup
        use dSFMT_interface, only:  genrand_real2
        use fciqmc_data, only: dmqmc_npsips
        use parallel
        use utils, only: binom_r
      
        integer(i0) :: total_hilbert_space, subspace_size
        integer :: i, up_spins, rand_basis, bits_set, total_bits_set
        integer :: bit_element, bit_position, npsips, basis_find, ipos
        integer(i0) :: f(basis_length)
        real(dp) prob_of_acceptance
        real :: rand_num
  
        total_hilbert_space = 2**(nbasis)
        npsips = int(dmqmc_npsips/nprocs)

        do i = 1, npsips

            ! First we need to decide how many bits will be set, total_bits_set, with
            ! probability equal to prob = fraction of configs with i bits set
            choose_bits_set: do
                rand_num = genrand_real2()
                total_bits_set = int(rand_num*(nbasis+1))
                subspace_size = binom_r(nbasis, total_bits_set)
                prob_of_acceptance = subspace_size/total_hilbert_space
                rand_num = genrand_real2()
                if (prob_of_acceptance < rand_num) exit choose_bits_set
            end do choose_bits_set

            ! Now form the initial bitstring, f, which will be modified to create the
            ! final desired bitstring
            if (total_bits_set <= int(nbasis/2)) then
                ! If less than half bits are to be set, start from all down.
                f = 0
                bits_set = 0
            else
                ! If more than half bits are to be set, start from all up,
                ! remembering to set the bits on the end down, which don't
                ! correspond to an orbital, so must always be 0.
                f = 0
                do ipos = 0, i0_end
                    basis_find = basis_lookup(ipos, basis_length)
                    if (basis_find == 0) then
                        f(basis_length) = ibset(f(basis_length),ipos)
                    end if
                end do
                f = not(f)
                bits_set = nbasis
            end if               

            ! Now create the bitstring.
            loop1: do
                ! If the correct number of bits are set, we have the
                ! bitstring that we want, so exit the loop.
                if (bits_set==total_bits_set) exit loop1
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
            end do loop1

            ! Now call a routine to add the corresponding diagonal element to
            ! the spawned walkers list.
            call create_diagonal_particle(f)

        end do

    end subroutine random_distribution_entire_space

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
                ! If half the spins are now flipped up, we have our basis
                ! function fully created, so exit the loop. 
                if (bits_set==up_spins) exit loop1
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

    subroutine decode_dm_bitstring(f, index1, index2)

        use basis, only: total_basis_length, basis_length
        use fciqmc_data, only: subsystem_A_bit_positions, subsystem_A_bit_positions
        use fciqmc_data, only: subsystem_A_size

        integer(i0), intent(in) :: f(total_basis_length)
        integer(i0), intent(out) :: index1, index2
        integer :: i

        index1 = 0
        index2 = 0
        
        do i = 1, subsystem_A_size
            if (btest(f(subsystem_A_bit_positions(i,2)),subsystem_A_bit_positions(i,1))) &
                index1 = ibset(index1,i-1)
            if (btest(f(subsystem_A_bit_positions(i,2)+basis_length),subsystem_A_bit_positions(i,1))) &
                index2 = ibset(index2,i-1)
        end do

        index1 = index1+1
        index2 = index2+1

    end subroutine decode_dm_bitstring

end module dmqmc_procedures
