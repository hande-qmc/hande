module gutzwiller_energy    

! For the Heisenberg model.
        
use const
implicit none

contains

    subroutine plot_gutzwiller_energy()
        
        ! For testing, calculate <psi_G|H|psi_G> where |psi_G> represents the Gutzwiller
        ! wavefunction, which is a function of b. We want to find b which will
        ! minimise this expectation energy.
        
        ! psi_G(S) = exp(-b*\sum_{i,j} S_i*S_j) is the (unnormalised) amplitude of a spin
        ! configuration, with S_i and S_j denoting either +1 or -1 for up and down spins
        ! on neighbouring sites.
        
        ! By minimising this energy we can find the best guiding wavefunction of the
        ! Gutzwiller form for the Heisenberg model.
        
        ! This isn't efficient but should atleast work. For a 4 by 4 lattice this shouldn't
        ! take long to run. You can just run it for a 4 by 4 lattice and get an idea of
        ! what value of b is needed.
        
        use basis, only: bit_lookup, basis_lookup, basis_length
        use bit_utils, only: count_set_bits
        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use hamiltonian, only: diagonal_element_heisenberg
        use hubbard_real, only: connected_orbs
        use system, only: nel, nsites, J_coupling
        use utils, only: binom_r
        
        integer :: ierr
        integer :: i, j, k, ipos, iel, basis_find, big_loop, counter
        integer :: spin_sum_1, spin_sum_2, iel_ket, ipos_ket, iel_bra, ipos_bra
        integer(i0) :: f0_not(basis_length), bra(basis_length), ket(basis_length)
        integer(i0) :: bonds_0_1(basis_length), ket_not(basis_length)
        integer(dp) :: hilbert_space_size
        real(dp) :: exp_energy, normalisation
        real(p) :: b_parameter
        integer, allocatable :: configuration(:)
        
        ! Print header for output data.
        print *, 'Expectation energy of the Gutwiller wavefunction against b:'
        print *, ''
        write (6,'(3X,a12,13X,a15)') 'Parameter: b', '<psi_G|H|psi_G>'
        
        ! Find hilbert space size and set initial values.
        hilbert_space_size = binom_r(nsites, nel)
        exp_energy = 0
        normalisation = 0
        
        ! Configuration will store the positions of the up spins on the lattice.
        ! Then we will be able to vary the components to find all possible spin
        ! configuartions, and sum over this.
        allocate(configuration(nsites/2), stat = ierr)
        call check_allocate('configuration',nsites/2,ierr)
        
        ! Initial configuration
        do i=1,nsites/2
            configuration(i) = i
        end do
        
        ! Initial value of the parameter b.
        b_parameter = -5.0
        
        ! Do the entire calculation several times for different values of b,
        ! from b = -5 to +5 in steps of 0.1.
        do big_loop = 0,100
        
            ! Loop over every configuration, and over every connected configuration.
            ! Find the diagonal matrix elements for each of these configurations,
            ! and use these to compute the total expectation value.
               
            ! Loop over every possible spin configuration.
            do i=1,hilbert_space_size
                ! Create the ket corresponding to this current configuration.
                ket = 0
                do j=1,nsites/2
                    ipos = bit_lookup(1,configuration(j))
                    iel = bit_lookup(2,configuration(j))
                    ket(iel) = ibset(ket(iel),ipos)
                end do
                
                ! Find \sum_{i,j} S_i S_j for this configuration, where
                ! S is +1 for an up spin, or -1 for a down spin.
                ! I just reuse the diagonal matrix element calculation
                ! for simplicity, since it contains \sum_{i,j} S_i S_j.
                spin_sum_1 = diagonal_element_heisenberg(ket)/J_coupling
            
                ! Now need to find connected (non-identical) determinants
                ket_not = not(ket)
                do iel_ket = 1, basis_length
                    do ipos_ket = 0, i0_end
                        if (btest(ket(iel_ket), ipos_ket)) then
                            basis_find = basis_lookup(ipos_ket, iel_ket)
                            bonds_0_1 = iand(ket_not, connected_orbs(:,basis_find))
                            ! Once we have found all the 0-1 bonds on the lattice,
                            ! we next want to create these connected basis functions.
                            counter = sum(count_set_bits(bonds_0_1))
                            ! Look for the connected sites so that we can flip them to
                            ! create the connected basis functions.
                            ! This is quite a lot of do loops now...!
                            find: do iel_bra = 1, basis_length
                                do ipos_bra = 0, i0_end
                                    if (btest(bonds_0_1(iel_bra), ipos_bra)) then
                                        bra = ibclr(ket(iel_ket), ipos_ket)
                                        bra = ibset(bra(iel_bra), ipos_bra)
                                        spin_sum_2 = diagonal_element_heisenberg(bra)/J_coupling
                                        exp_energy = exp_energy - 2.0_dp*J_coupling*&
                                              exp(b_parameter*(spin_sum_1+spin_sum_2))
                                        counter = counter - 1
                                        if (counter == 0) then
                                            ! Once we have all connected basis functions
                                            ! included in the sum, quit.
                                            exit find
                                        end if
                                    end if
                                end do
                            end do find
                            ! If we haven't found all connected basis functions, something
                            ! has gone wrong somewhere, so stop.
                            if (counter /= 0) call stop_all('init_fciqmc','Not all &
                                                               &connected sites found')
                        end if
                    end do
                end do
                ! Now need to add the term due to the diagonal element:
                exp_energy = exp_energy + diagonal_element_heisenberg(ket)*&
                                       (exp(b_parameter*spin_sum_1))**2
                normalisation = normalisation + (exp(b_parameter*spin_sum_1))**2
                ! This has given all the contributions from a single basis function
                
                ! Now repeat for the others!
            
                ! Find the spin up positions for the new basis function, ready
                ! for the next cycle. Basically we start with the state where
                ! all spins are up on the first nsites/2 sites, and then cycle
                ! through all other configurations. First we move the last up spin
                ! along through the sites. Then we move the second to last spin
                ! to next position and move the last spin through all remainig sites.
                ! Then we repeat so that the second to last spin has visited all
                ! sites, and then similarly with all other sites, etc... 
                find_basis: do j = nsites/2,1,-1
                    if (configuration(j) /= (nsites/2)+j) then
                        configuration(j) = configuration(j) + 1
                        do k = j+1,nsites/2
                            configuration(k) = configuration(j)+k-j
                        end do
                        exit find_basis
                    end if
                end do find_basis
            end do
        
            ! Print out the results for this b value.
            print *, b_parameter, exp_energy/normalisation
            b_parameter = b_parameter+0.1
            ! Reset parameters.
            exp_energy = 0
            normalisation = 0
            do i=1,nsites/2
                configuration(i) = i
            end do
        
        end do
        
        ! We are done. We can deallocate the array and quit.
        
        print *, ''
        
         if (allocated(configuration)) then
            deallocate(configuration, stat=ierr)
            call check_deallocate('configuration',ierr)
        end if
        
    end subroutine plot_gutzwiller_energy
    
end module gutzwiller_energy
