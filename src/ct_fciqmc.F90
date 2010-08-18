module ct_fciqmc

use fciqmc_data
use const, only: p
    
implicit none

contains
    subroutine do_ct_fciqmc(decoder, update_proj_energy, enumerator, sc0, matel)

        use annihilation, only: direct_annihilation
        use basis, only: basis_length
        use determinants, only: det_info, alloc_det_info
        use energy_evaluation, only: update_energy_estimators
        use excitations, only: excit
        use fciqmc_common, only: load_balancing_report, initial_fciqmc_status
        use fciqmc_restart
        use interact
        use system, only: ndim, nsites, nalpha, nbeta, system_type, hub_k, hub_real
        use parallel

        interface  
            subroutine enumerator(det, nexcits, connection_list)
                use determinants, only: det_info
                use const, only: p
                use excitations, only: excit,&
                enumerate_all_excitations_hub_real,&
                enumerate_all_excitations_hub_k 
                implicit none
                type(det_info), intent(in) :: det
                integer, intent(out) :: nexcits
                type(excit), intent(out) :: connection_list(:)
            end subroutine enumerator
            subroutine decoder(f,d)
                use basis, only: basis_length
                use const, only: i0
                use determinants, only: det_info
                implicit none
                integer(i0), intent(in) :: f(basis_length)
                type(det_info), intent(inout) :: d
            end subroutine decoder
            subroutine update_proj_energy(idet)
                use const, only: p
                implicit none
                integer, intent(in) :: idet
            end subroutine update_proj_energy
            function sc0(f) result(hmatel)
                use basis, only: basis_length
                use const, only: i0, p
                implicit none
                real(p) :: hmatel
                integer(i0), intent(in) :: f(basis_length)
            end function sc0
        end interface

        real(p), intent(in) :: matel ! either U or t, depending whether we are working in the real or k-space

        integer :: nspawned, nexcitations, tot_spawned, nparticles_old(sampling_size), ireport, idet, iparticle, tmp_pop
        integer, allocatable :: current_pos(:) ! (0:max(1,nprocs-1))
        real(p) :: time, t_barrier, K_ii, R
        real :: t1, t2
        type(det_info) :: cdet
        type(excit) :: connections
        type(excit), allocatable :: connection_list(:)
        logical :: soft_exit

        if (system_type == hub_k) then
            allocate(connection_list(nalpha*nbeta*min(nsites-nalpha,nsites-nbeta)))  
        else if (system_type == hub_real) then
            allocate(connection_list(2*ndim*nel))
        end if

        allocate(current_pos(0:max(1,nprocs-1)))
        call alloc_det_info(cdet)

        nparticles_old = nparticles_old_restart

        t_barrier = tau ! or we could just not bother with the t_barrier var...

        if (parent) call write_fciqmc_report_header()
        call initial_fciqmc_status(update_proj_energy)

        ! Main fciqmc loop
        do ireport = 1, nreport
    
            ! time the report loop
            call cpu_time(t1)
            
            ! Zero cycle quantities
            rspawn = 0.0_p
            proj_energy = 0.0_p
            D0_population = 0.0_p
            tot_spawned = 0
            
            ! Reset the pointer to the current position in the spawning array to 
            ! be the slot preceding the first
            spawning_head = spawning_block_start

            do idet = 1, tot_walkers ! loop over determinants in the walker list
            
                ! get the determinant bitstring once so we do not need to keep
                ! doing it. Then find list of occupied orbitals
                cdet%f = walker_dets(:,idet) 
                call decoder(cdet%f, cdet)
                call enumerator(cdet,nexcitations,connection_list)
                tmp_pop = walker_population(1,idet)

                !evaluate the projected energy
                call update_proj_energy(idet)
                 
                ! loop over each walker on the determinant
                do iparticle = 1, abs(walker_population(1,idet))

                    time = 0.0_p
                    do
                        time = time + timestep(abs(walker_energies(1,idet)) + nexcitations*abs(matel))
                        if ( time > t_barrier ) exit
                        
                        call ct_spawn(nexcitations, connection_list, walker_energies(1,idet), walker_population(1,idet), matel, nspawned, connections)

                        ! If death then kill the walker immediately and move
                        ! onto the next one
                        ! if the spawned walker and the parent (all the
                        ! walkers on a perticular det. have the same sgn due
                        ! to annihilation) are of opposite sgn we get death
                        if (connections%nexcit == 0 .and. &
                        walker_population(1,idet)*nspawned < 0.0_p) then
                            tmp_pop = tmp_pop + nspawned 
                            nparticles(1) = nparticles(1) - abs(nspawned) 
                            exit ! the walker is dead
                        end if
                         
                        ! used for calculating rspawn
                        tot_spawned = tot_spawned + abs(nspawned)
                        
                        ! If there were some walkers spawned, append them to the
                        ! spawned array - maintaining processor blocks if going in
                        ! parallel. We now also have an extra "time" array giving
                        ! the birth time of the walker
                        if (nspawned /= 0) call create_spawned_particle_ct(cdet, connections, nspawned, spawned_pop, time)
                    end do

                end do
                
                walker_population(1,idet) = tmp_pop
            
            end do


            ! now we advance all the spawned walkers to the barrier from their
            ! respective birth times - Any walkers spawned as a consequence of
            ! this  must be appened to the spawned array and themselves advanced
            ! to the barrier.

            !start @ the start of the block
            current_pos = spawning_block_start
            do
                do iproc = 0, nprocs-1

                    if (current_pos(iproc) /= spawning_head(iproc) + 1 .and. current_pos(iproc) /= spawning_block_start(iproc)) then

                        ! decode the spawned walker bitstring
                        cdet%f = spawned_walkers(:basis_length,current_pos(iproc))
                        K_ii = sc0(cdet%f)
                        call decoder(cdet%f,cdet)
                        call enumerator(cdet,nexcitations,connection_list)
!                        R = calc_R(cdet)

                        ! Spawn from this walker & append to the spawned array until
                        ! we hit the barrier
                        time = spawn_times(current_pos(iproc))
                        do

                            time = time + timestep(abs(sc0(cdet%f)) + nexcitations*abs(matel))
                            if ( time > t_barrier ) exit

                            call ct_spawn(nexcitations, connection_list, sc0(cdet%f), spawned_walkers(spawned_pop,current_pos(iproc)), matel, nspawned, connections)
                           
                            ! Handle walker death
                            if(connections%nexcit == 0 .and. &
                            spawned_walkers(spawned_pop,current_pos(iproc))*nspawned < 0) then
                                spawned_walkers(spawned_pop,current_pos(iproc)) = &
                                    spawned_walkers(spawned_pop,current_pos(iproc)) + nspawned 
                                exit ! the walker is dead - do not continue
                            end if

                            ! add a walker to the end of the spawned walker list in the
                            ! appropriate block - this will increment the appropriate
                            ! spawning heads for the processors which were spawned on
                            call create_spawned_particle_ct(cdet, connections, nspawned, spawned_pop, time)
                        end do

                        ! go on to the next element
                        current_pos(iproc) = current_pos(iproc) + 1

                    end if

                end do

                if(all(current_pos == spawning_head+1 .or. current_pos == spawning_block_start)) exit
                
            end do


            ! calculate rspawn
            ! JSS: don't understand why this is a problem. Will fix later.
!            rspawn = tot_spawned/nparticles_old(1)

            !update spawn rate 
            call direct_annihilation(sc0)

            !update projected energy and shift
            call update_energy_estimators(ireport, nparticles_old)

            call cpu_time(t2)

            if(parent) call write_fciqmc_report(ireport, nparticles_old(1), t2-t1)
            
            t1 = t2

            call fciqmc_interact(ireport, soft_exit)
            if(soft_exit) exit

        end do

        if (parent) then
            call write_fciqmc_final(ireport)
            write(6,'()')
        end if

        call load_balancing_report()

        if (dump_restart_file) call dump_restart(mc_cycles_done+ncycles*nreport, nparticles_old(1))
        
        deallocate(current_pos)
        deallocate(connection_list)
    end subroutine do_ct_fciqmc

    
    subroutine ct_spawn(num_excitations, connection_list, K_ii, parent_sgn, matel, nspawned, connection)
    
        ! randomly select a (valid) excitation 
        !
        ! In: 
        !    det: info on current determinant that we will spawn from.
        !    R_ii: the diagonal Hamiltonian matrix element for the determinant d
        !    parent_sgn: sgn on the parent determinant (i.e. +ve or -ve integer)
        !
        ! Out:
        !    nspawn: +/- 1 as @ the end of each time "jump" we only spawn
        !            1 walker.
        !    connection: the excitation connection between the parent and child
        !                determinants

        use excitations, only: enumerate_all_excitations_hub_real, excit
        use determinants, only: det_info
        use dSFMT_interface, only: genrand_real2
        use system, only: ndim, nel, system_type, hub_real, hub_k
        use hamiltonian, only: slater_condon1_hub_real, slater_condon2_hub_k

        integer, intent(in) :: parent_sgn, num_excitations
        real(p), intent(in) :: K_ii, matel
        type(excit), intent(in) :: connection_list(:)
        
        integer, intent(out) :: nspawned
        type(excit), intent(out) :: connection
        
        real(p) :: rand, test, R_ii, R, abs_matel, K_ij
        integer :: j

        R_ii = abs(K_ii-shift)
        abs_matel = abs(matel)
        R = R_ii + abs_matel*num_excitations
        rand = genrand_real2()*R



        if (rand < R_ii) then
            connection%nexcit = 0 ! spawn onto the same determinant (death/cloning)
            return
        else
            test = R_ii
            do j = 2, num_excitations
                test = test + abs_matel
                if (rand < test) then
                    connection = connection_list(j)
                    exit
                end if
            end do
        end if

        if (system_type == hub_k) then
            K_ij = slater_condon2_hub_k(connection%from_orb(1), connection%from_orb(2),&
                                        connection%to_orb(1), connection%to_orb(2), connection%perm)
        else if (system_type == hub_real) then
            K_ij = slater_condon1_hub_real(connection%from_orb(1), connection%to_orb(1), connection%perm)
        end if

            if (K_ij < 0.0_p) then    ! child is same sign as parent
                nspawned = -sign(1,parent_sgn)
            else if (K_ij > 0.0_p) then
                nspawned = sign(1,parent_sgn)
            else
                nspawned = 0
            end if
    end subroutine ct_spawn


!    subroutine ct_spawn_kspace(connection_list,num_excitations, K_ii, parent_sgn, matel, nspawned, connection)
!    
!        ! randomly select a (valid) excitation from the current determinant
!        ! stored in "det" for the Hubbard model in realspace
!        !
!        ! In: 
!        !    det: info on current determinant that we will spawn from.
!        !    R_ii: the diagonal Hamiltonian matrix element for the determinant d
!        !    parent_sgn: sgn on the parent determinant (i.e. +ve or -ve integer)
!        !
!        ! Out:
!        !    nspawn: +/- 1 as @ the end of each time "jump" we only spawn
!        !            1 walker.
!        !    connection: the excitation connection between the parent and child
!        !                determinants
!
!        use determinants, only: det_info
!        use excitations, only: enumerate_all_excitations_hub_k, excit
!        use dSFMT_interface, only: genrand_real2
!        use system, only: ndim, nel, nalpha, nbeta, nsites
!
!        type(det_info), intent(in) :: det
!        integer, intent(in) :: parent_sgn
!        real(p), intent(in) :: K_ii, matel
!        
!        integer, intent(out) :: nspawned
!        type(excit), intent(out) :: connection
!        
!        real(p) :: rand, test, R_ii, R, abs_matel
!        integer :: num_excitations, j
!        type(excit) :: connection_list(nalpha*nbeta*min(nsites-nalpha,nsites-nbeta))
!
!        rand = genrand_real2()*R
!
!        if (K_ii < 0) then    ! child is same sign as parent
!            nspawned = sign(1,parent_sgn)
!        else if (K_ii > 0) then
!            nspawned = -sign(1,parent_sgn)
!        else
!            nspawned = 0
!        end if
!
!        if (rand < R_ii) then
!            connection%nexcit = 0 ! spawn onto the same determinant (death/cloning)
!        else
!            test = R_ii
!            do j = 2, size(connection_list) ! cycle over connections and test for each one
!                test = test + abs_matel
!                if (rand < test) then
!                    connection = connection_list(j)
!                    exit
!                end if
!            end do
!        end if
!
!    end subroutine ct_spawn_kspace

    
    subroutine create_spawned_particle_ct(cdet, connection, nspawn, particle_type, spawn_time)

        ! Create a spawned walker in the spawned walkers lists.
        ! The current position in the spawning array is updated.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connecting the current determinant to its
        !        offspring.  Note that the perm field is not used.
        !    nspawn: the (signed) number of particles to create on the
        !        spawned determinant.
        !    particle_type: the type of particle created.  Must correspond to
        !        the desired element in the spawning array (i.e. be spawned_pop
        !        for Hamiltonian particles and spawned_hf_pop for
        !        Hellmann--Feynman particles).

        use hashing
        use parallel, only: iproc, nprocs

        use basis, only: basis_length
        use determinants, only: det_info
        use excitations, only: excit, create_excited_det
        use fciqmc_data, only: spawned_walkers, spawning_head, spawned_pop

        type(det_info), intent(in) :: cdet
        type(excit), intent(in) :: connection
        integer, intent(in) :: nspawn
        integer, intent(in) :: particle_type
        real(p), intent(in) :: spawn_time

        integer(i0) :: f_new(basis_length)
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn 
#endif

        ! Create bit string of new determinant.
        call create_excited_det(cdet%f, connection, f_new)

#ifdef PARALLEL
        ! (Extra credit for parallel calculations)
        ! Need to determine which processor the spawned walker should be sent
        ! to.  This communication is done during the annihilation process, after
        ! all spawning and death has occured..
        iproc_spawn = modulo(murmurhash_bit_string(f_new, basis_length), nprocs)
#endif

        ! Move to the next position in the spawning array.
        spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawned_walkers(:,spawning_head(iproc_spawn)) = 0
        spawned_walkers(:basis_length,spawning_head(iproc_spawn)) = f_new
        spawned_walkers(particle_type,spawning_head(iproc_spawn)) = nspawn
        spawn_times(spawning_head(iproc_spawn)) = spawn_time

    end subroutine create_spawned_particle_ct

    function timestep(R) result(dt)

        ! Returns a random timestep to advance a walker by for the continuous
        ! time algorithm.

        use dSFMT_interface, only: genrand_real2

        real(p) :: dt
        real(p), intent(in)  :: R

        dt = -R*log(genrand_real2())
    
    end function timestep


!    pure function calc_R(d,matel,K_ii) result(R)
!
!        ! Returns \sum_i R_ij where R_ij is the unsigned matrix element
!        ! connecting |D_i> to |D_j>. Used for selecting a time to jump to and
!        ! also which excitation to choose when spawning.
!
!        use determinants, only: det_info
!
!        real(p) :: R
!        type(det_info), intent(in) :: d
!        real(p), intent(in) :: matel, K_ii
!        
!        call enumerate_all_excitations_hub_k(d, num_excitations, connection_list)
!        R = abs(K_ii) + abs(matel)*num_excitations
!
!    end function calc_R

end module ct_fciqmc
