module tc_fciqmc

use fciqmc_data
    
implicit none

    
    ! Array to hold the spawned times of the walkers
real(p), allocatable :: spawn_times(:) ! (spawned_walker_length)

contains
    subroutine do_tc_fciqmc()
        
        integer :: nspawned, ireport, idet, iparticle
        integer, allocatable :: current_pos(:) ! (0:max(1,nprocs-1))
        real(p) :: time, t_barrier, t1, t2

        allocate(current_pos(0:max(1,nprocs-1)))

        ! time the report loop
        call cpu_time(t1)

        do ireport = 1, nreports
    
            ! Zero cycle quantities
            rspawn = 0.0_p
            proj_energy = 0.0_p
            D0_population = 0.0_p
            
            ! Reset the pointer to the current position in the spawning array to 
            ! be the slot preceding the first
            spawning_head = spawning_block_start

            do idet = 1, tot_walkers ! loop over determinants in the walker list
            
                ! get the determinant bitstring once so we do not need to keep
                ! doing it. Then find list of occupied orbitals
                cdet%f = walker_dets(:,idet) 
                call decoder(cdet%f, cdet)
                R = calc_R(cdet)
                tmp_pop = walker_population(1,idet) ! not sure if this should be a "2" here?

                !evaluate the projected energy
                call update_proj_energy(idet)
                 
                ! loop over each walker on the determinant
                do iparticle = 1, abs(walker_population(1,idet))

                    time = 0
                    do
                        time = time + timestep(R)
                        if ( time > t_barrier ) exit

                        call ct_spawn(cdet, walker_population(1,idet), nspawned, connections)
                        
                        ! If death then kill the walker immediately and move
                        ! onto the next one
                        if (connections%nexcit == 0) then
                            !if the spawned walker and the parent (all the
                            !walkers on a perticular det. have the same sgn due
                            !to annihilation) are of opposite sgn we get death
                            if(walker_population(1,idet)*nspawned < 0) then
                                if (sgn(nspawned) == 1) then
                                    tmp_pop = tmp_pop + 1 ! if nspawned +ve then add it
                                else
                                    tmp_pop = tmp_pop - 1 ! if nspawned -ve then subtract it
                                end if
                                exit ! the walker is dead
                            end if
                        end if

                        ! If there were some walkers spawned, append them to the
                        ! spawned array - maintaining processor blocks if going in
                        ! parallel. We now also have an extra "time" array giving
                        ! the birth time of the walker
                        if (nspawned /= 0) call create_spawned_particle_ct(cdet, connections, time, nspawned, spawned_pop)
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

                    if (current_pos(iproc) /= spawning_head(iproc) + 1) then

                        ! decode the spawned walker bitstring
                        cdet%f = spawned_walkers(:basis_length,iproc)
                        call decoder(cdet%f,cdet)
                        R = calc_R(cdet)

                        ! Spawn from this walker & append to the spawned array until
                        ! we hit the barrier
                        time = spawn_times(current_pos(iproc))
                        do

                            time = time - timestep(R)
                            if ( time > t_barrier ) exit

                            call ct_spawn(cdet, spawned_walkers(spawned_pop,idet), nspawned, connections)
                            
                            ! add a walker to the end of the spawned walker list in the
                            ! appropriate block - this will increment the appropriate
                            ! spawning heads for the processors which were spawned on
                            call create_spawned_particle_ct(cdet, connections, time, nspawned, spawned_pop)
                        end do

                        ! go on to the next element
                        current_pos(iproc) = current_pos(iproc) + 1

                    end if

                end do

                if(all(current_pos == spawning_head+1)) exit
                
            end do
        end do

    end subroutine do_tc_fciqmc

    
    subroutine ct_spawn(d, parent_sgn, nspawned, connection)
    
        ! randomly select a (valid) excitation from the current determinant
        ! stored in "d"
        !
        ! In: 
        !    d: info on current determinant that we will spawn from.
        !    parent_sgn: sgn on the parent determinant (i.e. +ve or -ve integer)
        !
        ! Out:
        !    nspawn: +/- 1 as @ the end of each time "jump" we only spawn
        !            1 walker.
        !    connection: the excitation connection between the parent and child
        !                determinants

    end subroutine ct_spawn
    

    subroutine create_spawned_particle_ct(d, connection, time, nspawned, spawned_pop)

        ! Create a spawned walker in the spawned walkers lists.
        ! The spawning_head is updated
        !
        ! In:
        !    d: info on the current determinant (which we have spawned from).
        !    connection: excitation connecting the current determinant to its
        !        progeny.
        !    time: the time (from the last barrier) @ which the walker was
        !        spawned. This is needed so that children can themselves be
        !        advanced and allowed to spawn chidlren of their own but with
        !        the correct probability before they hit the next barrier.
        !    nspawned: the signed number of particles to create on the spawned
        !        determinant.
        !    particle_type: the type of particle created ( i.e. spawned_pop is
        !        passed for Hamiltonian particles and spawned_hf_pop for
        !        Hellman--Feynman particles).



    end subroutine create_spawned_particle

    
    pure function timestep(R)

        ! Returns a random timestep to advance a walker by for the continuous
        ! time algorithm.
        use dSFMT_interface, only: genrand_real2

        real(p), intent(out) :: timestep
        real(p), intent(in)  :: R

        timestep = -R*log(genrand_real2())
    
    end function timestep


    pure function calc_R(d) result(R)

        ! Returns \sum_i R_ij where R_ij is the unsigned matrix element
        ! connecting |D_i> to |D_j>. Used for selecting a time to jump to and
        ! also which excitation to choose when spawning.

    end function calc_R

end module tc_fciqmc
