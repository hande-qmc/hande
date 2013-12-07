module load_balancing

! Module for performing load balancing

use loadbal_data
implicit none

contains

    subroutine do_load_balancing(proc_map)

        ! Main subroutine in module, carries out load balancing as follows:
        ! 1. Work out if we want to perform load balancing. 
        ! 2. If we do then:
        !   * Find processors which have above/below average population
        !   * For processors with above average populations enter slots from slot_list
        !   * into smaller array which is then ranked according to population.
        !   * Attempt to move these slots to processors with below average population.
        ! 3. Once proc_map is modified so that its entries contain the new locations
        !   of of donor slots, we then add these determinants to spawned walker list so
        !   that they can be moved to their new processor.
        ! 4. Set load_tag to be other than one to prevent another call this report loop

        ! In/Out:
        ! proc_map: array which maps determinants to processors
        !       proc_map(modulo(hash(d),load_balancing_slots*nprocs)=processor
        
        use parallel
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        use spawn_data, only: spawn_t
        use fciqmc_data, only: qmc_spawn, walker_dets, walker_population, tot_walkers, &
                               nparticles, load_balancing_slots, nparticles_proc, sampling_size,  &
                               load_tag, pop_av, max_tries
        use ccmc, only: redistribute_excips
        use ranking, only: insertion_rank_int

        integer, intent(inout) :: proc_map(:)

        integer(lint) :: slot_pop(0:size(proc_map)-1)
        integer(lint) :: slot_list(0:size(proc_map)-1) 
        integer(lint) :: procs_pop(0:nprocs-1)

        integer, allocatable :: donors(:), receivers(:)
        integer, allocatable :: d_rank(:), d_index(:)
        integer(lint), allocatable :: d_map(:)
        
        integer :: ierr, i
        integer(lint) :: p_av
        integer ::  d_siz, r_siz, d_map_size
        integer :: up_thresh, low_thresh
        real(dp) :: perc_diff=0.01
        integer ::  imbal=0

        slot_list=0
        ! Average population across processors.
        
        ! Find slot populations.
        call initialise_slot_pop(proc_map, load_balancing_slots, qmc_spawn, slot_pop)
#ifdef PARALLEL
        ! Gather these from every process into slot_list.
        call MPI_AllReduce(slot_pop, slot_list, p_map_size, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        ! Find population per processor, store in procs_pop.
        procs_pop(:nprocs-1)=nparticles_proc(1,:nprocs)

        ! Upper threshold.
        up_thresh=pop_av+int(real(pop_av*perc_diff))
        ! Lower threshold.
        low_thresh=pop_av-int(real(pop_av*perc_diff))

        ! Find donor/receiver processors.
        call find_processors(procs_pop, up_thresh, low_thresh, receivers, donors, d_map_size, proc_map)
        ! Number of processors we can donate from.
        d_siz=size(donors)
        ! Number of processors we can donate to.
        r_siz=size(receivers)
        ! Smaller list of donor slot populations.
        allocate(d_map(d_map_size))
        ! Contains ranked version of d_map.
        allocate(d_rank(d_map_size))
        ! Contains index in proc_map of donor slots.
        allocate(d_index(d_map_size))

        ! Put donor slots into array so we can sort them. 
        call reduce_slots(donors, slot_list, d_index, d_map, proc_map)
        ! Rank d_map.
        call insertion_rank_int(d_map, d_rank, 0)

        ! Attempt to modify proc map to get more even population distribution.
        call redistribute_slots(proc_map, d_map, d_index, d_rank, procs_pop, donors, receivers, up_thresh, low_thresh)

        ! Send slots of determinants to their new processor.
        do i=1, d_siz
            if(iproc==donors(i)) then
                call redistribute_excips(walker_dets, walker_population, tot_walkers, nparticles, qmc_spawn)
            end if 
        end do
        max_tries=max_tries+1
        load_tag=load_tag+1

    end subroutine do_load_balancing
    
    subroutine redistribute_particles(walker_dets, walker_populations, tot_walkers, nparticles, spawn)

        ! Due to the cooperative spawning (ie from multiple excitors at once) in
        ! CCMC, we need to give each excitor the chance to be on the same
        ! processor with all combinations of excitors, unlike in FCIQMC where
        ! the spawning events are independent.  We satisfy this by periodically
        ! moving an excitor to a different processor (MPI rank). 

        ! WARNING: if the number of processors is large or the system small,
        ! this introduces a bias as load balancing prevents all possible
        ! clusters from being on the same processor at the same time.
 
        ! In:
        !    walker_dets: list of occupied excitors on the current processor.
        !    total_walkers: number of occupied excitors on the current processor.
        ! In/Out:
        !    nparticles: number of excips on the current processor.
        !    walker_populations: Population on occupied excitors.  On output the
        !        populations of excitors which are sent to other processors are
        !        set to zero.
        !    spawn: spawn_t object.  On output particles which need to be sent
        !        to another processor have been added to the correct position in
        !        the spawned store.

        use basis, only: basis_length
        use const, only: i0, lint
        use spawn_data, only: spawn_t
        use spawning, only: assign_particle_processor, add_spawned_particles
        use parallel, only: iproc, nprocs

        integer(i0), intent(in) :: walker_dets(:,:)
        integer, intent(inout) :: walker_populations(:,:)
        integer, intent(inout) :: tot_walkers
        integer(lint), intent(inout) :: nparticles(:)
        type(spawn_t), intent(inout) :: spawn

        integer :: iexcitor, pproc, slot=0

        !$omp parallel do default(none) &
        !$omp shared(tot_walkers, walker_dets, walker_populations, basis_length, spawn, iproc, nprocs) &
        !$omp private(pproc)
        do iexcitor = 1, tot_walkers
            !  - set hash_shift and move_freq
            call assign_particle_processor(walker_dets(:,iexcitor), basis_length, spawn%hash_seed, &
                                              spawn%hash_shift, spawn%move_freq, nprocs, pproc, slot)
            if (pproc /= iproc) then
                ! Need to move.
                ! Add to spawned array so it will be sent to the correct
                ! processor during annihilation.
                ! NOTE: for initiator calculations we need to keep this
                ! population no matter what.  This relies upon the
                ! (undocumented) 'feature' that a flag of 0 indicates the parent
                ! was an initiator...
                call add_spawned_particles(walker_dets(:,iexcitor), walker_populations(:,iexcitor), pproc, spawn)
                ! Update population on the sending processor.
                nparticles = nparticles - abs(walker_populations(:,iexcitor))
                ! Zero population here.  Will be pruned on this determinant
                ! automatically during annihilation (which will also update tot_walkers).
                walker_populations(:,iexcitor) = 0
            end if
        end do
        !$omp end parallel do

    end subroutine redistribute_particles
    
    pure function check_imbalance() result(imbal) 

        use parallel
        use fciqmc_data, only: nparticles_proc, pop_av

        integer :: i, imbal
        integer(lint) :: max_i
        integer (lint) ::  procs_pop(nprocs)

        procs_pop(:nprocs)=nparticles_proc(1,:nprocs)

        max_i=pop_av+5000
        imbal=0 
        
        do i=1, nprocs
            if(procs_pop(i)>max_i) then
                imbal=1
                exit
            end if 
        end do

    end function check_imbalance

    subroutine redistribute_slots(proc_map, d_map, d_index, d_rank, procs_pop, donors, receivers, up_thresh, low_thresh)

        ! Attempt to modify entries in proc_map to get a more even population distribution across processors.
        ! Slots from d_map are currently donated in increasing slot population.
        ! This is carried out while the donor processor's population is above a specified threshold
        ! or the receiver processor's population is below a certain threshold.

        ! In/Out:
        !   proc_map: array which maps determinants to processors.
        !       proc_map(modulo(hash(d),load_balancing_slots*nprocs)=processor.
        ! In: 
        !   d_map: array containing populations of donor slots which we try and redistribute/.
        !   d_index: array containing index of entries in d_map in proc_map.
        !   d_rank: array containing indices of d_map ranked in increasing population.
        !   procs_pop: array containing populations on each processor.
        !   donors/receivers: array containing donor/receiver processors
        !       (ones with above/below average population).
        !   up_thresh: Upper population threshold for load imbalance.
        !   low_ thresh: lower population threshold for load imbalance.

        use parallel, only: nprocs
        use fciqmc_data, only: load_balancing_slots

        integer(lint), intent(in) :: d_map(:)
        integer, intent(in) ::  d_index(:), d_rank(:)
        integer, intent(in) :: donors(:), receivers(:)
        integer, intent(in) :: up_thresh, low_thresh
        integer, intent(inout) :: proc_map(0:p_map_size-1)
        integer(lint), intent(inout) :: procs_pop(0:nprocs-1)

        integer :: pos
        integer :: i, j, total, donor_pop, new_pop 

        donor_pop=0
        new_pop=0

        do i=1, size(d_map)
            ! Loop over receivers.
            pos=d_rank(i)
            do j=1, size(receivers)
                ! Try to add this to below average population.
                new_pop=d_map(pos) + procs_pop(receivers(j))             
                ! Modify donor population. 
                donor_pop=procs_pop(proc_map(d_index(pos)))-d_map(pos)
                ! If adding subtracting slot doesn't move processor pop past a bound.
                if (new_pop .le. up_thresh .and. donor_pop .ge. low_thresh ) then
                    ! Changing processor population
                    procs_pop(proc_map(d_index(pos)))=donor_pop
                    procs_pop(receivers(j))=new_pop
                    ! Updating proc_map.
                    proc_map(d_index(pos))=receivers(j)
                    ! Leave the j loop, could be more than one receiver.
                    exit
	            end if
	        end do
	    end do

    end subroutine redistribute_slots

    subroutine reduce_slots(donors, slot_list, d_index, d_map, proc_map)

        ! Reduce the size of array we have to search when finding large/small slots to redistribute.

        ! In: 
        !   donors: array containing donor processors.
        !   slot_list: array containing populations of slots across all processors.
        ! In/Out:
        !   d_index: array containing index of entries in d_map in proc_map.
        !   d_map: array containing populations of donor slots which we try and redistribute
        !   proc_map: array which maps determinants to processors.
        !       proc_map(modulo(hash(d),load_balancing_slots*nprocs)=processor.

        use parallel, only: iproc

        integer, intent (in) :: donors(:)
        integer(lint), intent(in) :: slot_list(0:p_map_size-1) 
        integer, intent (in) :: proc_map(0:p_map_size-1)
        integer(lint), intent(inout) :: d_map(:)
        integer, intent(inout) :: d_index(:)

        integer :: i, j, k

        k=1
        
        do i=1, size(donors)
            do j=0, size(slot_list)-1
                ! Putting appropriate blocks of slots in d_map.
                if(proc_map(j)==donors(i)) then
                    d_map(k)=slot_list(j)
                    ! Index is important as well.
                    d_index(k)=j
                    k=k+1
                end if 
            end do
        end do

    end subroutine reduce_slots
      
    subroutine find_processors(procs_pop, up_thresh, low_thresh, rec_dummy, don_dummy, donor_slots, proc_map)
        
        ! Find donor/receiver processors.
        ! Put these into varying size array receivers/donors.

        ! In:
        !   procs_pop: number particles on each processor.
        !   upper/lower_thresh: upper/lower thresholds for load imblance i.e. how close to the average population
        !       we aspire to.
        ! In/Out:
        !   rec_dummy/don_dummy: arrays which contain donor/receivers processors.
        !   donor_slots: number of slots which we can donate, this varies as more entries in proc_map are 
        !       modified.
        !   proc_map: array which maps determinants to processors.
        !       proc_map(modulo(hash(d),load_balancing_slots*nprocs)=processor.

        use parallel, only: nprocs
        use ranking, only: insertion_rank_int

        integer(lint), intent(in) :: procs_pop(0:nprocs-1)
        integer, intent(in) :: proc_map(0:p_map_size-1)
        integer, intent(in) :: up_thresh, low_thresh
        integer, intent(inout) :: donor_slots

        integer ::  i, j, k, upper, lower
        integer, allocatable, dimension(:) ::  tmp_rec, tmp_don, rec_sort
        integer, allocatable :: rec_dummy(:), don_dummy(:)
        integer :: rank_nparticles(nprocs)

        allocate(tmp_rec(nprocs))
        allocate(tmp_don(nprocs))
        k=1
        j=1

        ! Find donor/receiver processors.
        
        do i=0, size(procs_pop)-1
            if(procs_pop(i) .lt. low_thresh) then
                tmp_rec(j)=i
                j=j+1
            else if (procs_pop(i) .gt. up_thresh) then
                tmp_don(k)=i
                k=k+1
            end if 
        end do

        ! Put processor ID into smaller array.

        allocate(rec_dummy(j-1))
        allocate(rec_sort(j-1))
        allocate(don_dummy(k-1))

        don_dummy=tmp_don(:k-1)
        rec_dummy=tmp_rec(:j-1)

        ! Sort receiver processers.

        call insertion_rank_int(procs_pop, rank_nparticles, 0) 
        do i=1, size(rec_dummy)
            rec_sort(i)=rank_nparticles(i)-1
        end do
        rec_dummy=rec_sort

        ! Calculate number of donor slots which we can move.

        donor_slots=0
        do i=0, size(proc_map)-1
            do j=1, size(don_dummy)
                if(proc_map(i)==don_dummy(j)) then
                    donor_slots=donor_slots+1
                end if 
            end do 
        end do 

    end subroutine find_processors

    subroutine particles_per_proc(proc_map, slot_list, procs_pop)

        ! Find number of particles per processor and store in array procs_pop.

        ! In:
        !   proc_map: array which maps determinants to processors.
        !       proc_map(modulo(hash(d),load_balancing_slots*nprocs)=processor.
        !   slot_list: array containing populations of slots in proc_map across all
        !       processors.
        ! Out:
        !   procs_pop(nprocs): array containing population on each processor.

        use parallel, only : nprocs

        integer, intent(in) :: proc_map(0:p_map_size-1)
        integer(lint), intent(in) :: slot_list(0:p_map_size-1) 
        integer(lint), intent(out) :: procs_pop(0:nprocs-1)

        integer :: i

        procs_pop=0  

        do i=0, size(proc_map)-1
            procs_pop(proc_map(i))=procs_pop(proc_map(i))+slot_list(i)
        end do

    end subroutine particles_per_proc

    subroutine initialise_slot_pop(proc_map, load_balancing_slots, spawn, slot_pop)

        ! In: 
        !   proc_map(p_map_size): array which maps determinants to processors.
        !       proc_map(modulo(hash(d),load_balancing_slots*nprocs)=processor
        !   load_balancing_slots: number of slots which we divide slot_pop (and similar arrays) into.
        !   spawn: spawn_t object.     
        ! In/Out:
        !   slot_pop(p_map_size): array containing population of slots in proc_map
        !       processor dependendent.

        use parallel, only: nprocs, iproc
        use basis, only: basis_length
        use fciqmc_data, only: tot_walkers, walker_dets, walker_population
        use spawning, only: assign_particle_processor
        use spawn_data, only: spawn_t

        integer, intent(in):: load_balancing_slots
        type(spawn_t), intent(in) :: spawn
        integer, intent(in) :: proc_map(0:p_map_size-1)
        integer(lint), intent(out) :: slot_pop(0:size(proc_map)-1) 

        integer :: i, det_pos, iproc_slot

        slot_pop=0
        do i=1, tot_walkers
            call assign_particle_processor(walker_dets(:,i), basis_length, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, &
                                           nprocs, iproc_slot, det_pos)
            slot_pop(det_pos)=slot_pop(det_pos)+abs(walker_population(1,i))
        end do

   end subroutine initialise_slot_pop

end module load_balancing
