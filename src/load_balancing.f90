module load_balancing

! Module for performing load balancing
! Determinants are assigned to processors using proc_map(num_slots*nprocs) which 
! is initialised so that its entries cyclically contain processors 0,..,nprocs-1
! Population for each entry in proc_map is calculated 
! Attempt to donate slots from processors with above average populations to those
! with below average populations

    use loadbal_data
    implicit none
    ! Upper and lower population thresholds
    integer :: up_thresh, low_thresh
    ! How much of a load imbalance is acceptable
    real(dp) :: perc_diff=0.01
    
contains 
    
    subroutine redistribute_processors(proc_map)
        
        ! Main subroutine in module, carries out load balancing is relatively logical order
        ! In/Out:
        ! proc_map: array which maps determinants to processors
        !       proc_map(modulo(hash(d),num_slots*nprocs)=processor
        
        use parallel
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        use spawn_data, only: spawn_t
        use fciqmc_data, only: qmc_spawn
        
        integer ::  d_siz, r_siz, p_total, d_map_size
        integer, intent(inout) :: proc_map(p_map_size)
        integer(lint) :: slot_pop(p_map_size)
        integer(lint) :: slot_list(p_map_size)  
        integer(lint) :: nparticles_proc(nprocs)
        integer, allocatable :: donors(:), receivers(:)
        integer, allocatable :: d_rank(:), d_index(:)
        integer(lint), allocatable :: d_map(:)
        integer :: ierr,i,icycle
        slot_list=0
        nparticles_proc=0
        ! average population across processors
        p_total=0
        ! find slot populations
        call initialise_slot_pop(slot_pop, proc_map, p_map_size, num_slots) 
        ! gather these from every process into slot_list
        call MPI_AllReduce(slot_pop, slot_list,p_map_size, MPI_INTEGER , MPI_SUM, MPI_COMM_WORLD, ierr)   
        ! find population per processor, store in nparticles_proc
        call particles_per_proc(proc_map, slot_list, nparticles_proc, p_map_size)
        ! average pop across proccessors
        p_total = int(real(sum(nparticles_proc)/nprocs))
        ! upper threshold
        up_thresh=p_total+int(real(p_total*perc_diff))
        ! lower threshold
        low_thresh=p_total-int(real(p_total*perc_diff))
        ! find donor/receiver processors
        call find_processors(receivers, donors, nparticles_proc, d_map_size)
        ! donor size
        d_siz=size(donors)
        ! receiver size
        r_siz=size(receivers)
        ! smaller list of donor slot populations
        allocate(d_map(d_map_size))
        ! contains indices ranked d_map
        allocate(d_rank(d_map_size))
        ! contains index in proc_map of donor slots
        allocate(d_index(d_map_size))
        ! put donor slots into array so we can sort them 
        call reduce_slots(donors, d_map, d_index, d_map_size, slot_list, d_siz)
        ! rank d_map
        call insertion_rank_int(d_map, d_rank, 0) 
        ! find number of particles per processor
        call particles_per_proc(proc_map, slot_list, nparticles_proc, p_map_size)
        print *, "before", iproc, nparticles_proc
        ! attempt to modify proc map to get more even population distribution 
        call redistribute(proc_map,d_map,d_index,d_rank, nparticles_proc,d_map_size,donors,d_siz,receivers,r_siz)
        ! only have to call for donors
        call move_determinants(walker_dets, walker_population, tot_walkers, nparticles, qmc_spawn)
        ! this is just for counting
        call initialise_slot_pop(slot_pop, proc_map, p_map_size, num_slots)
        call particles_per_proc(proc_map, slot_list, nparticles_proc, p_map_size)
        print *,"after ", iproc, nparticles_proc

    end subroutine redistribute_processors
   
    subroutine move_determinants(walker_dets, walker_population, tot_walkers, nparticles, spawn)
        
        ! Loops through determinants, D, and checks if proc_map(mod(hash(D),nprocs*num_slots))
        ! is the current processor. Sends determinants if not. Copied entirely from ccmc.F90 
        
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
        use hashing
        use spawn_data, only: spawn_t
        use spawning, only:  add_spawned_particles_load
        use parallel, only: iproc, nprocs
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        
        integer(i0), intent(in) :: walker_dets(:,:)
        integer, intent(inout) :: walker_population(:,:)
        integer, intent(inout) :: tot_walkers
        integer(lint), intent(inout) :: nparticles(:)
        type(spawn_t), intent(inout) :: spawn
        integer :: idet, new_proc, det_pos

        do idet = 1, tot_walkers
            det_pos=modulo(murmurhash_bit_string(walker_dets(:,idet), basis_length, 7),num_slots*nprocs)+1
            new_proc=proc_map(det_pos) 
            if (new_proc /= iproc) then
                ! Send determinant to new processor                
                call add_spawned_particles_load(walker_dets(:,idet), walker_population(:,idet), new_proc, spawn)
                ! Update population on the sending processor.
                nparticles = nparticles - abs(walker_population(:,idet))
                ! Zero population here.  Will be pruned on this determinant
                ! automatically during annihilation (which will also update tot_walkers).
                walker_population(:,idet) = 0
            end if
        end do

    end subroutine move_determinants
 
    subroutine redistribute(proc_map, d_map, d_index, d_rank, nparticles_proc, d_map_size, donors, d_siz, receivers, r_siz)
        
        ! Attempt to modify entries in proc_map to get a more even population distribution across processors
        ! In: 
        !   d_map: array containing populations of donor slots which we try and redistribute
        !   d_index: array containing index of entries in d_map in proc_map
        !   d_rank: array containing indices of d_map ranked in increasing population
        !   d_map_size: length of above arrays
        !   nparticles_proc: array containing populations on each processor
        !   donors/receivers: array containing donor/receiver processors (ones with above/below average population)
        !   d_siz/r_siz: length of donors/receivers arrays
        ! In/Out:
        !   proc_map: array which maps determinants to processors
        !       proc_map(modulo(hash(d),num_slots*nprocs)=processor
        
        use parallel, only : nprocs
        
        integer :: i, j, total, donor_pop, new_pop 
        integer, intent(in) ::d_map_size, d_siz, r_siz
        integer(lint), intent(in) :: d_map(d_map_size)
        integer, intent(in) ::  d_index(d_map_size), d_rank(d_map_size)
        integer, intent(in) :: donors(d_siz), receivers(r_siz)
        integer, intent(inout) :: proc_map(p_map_size)
        integer(lint), intent(inout) :: nparticles_proc(nprocs)
        integer :: pos
        donor_pop=0
        new_pop=0
        do i=1, d_map_size
            ! loop over receivers
            pos=d_rank(i)
            do j=1, r_siz
                ! try to add this to below average population
                new_pop=d_map(pos) + nparticles_proc(receivers(j)+1)             
                ! modify donor population 
                donor_pop=nparticles_proc(proc_map(d_index(pos))+1)-d_map(pos)
                ! if adding subtracting slot doesn't move processor pop past a bound
                if (new_pop .le. up_thresh .and. donor_pop .ge. low_thresh ) then
                    ! changing processor population
                    nparticles_proc(proc_map(d_index(pos))+1)=donor_pop
                    nparticles_proc(receivers(j)+1)=new_pop
                    ! changing proc_map
                    proc_map(d_index(pos))=receivers(j)
                    ! leave the j loop, could be more than one receiver
                    exit
	            end if
	        end do
	    end do
    
    end subroutine redistribute

    subroutine reduce_slots(donors, d_map, d_index, d_map_size, slot_list, d_siz)
        
        ! reduce the size of array we have to search when finding large/small slots to redistribute
        ! In: 
        !   donors: array containing donor processors
        !   d_siz: length of donors
        !   d_map: array containing populations of donor slots which we try and redistribute
        !   d_index: array containing index of entries in d_map in proc_map
        !   d_map_size: length of d_map/d_index
        !   slot_list: array containing populations of slots across all processors

        use parallel, only: iproc
        integer, intent(in) :: d_map_size, d_siz
        integer, intent (in) :: donors(d_siz)
        integer(lint), intent(in) :: slot_list(p_map_size) 
        integer(lint), intent(inout) :: d_map(d_map_size)
        integer, intent(inout) :: d_index(d_map_size)
        integer :: i, j, k
        k=1
        do i=1, d_siz
            do j=1, p_map_size
                ! putting appropriate blocks of slots in d_map 
                if(proc_map(j)==donors(i)) then
                    d_map(k)=slot_list(j)
                    ! index is important as well
                    d_index(k)=j
                    k=k+1
                end if 
            end do
        end do

    end subroutine reduce_slots
     
    subroutine insertion_rank_int(arr, rank, tolerance)
        
        ! ranking algorithm modified for integer valued arrays

        ! Rank a real(p) array in increasing order using the insertion sort
        ! algorithm.
        !
        ! Resultant ranking is *stable* and insertion sort is really pretty
        ! decent, especially when the array is small.  Naturally one should
        ! investigate quicksort and the like for large arrays, but insertion
        ! sort is a good compromise and easy to code.
        !
        ! Based upon the F90 insertion sort implementation in Rosetta Stone
        ! http://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran.
        !
        ! ***WARNING***
        ! We assume that the arrays are 1-indexed due to a feature of how array
        ! bounds are handled in procedure arguments in Fortran.
        !
        ! In:
        !   arr: array of real values.
        !   tolerance (optional, default 0.0): tolerance for comparing values.
        !   If present, then two values which differ by less than tolerance are
        !   treated as equal.  This allows one to avoid changing the order of
        !   items that have identical values to within a desired precision.
        !   Naturally, one should have this as small as possible.
        ! In/Out:
        !    rank: on output rank(i) contains the ranked index (in increasing
        !    order) of the value in arr(i), that is arr(rank(i)) returns the
        !    i-th element in the sorted list of values of arr.
        !    NOTE: rank must have at least the dimensions of arr on input.

        integer(lint), intent(in) :: arr(:)
        integer, intent(inout) :: rank(:) ! inout to avoid automatic deallocation
                                          ! of an allocatable array on entry
        integer :: i, j, tmp, tolerance,tol

        forall (i=1:size(arr)) rank(i) = i
        tol=tolerance
        do i = 2, size(arr)
            j = i - 1
            tmp = rank(i)
            do while ( j >= 1 )
                if (arr(rank(j)) - arr(tmp) < tol) exit
                rank(j+1) = rank(j)
                j = j - 1
            end do
            rank(j+1) = tmp
        end do

    end subroutine insertion_rank_int
    
    subroutine find_processors(rec_dummy, don_dummy, nparticles_proc, donor_slots)
        
        ! Find donor/receiver processors
        ! put these into varying size array receivers/donors
        ! In:
        !   nparticles_proc: num particles on each processor
        ! In/Out:
        !   rec_dummy/don_dummy: arrays which contain donor/receivers, entries are then kept in 
        !       smaller arrays 
        !   donor_slots: number of slots which we can donate, this varies as more entries in proc_map are 
        !       modified.

        use parallel, only: nprocs
        
        integer ::  i, j, k, upper, lower
        integer, allocatable, dimension(:) ::  tmp_rec, tmp_don
        integer, allocatable :: rec_dummy(:), don_dummy(:)
        integer(lint), intent(in) :: nparticles_proc(nprocs)
        integer, intent(inout) :: donor_slots
        allocate(tmp_rec(nprocs))
        allocate(tmp_don(nprocs))

        j=1
        k=1
        ! find donors/receivers
        do i=1, nprocs
            if(nparticles_proc(i) .lt. low_thresh) then
                ! rank = 0 is a problem
                tmp_rec(j)=i-1
                j=j+1
            else if (nparticles_proc(i) .gt. up_thresh) then
                tmp_don(k)=i-1
                k=k+1
            end if 
        end do

        ! put index of processors into smaller array
        
        allocate(rec_dummy(j-1))
        allocate(don_dummy(k-1))

        do i=1, j-1
            rec_dummy(i)=tmp_rec(i)
        end do
        
        do i=1, k-1
            don_dummy(i)=tmp_don(i)
        end do

        ! calculate number of donor slots which we can move
        
        donor_slots=0
        do i=1, p_map_size
            do j=1, size(don_dummy)
                if(proc_map(i)==don_dummy(j)) then
                    donor_slots=donor_slots+1
                end if 
            end do 
        end do 

    end subroutine find_processors
    
    subroutine particles_per_proc(proc_map, slot_list, nparticles_proc, p_map_size)

        ! Find number of particles per processor and store in array nparticles_proc
        ! In:
        !   proc_map: array which maps determinants to processors
        !       proc_map(modulo(hash(d),num_slots*nprocs)=processor
        !   p_map_size: length of proc_map
        !   slot_list: array containing populations of slots in proc_map across all
        !       processors
        ! Out:
        !   nparticles_proc(nprocs): array containing population on each processor

        use parallel, only : nprocs

        integer :: i
        integer, intent(in) :: p_map_size
        integer(lint), intent(out) :: nparticles_proc(nprocs)
        integer, intent(in) :: proc_map(p_map_size)
        integer(lint), intent(in) :: slot_list(p_map_size) 
        nparticles_proc=0  
        
        do i=1, p_map_size
            nparticles_proc(proc_map(i)+1)=nparticles_proc(proc_map(i)+1)+slot_list(i)
        end do

    end subroutine particles_per_proc

    subroutine initialise_slot_pop(slot_pop, proc_map, p_map_size, num_slots)

        ! In/Out:
        !   slot_pop(p_map_size): array containing population of slots in proc_map
        !       processor dependendent
        ! In: 
        !   proc_map(p_map_size): array which maps determinants to processors
        !       proc_map(modulo(hash(d),num_slots*nprocs)=processor
        !   num_slots: number of slots which we divide slot_pop into

        use parallel, only: nprocs, iproc
        use hashing 
        use basis, only: basis_length
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
       
        integer :: i, det_pos
        integer, intent(in):: p_map_size, num_slots
        integer(lint), intent(out) :: slot_pop(p_map_size)        
        type(det_info) :: cdet
        integer, intent(in) :: proc_map(p_map_size)
        slot_pop=0
        do i=1, tot_walkers
            cdet%f => walker_dets(:,i)
            det_pos=modulo(murmurhash_bit_string(cdet%f, basis_length, 7),num_slots*nprocs)+1
            slot_pop(det_pos)=slot_pop(det_pos)+abs(walker_population(1,i))
        end do
   
   end subroutine initialise_slot_pop
    
end module load_balancing
