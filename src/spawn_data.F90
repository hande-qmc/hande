module spawn_data

! Generic structures for storing and manipulating spawned particles (ie
! psips/excips in FCIQMC/CCMC/DMQMC calculations).

use const, only: p, dp, i0
use, intrinsic :: iso_c_binding, only: c_int

implicit none

! Holder for spawned objects.
! Each object consists of a single entry, containing:
! * bit string identifying location of the spawned particle(s).
! * population (number) of spawned particles in one or more 'spaces' (e.g. real
!   and imaginary, sampling different operators, etc).
! * an optional flag containing various additional information (e.g. if it was
!   produced by an initiator in iFCIQMC).
! By combining the sets of info into one array, we can reduce the number of MPI
! communication calls during annihilation.
type spawn_t
    ! Maximum number of spawned objects.
    integer :: array_len
    ! Number of integers which make up the identifying bit string (ie some kind
    ! of location of the particle).
    integer :: bit_str_len
    ! Number of types of of different particles which are spawned at the same
    ! time.
    integer :: ntypes
    ! sdata(flag_indx,i) is the element containing flags (ie additional info)
    ! for the i-th element.
    integer :: flag_indx
    ! Total number of elements in each spawned object/element (ie len(sdata(:,1)).
    integer :: element_len
    ! Not actually allocated but actually points to an internal store for
    ! efficient communication.
    integer(i0), pointer :: sdata(:,:)
    ! When creating spawned particles, a particle residing on one processor can
    ! create a particle which needs to reside on any other processor.  To make
    ! the communication easy, we assign different sections of sdata to particles
    ! which should be sent to a given processor.  For thread-safety without
    ! requiring critical/atomic blocks, we assign different indices (using
    ! modulo arithmetic) of each block to a given thread of the processor
    ! producing the child particles.  This is merely for convenience.
    ! NOTE: this assumes that each thread will produce a roughly equal number of
    ! spawned particles and each processor will be sent a roughly equal number
    ! of spawned particles.
    ! block_size is the size of each block of sdata (in terms of the outer-most index)
    ! which is associated with a single processor.
    integer :: block_size
    ! head(j,i) gives the position in sdata of the last spawned particle created
    ! by thread j to be sent to processor i.  The meaning of head is changed by
    ! the compression and communication routines (see below).  head must be set
    ! to head_start at the start of each set of spawning events.
    integer, allocatable :: head(:,:)  ! (0:nthreads-1,0:nprocs-1)
    ! head_start(j,i)+nthreads gives the position in sdata of the first spawned
    ! particle created by thread j to be sent to processor i.
    ! (Note that the first thing the particle creation routines do is add nthread to the
    ! appropriate element of head to find the next empty slot to spawn into.)
    integer, allocatable :: head_start(:,:) ! (0:nthreads-1,0:nprocs-1)
    ! seed for hash routine used to determine which processor a given bit string
    ! should be sent to.
    integer(c_int) :: hash_seed
    ! if non-zero, then hash_shift is mixed with a hash of the bit string and
    ! then rehashed to obtain the process to which a bit string is assigned.
    ! See assign_particle_processor.
    integer :: hash_shift = 0
    ! if hash_shift is non-zero, then the processor a bit string is assigned to
    ! will change once in every 2^move_freq blocks of values for hash_shift.
    ! See assign_particle_processor.
    integer :: move_freq = 0
    ! Time spent doing MPI communication.
    real(dp) :: comm_time = 0.0_dp
    ! Private storage arrays for communication.
    integer(i0), pointer, private :: sdata_recvd(:,:)
    integer(i0), pointer, private :: store1(:,:), store2(:,:) ! (element_len,array_len)
end type spawn_t

interface annihilate_wrapper_spawn_t
    module procedure annihilate_wrapper_spawn_t_single
    module procedure annihilate_wrapper_spawn_t_arr
end interface

contains

!--- Initialisation/finalisation ---

    subroutine alloc_spawn_t(bit_str_len, ntypes, flag, array_len, hash_seed, spawn)

        ! Allocate and initialise a spawn_t object.

        ! In:
        !    bit_str_len, ntypes, array_len, hash_seed: see description of
        !       matching components in the spawn_t definition.
        !    flag: whether or not to append an element for storing flags (ie
        !       additional information in bit string format) for the data stored
        !       for each entry.
        ! Out:
        !    spawn: initialised object for storing, communicating and
        !       annihilation spawned particles.

        use parallel, only: nthreads, nprocs
        use checking, only: check_allocate

        integer, intent(in) :: bit_str_len, ntypes, array_len, hash_seed
        logical, intent(in) :: flag
        type(spawn_t), intent(out) :: spawn

        integer :: ierr, block_size, i, j

        spawn%bit_str_len = bit_str_len
        spawn%ntypes = ntypes
        spawn%element_len = spawn%bit_str_len + spawn%ntypes
        if (flag) then
            spawn%element_len = spawn%element_len + 1
            spawn%flag_indx = spawn%element_len
        else
            spawn%flag_indx = -1
        end if
        spawn%array_len = array_len
        spawn%hash_seed = hash_seed
        spawn%comm_time = 0.0_dp

        allocate(spawn%store1(spawn%element_len, spawn%array_len), stat=ierr)
        call check_allocate('spawn%store1', size(spawn%store1), ierr)
        allocate(spawn%store2(spawn%element_len, spawn%array_len), stat=ierr)
        call check_allocate('spawn%store2', size(spawn%store2), ierr)

        allocate(spawn%head(0:nthreads-1,0:nprocs-1), stat=ierr)
        call check_allocate('spawn%head', size(spawn%head), ierr)
        allocate(spawn%head_start(0:nthreads-1,0:nprocs-1), stat=ierr)
        call check_allocate('spawn%head_start', size(spawn%head_start), ierr)

        spawn%sdata => spawn%store1
        spawn%sdata_recvd => spawn%store2

        ! Allocate uniform chunks of the spawned list to each processor.
        ! Each thread should spawn to once per nthread elements within its
        ! processor block (ie modulo arithmetic).
        ! We manage this by keeping track of the 'head' of the data array
        ! for each thread on each processor (ie the last position spawned
        ! to).  We need to know where to start though...
        spawn%block_size = array_len / nprocs
        forall (i=0:nprocs-1)
            forall (j=0:nthreads-1) spawn%head_start(j,i) = i*spawn%block_size + j - nthreads + 1
        end forall

        spawn%head = spawn%head_start

    end subroutine alloc_spawn_t

    subroutine dealloc_spawn_t(spawn)

        ! In/Out:
        !    spawn: On output the spawn_t object is completely deallocated and
        !       pointers nullified.

        use checking, only: check_deallocate

        type(spawn_t), intent(inout) :: spawn

        integer :: ierr

        if (associated(spawn%store1)) then
            deallocate(spawn%store1, stat=ierr)
            call check_deallocate('spawn%store1', ierr)
        end if
        if (associated(spawn%store2)) then
            deallocate(spawn%store2, stat=ierr)
            call check_deallocate('spawn%store2', ierr)
        end if
        if (allocated(spawn%head_start)) then
            deallocate(spawn%head_start, stat=ierr)
            call check_deallocate('spawn%head_start', ierr)
        end if
        if (allocated(spawn%head)) then
            deallocate(spawn%head, stat=ierr)
            call check_deallocate('spawn%head', ierr)
        end if
        nullify(spawn%sdata)
        nullify(spawn%sdata_recvd)

    end subroutine dealloc_spawn_t

!--- Helper procedures ---

    subroutine annihilate_wrapper_spawn_t_single(spawn, tinitiator)

        ! Helper procedure for performing annihilation within a spawn_t object.

        ! In:
        !    tinitiator: true if the initiator approximation is being used.
        ! In/Out:
        !    spawn: spawn_t object containing spawned particles.  On output, the
        !        spawned particles are sent to the processor which 'owns' the
        !        determinant they are located on and annihilation is performed
        !        internally, so each determinant appears (at most) once in the
        !        spawn%sdata array.

        use parallel, only: nthreads, nprocs
        use sort, only: qsort

        type(spawn_t), intent(inout) :: spawn
        logical, intent(in) :: tinitiator

        integer, parameter :: thread_id = 0

        ! Compress the successful spawning events from each thread so the
        ! spawned list being sent to each processor contains no gaps.
        if (nthreads > 1) call compress_threaded_spawn_t(spawn)

        ! Send spawned walkers to the processor which "owns" them and receive
        ! the walkers "owned" by this processor.
        if (nprocs > 1) call comm_spawn_t(spawn)

        if (spawn%head(thread_id,0) > 0) then

            ! Have spawned walkers on this processor.

            call qsort(spawn%sdata, spawn%head(thread_id,0), spawn%bit_str_len)

            ! Annihilate within spawned walkers list.
            ! Compress the remaining spawned walkers list.
            if (tinitiator) then
                call annihilate_spawn_t_initiator(spawn)
            else
                call annihilate_spawn_t(spawn)
            end if

        end if

    end subroutine annihilate_wrapper_spawn_t_single

    subroutine annihilate_wrapper_spawn_t_arr(spawn_arr, tinitiator)

        ! Helper procedure for performing annihilation for an array of
        ! spawn_t objects.

        ! In:
        !    tinitiator: true if the initiator approximation is being used.
        ! In/Out:
        !    spawn_arr: array of spawn_t objects.  See
        !        annihilate_wrapper_spawn_t_single, which is called for each
        !        element..

        type(spawn_t), intent(inout) :: spawn_arr(:)
        logical, intent(in) :: tinitiator

        integer :: i

        do i = 1, ubound(spawn_arr,dim=1)
            call annihilate_wrapper_spawn_t_single(spawn_arr(i), tinitiator)
        end do

    end subroutine annihilate_wrapper_spawn_t_arr

    subroutine annihilate_wrapper_spawned_list(spawn, tinitiator)

        ! Helper procedure for annihiliating within a portion of a spawn_t object.
        ! For non-blocking communications only particles spawned into chunk of
        ! spawned list corresponding to current processor are annihilated with main
        ! list.

        ! In:
        !    tinitiator: true if the initiator approximation is being used.
        ! In/Out:
        !    spawn: spawn_t object containing spawned particles.
        !       On output subsection of spawn_t object which contains information
        !       of walkers spawned onto current processor is compressed so that
        !       each determinant appears at most once.

        use parallel, only: nthreads, iproc
        use sort, only: qsort

        type(spawn_t), intent(inout) :: spawn
        logical, intent(in) :: tinitiator

        integer, parameter :: thread_id = 0
        integer :: start, endp

        ! Compress the sucessful spawning events from each thread so the spawned list
        ! contains no gaps.
        if (nthreads > 1) call compress_threaded_spawn_t(spawn)

        ! Want to perform annihilation on spawn%sdata(:,start:endp).
        start = spawn%head_start(thread_id, iproc)
        endp = spawn%head(thread_id, iproc)
        if (endp > start) then
            start = start + 1
            call qsort(spawn%sdata, endp, spawn%bit_str_len, start)
            ! Annihilate within spawned walkers list.
            ! Compress the remaining spawned walkers list.
            if (tinitiator) then
                call annihilate_spawn_t_initiator(spawn, start, endp)
            else
                call annihilate_spawn_t(spawn, start, endp)
            end if
        end if

    end subroutine annihilate_wrapper_spawned_list

    subroutine annihilate_wrapper_received_list(spawn, tinitiator)

        ! Helper procedure for preparing list of walkers spawned onto
        ! current processor from last iteration for annihilation with the main list.
        ! This is slightly different from the annihilation of the spawned
        ! list, as walkers here have already undergone compression and the whole
        ! list needs to be sorted and annihilated.

        ! In:
        !    tinitiator: true if the initiator approximation is being used.
        ! In/Out:
        !    spawn: spawn_t object containing spawned particles.
        !       On output spawn_t object will have undergone annihilation within
        !       the list itself and will have been sorted for annihilation with
        !       main list.

        use parallel, only: iproc
        use sort, only: qsort

        logical, intent(in) :: tinitiator
        type(spawn_t), intent(inout) :: spawn

        integer, parameter :: thread_id = 0
        integer :: start, endp

        ! Here we annihilate among the whole list.
        start = 1
        endp = spawn%head(thread_id, 0)
        if (endp > 0) then
            call qsort(spawn%sdata, endp, spawn%bit_str_len, start)

            ! Annihilate within spawned walkers list.
            ! Compress the remaining spawned walkers list.
            if (tinitiator) then
                call annihilate_spawn_t_initiator(spawn, start, endp)
            else
                call annihilate_spawn_t(spawn, start, endp)
            end if
        end if

    end subroutine annihilate_wrapper_received_list

    subroutine calculate_displacements(spawn, send_disp)

        ! Work out how many particles we are sending from the current processor
        ! to all other processors. Necessary for non-blocking communications.

        ! In:
        !    spawn: spawn_t object containing spawned particles in blocks (one
        !      per processor).
        ! Out:
        !    send_disp: Each element of this array contains how many walker will be sent
        !      from current processor to every other processor.

        use parallel, only: nprocs, iproc

        type(spawn_t), intent(in) :: spawn
        integer, intent(out) :: send_disp(0:)

        integer :: i
        integer, parameter :: thread_id = 0

        do i = 0, nprocs-1
            send_disp(i) = spawn%head(thread_id,i) - spawn%head_start(thread_id,i)
        end do
        send_disp(iproc) = 0

    end subroutine calculate_displacements

!--- Thread handling and communication ---

    elemental subroutine compress_threaded_spawn_t(spawn)

        ! Compress the list of spawned particles when created by multiple
        ! threads.

        ! As each thread has its own set of indices (using modulo arithmetic) in
        ! which to place new particles and each thread can create a different
        ! number of particles, we need to remove any gaps in each processor
        ! section in the list of spawned particles.

        ! NOTE: all other annihilation routines assume a contiguous block of
        ! spawned particles in each processor section.

        ! In/Out:
        !    spawn: spawn_t object containing spawned particles.
        !        On output, spawn%head(0,i) contains the number of elements to
        !        be sent to processor i (all other elements of spawn%head are
        !        meaningless) and there are no gaps in each processor-specific
        !        section of spawn%sdata.

        ! (See comments at for spawn_t for more information about how threaded
        ! spawning works.)

        use parallel, only: nthreads, nprocs

        type(spawn_t), intent(inout) :: spawn

        integer :: ip, i, offset, istart, iend, it

        do ip = 0, nprocs-1
            offset = 0
            ! Assuming a large enough number of spawning events, each thread
            ! should have had roughly the same number of successful spawning
            ! events, so this loop should be fast.
            ! First element *not* spawned to:
            istart = minval(spawn%head(:,ip)+nthreads)
            ! Last element spawned to:
            iend = maxval(spawn%head(:,ip))
            do i = istart, iend
                ! Find which thread (should have) created this element.
                ! block_size*ip is the first element in the processor block.
                ! thread indices are from 0 to nthreads-1.
                it = mod(i-spawn%block_size*ip-1,nthreads)
                if (spawn%head(it,ip) < i) then
                    ! This element was *not* spawned into.  Filling in this gap...
                    offset = offset + 1
                else
                    spawn%sdata(:,i-offset) = spawn%sdata(:,i)
                end if
            end do
            spawn%head(0,ip) = i - 1 - offset
        end do

    end subroutine compress_threaded_spawn_t

    subroutine comm_spawn_t(spawn)

        ! Send spawned particles to the pre-designated processor which hosts the
        ! determinant upon which the particle has been spawned.

        ! NOTE: all other annihilation routines assume a contiguous block of
        ! spawned particles (with the exception of compress_spawn_t, as
        ! described above).

        ! In/Out:
        !    spawn: spawn_t object containing spawned particles in blocks (one
        !       per processor).  All particles are sent to the
        !       processor on which they reside and particles are received from
        !       all other processors.  On output, spawn%sdata contains
        !       a contiguous block of spawned objects which reside on this
        !       processor and spawn%head(0,0) contains the number of such
        !       spawned objects (all other elements of spawn%head are
        !       meaningless).

#ifdef PARALLEL

        use parallel

        type(spawn_t), intent(inout) :: spawn

        integer :: send_counts(0:nprocs-1), send_displacements(0:nprocs-1)
        integer :: receive_counts(0:nprocs-1), receive_displacements(0:nprocs-1)
        integer :: i, ierr
        integer(i0), pointer :: tmp_data(:,:)

        real(dp) :: t1
        ! Must be compressed by now, if using threads.
        integer, parameter :: thread_id = 0

        t1 = MPI_WTIME()

        ! Send spawned particles to the processor which "owns" them and receive
        ! the particles "owned" by this processor.

        ! The particles are already stored in the sdata array in blocks,
        ! where each block corresponds to determinants owned by a given
        ! processor.

        ! Tests on cx2 indicate that there is not much difference between
        ! sending messages of the same size using MPI_AlltoAll and
        ! MPI_AlltoAllv (though MPI_AlltoAllv is very slightly slower, by a few
        ! percent).  Therefore it is likely that using MPI_AlltoAllv will be
        ! more efficient as it allows us to only send spawned paticles rather
        ! than the entire spawned lists.  It does require an additional
        ! communication to set up however, so for calculations with large
        ! numbers of particles maybe MPI_AlltoAll would be more efficient?

        ! Find out how many particles we are going to send and receive.
        forall (i=0:nprocs-1)
            send_counts(i) = spawn%head(thread_id,i) - spawn%head_start(thread_id,i) + nthreads - 1
        end forall

        call MPI_AlltoAll(send_counts, 1, MPI_INTEGER, receive_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

        ! Want spawning data to be continuous after move, hence need to find the
        ! receive displacements.
        receive_displacements(0) = 0
        do i=1, nprocs-1
            receive_displacements(i) = receive_displacements(i-1) + receive_counts(i-1)
        end do

        ! Set spawn%head(thread_id,0) to be the number of particles now on this
        ! processor.
        ! This is given by the number of items sent by the last processor plus
        ! the displacement used by the last processor (which is the number of
        ! items sent by all other processors).
        spawn%head(thread_id,0) = receive_displacements(nprocs-1) + receive_counts(nprocs-1)

        ! Send data.
        ! Each element contains element_len integers (of type
        ! i0/mpi_det_integer) so we need to change the counts and
        ! displacements accordingly:
        send_counts = send_counts*spawn%element_len
        receive_counts = receive_counts*spawn%element_len
        ! displacement is the number of elements in the preceding processor blocks.
        send_displacements = (spawn%head_start(thread_id,:nprocs-1)+nthreads-1)*spawn%element_len
        receive_displacements = receive_displacements*spawn%element_len

        call MPI_AlltoAllv(spawn%sdata, send_counts, send_displacements, mpi_det_integer, &
                           spawn%sdata_recvd, receive_counts, receive_displacements, mpi_det_integer, &
                           MPI_COMM_WORLD, ierr)

        spawn%comm_time = spawn%comm_time + MPI_WTIME() - t1

        ! Swap pointers so that spawn%sdata points to the received data.
        tmp_data => spawn%sdata
        spawn%sdata => spawn%sdata_recvd
        spawn%sdata_recvd => tmp_data

#else
        type(spawn_t), intent(inout) :: spawn
#endif

    end subroutine comm_spawn_t

    subroutine receive_spawned_walkers(received_list, req_size_s, req_data_s)

        ! Receive walkers spawned onto this processor from previous iteration.

        ! To avoid having to deal with non-contiguous received list array
        ! we received messages in order of message size so any zero messages
        ! appear at the end of the array and are ignored.

        ! In/Out:
        !   received_list: spawn_t list we receive spawned walkers into.
        !   req_size_s: array of requests initialised from previous iterations's send of message sizes.
        !   req_data_s: array of requests initialised from previous iteration's send of walker list.

#ifdef PARALLEL

        use parallel
        use ranking, only: insertion_rank

        type(spawn_t), intent(inout):: received_list
        integer, intent(inout) :: req_size_s(0:), req_data_s(0:)

        integer :: i, ierr, start, endp
        integer, parameter :: thread_id = 0
        ! These can probably be combined into a type.
        integer :: message_size(0:nprocs-1), message_rank(0:nprocs-1), displacements(0:nprocs-1)
        integer :: req_size_r(0:nprocs-1), req_data_r(0:nprocs-1)
        integer :: stat_size_r(MPI_STATUS_SIZE, nprocs), stat_size_s(MPI_STATUS_SIZE, nprocs)
        integer :: stat_data_s(MPI_STATUS_SIZE, nprocs), stat_data_r(MPI_STATUS_SIZE, nprocs)

        ! Find out the number of walkers we are about to receive onto this processor from all others.
        do i = 0, nprocs-1
            call MPI_IRecv(message_size(i), 1, mpi_integer, i, 123, mpi_comm_world, req_size_r(i), ierr)
        end do
        ! Need to complete the receive and send so that the receive/send buffer and array of
        ! requests can safely be used again.
        call MPI_Waitall(nprocs, req_size_r, stat_size_r, ierr)
        call MPI_Waitall(nprocs, req_size_s, stat_size_s, ierr)

        ! Sort the messages array in descending order so that any messages of size zero
        ! will be at the end of the received list i.e. they will appear after received_list%head(thread_id, 0).
        call insertion_rank(message_size, message_rank)
        ! Fortran indexing.
        message_rank = message_rank - 1
        ! Insertion rank ranks in increasing order, we want decreasing order.
        ! Maybe copying array is a bad idea and could use different iterator in the two loops below.
        message_rank(0:nprocs-1) = message_rank(nprocs-1:0:-1)
        do i = 0, nprocs-1
            ! Messages of size 1 have a zero displacement in an array i.e. they go from
            ! arr(pos:pos). Don't want any negative displacements.
            ! Also message_size = num_dets*element_len.
            if (message_size(message_rank(i)) > 0) then
                displacements(i) = message_size(message_rank(i))/received_list%element_len - 1
            else
                displacements(i) = 0
            end if
        end do

        ! Receive the spawned walkers from other processors.
        ! We collect the messages in descending order of message size. I think this is well behaved
        ! as uniquely size messages will be received from correct processor while for messages of same
        ! size there is nothing to worry about.
        start = 1
        endp = 1
        do i = 0, nprocs-1
            endp = start + displacements(i)
            call MPI_IRecv(received_list%sdata(:,start:endp), message_size(message_rank(i)), &
                           mpi_det_integer, message_rank(i), 456, mpi_comm_world, req_data_r(i), ierr)
            start = endp + 1
        end do
        call MPI_Waitall(nprocs, req_data_r, stat_data_r, ierr)
        call MPI_Waitall(nprocs, req_data_s, stat_data_s, ierr)
        ! Define the following as the number of determinants in received list.
        received_list%head(thread_id, 0) = sum(message_size)/received_list%element_len

#else
        type(spawn_t), intent(inout):: received_list
        integer, intent(inout) :: req_size_s(0:), req_data_s(0:)
#endif

    end subroutine receive_spawned_walkers

    subroutine non_blocking_send(spawn, send_counts, req_size_s, req_data_s)

        ! Send remaining walkers from spawned walker list to their
        ! new processors.

        ! In:
        !   send_counts: array containing the number of spawned walkers we wish to send
        !       to each processor.
        ! In/Out:
        !   spawn: spawn_t object containing spawned particles which
        !       need to be sent to their processor.
        !   req_size_s: array of requests used when sending send_counts array.
        !   req_data_s: array of requests used when sending spawned list.

#ifdef PARALLEL

        use parallel

        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout) :: send_counts(0:)
        integer, intent(inout) :: req_size_s(0:), req_data_s(0:)

        integer, parameter :: thread_id = 0
        integer :: i, start_point, end_point
        integer(i0), pointer :: tmp_data(:,:)
        integer :: ierr, send_copy(0:nprocs-1)

        ! We want to copy the spawned walker list to another store for communication.
        ! This is to avoid any potential accesses of the send buffer. Potentially need
        ! different spawn_t type. I'm not sure if there is an issue because they
        ! are both elements of the same type?
        tmp_data => spawn%sdata_recvd
        spawn%sdata_recvd => spawn%sdata
        spawn%sdata => tmp_data

        send_counts = send_counts*spawn%element_len
        send_copy = send_counts

        ! To receive the walkers we need to know the size of the message.
        ! I think it's easier to send the sizes rather than use MPI_Probe
        ! as these can block, although the receive is essentially blocking anyway.
        ! These are single integers so they should send quickly.
        do i = 0, nprocs-1
            call MPI_ISend(send_counts(i), 1, mpi_integer, i, 123, MPI_COMM_WORLD, req_size_s(i), ierr)
        end do
        ! Send the walkers to their processors. The information may not send immediately
        ! due to potentially large messages. So we don't want to access the array.
        ! In the worse case scenario the send will only complete when there is a matching receive posted.
        ! As a result the MPI_Waitall, which is required to ensure message completion, is called in the receive
        ! subroutine to prevent blocking at this point.
        do i = 0, nprocs-1
            start_point = spawn%head_start(thread_id, i) + 1
            end_point = start_point + send_copy(i)/spawn%element_len - 1
            call MPI_ISend(spawn%sdata_recvd(:,start_point:end_point), send_copy(i), mpi_det_integer, i, 456, MPI_COMM_WORLD, req_data_s(i), ierr)
        end do

#else
        type(spawn_t), intent(inout) :: spawn
        integer, intent(in) :: send_counts(0:)
        integer, intent(inout) :: req_size_s(0:), req_data_s(0:)
#endif

    end subroutine non_blocking_send

!--- Annihilation within spawned data (ie combine particles on same element) ---

    elemental subroutine annihilate_spawn_t(spawn, start, endp)

        ! Annihilate the list of spawned particles: ie sum the populations of
        ! particles residing on each location and remove any locations which
        ! have a zero total.

        ! NOTE: assume the list of spawned particles is sorted by the
        ! identifying bit string.

        ! In/Out:
        !    spawn: spawn_t object containing spawned particles.  On output,
        !        this is compressed so that each location occurs (at most) once
        !        and locations with a total of zero spawned particles (ie all
        !        particles cancel out) are removed.
        ! In:
        !   start (optional): Starting point in spawn_t object we search from.
        !   endp (optional): Search spawn_t object until we reach endp.

        type(spawn_t), intent(inout) :: spawn
        integer, intent(in), optional :: start, endp

        integer :: islot, k, upper_bound
        integer, parameter :: thread_id = 0

        ! The spawned list is already sorted, so annihilation amounts to
        ! looping through the list and adding consective populations together if
        ! they're on the same location.

        ! islot is the current element in the spawned lists.
        ! k is the current element which is being compressed into islot (if
        ! k and islot refer to the same determinants).

        if (present(start)) then
            islot = start
            k = start
        else
            islot = 1
            k = 1
        end if

        if (present(endp)) then
            upper_bound = endp
        else
            upper_bound = spawn%head(thread_id,0)
        end if

        self_annihilate: do
            ! Set the current free slot to be the next unique spawned location.
            spawn%sdata(:,islot) = spawn%sdata(:,k)
            compress: do
                k = k + 1
                if (k > upper_bound) exit self_annihilate
                if (all(spawn%sdata(:spawn%bit_str_len,k) == spawn%sdata(:spawn%bit_str_len,islot))) then
                    ! Add the populations of the subsequent identical particles.
                    spawn%sdata(spawn%bit_str_len+1:,islot) =    &
                         spawn%sdata(spawn%bit_str_len+1:,islot) &
                       + spawn%sdata(spawn%bit_str_len+1:,k)
                else
                    ! Found the next unique spawned particle.
                    exit compress
                end if
            end do compress
            ! All done?
            if (islot == upper_bound) exit self_annihilate
            ! go to the next slot if the current determinant wasn't completed
            ! annihilated.
            if (any(spawn%sdata(spawn%bit_str_len+1:,islot) /= 0)) islot = islot + 1
        end do self_annihilate

        ! We didn't check if the population on the last determinant is
        ! completely annihilated or not.
        if (all(spawn%sdata(spawn%bit_str_len+1:,islot) == 0)) islot = islot - 1

        ! update spawn%head(thread_id,0)
        spawn%head(thread_id,0) = islot

    end subroutine annihilate_spawn_t

    elemental subroutine annihilate_spawn_t_initiator(spawn, start, endp)

        ! Annihilate the list of spawned particles: ie sum the populations of
        ! particles residing on each location and remove any locations which
        ! have a zero total.

        ! NOTE: assume the list of spawned particles is sorted by the
        ! identifying bit string.

        ! In/Out:
        !    spawn: spawn_t object containing spawned particles.  On output,
        !        this is compressed so that each location occurs (at most) once
        !        and locations with a total of zero spawned particles (ie all
        !        particles cancel out) are removed.
        ! In:
        !   start (optional): Starting point in spawn_t object we search from.
        !   end (optional): Search spawn_t object until we reach head.

        ! This version is for the initiator algorithm, whereby we also need to
        ! take care of the parent flag (ie handle the origin of the spawned
        ! particles).

        type(spawn_t), intent(inout) :: spawn
        integer, intent(in), optional :: start, endp

        integer :: islot, ipart, k, pop_sign, upper_bound
        integer, allocatable :: events(:)
        integer(i0), allocatable :: initiator_pop(:)
        integer, parameter :: thread_id = 0
        logical :: new_slot

        allocate(events(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes))
        allocate(initiator_pop(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes))

        ! islot is the current element in the spawned list.
        ! k is the current element which is being compressed into islot (if
        ! k and islot refer to the same determinants).
        if (present(start)) then
            islot = start
            k = start
        else
            islot = 1
            k = 1
        end if

        if (present(endp)) then
            upper_bound = endp
        else
            upper_bound = spawn%head(thread_id,0)
        end if

        self_annihilate: do
            ! Set the current free slot to be the next unique spawned location.
            spawn%sdata(:,islot) = spawn%sdata(:,k)
            do ipart = spawn%bit_str_len+1, spawn%bit_str_len+spawn%ntypes
                if (spawn%sdata(spawn%flag_indx,k) == 0) then
                    ! from an initiator
                    initiator_pop(ipart) = spawn%sdata(ipart,k)
                    events(ipart) = 0
                else
                    initiator_pop(ipart) = 0
                    events(ipart) = sign(1_i0,spawn%sdata(ipart,k))
                end if
            end do
            compress: do
                k = k + 1
                new_slot = k > spawn%head(thread_id,0)
                if (.not. new_slot) new_slot = any(spawn%sdata(:spawn%bit_str_len,k) /= spawn%sdata(:spawn%bit_str_len,islot))
                if (new_slot) then
                    ! Found the next unique spawned location.
                    ! Set the overall parent flag of the population on the
                    ! current determinant and move on.
                    ! Rules:
                    !   * keep psips from initiator determinants
                    !   * keep psips from multiple coherent events
                    ! These are indicated by a 0 flag.
                    !   * keep other psips only if determinant is already
                    !     occupied.
                    ! We also want the outcome to be independent of the order
                    ! the spawning events occured in, hence accumulating the
                    ! signed number of events and number of initiator particles
                    ! separately.
                    ! Corner cases are tricky!
                    ! * If multiple initiator events exactly cancel out, then the
                    !   flag is determined by the number of coherent events from
                    !   non-initiator parent determinants.
                    ! * If more psips from non-initiators are spawned than
                    !   initiators and the two sets have opposite sign, the flag
                    !   is determined by number of coherent events from
                    !   non-initiator parents.
                    spawn%sdata(spawn%flag_indx,islot) = 0
                    do ipart = spawn%bit_str_len+1, spawn%bit_str_len+spawn%ntypes
                        if (initiator_pop(ipart) /= 0 .and.  &
                                sign(1_i0,spawn%sdata(ipart,islot)) == sign(1_i0,initiator_pop(ipart)) ) then
                            ! Keep all.  We should still annihilate psips of
                            ! opposite sign from non-initiator events(spawn%bit_str_len+1).
                        else if (abs(events(spawn%bit_str_len+1)) > 1) then
                            ! Multiple coherent spawning events(spawn%bit_str_len+1) after removing pairs
                            ! of spawning events(spawn%bit_str_len+1) of the opposite sign.
                            ! Keep.
                        else
                            ! Should only keep if determinant is already occupied.
                            spawn%sdata(spawn%flag_indx,islot) = spawn%sdata(spawn%flag_indx,islot) + 2**ipart
                        end if
                    end do
                    exit compress
                else
                    ! Accumulate the population on this determinant, how much of the population came
                    ! from an initiator and the sign of the event.
                    if (spawn%sdata(spawn%flag_indx,k) == 0) then
                        initiator_pop = initiator_pop + spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes,k)
                    else
                        do ipart = spawn%bit_str_len+1, spawn%bit_str_len+spawn%ntypes
                            if (spawn%sdata(ipart,k) < 0) then
                                events(ipart) = events(ipart) - 1
                            else
                                events(ipart) = events(ipart) + 1
                            end if
                        end do
                    end if
                    spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes,islot) = &
                         spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes,islot) + &
                         spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes,k)
                end if
            end do compress
            ! All done?
            if (islot == upper_bound .or. k > upper_bound) exit self_annihilate
            ! go to the next slot if the current determinant wasn't completed
            ! annihilated.
            if (any(spawn%sdata(spawn%bit_str_len+1:,islot) /= 0)) islot = islot + 1
        end do self_annihilate

        ! We didn't check if the population on the last determinant is
        ! completely annihilated or not.
        if (all(spawn%sdata(spawn%bit_str_len+1:,islot) == 0)) islot = islot - 1

        ! update spawn%head(thread_id,0)
        spawn%head(thread_id,0) = islot

    end subroutine annihilate_spawn_t_initiator

end module spawn_data
