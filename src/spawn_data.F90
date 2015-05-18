module spawn_data

! Generic structures for storing and manipulating spawned particles (ie
! psips/excips in FCIQMC/CCMC/DMQMC calculations).

use const, only: p, dp, int_s, int_p
use parallel, only: parallel_timing_t
use, intrinsic :: iso_c_binding, only: c_int

implicit none

! Holder for mapping the hash of a state to a processor.
type proc_map_t
    ! Number of slots (initially) assigned to each processor.
    integer :: nslots
    ! Array which maps particles to processors. If attempting load balancing
    ! then proc_map is initially subdivided into load_balancing_slots number of
    ! slots which cyclically map particles to processors using modular
    ! arithmetic. Otherwise it's entries are
    ! 0,1,..,nprocs-1.
    ! map(modulo(hash(state),nslots*nprocs)) gives the processor to which the
    ! state is assigned.
    integer, allocatable :: map(:) ! nslots*nprocs
end type proc_map_t

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
    ! Total number of bits in the bit string used to represent the label.
    ! This allows us to assign a processor based upon hashing a given number of bits and
    ! hence make the processor assignment independent of the size of the integer used in
    ! the bit string (which naturally includes some padding if bit_str_nbits is not
    ! a multiple of i0_length).
    integer :: bit_str_nbits
    ! Number of types of of different particles which are spawned at the same
    ! time.
    integer :: ntypes
    ! sdata(flag_indx,i) is the element containing flags (ie additional info)
    ! for the i-th element.
    integer :: flag_indx
    ! Total number of elements in each spawned object/element (ie len(sdata(:,1))).
    integer :: element_len
    ! The minimum allowed population of a spawning event. Any events with
    ! with populations below this threshold should be stochastically rounded
    ! up to it or down to zero.
    ! During creation of a spawn_t instance, cutoff is multiplied by
    ! 2^(real_bit_shift) and rounded up to nearest integer, and then stored in
    ! this format.
    integer(int_p) :: cutoff

    ! sdata holds the spawned particles and associated information.
    ! The structure of sdata is:
    !   sdata(1:bit_str_len,i)
    !       the bit string of the location (ie determinant) of the i-th spawned particle(s).
    !   sdata(bit_str_len+j,i), 1<=j<=ntypes
    !       population of a particle of type j on the location above.
    !   sdata(flag_indx,i)
    !       any flags (by setting bits) associated with the spawned particles.
    ! element_len = bit_str_len + ntypes + 1.

    ! With MPI parallelisation, a particle on one processor can create a spawned
    ! particle on any other processor.  To make communication simple, we assign
    ! the first block_size entries to particles being sent to processor 0, the next
    ! block_size entries to those being sent to processor 1 and so on.

    ! OpenMP parallelisation complicates matters slightly as we wish to avoid
    ! critical/atomic blocks when spawning.  sdata is therefore an *interleaved* array.
    ! If sdata(:,i) is created by thread 0, then sdata(:,i+1) is created by thread 1,
    ! sdata(:,i+2) by thread 2, ..., sdata(:,i+nthreads-1) by thread nthreads-1 and
    ! sdata(:,i+nthreads) by thread 0 again.  The location to use in sdata is controlled
    ! by head.

    ! NOTE: the above assumes that each thread will produce a roughly equal number of
    ! spawned particles and each processor will be sent a roughly equal number
    ! of spawned particles from a given processor.

    ! Note sdata is actually allocated but actually points to an internal store (store1 or
    ! store2) for efficient communication and simpler code.
    integer(int_s), pointer :: sdata(:,:) ! (element_len,array_len)

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
    ! (Note:
    !       1) the first thing the particle creation routines do is add nthread to
    !          the appropriate element of head to find the next empty slot to spawn
    !          into.
    !       2) head_start(nthreads-1,i) gives the index before the first element of
    !          a particle to be sent to processor i.
    ! )
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
    ! Processor map for assigning states to a given processor.
    ! WARNING: it is the programmer's responsibility to ensure this is consistent
    ! across all spawn_t objects being used with the same set of particles.
    type(proc_map_t) :: proc_map
    ! Information on timings related to MPI communication.
    type(parallel_timing_t) :: mpi_time
    ! Private storage arrays for communication.
    ! Array for receiving spawned particles from other processes (like sdata points to
    ! store1/store2).
    integer(int_s), pointer, private :: sdata_recvd(:,:) ! (element_len,array_len)
    ! Memory stores for spawned particles.
    ! In order to avoid an additional copy, we receive the spawned particles from sdata
    ! to sdata_recvd and then point sdata to the memory previously associated with
    ! sdata_recvd and vice versa.  This allows for the spawning and annihilation
    ! procedures to be identical with and without MPI parallelisation.
    integer(int_s), pointer, private :: store1(:,:), store2(:,:) ! (element_len,array_len)
end type spawn_t

interface annihilate_wrapper_spawn_t
    module procedure annihilate_wrapper_spawn_t_single
    module procedure annihilate_wrapper_spawn_t_arr
end interface

contains

!--- Initialisation/finalisation ---

    subroutine alloc_spawn_t(bit_str_len, bit_str_nbits, ntypes, flag, array_len, cutoff, bit_shift, proc_map, hash_seed, &
                             mpi_barriers, spawn)

        ! Allocate and initialise a spawn_t object.

        ! In:
        !    bit_str_len, bit_str_nbits, ntypes, array_len, hash_seed: see description of
        !       matching components in the spawn_t definition.
        !    flag: whether or not to append an element for storing flags (ie
        !       additional information in bit string format) for the data stored
        !       for each entry.
        !    cutoff: The size of the minimum spawning event allowed.
        !    bit_shift: The number of bits to shift the cutoff by when encoding it.
        !    mpi_barriers: If true then use an mpi_barrier call to measure
        !        load balancing before semi-stochastic communication.
        ! Out:
        !    spawn: initialised object for storing, communicating and
        !       annihilation spawned particles.

        use const, only: int_s_length

        use parallel, only: nthreads, nprocs
        use checking, only: check_allocate
        use errors, only: stop_all

        integer, intent(in) :: bit_str_len, ntypes, array_len, hash_seed, bit_shift, bit_str_nbits
        real(p) :: cutoff
        logical, intent(in) :: flag, mpi_barriers
        type(proc_map_t), intent(in) :: proc_map
        type(spawn_t), intent(out) :: spawn

        integer :: ierr, block_size, i, j

        spawn%bit_str_len = bit_str_len
        spawn%bit_str_nbits = bit_str_nbits
        spawn%ntypes = ntypes
        spawn%element_len = spawn%bit_str_len + spawn%ntypes
        if (flag) then
            spawn%element_len = spawn%element_len + 1
            spawn%flag_indx = spawn%element_len
            if (spawn%ntypes > int_s_length) then
                ! A single int_s integer cannot hold the flags for all spaces.
                call stop_all('alloc_spawn_t', 'Cannot support more than int_s_length particle types.')
            end if
        else
            spawn%flag_indx = -1
        end if
        spawn%array_len = array_len
        spawn%hash_seed = hash_seed
        if (mpi_barriers) spawn%mpi_time%check_barrier_time = .true.
        spawn%mpi_time%barrier_time = 0.0_p
        spawn%mpi_time%comm_time = 0.0_p
        ! Convert the spawning cutoff to the encoded representation for walker
        ! populations (see comments for particle_t%pops) and round up to
        ! nearest integer. (It may not be possible to use the exact cutoff
        ! requested. This will be the case if rounding is required).
        spawn%cutoff = ceiling(cutoff*(2_int_p**int(bit_shift,int_p)), int_p)

        ! Dare not risk allocate(A, mold=B) from F2008 yet...
        allocate(spawn%proc_map%map(lbound(proc_map%map,dim=1):ubound(proc_map%map,dim=1)), stat=ierr)
        call check_allocate('spawn%proc_map%map', size(spawn%proc_map%map), ierr)
        spawn%proc_map = proc_map

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

    subroutine memcheck_spawn_t(spawn, warn_level, dont_warn)

        ! In:
        !    spawn: spawn_t object containing spawned particles.
        !    warn_level (optional, default: 0.95): the fraction of slots in the spawn_t
        !       object which, if filled beyond, triggers a warning.  Note that each
        !       processor block/section is considered separately, so this should be called
        !       before any of the compression/annihilation procedures.
        !    dont_warn (optional): if present and true then don't actually print any
        !        warning. This can prevent spamming the user with messages.

        use parallel, only: nthreads, nprocs, iproc
        use utils, only: int_fmt
        use errors, only: stop_all

        type(spawn_t), intent(in) :: spawn
        real, intent(in), optional :: warn_level
        logical, intent(in), optional :: dont_warn

        real :: fill(nprocs), level
        integer :: it
        logical :: warn

        ! Do we actually want to print any warning?
        warn = .true.
        if (present(dont_warn)) warn = .not. dont_warn

        if (warn) then
            if (present(warn_level)) then
                level = warn_level
            else
                level = 0.95
            end if

            it = nthreads - 1
            fill = real(maxval(spawn%head(:,:),dim=1) - spawn%head_start(it,:)) / spawn%block_size
            if (any(fill - level > 0.0)) then
                write (6,'(1X,"# Warning: filled over 95% of spawning array on processor",'//int_fmt(iproc,1)//',".")') iproc
            end if
        end if

    end subroutine memcheck_spawn_t

!--- Helper procedures ---

    subroutine annihilate_wrapper_spawn_t_single(spawn, tinitiator, determ_size)

        ! Helper procedure for performing annihilation within a spawn_t object.

        ! In:
        !    tinitiator: true if the initiator approximation is being used.
        !    determ_size (optional): The size of the deterministic space in
        !       use, on this process. If input then the deterministic states
        !       received from the various processes will be combined in a
        !       separate call to compress_determ_repeats.
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
        integer, intent(in), optional :: determ_size

        integer :: nstates_received(0:nprocs-1)
        integer, parameter :: thread_id = 0

        ! Compress the successful spawning events from each thread so the
        ! spawned list being sent to each processor contains no gaps.
        if (nthreads > 1) call compress_threaded_spawn_t(spawn)

        if (nprocs > 1) then
            ! Send spawned walkers to the processor which "owns" them and
            ! receive the walkers "owned" by this processor.
            call comm_spawn_t(spawn, nstates_received)

            ! Compress the repeats of the various deterministic states, each of
            ! which is received once from each process.
            if (present(determ_size)) call compress_determ_repeats(spawn, nstates_received, determ_size)
        end if

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

    pure function calc_events_spawn_t(spawn) result(nevents)

        ! In:
        !    spawn: spawn_t object containing spawned particles *before* any
        !        communication/thread compression/etc.

        ! Returns:
        !     number of spawning events which occurred on the processor.

        use parallel, only: nthreads

        integer :: nevents
        type(spawn_t), intent(in) :: spawn

        ! Each event in a given thread is nthreads elements away from the
        ! previous event.
        nevents = sum(spawn%head - spawn%head_start) / nthreads

    end function calc_events_spawn_t
    
    subroutine annihilate_wrapper_non_blocking_spawn(spawn, tinitiator, proc)

        ! Helper procedure for performing annihilation within the two spawned
        ! lists required for non blocking communications.

        ! In:
        !    tinitiator: true if the initiator approximation is being used.
        ! In/Out:
        !    spawn: spawn_t object containing spawned particles.
        !        On output subsection of spawn_t object which contains information
        !        of walkers spawned onto current processor is compressed so that
        !        each determinant appears at most once.
        ! In (optional):
        !    proc: if present only annihilate section of spawned list relating
        !        to proc. Also if proc is present we are annihilating the
        !        actual spawned list rather than the received list so need to
        !        compress across threads. Otherwise we annihilate the entirety
        !        of the received list, but don't compress as it already has
        !        been.

        use parallel, only: nthreads, iproc
        use sort, only: qsort

        type(spawn_t), intent(inout) :: spawn
        logical, intent(in) :: tinitiator
        integer, optional, intent(in) :: proc

        integer, parameter :: thread_id = 0
        integer :: start, endp, spawn_zero

        if (present(proc)) then
            ! Compress the sucessful spawning events from each thread so the spawned list
            ! contains no gaps.
            if (nthreads > 1) call compress_threaded_spawn_t(spawn)
            ! Annihilate within spawned list between start and endp.
            ! spawn_zero is the index before the first entry in the spawned
            ! walker list.
            spawn_zero = spawn%head_start(thread_id, proc) + nthreads - 1
            start = spawn_zero + 1
            endp = spawn%head(thread_id, proc)
        else
            spawn_zero = 0
            start = 1
            endp = spawn%head(thread_id, 0)
        end if

        if (endp > spawn_zero) then
            call qsort(spawn%sdata(:,start:endp), endp-spawn_zero, spawn%bit_str_len)
            ! Annihilate within spawned walkers list.
            ! Compress the remaining spawned walkers list.
            if (tinitiator) then
                call annihilate_spawn_t_initiator(spawn, start, endp)
            else
                call annihilate_spawn_t(spawn, start, endp)
            end if
        end if

    end subroutine annihilate_wrapper_non_blocking_spawn

    subroutine calculate_displacements(spawn, send_disp, non_block_spawn)

        ! Work out how many particles we are sending from the current processor
        ! to all other processors. Necessary for non-blocking communications.

        ! In:
        !    spawn: spawn_t object containing spawned particles in blocks (one
        !      per processor).
        ! Out:
        !    send_disp: Each element of this array contains how many walker will be sent
        !      from current processor to every other processor.
        !    non_block_spawn: number of spawned particles on current processor
        !      during current MC cycle.

        use parallel, only: nprocs, iproc, nthreads

        type(spawn_t), intent(in) :: spawn
        integer, intent(out) :: send_disp(0:)
        integer, intent(out) :: non_block_spawn(:)

        integer :: i
        integer, parameter :: thread_id = 0

        do i = 0, nprocs-1
            send_disp(i) = spawn%head(thread_id,i) - spawn%head_start(thread_id,i) + nthreads - 1
        end do
        ! Need to copy the number of walkers we've spawned during current time
        ! step for estimating spawning rate.
        non_block_spawn(1) = sum(send_disp)
        ! After each time step we have main_list + walkers spawned onto current
        ! processor + walkers spawned to other processors. Walkers on current
        ! processor will annihilate with main list and update nparticles
        ! accordingly, so the total number of walkers across processors is
        ! \sum_i nparticles(i) + non_block_spawn(1)_i - send_disp(i)_i.
        non_block_spawn(2) = non_block_spawn(1) - send_disp(iproc)
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

    subroutine comm_spawn_t(spawn, nstates_received)

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
        ! Out:
        !    nstates_received: Array holding the number of spawned states
        !       received from each process.

#ifdef PARALLEL

        use parallel

        type(spawn_t), intent(inout) :: spawn
        integer, intent(out) :: nstates_received(0:nprocs-1)

        integer :: send_counts(0:nprocs-1), send_displacements(0:nprocs-1)
        integer :: receive_counts(0:nprocs-1), receive_displacements(0:nprocs-1)
        integer :: i, ierr
        integer(int_s), pointer :: tmp_data(:,:)

        real(p) :: t1
        ! Must be compressed by now, if using threads.
        integer, parameter :: thread_id = 0

        ! Start by timing an MPI_Barrier call, which can indicate potential
        ! load balancing issues.
        if (spawn%mpi_time%check_barrier_time) then
            t1 = real(MPI_WTIME(), p)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            spawn%mpi_time%barrier_time = spawn%mpi_time%barrier_time + real(MPI_WTIME(), p) - t1
        end if

        t1 = real(MPI_WTIME(), p)

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
            send_counts(i) = spawn%head(thread_id,i) - spawn%head_start(nthreads-1,i)
        end forall

        call MPI_AlltoAll(send_counts, 1, MPI_INTEGER, receive_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

        nstates_received = receive_counts

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
        ! int_s/mpi_sdata_integer) so we need to change the counts and
        ! displacements accordingly:
        send_counts = send_counts*spawn%element_len
        receive_counts = receive_counts*spawn%element_len
        ! displacement is the number of elements in the preceding processor blocks.
        send_displacements = spawn%head_start(nthreads-1,:)*spawn%element_len
        receive_displacements = receive_displacements*spawn%element_len

        call MPI_AlltoAllv(spawn%sdata, send_counts, send_displacements, mpi_sdata_integer, &
                           spawn%sdata_recvd, receive_counts, receive_displacements, mpi_sdata_integer, &
                           MPI_COMM_WORLD, ierr)

        spawn%mpi_time%comm_time = spawn%mpi_time%comm_time + real(MPI_WTIME(), p) - t1

        ! Swap pointers so that spawn%sdata points to the received data.
        tmp_data => spawn%sdata
        spawn%sdata => spawn%sdata_recvd
        spawn%sdata_recvd => tmp_data

#else
        use parallel, only: nprocs

        type(spawn_t), intent(inout) :: spawn
        integer, intent(out) :: nstates_received(0:nprocs-1)
#endif

    end subroutine comm_spawn_t

    subroutine receive_spawned_walkers(received_list, req_data_s)

        ! Receive walkers spawned onto this processor from previous iteration.

        ! In/Out:
        !   received_list: spawn_t list we receive spawned walkers into.
        !       Upon output received_list%head(0,0) contains the number of
        !       walkers spawned onto current processor.
        !   req_data_s: array of requests initialised from previous iteration's
        !       send of walker list.

#ifdef PARALLEL

        use parallel

        type(spawn_t), intent(inout):: received_list
        integer, intent(inout) :: req_data_s(0:)

        integer :: i, ierr, start, empty_message, counter
        integer, parameter :: thread_id = 0
        integer :: stat_probe_r(MPI_STATUS_SIZE), stat_data_r(MPI_STATUS_SIZE)
        integer :: stat_data_s(MPI_STATUS_SIZE, nprocs)

        start = 1
        do i = 0, nprocs-1
            ! Probe incoming messages. We receive the messages as they arrive
            ! and decide how to insert them into received_list.
            call MPI_Probe(i, iproc, MPI_COMM_WORLD, stat_probe_r, ierr)
            ! Find the size of the message size.
            call MPI_Get_Count(stat_probe_r, mpi_sdata_integer, counter, ierr)
            counter = counter / received_list%element_len
            if (counter > 0) then
                ! Actually receiving walkers from another processor so enter
                ! them in the received_list one after the other
                call MPI_Recv(received_list%sdata(:,start:start+counter-1), counter*received_list%element_len, &
                              mpi_sdata_integer, stat_probe_r(MPI_SOURCE), stat_probe_r(MPI_TAG), MPI_COMM_WORLD, &
                              stat_data_r, ierr)
            else
                ! We've been sent a zero sized message. Do nothing with it, but still
                ! need to receive to complete the send.
                call MPI_Recv(empty_message, 0, mpi_sdata_integer, stat_probe_r(MPI_SOURCE), &
                              stat_probe_r(MPI_TAG), MPI_COMM_WORLD, stat_data_r, ierr)
            end if
            start = start + counter
        end do

        ! Complete the non-blocking send.
        call MPI_Waitall(nprocs, req_data_s, stat_data_s, ierr)
        ! Define the following as the number of determinants in received list.
        received_list%head(thread_id, 0) = start - 1

#else
        type(spawn_t), intent(inout):: received_list
        integer, intent(inout) :: req_data_s(0:)
#endif

    end subroutine receive_spawned_walkers

    subroutine non_blocking_send(spawn, send_counts, req_data_s)

        ! Send remaining walkers from spawned walker list to their
        ! new processors.

        ! In:
        !   send_counts: array containing the number of spawned walkers
        !       we wish to send to each processor.
        ! In/Out:
        !   spawn: spawn_t object containing spawned particles which
        !       need to be sent to the appropriate processor.
        !   req_data_s: array of requests used when sending spawned list.

#ifdef PARALLEL

        use parallel

        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout) :: send_counts(0:)
        integer, intent(inout) :: req_data_s(0:)

        integer, parameter :: thread_id = 0
        integer :: i, start_point, end_point
        integer(int_s), pointer :: tmp_data(:,:)
        integer :: ierr

        ! We want to copy the spawned walker list to another store for communication.
        ! This is to avoid any potential accesses of the send buffer. Potentially need
        ! different spawn_t type. I'm not sure if there is an issue because they
        ! are both elements of the same type?
        tmp_data => spawn%sdata_recvd
        spawn%sdata_recvd => spawn%sdata
        spawn%sdata => tmp_data

        ! Each element contains element_len integers (of type
        ! i0/mpi_det_integer) so we need to change the counts and
        send_counts = send_counts*spawn%element_len

        ! Send the walkers to their processors. The information may not send immediately
        ! due to potentially large messages. So we don't want to access the array.
        ! In the worse case scenario the send will only complete when there is a matching receive posted.
        ! As a result the MPI_Waitall, which is required to ensure message completion, is called in the receive
        ! subroutine to prevent blocking at this point.
        do i = 0, nprocs-1
            start_point = spawn%head_start(thread_id, i) + nthreads
            end_point = start_point + send_counts(i)/spawn%element_len - 1
            call MPI_ISend(spawn%sdata_recvd(:,start_point:end_point), send_counts(i), mpi_sdata_integer, i, &
                           i, MPI_COMM_WORLD, req_data_s(i), ierr)
        end do

#else
        type(spawn_t), intent(inout) :: spawn
        integer, intent(in) :: send_counts(0:)
        integer, intent(inout) :: req_data_s(0:)
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
        !    start (optional): Starting point in spawn_t object we search from.
        !    endp (optional): Search spawn_t object until we reach endp.

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
            if (any(spawn%sdata(spawn%bit_str_len+1:,islot) /= 0_int_s)) islot = islot + 1
        end do self_annihilate

        ! We didn't check if the population on the last determinant is
        ! completely annihilated or not.
        if (all(spawn%sdata(spawn%bit_str_len+1:,islot) == 0_int_s)) islot = islot - 1

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
        !    start (optional): Starting point in spawn_t object we search from.
        !    end (optional): Search spawn_t object until we reach head.

        ! This version is for the initiator algorithm, whereby we also need to
        ! take care of the parent flag (ie handle the origin of the spawned
        ! particles).

        type(spawn_t), intent(inout) :: spawn
        integer, intent(in), optional :: start, endp

        integer :: islot, ipart, k, pop_sign, upper_bound
        integer, allocatable :: events(:)
        integer(int_s), allocatable :: initiator_pop(:)
        ! thread_id is a convention from when OpenMP threading support was added to spawn_t.
        ! Here we're single-threaded, so just use the first element, thread_id=0
        integer, parameter :: thread_id = 0
        logical :: same_slot

        allocate(events(spawn%ntypes))
        allocate(initiator_pop(spawn%ntypes))

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
            ! After compression and communication, spawn%sdata is contiguous and sorted
            ! and spawn%head(0,0) contains the last index in spawn%sdata which has
            ! a particle spawned in it from this iteration.
            upper_bound = spawn%head(thread_id,0)
        end if

        ! initiator_pop(:) and events(:) is used to store the data (for all ntypes of
        ! particle) at location which is being compressed into islot.
        associate(spawn_parts => spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes,:), &
                  spawn_flag => spawn%sdata(spawn%flag_indx,:), spawn_dets => spawn%sdata(:spawn%bit_str_len,:))
            self_annihilate: do
                ! Set the current free slot to be the next unique spawned location.
                spawn%sdata(:,islot) = spawn%sdata(:,k)
                do ipart = 1, spawn%ntypes
                    if (spawn_flag(k) == 0_int_s) then
                        ! from an initiator
                        initiator_pop(ipart) = spawn_parts(ipart,k)
                        events(ipart) = 0
                    else
                        initiator_pop(ipart) = 0
                        events(ipart) = sign(1_int_s,spawn_parts(ipart,k))
                    end if
                end do
                compress: do
                    k = k + 1
                    ! Are we still on the same determinant?  2 conditions:
                    ! 1. Have not yet reached the end of the list.
                    same_slot = k <= spawn%head(thread_id,0)
                    ! 2. We've found a particle on a different determinant to the particles we're compressing to islot.
                    if (same_slot) same_slot = all(spawn_dets(:,k) == spawn_dets(:,islot))

                    if (same_slot) then
                        ! Accumulate the population on this determinant
                        if (spawn_flag(k) == 0_int_s) then
                            ! This slot (k) was spawned from an initiator, so we accumulate
                            ! the population from each type (separately).
                            initiator_pop = initiator_pop + spawn_parts(:,k)
                        else
                            do ipart = 1, spawn%ntypes
                                ! Not an initiator, so depending on the sign of the spawning
                                ! we accumulate (signed) events.  Take care not to accumulate
                                ! events if no ipart particle was spawned.
                                if (spawn_parts(ipart,k) < 0_int_s) then
                                    events(ipart) = events(ipart) - 1
                                else if (spawn_parts(ipart,k) > 0_int_s) then
                                    events(ipart) = events(ipart) + 1
                                end if
                            end do
                        end if
                        ! For each type of particle, add the spawned particles into this slot.
                        spawn_parts(:,islot) = spawn_parts(:,islot) + spawn_parts(:,k)
                    else
                        ! Found the next unique spawned location.
                        ! Finalise the data we've been compressing to islot.
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
                        ! Assume an initiator to start with.
                        spawn_flag(islot) = 0_int_s
                        do ipart = 1, spawn%ntypes
                            ! For each type of particle:
                            if (initiator_pop(ipart) /= 0_int_s .and.  &
                                    sign(1_int_s,spawn_parts(ipart,islot)) == sign(1_int_s,initiator_pop(ipart)) ) then
                                ! There were some particles spawned from initiators.  If any
                                ! particles were spawned from non-initiators, they (in total) are
                                ! either of the same sign as the total from initiators (i.e. multiple
                                ! coherent events) or do not entirely annihilate the particles
                                ! from initiators (i.e. the remaining particles are regarded
                                ! as coming from initiators).
                                ! Keep.
                            else if (abs(events(ipart)) > 1) then
                                ! If the net number of events is over 1 at this particle type,
                                ! then there are multiple coherent spawning events remaining
                                ! after removing pairs of spawning events of the opposite sign.
                                ! Keep.
                            else
                                ! Set flag to only keep if determinant is already occupied
                                ! with particles of this type.
                                spawn_flag(islot) = spawn_flag(islot) + 2**(ipart-1)
                            end if
                        end do
                        ! Now move onto the next determinant (which we just found!).
                        exit compress
                    end if
                end do compress
                ! All done?
                if (islot == upper_bound .or. k > upper_bound) exit self_annihilate
                ! go to the next slot if the current determinant wasn't completed
                ! annihilated.
                if (any(spawn_parts(:,islot) /= 0_int_s)) islot = islot + 1
            end do self_annihilate

            ! We didn't check if the population on the last determinant is
            ! completely annihilated or not.
            if (all(spawn_parts(:,islot) == 0_int_s)) islot = islot - 1
        end associate

        ! update spawn%head(thread_id,0)
        spawn%head(thread_id,0) = islot

    end subroutine annihilate_spawn_t_initiator

    subroutine compress_determ_repeats(spawn, nstates_received, determ_size)

        ! In/Out:
        !    spawn: spawn_t object containing spawned particles.  On output,
        !        this is compressed so that each deterministic state appears in
        !        the spawned list only once.
        ! In:
        !    nstates_received: Array holding the number of spawned states
        !       received from each process.
        !    determ_size: The number of deterministic states belonging to this
        !       process.

        use parallel, only: nprocs

        type(spawn_t), intent(inout) :: spawn
        integer, intent(in) :: nstates_received(0:nprocs-1)
        integer, intent(in) :: determ_size

        integer :: i, j, displacement, ierr
        integer :: min_ind_0, max_ind_0, min_ind_i, max_ind_i
        integer :: nstates_left(0:nprocs-1)
        integer, parameter :: thread_id = 0

        ! The minimum and maximum indices of deterministic states received from
        ! process 0 in the spawned list.
        min_ind_0 = nstates_received(0) - determ_size + 1
        max_ind_0 = nstates_received(0)

        do i = 1, nprocs-1
            displacement = sum(nstates_received(0:i-1))
            ! The minimum and maximum indices of deterministic states received
            ! from process i in the spawned list.
            min_ind_i = displacement + nstates_received(i) - determ_size + 1
            max_ind_i = displacement + nstates_received(i)
            associate(bsl => spawn%bit_str_len)
                ! Add the deterministic spawned amplitudes from process i to
                ! those on process 0, to combine them.
                ! The sign of the spawned state is held indices from bsl+1 to
                ! bsl+spawn%ntypes.
                spawn%sdata(bsl+1:bsl+spawn%ntypes, min_ind_0:max_ind_0) = &
                    spawn%sdata(bsl+1:bsl+spawn%ntypes, min_ind_0:max_ind_0) + &
                    spawn%sdata(bsl+1:bsl+spawn%ntypes, min_ind_i:max_ind_i)
            end associate
        end do

        ! The number of states left in the spawning array, for each process'
        ! section. Process 0 still holds all its states, whilst all other
        ! processes no longer hold a determinstic state.
        nstates_left(0) = nstates_received(0)
        nstates_left(1:nprocs-1) = nstates_received(1:nprocs-1) - determ_size

        ! Now shuffle all non-deterministic states down to fill the gaps.
        ! No need to the move the states from the process numbers 0 and 1.
        do i = 2, nprocs-1
            ! The minimum and maximum indices to move the states to.
            displacement = sum(nstates_left(0:i-1))
            min_ind_0 = displacement + 1
            max_ind_0 = displacement + nstates_left(i)
            ! The minimum and maximum indices to move the states from.
            displacement = sum(nstates_received(0:i-1))
            min_ind_i = displacement + 1
            max_ind_i = displacement + nstates_left(i)

            ! Perform the shuffle.
            spawn%sdata(:, min_ind_0:max_ind_0) = spawn%sdata(:, min_ind_i:max_ind_i)
        end do

        ! Update the number of states in the spawning list.
        spawn%head(thread_id,0) = sum(nstates_left)

    end subroutine compress_determ_repeats

end module spawn_data
