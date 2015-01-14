module semi_stoch

! Code relating to the semi-stochastic algorithm. This includes initialisation
! code and also routines used throughout the simulation.

! Semi-stochastic
! ===============
!
! The semi-stochastic algorithm (PRL 109, 230201) is an adaptation to the
! FCIQMC algorithm whereby the FCIQMC projection is performed exactly within a
! small region of the space. If this region is chosen so as to contain the
! most dominant determinants in the ground state wave function then a
! significant reduction in stochastic noise can be obtained.
!
! Split all the states in the FCI space into 'deterministic' states and
! 'stochastic' states. Let D refer to the deterministic space as a whole and S
! to the stochastic space as a whole. Then, in the semi-stochastic method we
! write the full projection matrix as
!
! P = P_{DD} + P_{SD} + P_{DS} + P_{SS}.
!
! In normal FCIQMC all the elements of this matrix are sampled stochastically.
! In the semi-stochastic algorithm, while P_{SD} (deterministic to stochastic
! spawning), P_{DS} (stochastic to deterministic spawning) and P_{SS}
! stochastic to stochastic spawning) are calculated stochastically, P_{DD} is
! calculated exactly via matrix-vector multiplication.
!
! Key implementation points
! -------------------------
!
! Thus, in this semi-stochastic code, we generate a set of states which will be
! designated 'deterministic states'. These are stored in semi_stocht%dets.
! We then calculate and store the entire Hamiltonian between all pairs of
! deterministic states (stored in a sparse format in semi_stocht%hamil). Note
! that we store the Hamiltonian, not the projection operator
! (1 + \tau S) - \tau \hat{H}, but the latter is easily calculated from the
! former.
!
! Only part of the deterministic Hamiltonian is stored on each processor.
! Each processor will store all columns of the Hamiltonian corresponding to
! deterministic states on that processor only (but the entire column is
! stored). That is, if I is a deterministic state on this processor and J is a
! deterministic state belonging to another processor, then H_{II} and H_{JI}
! will be stored, but H_{IJ} and H_{JJ} will not.
!
! Deterministic states still reside in the main walker arrays, walker_dets,
! walker_populations and walker_data. However, we want to perform spawning from
! deterministic to deterministic states exactly. Thus, whenever a stochastic
! spawning event of this type is generated, it is cancelled. We still have to
! attempt stochastic spawning from deterministic states because deterministic to
! stochastic spawning *is* allowed, and must be treated using the excitation
! generators as usual.
!
! Because all deterministic states are still stored in walker_dets, they do not
! have to be treated differently for the most part. For example, energy
! evaluation is performed exactly as it is without semi-stochastic.
!
! We store an array of flags (semi_stoch_t%flags) which specify whether or not
! states in walker_dets belong to the deterministic space or not. The status
! of newly spawned states is checked on-the-fly.
!
! To check if a newly spawned determinant belongs to the deterministic space
! or not, we use a hash table look-up. This is possible because
! semi_stocht%dets stores *all deterministic states on all processors*
! (a spawned state can of course belong to any processor).
!
! Note that all deterministic states are stored in walker_dets only for that
! one processor (just as for all other states). Deterministic states are never
! removed from walker_dets, even if they have an amplitude of exactly zero
! (which is very unlikely in practice). This simplifies the implemenation in
! some places. No stochastic rounding is ever performed for the populations of
! deterministic states.
!
! As we run through all states in the main algorithm, the populations on
! deterministic states are copied across to a vector (semi_stoch_t%vector).
! This vector is what is multiplied in the deterministic projection itself.
!
! There are two possible methods for adding deterministic spawnings back into
! the main psip array. If determ%separate_annihilation is true, then each
! processor will receive the full list of deterministic amplitudes from all
! processors via an extra MPI call per iteration. Using this complete list, the
! deterministic spawning amplitudes for each processor can be calculated by
! a single matrix multiplication without any further communication. These
! deterministic spawnings are then added back into the main psip array in
! the annihilation routines in deterministic_annihilation, which treats these
! spawnings separately.
!
! If determ%separate_annihilation is false then the deterministic spawning
! amplitudes from a given processor will be calculated by performing the exact
! projection, and the resulting 'spawnings' will then be added to the spawned
! list, just like non-deterministic spawnings. This prevents the need for an
! extra MPI call each iteration, but the spawned list which is communicated
! becomes larger.

use const
use fciqmc_data, only: semi_stoch_t, determ_hash_t

implicit none

contains

    subroutine init_semi_stoch_t(determ, sys, spawn, space_type, target_size, separate_annihilation, write_determ_in)

        ! Create a semi_stoch_t object which holds all of the necessary
        ! information to perform a semi-stochastic calculation. The type of
        ! deterministic space is determined by space_type.

        ! In/Out:
        !    determ: Deterministic space being used.
        ! In:
        !    sys: system being studied
        !    spawn: spawn_t object to which deterministic spawning will occur.
        !    space_type: Integer parameter specifying which type of
        !        deterministic space to use.
        !    target_size: A size of deterministic space to aim for. This
        !        is only necessary for particular deterministic spaces.
        !    separate_annihilation: If true then routines which use the created
        !        determ object will not use the standard annihilation routine
        !        to treat deterministic spawnings, but will handle them
        !        separately.
        !    write_determ_in: If true then write out the deterministic states to
        !        a file.

        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: walker_length, empty_determ_space, high_pop_determ_space
        use fciqmc_data, only: read_determ_space
        use parallel
        use sort, only: qsort
        use spawn_data, only: spawn_t
        use system, only: sys_t
        use utils, only: int_fmt

        type(semi_stoch_t), intent(inout) :: determ
        type(sys_t), intent(in) :: sys
        type(spawn_t), intent(in) :: spawn
        integer, intent(in) :: space_type
        integer, intent(in) :: target_size
        logical, intent(in) :: separate_annihilation
        logical, intent(in) :: write_determ_in

        integer :: i, ierr, determ_dets_mem
        integer :: displs(0:nprocs-1)
        ! dtes_this_proc will hold deterministic states on this processor only.
        ! This is only needed during initialisation.
        integer(i0), allocatable :: dets_this_proc(:,:)
        logical :: print_info, write_determ

        ! Only print information if the parent processor and if we are using a
        ! non-trivial deterministic space.
        print_info = parent .and. space_type /= empty_determ_space

        ! Copy across this input option to the derived type instance.
        determ%separate_annihilation = separate_annihilation

        ! If an empty space is being used then don't dump a semi-stoch file.
        write_determ = write_determ_in .and. space_type /= empty_determ_space

        if (print_info) write(6,'(1X,a43)') '# Beginning semi-stochastic initialisation.'

        allocate(determ%flags(walker_length), stat=ierr)
        call check_allocate('determ%flags', walker_length, ierr)
        allocate(determ%sizes(0:nprocs-1), stat=ierr)
        call check_allocate('determ%sizes', nprocs, ierr)

        determ%sizes = 0
        determ%space_type = space_type

        ! Create temporary space for enumerating the deterministic space
        ! belonging to this processor only.
        allocate(dets_this_proc(sys%basis%tensor_label_len, walker_length), stat=ierr)
        call check_allocate('dets_this_proc', size(dets_this_proc), ierr)
        dets_this_proc = 0_i0

        if (print_info) write(6,'(1X,a31)') '# Creating deterministic space.'

        ! If space_type does not take one of the below values then an empty
        ! deterministic space will be used. This is the default behaviour
        ! (space_type = empty_determ_space).
        if (space_type == high_pop_determ_space) then
            call create_high_pop_space(dets_this_proc, spawn, target_size, determ%sizes(iproc))
        else if (space_type == read_determ_space) then
            call read_determ_from_file(dets_this_proc, determ, spawn, sys, print_info)
        end if

        ! Let each process hold the number of deterministic states on each process.
#ifdef PARALLEL
        ! Note: just using displs as scratch space as experience has taught us
        !to avoid MPI_IN_PLACE (perhaps things are better now though?).
        call mpi_allgather(determ%sizes(iproc), 1, mpi_integer, displs, 1, mpi_integer, MPI_COMM_WORLD, ierr)
        determ%sizes = displs
#endif
        determ%tot_size = sum(determ%sizes)

        if (print_info .and. nprocs > 1) then
            write(6,'(1X,a46,'//int_fmt(minval(determ%sizes),1)//')') &
                '# Min deterministic space size on a processor:', minval(determ%sizes)
            write(6,'(1X,a46,'//int_fmt(maxval(determ%sizes),1)//')') &
                '# Max deterministic space size on a processor:', maxval(determ%sizes)
            write(6,'(1X,a51,'//int_fmt(determ%tot_size,1)//')') &
                '# Total deterministic space size on all processors:', determ%tot_size
        else if (print_info) then
            write(6,'(1X,a27,'//int_fmt(determ%tot_size,1)//')') '# Deterministic space size:', determ%tot_size
        end if

        ! Displacements used for MPI communication.
        displs(0) = 0
        do i = 1, nprocs-1
            displs(i) = displs(i-1) + determ%sizes(i-1)
        end do

        ! Vector to hold deterministic amplitudes from this process.
        allocate(determ%vector(determ%sizes(iproc)), stat=ierr)
        call check_allocate('determ%vector', determ%sizes(iproc), ierr)
        determ%vector = 0.0_p

        ! Vector to hold deterministic amplitudes from all processes.
        if (determ%separate_annihilation) then
            allocate(determ%full_vector(determ%tot_size), stat=ierr)
            call check_allocate('determ%full_vector', determ%tot_size, ierr)
            determ%full_vector = 0.0_p

            allocate(determ%indices(determ%sizes(iproc)), stat=ierr)
            call check_allocate('determ%indices', determ%sizes(iproc), ierr)
            determ%indices = 0
        end if

        call qsort(dets_this_proc, determ%sizes(iproc)) 

        ! Array to hold all deterministic states from all processes.
        ! The memory required in MB.
        determ_dets_mem = sys%basis%tensor_label_len*determ%tot_size*i0_length/(8*10**6)
        if (print_info) write(6,'(1X,a60,'//int_fmt(determ_dets_mem,1)//')') &
            '# Memory required per core to store deterministic dets (MB):', determ_dets_mem
        allocate(determ%dets(sys%basis%tensor_label_len, determ%tot_size), stat=ierr)
        call check_allocate('determ%dets', size(determ%dets), ierr)

        ! Join and store all deterministic states from all processes.
#ifdef PARALLEL
        associate(tbl=>sys%basis%tensor_label_len)
            call mpi_allgatherv(dets_this_proc(:,1:determ%sizes(iproc)), tbl*determ%sizes(iproc), &
                            mpi_det_integer, determ%dets, tbl*determ%sizes, tbl*displs, &
                            mpi_det_integer, MPI_COMM_WORLD, ierr)
        end associate
#else
        determ%dets = dets_this_proc(:,1:determ%sizes(iproc))
#endif

        call create_determ_hash_table(determ, print_info)

        call create_determ_hamil(determ, sys, displs, dets_this_proc, print_info)

        ! All deterministic states on this processor are always stored in
        ! walker_dets, even if they have a population of zero, so they are
        ! added in here.
        call add_determ_dets_to_walker_dets(determ, sys, dets_this_proc)

        ! We don't need this temporary space anymore. All deterministic states
        ! from all processors are stored in determ%dets.
        deallocate(dets_this_proc, stat=ierr)
        call check_deallocate('dets_this_proc', ierr)

        if (write_determ .and. parent) call write_determ_to_file(determ, print_info)

        if (print_info) write(6,'(1X,a42)') '# Semi-stochastic initialisation complete.'

    end subroutine init_semi_stoch_t

    subroutine dealloc_semi_stoch_t(determ)

        ! In/Out:
        !    determ: Deterministic space object to be deallocated.

        use checking, only: check_deallocate
        use csr, only: end_csrp
        use fciqmc_data, only: empty_determ_space

        type(semi_stoch_t), intent(inout) :: determ
        integer :: ierr

        ! Reset the space type to an empty space, the default value, in case
        ! this semi_stoch_t object is reused.
        determ%space_type = empty_determ_space

        if (allocated(determ%sizes)) then
            deallocate(determ%sizes, stat=ierr)
            call check_deallocate('determ%sizes', ierr)
        end if
        if (allocated(determ%vector)) then
            deallocate(determ%vector, stat=ierr)
            call check_deallocate('determ%vector', ierr)
        end if
        if (allocated(determ%full_vector)) then
            deallocate(determ%full_vector, stat=ierr)
            call check_deallocate('determ%full_vector', ierr)
        end if
        if (allocated(determ%indices)) then
            deallocate(determ%indices, stat=ierr)
            call check_deallocate('determ%indices', ierr)
        end if
        if (allocated(determ%dets)) then
            deallocate(determ%dets, stat=ierr)
            call check_deallocate('determ%dets', ierr)
        end if
        if (allocated(determ%flags)) then
            deallocate(determ%flags, stat=ierr)
            call check_deallocate('determ%flags', ierr)
        end if

        call dealloc_determ_hash_t(determ%hash_table)
        call end_csrp(determ%hamil)

    end subroutine dealloc_semi_stoch_t

    subroutine dealloc_determ_hash_t(ht)

        ! In/Out:
        !    ht: Deterministic hash table to be deallocated.

        use checking, only: check_deallocate

        type(determ_hash_t), intent(inout) :: ht
        integer :: ierr
        
        if (allocated(ht%ind)) then
            deallocate(ht%ind, stat=ierr)
            call check_deallocate('ht%ind', ierr)
        end if
        if (allocated(ht%hash_ptr)) then
            deallocate(ht%hash_ptr, stat=ierr)
            call check_deallocate('ht%hash_ptr', ierr)
        end if

    end subroutine dealloc_determ_hash_t

    subroutine create_determ_hash_table(determ, print_info)

        ! Create the hash table for the deterministic space.

        ! In/Out:
        !    determ: Deterministic space being used. On input, determ%tot_size
        !        and determ%dets should be created and set. On output the
        !        determ%hash_table object will be created.
        ! In:
        !    print_info: Should we print information to the screen?

        use checking, only: check_allocate, check_deallocate
        use hashing, only: murmurhash_bit_string
        use utils, only: int_fmt

        type(semi_stoch_t), intent(inout) :: determ
        logical, intent(in) :: print_info

        integer :: i, iz, hash, mem_reqd, ierr, tensor_label_len
        integer, allocatable :: nclash(:)

        tensor_label_len = size(determ%dets, dim=1)

        associate(ht => determ%hash_table)

            ht%seed = 37
            ! For now just let there be as many hash values as deterministic states.
            ht%nhash = determ%tot_size

            if (print_info) then 
                ! The memory required in MB.
                ! Two lists with length ht%nhash taking 4 bytes for each element.
                mem_reqd = 2*ht%nhash*4/10**6
                write(6,'(1X,a52,'//int_fmt(mem_reqd,1)//')') &
                    '# Memory required per core to store hash table (MB):', mem_reqd
            end if

            allocate(ht%ind(determ%tot_size), stat=ierr) 
            call check_allocate('determ%hash_table%ind', determ%tot_size, ierr)
            allocate(ht%hash_ptr(ht%nhash+1), stat=ierr) 
            call check_allocate('determ%hash_table%hash_ptr', ht%nhash, ierr)

            ! Array to count the number of deterministic states with each hash value.
            allocate(nclash(ht%nhash), stat=ierr)
            call check_allocate('nclash', ht%nhash, ierr)
            nclash = 0
                        
            ! Count the number of deterministic states with each hash value.
            do i = 1, determ%tot_size
                hash = modulo(murmurhash_bit_string(determ%dets(:,i), tensor_label_len, ht%seed), ht%nhash) + 1
                nclash(hash) = nclash(hash) + 1
            end do

            ! Fill in the hash_ptr array (see comments in the determ_hash_t type).
            ht%hash_ptr(1) = 1
            do i = 2, ht%nhash
                ht%hash_ptr(i) = ht%hash_ptr(i-1) + nclash(i-1)
            end do
            ht%hash_ptr(ht%nhash+1) = determ%tot_size + 1

            ! Now loop over all states again and this time fill in the hash table.
            nclash = 0
            do i = 1, determ%tot_size
                hash = modulo(murmurhash_bit_string(determ%dets(:,i), tensor_label_len, ht%seed), ht%nhash) + 1
                iz = ht%hash_ptr(hash) + nclash(hash)
                ht%ind(iz) = i
                nclash(hash) = nclash(hash) + 1
            end do

            deallocate(nclash, stat=ierr)
            call check_deallocate('nclash', ierr)

        end associate

    end subroutine create_determ_hash_table

    subroutine create_determ_hamil(determ, sys, displs, dets_this_proc, print_info)

        ! In/Out:
        !    determ: Deterministic space being used. On input, determ%sizes,
        !        determ%tot_size and determ%dets should be created and set. On
        !        output, determ%hamil will have been created.
        ! In:
        !    sys: system being studied
        !    displs: displs(i) holds the cumulative sum of the number of
        !        deterministic states belonging to processor numbers 0 to i-1.
        !    dets_this_proc: The deterministic states belonging to this
        !        processor.
        !    print_info: Should we print information to the screen?

        use checking, only: check_allocate
        use csr, only: init_csrp
        use fciqmc_data, only: H00
        use hamiltonian, only: get_hmatel
        use parallel
        use system, only: sys_t
        use utils, only: int_fmt

        type(semi_stoch_t), intent(inout) :: determ
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: displs(0:nprocs-1)
        integer(i0), intent(in) :: dets_this_proc(:,:)
        logical, intent(in) :: print_info

        integer :: i, j, nnz, imode, ierr
        integer :: mem_reqd, max_mem_reqd
        real(p) :: hmatel
        real :: t1, t2
        logical :: diag_elem

        if (print_info) write(6,'(1X,a74)') '# Counting number of non-zero deterministic Hamiltonian elements to store.'

        associate(hamil => determ%hamil)

            ! For imode = 1 count the number of non-zero elements, for imode = 2 store them.
            do imode = 1, 2
                if (imode == 1) call cpu_time(t1)
                nnz = 0
                ! Over all deterministic states on all processes (all rows).
                do i = 1, determ%tot_size
                    ! Over all deterministic states on this process (all columns).
                    do j = 1, determ%sizes(iproc)
                        hmatel = get_hmatel(sys, determ%dets(:,i), dets_this_proc(:,j))
                        diag_elem = i == j + displs(iproc)
                        ! Take the Hartree-Fock energy off the diagonal elements.
                        if (diag_elem) hmatel = hmatel - H00
                        if (abs(hmatel) > depsilon) then
                            nnz = nnz + 1
                            if (imode == 2) then
                                hamil%mat(nnz) = hmatel
                                hamil%col_ind(nnz) = j
                                if (hamil%row_ptr(i) == 0) hamil%row_ptr(i) = nnz
                            end if
                        end if
                    end do
                    if (imode == 2) then
                        if (hamil%row_ptr(i) == 0) hamil%row_ptr(i) = nnz + 1
                    end if
                end do

                ! Allocate the CSR type components and print information.
                if (imode == 1) then
                    call cpu_time(t2)
                    if (print_info) then 
                        write(6,'(1X,a41,1X,f10.2,a1)') '# Time taken to generate the Hamiltonian:', t2-t1, "s"
                        ! The memory required in MB.
#ifdef SINGLE_PRECISION
                        mem_reqd = ((determ%tot_size+1)*4 + nnz*(4 + 4))/10**6
#else
                        mem_reqd = ((determ%tot_size+1)*4 + nnz*(4 + 8))/10**6
#endif
                    end if

#ifdef PARALLEL
                    call mpi_allreduce(mem_reqd, max_mem_reqd, 1, mpi_integer, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
                    max_mem_reqd = mem_reqd
#endif
                    if (print_info) then 
                        write(6,'(1X,a75,'//int_fmt(mem_reqd,1)//')') &
                            '# Maximum memory required by a core for the deterministic Hamiltonian (MB):', mem_reqd
                        write(6,'(1X,a54)') '# The Hamiltonian will now be recalculated and stored.'
                    end if

                    ! Allocate the CSR Hamiltonian arrays.
                    call init_csrp(hamil, determ%tot_size, nnz)
                    hamil%row_ptr(1:determ%tot_size) = 0
                end if
            end do

        end associate

    end subroutine create_determ_hamil

    subroutine add_determ_dets_to_walker_dets(determ, sys, dets_this_proc)

        ! Also set the deterministic flags of any deterministic states already
        ! in walker_dets, and add deterministic data to walker_populations and
        ! walker_data. All deterministic states not already in walker_dets are
        ! set to have an initial sign of zero.

        ! In/Out:
        !    determ: Deterministic space being used.
        ! In:
        !    sys: system being studied
        !    dets_this_proc: The deterministic states belonging to this
        !        processor.

        use annihilation, only: insert_new_walker
        use fciqmc_data, only: walker_dets, walker_population, walker_data
        use fciqmc_data, only: tot_walkers, sampling_size
        use parallel, only: iproc
        use search, only: binary_search
        use system, only: sys_t

        type(semi_stoch_t), intent(inout) :: determ
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: dets_this_proc(:,:)

        integer :: i, istart, iend, pos
        integer(int_p) :: zero_population(sampling_size)
        logical :: hit

        zero_population = 0_int_p
        determ%flags = 1

        istart = 1
        iend = tot_walkers
        do i = 1, determ%sizes(iproc)
            call binary_search(walker_dets, dets_this_proc(:,i), istart, iend, hit, pos)
            if (.not. hit) then
                ! This deterministic state is not in walker_dets. Move all
                ! determinants with index pos or greater down one and insert
                ! this determinant with an initial sign of zero.
                walker_dets(:,pos+1:tot_walkers+1) = walker_dets(:,pos:tot_walkers)
                walker_population(:,pos+1:tot_walkers+1) = walker_population(:,pos:tot_walkers)
                walker_data(:,pos+1:tot_walkers+1) = walker_data(:,pos:tot_walkers)

                ! Insert a determinant with population zero into the walker arrays.
                call insert_new_walker(sys, pos, dets_this_proc(:,i), zero_population)

                tot_walkers = tot_walkers + 1
            end if

            ! Set this flag to specify a deterministic state.
            determ%flags(pos) = 0

            istart = pos + 1
            iend = tot_walkers
        end do

    end subroutine add_determ_dets_to_walker_dets

    function check_if_determ(ht, dets, f) result(is_determ)

        ! If f is a deterministic state, return is_determ as true.

        ! In:
        !    ht: Hash table used to index deterministic states.
        !    dets: Array containing all deterministic states from all processors.
        !    f: Determinant whose status is to be returned.

        use hashing, only: murmurhash_bit_string

        type(determ_hash_t), intent(in) :: ht
        integer(i0), intent(in) :: dets(:,:)
        integer(i0), intent(in) :: f(:)
        integer :: iz, hash
        logical :: is_determ

        is_determ = .false.

        hash = modulo(murmurhash_bit_string(f, size(f), ht%seed), ht%nhash) + 1
        ! Search the region of the hash table corresponding to this hash value.
        do iz = ht%hash_ptr(hash), ht%hash_ptr(hash+1)-1
            if (all(f == dets(:,ht%ind(iz)))) then
                is_determ = .true.
                return
            end if
        end do

    end function check_if_determ

    subroutine determ_projection(rng, spawn, determ)

        ! Apply the deterministic part of the FCIQMC projector to the
        ! amplitudes in the deterministic space. The corresponding spawned
        ! amplitudes are then added to the spawning array.

        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object to which deterministic spawning will occur.
        ! In:
        !    determ: Deterministic space being used.

        use calc, only: initiator_approximation
        use csr, only: csrpgemv_single_row
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use fciqmc_data, only: shift, tau, real_factor
        use omp_lib
        use parallel, only: nprocs, iproc, nthreads
        use spawn_data, only: spawn_t
        use spawning, only: add_spawned_particle, add_flagged_spawned_particle

        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn
        type(semi_stoch_t), intent(in) :: determ
        integer :: i, proc, row
        real(p) :: out_vec

        row = 0

        do proc = 0, nprocs-1
            ! To avoid an if statement which could be called a very large number
            ! of times, separate out states on this processor from those not.
            if (proc == iproc) then
                do i = 1, determ%sizes(proc)
                    row = row + 1
                    ! Perform the projetion.
                    call csrpgemv_single_row(determ%hamil, determ%vector, row, out_vec)
                    ! For states on this processor (proc == iproc), add the
                    ! contribution from the shift.
                    out_vec = -out_vec + shift(1)*determ%vector(i)
                    out_vec = out_vec*tau
                    call create_spawned_particle_determ(out_vec, row, proc)
                end do
            else
                do i = 1, determ%sizes(proc)
                    ! The same as above, but without the shift contribution.
                    row = row + 1
                    call csrpgemv_single_row(determ%hamil, determ%vector, row, out_vec)
                    out_vec = -out_vec*tau
                    call create_spawned_particle_determ(out_vec, row, proc)
                end do
            end if
        end do

        contains

            subroutine create_spawned_particle_determ(target_nspawn, idet, proc)

                ! Add a deterministic spawning to the spawning array. Before
                ! this can be done, the target population must be encoded as an
                ! integer using the reals encoding scheme.

                ! In:
                !    target_nspawn: The number of spawns we want to add to the
                !        spawning array, before the integer encoding has
                !        been performed.
                !    idet: The index of the determinant to be added in the
                !        determ%dets array.
                !    proc: The processor to which the determinant to be added
                !        belongs.

                real(p), intent(in) :: target_nspawn
                integer, intent(in) :: idet, proc

                integer(int_p) :: nspawn
                real(p) :: sgn, target_nspawn_scaled
#ifndef _OPENMP
                integer, parameter :: thread_id = 0
#else
                integer :: thread_id
                thread_id = omp_get_thread_num()
#endif

                ! Multiply target_nspawn by real_factor to allow it to be
                ! encoded as an integer.
                target_nspawn_scaled = target_nspawn*real_factor
                sgn = sign(1.0_p, target_nspawn)
                ! Stochastically round up or down, to encode as an integer.
                nspawn = int(abs(target_nspawn_scaled), int_p)
                if (abs(target_nspawn_scaled) - nspawn > get_rand_close_open(rng)) nspawn = nspawn + 1_int_p
                nspawn = nspawn*nint(sgn, int_p)
                if (initiator_approximation) then
                    call add_flagged_spawned_particle(determ%dets(:,idet), nspawn, 1, 0, proc, spawn)
                else
                    call add_spawned_particle(determ%dets(:,idet), nspawn, 1, proc, spawn)
                end if

            end subroutine create_spawned_particle_determ

    end subroutine determ_projection

    subroutine determ_projection_separate_annihil(determ)

        ! Perform the deterministic part of the projection. This is done here
        ! without adding deterministic spawnings to the spawned list, but
        ! rather by using an extra MPI call to perform the annihilation of
        ! these spawnings among themselves directly.

        ! In/Out:
        !    determ: Deterministic space being used. On input determ%vector
        !       should hold the amplitudes of deterministic states on this
        !       process. On output this will be overwritten by -tau*(Hv-Sv),
        !       where H is the Hamiltonian matrix, v is the vector of all
        !       deterministic amplitudes across all processes, tau is the
        !       timestep and S is the shift. It will hold the components of
        !       this vector belonging to this process only.
        !       determ%full_vector will hold v, the vector of all deterministic
        !       amplitudes across all processes, on output.

        use csr, only: csrpgemv
        use fciqmc_data, only: shift, tau
        use parallel

        type(semi_stoch_t), intent(inout) :: determ

        integer :: i, ierr
        integer :: send_counts(0:nprocs-1), receive_counts(0:nprocs-1)
        integer :: disps(0:nprocs-1)

        ! Create displacements used for MPI communication.
        disps(0) = 0
        do i = 1, nprocs-1
            disps(i) = disps(i-1) + determ%sizes(i-1)
        end do

#ifdef PARALLEL
        ! 'Stick together' the deterministic vectors from each process, on
        ! each process.
        call mpi_allgatherv(determ%vector, determ%sizes(iproc), mpi_preal, determ%full_vector, &
                             determ%sizes, disps, mpi_preal, MPI_COMM_WORLD, ierr)
#endif

        ! Perform the multiplication of the deterministic Hamiltonian on the
        ! deterministic vector.
        ! [review] - JSS: how is determ%full_vector set in serial?  I suspect ! it's not...
        call csrpgemv(.true., determ%hamil, determ%full_vector, determ%vector)

        ! We want the final vector to hold -tau*(Hv-Sv), where tau is the
        ! timestep, H is the determinstic Hamiltonian, v is the vector of
        ! deterministic amplitudes and S is the shift. We therefore need to
        ! apply a minus sign to the result of the above matrix multiplication,
        ! add the shift on and then multiply by tau.
        ! The vector v is still stored in determ%full_vector, even though
        ! determ%vector has been overwritten by the matrix multiplication.
        ! Explicitly do the full loop, to avoid creating a potentially large
        ! array on the stack.

        ! [review] - JSS: as an optimisation, why not do determ%vector = tau*shift(1)*determ%vector
        ! [review] - JSS: after the MPI_AllGatherV and beforethe csrpgemv call and add an
        ! [review] - JSS: optional argument to csrpgemv to avoid zeroing the input array?
        do i = 1, determ%sizes(iproc)
            determ%vector(i) = tau*(shift(1)*determ%full_vector(i+disps(iproc)) - determ%vector(i))
        end do

    end subroutine determ_projection_separate_annihil

    subroutine add_det_to_determ_space(determ_size_this_proc, dets_this_proc, spawn, f, check_proc)

        ! In/Out:
        !    determ_size_this_proc: Size of the deterministic space being
        !        created on this processor only.
        !    dets_this_proc: The deterministic states belonging to this
        !        processor.
        ! In:
        !    spawn: spawn_t object to which deterministic spawning will occur.
        !    f: det to be added to the determ object.
        !    check_proc: If true then first check if f belongs to this
        !        processor. If not then don't add it.

        use hashing, only: murmurhash_bit_string
        use parallel, only: iproc, nprocs
        use spawn_data, only: spawn_t
        use spawning, only: assign_particle_processor

        integer, intent(inout) :: determ_size_this_proc
        integer(i0), intent(inout) :: dets_this_proc(:,:)
        type(spawn_t), intent(in) :: spawn
        integer(i0), intent(in) :: f(:)
        logical, intent(in) :: check_proc

        integer :: proc, slot

        ! If check_proc is true then make sure that the determinant does belong
        ! to this processor. If it doesn't, don't add it and return.
        if (check_proc) then
            call assign_particle_processor(f, size(f), spawn%hash_seed, spawn%hash_shift, &
                                            spawn%move_freq, nprocs, proc, slot)
        else
            proc = iproc
        end if

        if (proc == iproc) then
            determ_size_this_proc = determ_size_this_proc + 1
            dets_this_proc(:, determ_size_this_proc) = f
        end if

    end subroutine add_det_to_determ_space

    subroutine create_high_pop_space(dets_this_proc, spawn, target_size, determ_size_this_proc)

        ! Find the most highly populated determinants in walker_dets and use
        ! these to define the deterministic space.

        ! In/Out:
        !    dets_this_proc: The deterministic states belonging to this
        !        processor in the final created deterministic space.
        ! In:
        !    spawn: spawn_t object to which deterministic spawning will occur.
        !    target_size: Size of deterministic space to use if possible. If
        !        not then use the largest space possible (all determinants in
        !        walker_dets).
        ! Out:
        !    determ_size_this_proc: Size of the deterministic space created,
        !        on this processor only.

        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: tot_walkers, walker_dets, walker_population
        use parallel
        use spawn_data, only: spawn_t

        integer(i0), intent(out) :: dets_this_proc(:,:)
        type(spawn_t), intent(in) :: spawn
        integer, intent(in) :: target_size
        integer, intent(out) :: determ_size_this_proc

        integer :: ndets, ndets_tot, determ_size
        integer :: all_ndets(0:nprocs-1), displs(0:nprocs)
        ! Temporary arrays used for finding the desired deterministic space.
        integer(i0), allocatable :: determ_dets(:,:)
        integer(int_p), allocatable :: determ_pops(:), all_determ_pops(:)
        integer, allocatable :: indices(:)
        integer :: i, ind_local, ierr
        logical :: on_this_proc

        ! If there are less determinants on this processor than the target
        ! number then obviously we need to consider a smaller number.
        ndets = min(target_size, tot_walkers)

        ! all_ndets will hold the values of ndets from each processor.
#ifdef PARALLEL
        call mpi_allgather(ndets, 1, mpi_integer, all_ndets, 1, mpi_integer, MPI_COMM_WORLD, ierr)
#else
        all_ndets = ndets
#endif
        ndets_tot = sum(all_ndets)

        ! Displacements used for MPI communication.
        displs(0) = 0
        do i = 1, nprocs-1
            displs(i) = displs(i-1) + all_ndets(i-1)
        end do
        displs(nprocs) = ndets_tot

        determ_size = target_size
        ! If there are fewer determinants on all processors than the target
        ! number then we will have to use a smaller deterministic space than
        ! requested.
        if (target_size > ndets_tot) determ_size = ndets_tot

        allocate(determ_dets(size(dets_this_proc, dim=1), ndets), stat=ierr)
        call check_allocate('determ_dets', size(determ_dets), ierr)
        allocate(determ_pops(ndets), stat=ierr)
        call check_allocate('determ_pops', ndets, ierr)
        allocate(all_determ_pops(ndets_tot), stat=ierr)
        call check_allocate('all_determ_pops', ndets_tot, ierr)
        allocate(indices(determ_size), stat=ierr)
        call check_allocate('indices', determ_size, ierr)

        ! In determ_dets and determ_pops, return the determinants and
        ! populations of the most populated determinants on this processor.
        call find_most_populated_dets(walker_dets, walker_population, tot_walkers, ndets, determ_dets, determ_pops)

        ! Create a joined list, all_determ_pops, of the most populated
        ! determinants from each processor.
#ifdef PARALLEL
        call mpi_allgatherv(determ_pops, ndets, mpi_pop_integer, all_determ_pops, all_ndets, &
                            displs(0:nprocs-1), mpi_pop_integer, MPI_COMM_WORLD, ierr)
#else
        ! In serial these arrays are of equal size.
        all_determ_pops = determ_pops
#endif

        ! In the array indices return a list of indices of the determ_size
        ! populations in all_determ_pops which are largest.
        call find_indices_of_most_populated_dets(ndets_tot, determ_size, all_determ_pops, indices)

        ! Zero before we fill in below.
        determ_size_this_proc = 0

        do i = 1, determ_size
            ! In determ_pops populations corresponding to determinants on
            ! processor 0 are all stored first, then processor 1, etc...
            ! We can use this to find which processor a particular determinant
            ! belongs to.
            on_this_proc = indices(i) > displs(iproc) .and. indices(i) < displs(iproc+1) + 1
            if (on_this_proc) then
                ! Find the index of this determinant in determ_dets, which only
                ! contains determinants on this processor.
                ind_local = indices(i) - displs(iproc)
                call add_det_to_determ_space(determ_size_this_proc, dets_this_proc, spawn, determ_dets(:,ind_local), .false.)
            end if
        end do

        deallocate(determ_dets, stat=ierr)
        call check_deallocate('determ_dets', ierr)
        deallocate(determ_pops, stat=ierr)
        call check_deallocate('determ_pops', ierr)
        deallocate(all_determ_pops, stat=ierr)
        call check_deallocate('all_determ_pops', ierr)
        deallocate(indices, stat=ierr)
        call check_deallocate('indices', ierr)

    end subroutine create_high_pop_space

    pure subroutine find_most_populated_dets(dets_in, pops_in, ndets_in, ndets_out, dets_out, pops_out)

        ! On output dets_out and pops_out hold the determinants and populations
        ! corresponding to the ndets_out most populated determinants, as
        ! specified by dets_in and pops_in.

        ! NOTE: It is assumed that size(dets_in,2) >= size(dets_out,2) and
        ! similarly for pops_in and pops_out.

        ! In:
        !    dets_in: List of determinants to consider.
        !    pops_in: Populations corresponding to the determinants in dets_in.
        !    ndets_in: The number of determinants in dets_in to consider.
        !    ndets_out: The number of determinants to keep when finding the
        !        most populated determinants. Anything beyond ndets_out in
        !        dets_out and pops_out may be garbage.
        ! Out:
        !    dets_out: List of most populated determinants.
        !    pops_out: Populations corresponding to the determinants in dets_out.

        ! dets_in(:,i) holds determinant i.
        integer(i0), intent(in) :: dets_in(:,:) 
        ! pops_in(j,i) holds the population of particle type j on determinant i.
        integer(int_p), intent(in) :: pops_in(:,:)
        integer, intent(in) :: ndets_in, ndets_out
        ! dets_out(:,i) holds determinant i.
        integer(i0), intent(out) :: dets_out(:,:)
        ! pops_in(j,i) holds the population of particle type j on determinant i.
        integer(int_p), intent(out) :: pops_out(:)

        integer :: i, j, min_ind
        integer(int_p) :: min_pop

        ! To start with add the first ndets_out determinants to the list.
        dets_out(:,1:ndets_out) = dets_in(:,1:ndets_out)
        pops_out(1) = sum(abs(pops_in(:,1)))
        min_pop = pops_out(1)
        min_ind = 1
        do i = 2, ndets_out
            pops_out(i) = sum(abs(pops_in(:,i)))
            if (pops_out(i) < min_pop) then
                min_pop = pops_out(i)
                min_ind = i
            end if
        end do

        ! Now loop over all remaining determinants and see if any have a larger
        ! amplitude than those in the list already.
        do i = ndets_out+1, ndets_in
            if (sum(abs(pops_in(:,i))) > min_pop) then
                ! Add this determinant to the list in the position of the
                ! previous smallest population.
                dets_out(:,min_ind) = dets_in(:,i)
                pops_out(min_ind) = sum(abs(pops_in(:,i)))
                ! Now find the position and value of the new smallest
                ! population.
                min_pop = pops_out(1)
                min_ind = 1
                do j = 2, ndets_out
                    if (pops_out(j) < min_pop) then
                        min_pop = pops_out(j)
                        min_ind = j
                    end if
                end do
            end if
        end do

    end subroutine find_most_populated_dets

    pure subroutine find_indices_of_most_populated_dets(npops_in, nind_out, pops, indices)

        ! On output indices will store the indices of the nind_out largest
        ! populations in pops.

        ! NOTE: If nind_out > npops_in then all indices beyond npops_in  will be
        ! returned as zero.

        ! In:
        !    npops_in: The number of populations in pops to consider.
        !    nind_out: The number of indices to keep when finding the
        !        most populated entries in pops.
        !    pops: List of populations to be maximised.
        ! Out:
        !    indices: List of the indices of the largest values in pops.

        integer, intent(in) :: npops_in
        integer, intent(in) :: nind_out
        integer(int_p), intent(in) :: pops(npops_in)
        integer, intent(out) :: indices(nind_out)

        integer :: i, j, min_ind
        integer(int_p) :: min_pop

        ! To start with just choose the first nind_out populations.
        indices(1) = 1
        min_pop = abs(pops(1))
        min_ind = 1
        do i = 2, min(nind_out, npops_in)
            indices(i) = i
            if (abs(pops(i)) < min_pop) then
                min_pop = abs(pops(i))
                min_ind = i
            end if
        end do

        ! If the number of requested indices is larger than the number of
        ! populations provided then return all further indices as zero.
        do i = npops_in+1, nind_out
            indices(i) = 0
        end do

        ! Now loop over all remaining populations and see if any are larger
        ! than those already in the list.
        do i = nind_out+1, npops_in
            if (abs(pops(i)) > min_pop) then
                ! Replace the old smallest index with this new index.
                indices(min_ind) = i
                ! Now find the position and value of the new smallest
                ! population.
                min_pop = abs(pops(indices(1)))
                min_ind = 1
                do j = 2, nind_out
                    if (abs(pops(indices(j))) < min_pop) then
                        min_pop = abs(pops(indices(j)))
                        min_ind = j
                    end if
                end do
            end if
        end do

    end subroutine find_indices_of_most_populated_dets

    subroutine read_determ_from_file(dets_this_proc, determ, spawn, sys, print_info)

        ! Use states read in from a HDF5 file to form the deterministic space.

        ! Out:
        !    dets_this_proc: The deterministic states belonging to this
        !        process. On entry it should be large enough to hold all states
        !        on this process.
        ! In/Out:
        !    determ: Deterministic space being used. On input determ%dets should
        !        not be allocated. On output determ%dets is deallocated.
        !        determ%sizes should be allocated on input, and on output
        !        contains the number of deterministic states on each process.
        ! In:
        !    spawn: spawn_t object to which deterministic spawning will occur.
        !    sys: system being studied.
        !    print_info: Should we print information to the screen?

#ifndef DISABLE_HDF5
        use checking, only: check_allocate, check_deallocate
        use hashing, only: murmurhash_bit_string
        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_read, hdf5_kinds_init
        use parallel
        use spawn_data, only: spawn_t
        use spawning, only: assign_particle_processor
        use system, only: sys_t
        use utils, only: get_unique_filename

        integer(i0), intent(out) :: dets_this_proc(:,:)
        type(semi_stoch_t), intent(inout) :: determ
        type(spawn_t), intent(in) :: spawn
        type(sys_t), intent(in) :: sys
        logical, intent(in) :: print_info

        type(hdf5_kinds_t) :: kinds
        integer(hid_t) :: file_id, dset_id, dspace_id
        character(255) :: filename
        integer :: i, proc, slot, ndeterm, ndeterm_this_proc, ierr
        integer :: displs(0:nprocs-1)
        integer(HSIZE_T) :: dims(2), maxdims(2)

        ! Read the deterministic states in on just the parent processor.
        if (parent) then
            call get_unique_filename("SEMI.STOCH", ".H5", .false., 0, filename)
            if (print_info) write(6,'(1X,"# Reading deterministic space states from",1X,a,".")') trim(filename)

            ! Initialise HDF5 and open file.
            call h5open_f(ierr)
            call hdf5_kinds_init(kinds)
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)

            ! Find how many determinants are in the file.
            call h5dopen_f(file_id, 'dets', dset_id, ierr)
            call h5dget_space_f(dset_id, dspace_id, ierr)
            call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, ierr)
            call h5dclose_f(dset_id, ierr)
            ! Number of determinants is the last index...
            ndeterm = dims(2)

            allocate(determ%dets(sys%basis%tensor_label_len, ndeterm), stat=ierr)
            call check_allocate('determ%dets', ndeterm*sys%basis%tensor_label_len, ierr)

            ! Perform the reading in of determinants to determ%dets.
            call hdf5_read(file_id, 'dets', kinds, shape(determ%dets), determ%dets)

            ! Close HDF5 file and HDF5.
            call h5fclose_f(file_id, ierr)
            call h5close_f(ierr)
        end if

#ifndef PARALLEL
        determ%sizes = ndeterm
        dets_this_proc(:,1:ndeterm) = determ%dets
#else
        if (parent) then
            ! Find how many determinants belong to each process.
            determ%sizes = 0
            do i = 1, ndeterm
                call assign_particle_processor(determ%dets(:,i), size(determ%dets,1), spawn%hash_seed, &
                                                spawn%hash_shift, spawn%move_freq, nprocs, proc, slot)
                determ%sizes(proc) = determ%sizes(proc) + 1
            end do

            ! Displacements used for MPI communication.
            displs(0) = 0
            do i = 1, nprocs-1
                displs(i) = displs(i-1) + determ%sizes(i-1)
            end do
        end if

        ! Send the number of determinants on a process to that process, from
        ! the root process.
        call mpi_scatter(determ%sizes, 1, mpi_integer, ndeterm_this_proc, 1, mpi_integer, root, &
                         MPI_COMM_WORLD, ierr)
        determ%sizes(iproc) = ndeterm_this_proc

        ! Send the determinants to their process.
        associate(tbl=>sys%basis%tensor_label_len)
            call mpi_scatterv(determ%dets, tbl*determ%sizes, tbl*displs, mpi_det_integer, &
                             dets_this_proc(:,1:ndeterm_this_proc), determ%sizes(iproc), &
                             mpi_det_integer, root, MPI_COMM_WORLD, ierr)
        end associate
#endif

        ! determ%dets is used to store the list of all deterministic states, but
        ! this is done again later for all processes, not just the parent. So
        ! for now deallocate determ%dets on the parent process.
        if (parent) then
            deallocate(determ%dets, stat=ierr)
            call check_deallocate('determ%dets', ierr)
        end if

#else
        use parallel
        use spawn_data, only: spawn_t
        use system, only: sys_t
        use errors, only: stop_all

        integer(i0), intent(out) :: dets_this_proc(:,:)
        type(semi_stoch_t), intent(inout) :: determ
        type(spawn_t), intent(in) :: spawn
        type(sys_t), intent(in) :: sys
        logical, intent(in) :: print_info

        call stop_all('read_determ_from_file', '# Not compiled with HDF5 support.  Cannot read semi-stochastic file.')
#endif

    end subroutine read_determ_from_file

    subroutine write_determ_to_file(determ, print_info)

        ! Write determinants stored in determ to a file.

        ! In:
        !    determ: Deterministic space being used.
        !    print_info: Should we print information to the screen?

#ifndef DISABLE_HDF5
        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_write, hdf5_kinds_init
        use calc, only: GLOBAL_META
        use utils, only: get_unique_filename

        type(semi_stoch_t), intent(in) :: determ
        logical, intent(in) :: print_info

        type(hdf5_kinds_t) :: kinds
        integer(hid_t) :: file_id
        character(255) :: filename
        integer :: ierr

        call get_unique_filename("SEMI.STOCH", ".H5", .true., 0, filename)
        if (print_info) write(6,'(1X,"# Writing deterministic space states to",1X,a,".")') trim(filename)

        ! Open HDF5 and create HDF5 kinds.
        call h5open_f(ierr)
        call hdf5_kinds_init(kinds)

        ! Open HDF5 file.
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)

        ! Write UUID so this can be linked to the main output, restart files, etc.
        call hdf5_write(file_id, 'uuid', GLOBAL_META%uuid)

        ! Write deterministic states to file.
        call hdf5_write(file_id, 'dets', kinds, shape(determ%dets), determ%dets)

        ! Close HDF5 file and HDF5.
        call h5fclose_f(file_id, ierr)
        call h5close_f(ierr)
#else
        use errors, only: warning

        type(semi_stoch_t), intent(in) :: determ
        logical, intent(in) :: print_info

        call warning('write_determ_to_file', '# Not compiled with HDF5 support.  Cannot write out semi-stochastic file.')
#endif

    end subroutine write_determ_to_file

end module semi_stoch
