module semi_stoch

! Code relating to the semi-stochastic algorithm. This includes initialisation
! code and also routines used throughout the simulation.

use const
use csr, only: csrp_t

implicit none

! Array to hold the indices of deterministic states in the dets array, accessed
! by calculating a hash value. This type is used by the semi_stoch_t type.
type determ_hash_t
    ! Seed used in the MurmurHash function, to calculate hash values.
    integer :: seed
    ! The number of unique hash values used, from 1 to nhash.
    integer :: nhash
    ! The indicies of the determinants in the dets array.
    integer, allocatable :: ind(:)
    ! hash_ptr(i) stores the index of the first index in the array ind which
    ! corresponds to a determinant with hash value i.
    ! This is similar to what is done int he CSR sparse matrix type (see
    ! csr.f90).
    integer, allocatable :: hash_ptr(:)
end type determ_hash_t

type semi_stoch_t
    ! The total number of deterministic states on all processes.
    integer :: tot_size
    ! sizes(i) holds the number of deterministic states belonging to process i.
    integer, allocatable :: sizes(:) ! (0:nproc-1)
    ! displs(i)+1 holds the index of the first deterministic state belonging
    ! to process i in the dets array. This is used in MPI communication.
    integer, allocatable :: displs(:) ! (0:nproc-1)
    ! The Hamiltonian in the deterministic space, stored in a sparse CSR form.
    type(csrp_t) :: hamil
    ! This array is used to store the values of amplitudes of deterministic
    ! states throughout a QMC calculation.
    real(p), allocatable :: vector(:) ! determ_sizes(iproc)
    ! dets stores the deterministic states across all processes.
    ! All states on process 0 are stored first, then process 1, etc...
    integer(i0), allocatable :: dets(:,:) ! (basis_length, tot_size)
    ! A hash table which allows the index of a determinant in dets to be found.
    ! This is done by calculating the hash value of the given determinant.
    type(determ_hash_t) :: hash_table
    ! Temporary space used for storing determinants during the creation of the
    ! deterministic space (before it is known how big the space is).
    integer(i0), allocatable :: temp_dets(:,:)
    ! Deterministic flags of states in the main list. If determ_flags(i) is
    ! equal to 0 then the corresponding state in position i of the main list is
    ! a deterministic state, else it is not.
    integer, allocatable :: flags(:)
end type semi_stoch_t

contains

    subroutine init_semi_stochastic(sys, spawn, determ, determ_type, determ_target_size)

        ! Create a semi_stoch_t object which holds all of the necessary
        ! information to perform a semi-stochastic calculation. The type of
        ! deterministic space is determined by determ_type.

        ! In/Out:
        !    determ: Deterministic space being used.
        ! In:
        !    sys: system being studied
        !    spawn: spawn_t object to which deterministic spawning will occur.
        !    determ_type: Integer parameter specifying which type of
        !        deterministic space to use.
        !    determ_target_size: A size of deterministic space to aim for. This
        !        is only necessary for particular deterministic spaces.

        use basis, only: total_basis_length
        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: walker_length
        use parallel
        use sort, only: qsort
        use spawn_data, only: spawn_t
        use system, only: sys_t
        use utils, only: int_fmt

        type(sys_t), intent(in) :: sys
        type(spawn_t), intent(in) :: spawn
        type(semi_stoch_t), intent(inout) :: determ
        integer, intent(in) :: determ_type
        integer, intent(in) :: determ_target_size

        integer :: i, ierr, determ_dets_mem
        logical :: print_info

        ! Only print information if the parent processor and if we are using a
        ! non-trivial deterministic space.
        print_info = parent .and. determ_type /= 0

        if (print_info) write(6,'(1X,a41)') 'Beginning semi-stochastic initialisation.'

        ! Create the temporary space for enumerating the deterinistic space and
        ! also the arrays to hold deterministic flags and space sizes and
        ! displacements.
        allocate(determ%temp_dets(total_basis_length, walker_length), stat=ierr)
        call check_allocate('determ%temp_dets', size(determ%temp_dets), ierr)
        allocate(determ%flags(walker_length), stat=ierr)
        call check_allocate('determ%flags', walker_length, ierr)
        allocate(determ%sizes(0:nprocs-1), stat=ierr)
        call check_allocate('determ%sizes', nprocs, ierr)
        allocate(determ%displs(0:nprocs-1), stat=ierr)
        call check_allocate('determ%displs', nprocs, ierr)

        determ%temp_dets = 0_i0
        determ%sizes = 0
        determ%displs = 0

        if (print_info) write(6,'(1X,a29)') 'Creating deterministic space.'

        ! If determ_type does not take one of the below values then an empty
        ! deterministic space will be used. This is the default behaviour
        ! (determ_type = 0).
        if (determ_type == 1) then
            call create_restart_space(determ, spawn, determ_target_size)
        end if

        ! Let each process hold the number of deterministic states on each process.
#ifdef PARALLEL
        call mpi_allgather(determ%sizes(iproc), 1, mpi_integer, determ%sizes, 1, mpi_integer, MPI_COMM_WORLD, ierr)
#endif
        determ%tot_size = sum(determ%sizes)

        if (print_info .and. nprocs > 1) then
            write(6,'(1X,a44,'//int_fmt(minval(determ%sizes),1)//')') &
                'Min deterministic space size on a processor:', minval(determ%sizes)
            write(6,'(1X,a44,'//int_fmt(maxval(determ%sizes),1)//')') &
                'Max deterministic space size on a processor:', maxval(determ%sizes)
            write(6,'(1X,a49,'//int_fmt(determ%tot_size,1)//')') &
                'Total deterministic space size on all processors:', determ%tot_size
        else if (print_info) then
            write(6,'(1X,a25,'//int_fmt(determ%tot_size,1)//')') 'Deterministic space size:', determ%tot_size
        end if

        ! Displacements used for MPI communication.
        determ%displs(0) = 0
        do i = 1, nprocs-1
            determ%displs(i) = determ%displs(i-1) + determ%sizes(i-1)
        end do

        ! Vector to hold deterministic amplitudes.
        allocate(determ%vector(determ%sizes(iproc)), stat=ierr)
        call check_allocate('determ%vector', determ%sizes(iproc), ierr)
        determ%vector = 0.0_p

        call qsort(determ%temp_dets, determ%sizes(iproc)) 

        ! Array to hold all deterministic states from all processes.
        ! The memory required in MB.
        determ_dets_mem = total_basis_length*determ%tot_size*i0_length/(8*10**6)
        if (print_info) write(6,'(1X,a58,'//int_fmt(determ_dets_mem,1)//')') &
            'Memory required per core to store deterministic dets (MB):', determ_dets_mem
        allocate(determ%dets(total_basis_length, determ%tot_size), stat=ierr)
        call check_allocate('determ%dets', size(determ%dets), ierr)

        ! Join and store all deterministic states from all processes.
#ifdef PARALLEL
        call mpi_allgatherv(determ%temp_dets(:,1:determ%sizes(iproc)), determ%sizes(iproc), mpi_det_integer, &
                            determ%dets, determ%sizes, determ%displs, mpi_det_integer, MPI_COMM_WORLD, ierr)
#else
        determ%dets = determ%temp_dets(:,1:determ%sizes(iproc))
#endif

        call create_determ_hash_table(determ, print_info)

        call create_determ_hamil(sys, determ, print_info)

        call add_determ_dets_to_walker_dets(sys, determ)

        ! We don't need this temporary space anymore.
        deallocate(determ%temp_dets, stat=ierr)
        call check_deallocate('determ%temp_dets', ierr)

        if (print_info) write(6,'(1X,a40,/)') 'Semi-stochastic initialisation complete.'

    end subroutine init_semi_stochastic

    subroutine create_determ_hash_table(determ, print_info)

        ! In/Out:
        !    determ: Deterministic space being used.
        ! In:
        !    print_info: Should we print information to the screen?

        use basis, only: total_basis_length
        use checking, only: check_allocate, check_deallocate
        use hashing, only: murmurhash_bit_string
        use utils, only: int_fmt

        type(semi_stoch_t), intent(inout) :: determ
        logical, intent(in) :: print_info

        integer :: i, iz, hash, mem_reqd, ierr
        integer, allocatable :: nclash(:)

        associate(ht => determ%hash_table)

            ht%seed = 37
            ! For now just let there be as many hash values as deterministic states.
            ht%nhash = determ%tot_size

            if (print_info) then 
                ! The memory required in MB.
                ! Two lists with length ht%nhash taking 4 bytes for each element.
                mem_reqd = 2*ht%nhash*4/10**6
                write(6,'(1X,a50,'//int_fmt(mem_reqd,1)//')') &
                    'Memory required per core to store hash table (MB):', mem_reqd
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
                hash = modulo(murmurhash_bit_string(determ%dets(:,i), total_basis_length, ht%seed), ht%nhash) + 1
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
                hash = modulo(murmurhash_bit_string(determ%dets(:,i), total_basis_length, ht%seed), ht%nhash) + 1
                iz = ht%hash_ptr(hash) + nclash(hash)
                ht%ind(iz) = i
                nclash(hash) = nclash(hash) + 1
            end do

            deallocate(nclash, stat=ierr)
            call check_deallocate('nclash', ierr)

        end associate

    end subroutine create_determ_hash_table

    subroutine create_determ_hamil(sys, determ, print_info)

        ! In/Out:
        !    determ: Deterministic space being used.
        ! In:
        !    sys: system being studied
        !    print_info: Should we print information to the screen?

        use checking, only: check_allocate
        use csr, only: init_csrp
        use fciqmc_data, only: H00
        use hamiltonian, only: get_hmatel
        use parallel
        use system, only: sys_t
        use utils, only: int_fmt

        type(sys_t), intent(in) :: sys
        type(semi_stoch_t), intent(inout) :: determ
        logical, intent(in) :: print_info

        integer :: i, j, nnz, imode, ierr
        integer :: mem_reqd, max_mem_reqd
        real(p) :: hmatel
        real :: t1, t2
        logical :: diag_elem

        if (print_info) write(6,'(1X,a72)') 'Counting number of non-zero deterministic Hamiltonian elements to store.'

        associate(hamil => determ%hamil)

            ! For imode = 1 count the number of non-zero elements, for imode = 2 store them.
            do imode = 1, 2
                if (imode == 1) call cpu_time(t1)
                nnz = 0
                ! Over all deterministic states on all processes (all rows).
                do i = 1, determ%tot_size
                    ! Over all deterministic states on this process (all columns).
                    do j = 1, determ%sizes(iproc)
                        hmatel = get_hmatel(sys, determ%dets(:,i), determ%temp_dets(:,j))
                        diag_elem = i == j + determ%displs(iproc)
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
                        write(6,'(1X,a39,1X,f10.2,a1)') 'Time taken to generate the Hamiltonian:', t2-t1, "s"
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
                        write(6,'(1X,a73,'//int_fmt(mem_reqd,1)//')') &
                            'Maximum memory required by a core for the deterministic Hamiltonian (MB):', mem_reqd
                        write(6,'(1X,a52)') 'The Hamiltonian will now be recalculated and stored.'
                    end if

                    ! Allocate the CSR Hamiltonian arrays.
                    call init_csrp(hamil, determ%tot_size, nnz)
                    hamil%row_ptr(1:determ%tot_size) = 0
                end if
            end do

        end associate

    end subroutine create_determ_hamil

    subroutine add_determ_dets_to_walker_dets(sys, determ)

        ! Also set the deterministic flags of any deterministic states already
        ! in walker_dets, and add deterministic data to walker_populations and
        ! walker_data. All deterministic states not already in walker_dets are
        ! set to have an initial sign of zero.

        ! In/Out:
        !    determ: Deterministic space being used.
        ! In:
        !    sys: system being studied

        use basis, only: basis_length
        use calc, only: dmqmc_calc, hfs_fciqmc_calc, trial_function
        use calc, only: neel_singlet, doing_calc
        use fciqmc_data, only: walker_dets, walker_population, walker_data
        use fciqmc_data, only: sampling_size, H00, tot_walkers, replica_tricks
        use heisenberg_estimators, only: neel_singlet_data
        use hfs_data, only: O00
        use parallel, only: iproc
        use proc_pointers, only: sc0_ptr, op0_ptr
        use search, only: binary_search
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(semi_stoch_t), intent(inout) :: determ

        integer :: i, j, k, istart, iend, pos
        logical :: hit

        determ%flags = 1

        istart = 1
        iend = tot_walkers
        do i = 1, determ%sizes(iproc)
            call binary_search(walker_dets, determ%temp_dets(:,i), istart, iend, hit, pos)
            if (hit) then
                ! This deterministic state is already in walker_dets. We simply
                ! need to set the deterministic flag.
                determ%flags(pos) = 0
            else
                ! This deterministic state is not in walker_dets. Move all
                ! determinants with index pos or greater down one and insert
                ! this determinant with an initial sign or zero.
                walker_dets(:,pos:tot_walkers) = walker_dets(:,pos+1:tot_walkers+1)
                walker_population(:,pos:tot_walkers) = walker_population(:,pos+1:tot_walkers+1)
                walker_data(:,pos:tot_walkers) = walker_data(:,pos+1:tot_walkers+1)

                walker_dets(:,pos) = determ%temp_dets(:,i)
                walker_population(:,pos) = 0_int_p
                if (.not. doing_calc(dmqmc_calc)) walker_data(1,pos) = sc0_ptr(sys, determ%temp_dets(:,i)) - H00
                if (trial_function == neel_singlet) walker_data(sampling_size+1:sampling_size+2,pos) = &
                    neel_singlet_data(walker_dets(:,pos))
                if (doing_calc(hfs_fciqmc_calc)) then
                    ! Set walker_data(2:,k) = <D_i|O|D_i> - <D_0|O|D_0>.
                    walker_data(2,pos) = op0_ptr(sys, walker_dets(:,pos)) - O00
                else if (doing_calc(dmqmc_calc)) then
                    ! Set the energy to be the average of the two induvidual energies.
                    walker_data(1,pos) = (walker_data(1,pos) + &
                        sc0_ptr(sys, walker_dets((basis_length+1):(2*basis_length),pos)) - H00)/2
                    if (replica_tricks) then
                        walker_data(2:sampling_size,pos) = walker_data(1,pos)
                    end if
                end if

                tot_walkers = tot_walkers + 1
            end if
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

        use basis, only: total_basis_length
        use hashing, only: murmurhash_bit_string

        type(determ_hash_t), intent(in) :: ht
        integer(i0), intent(in) :: dets(:,:)
        integer(i0), intent(in) :: f(total_basis_length)
        integer :: iz, hash
        logical :: is_determ

        is_determ = .false.

        hash = modulo(murmurhash_bit_string(f, total_basis_length, ht%seed), ht%nhash) + 1
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
        use csr, only: csrp_row
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use fciqmc_data, only: shift, tau, real_factor
        use omp_lib
        use parallel, only: nprocs, iproc, nthreads
        use spawn_data, only: spawn_t

        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn
        type(semi_stoch_t), intent(in) :: determ
        integer :: i, proc, row
        real(p) :: out_vec, sgn
        integer(int_p) :: nspawn
#ifndef _OPENMP
        integer, parameter :: thread_id = 0
#else
        integer :: thread_id
        thread_id = omp_get_thread_num()
#endif

        row = 0

        do proc = 0, nprocs-1
            ! To avoid an if statement which could be called a very large number
            ! of times, separate out states on this processor from those not.
            if (proc == iproc) then
                do i = 1, determ%sizes(proc)
                    row = row + 1
                    ! Perform the projetion.
                    call csrp_row(determ%hamil, determ%vector, out_vec, row)
                    ! For states on this processor (proc == iproc), add the
                    ! contribution from the shift.
                    out_vec = -out_vec + shift(1)*determ%vector(i)
                    ! Multiply by the timestep, and also by real_factor, to
                    ! allow the amplitude to be encoded as an integer.
                    out_vec = out_vec*tau*real_factor
                    sgn = sign(1.0_p, out_vec)
                    ! Stochastically round up or down, to encode as an integer.
                    nspawn = int(abs(out_vec), int_p)
                    if (abs(out_vec) - nspawn > get_rand_close_open(rng)) nspawn = nspawn + 1_int_p
                    nspawn = nspawn*nint(sgn, int_p)
                    ! Upate the spawning slot...
                    spawn%head(thread_id,proc) = spawn%head(thread_id,proc) + nthreads
                    ! ...and finally add the state to the spawning array.
                    spawn%sdata(:,spawn%head(thread_id,proc)) = 0_int_s
                    spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,proc)) = int(determ%dets(:,row), int_s)
                    spawn%sdata(spawn%bit_str_len+1,spawn%head(thread_id,proc)) = int(nspawn, int_s)
                    ! Spawning has occurred from a deterministic state, an
                    ! initiator by definition.
                    if (initiator_approximation) spawn%sdata(spawn%flag_indx,spawn%head(thread_id,proc)) = 0_int_s
                end do
            else
                do i = 1, determ%sizes(proc)
                    ! The same as above, but without the shift contribution.
                    row = row + 1
                    call csrp_row(determ%hamil, determ%vector, out_vec, row)
                    out_vec = -out_vec*tau*real_factor
                    sgn = sign(1.0_p, out_vec)
                    nspawn = int(abs(out_vec), int_p)
                    if (abs(out_vec) - nspawn > get_rand_close_open(rng)) nspawn = nspawn + 1
                    nspawn = nspawn*nint(sgn, int_p)
                    spawn%head(thread_id,proc) = spawn%head(thread_id,proc) + nthreads
                    spawn%sdata(:,spawn%head(thread_id,proc)) = 0_int_s
                    spawn%sdata(:spawn%bit_str_len,spawn%head(thread_id,proc)) = int(determ%dets(:,row), int_s)
                    spawn%sdata(spawn%bit_str_len+1,spawn%head(thread_id,proc)) = int(nspawn, int_s)
                    if (initiator_approximation) spawn%sdata(spawn%flag_indx,spawn%head(thread_id,proc)) = 0_int_s
                end do
            end if
        end do

    end subroutine determ_projection

    subroutine add_det_to_determ_space(determ, spawn, f, check_proc)

        ! In/Out:
        !    determ: Deterministic space being used.
        ! In:
        !    spawn: spawn_t object to which deterministic spawning will occur.
        !    f: det to be added to the determ object.
        !    check_proc: If true then first check if f belongs to this
        !        processor. If not then don't add it.

        use basis, only: total_basis_length
        use hashing, only: murmurhash_bit_string
        use parallel, only: iproc, nprocs
        use spawn_data, only: spawn_t

        type(semi_stoch_t), intent(inout) :: determ
        type(spawn_t), intent(in) :: spawn
        integer(i0), intent(in) :: f(total_basis_length)
        logical, intent(in) :: check_proc

        integer :: proc

        ! If check_proc is true then make sure that the determinant does belong
        ! to this processor. If it doesn't, don't add it and return.
        if (check_proc) then
            proc = modulo(murmurhash_bit_string(f, total_basis_length, spawn%hash_seed), nprocs)
            if (proc /= iproc) return
        end if

        determ%sizes(iproc) = determ%sizes(iproc) + 1

        determ%temp_dets(:, determ%sizes(iproc)) = f

    end subroutine add_det_to_determ_space

    subroutine create_restart_space(determ, spawn, target_size)

        ! Find the most highly populated determinants in walker_dets and use
        ! these to define the deterministic space. When using this routine
        ! the restart option should have been used, although it is not required.

        ! In/Out:
        !    determ: Deterministic space being used.
        ! In:
        !    spawn: spawn_t object to which deterministic spawning will occur.
        !    target_size: Size of deterministic space to use if possible. If
        !        not then use the largest space possible (all determinants from
        !        the restart file)).

        use basis, only: total_basis_length
        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: tot_walkers, walker_dets, walker_population
        use parallel
        use spawn_data, only: spawn_t

        type(semi_stoch_t), intent(inout) :: determ
        type(spawn_t), intent(in) :: spawn
        integer, intent(in) :: target_size

        integer :: ndets, ndets_tot, determ_size
        integer :: all_ndets(0:nprocs-1), displs(0:nprocs)
        integer(i0), allocatable :: determ_dets(:,:)
        integer(int_p), allocatable :: determ_pops(:), all_determ_pops(:)
        integer, allocatable :: indices(:)
        integer :: i, ind_local, ierr
        logical :: on_this_proc

        ! If there are less determinants on this processor than the target
        ! number then obviously we need to consider a smaller number.
        ndets = min(target_size, tot_walkers)

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

        allocate(determ_dets(total_basis_length, ndets), stat=ierr)
        call check_allocate('determ_dets', total_basis_length*ndets, ierr)
        allocate(determ_pops(ndets), stat=ierr)
        call check_allocate('determ_pops', ndets, ierr)
        allocate(all_determ_pops(ndets_tot), stat=ierr)
        call check_allocate('all_determ_pops', ndets_tot, ierr)
        allocate(indices(determ_size), stat=ierr)
        call check_allocate('indices', determ_size, ierr)

        ! In determ_dets and determ_pops, return the determinants and
        ! populations of the most populated determinants on this processor.
        call find_most_populated_dets(walker_dets, walker_population, tot_walkers, determ_dets, determ_pops, ndets)

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
        call find_indices_of_most_populated_dets(all_determ_pops, ndets_tot, indices, determ_size)

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
                call add_det_to_determ_space(determ, spawn, determ_dets(:,ind_local), .false.)
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

    end subroutine create_restart_space

    subroutine find_most_populated_dets(dets_in, pops_in, ndets_in, dets_out, pops_out, ndets_out)

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

        integer(i0), intent(in) :: dets_in(:,:) 
        integer(int_p), intent(in) :: pops_in(:,:)
        integer, intent(in) :: ndets_in, ndets_out
        integer(i0), intent(out) :: dets_out(:,:)
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

    subroutine find_indices_of_most_populated_dets(pops, npops_in, indices, nind_out)

        ! On output indices will store the indices of the nind_out largest
        ! populations in pops.

        ! NOTE 1: It is assumed that the populations in pops are all
        ! non-negative. The absolute values of populations are not taken in
        ! this routine.

        ! NOTE 2: It is assumed that npops_in >= nind_out. It is up to the
        ! programmer to ensure this is true and there will likely be errors if
        ! it isn't.

        ! In:
        !    pops: List of populations to be maximised.
        !    npops_in: The number of populations in pops to consider.
        !    nind_out: The number of indices to keep when finding the
        !        most populated entries in pops. Anything beyond nind_out in
        !        indices may be garbage.
        ! Out:
        !    indices: List of the indices of the largest values in pops.

        integer(int_p), intent(in) :: pops(:)
        integer, intent(in) :: npops_in
        integer, intent(out) :: indices(:)
        integer, intent(in) :: nind_out

        integer :: i, j, min_ind
        integer(int_p) :: min_pop

        ! To start with just choose the first nind_out populations.
        indices(1) = 1
        min_pop = pops(1)
        min_ind = 1
        do i = 2, nind_out
            indices(i) = i
            if (pops(i) < min_pop) then
                min_pop = pops(i)
                min_ind = i
            end if
        end do

        ! Now loop over all remaining populations and see if any are larger
        ! than those already in the list.
        do i = nind_out+1, npops_in
            if (pops(i) > min_pop) then
                ! Replace the old smallest index with this new index.
                indices(min_ind) = i
                ! Now find the position and value of the new smallest
                ! population.
                min_pop = pops(indices(1))
                min_ind = 1
                do j = 2, nind_out
                    if (pops(indices(j)) < min_pop) then
                        min_pop = pops(indices(j))
                        min_ind = j
                    end if
                end do
            end if
        end do

    end subroutine find_indices_of_most_populated_dets

end module semi_stoch
