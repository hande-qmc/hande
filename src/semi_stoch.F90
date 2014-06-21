module semi_stoch

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
    ! equal to 1 then the corresponding state in position i of the main list is
    ! a deterministic state, else it is not.
    integer, allocatable :: flags(:)
end type semi_stoch_t

contains

    subroutine init_semi_stochastic(determ_type, determ_target_size, sys, determ)

        ! Create a semi_stoch_t object which holds all of the necessary
        ! information to perform a semi-stochastic calculation. The type of
        ! deterministic space is determined by determ_type.

        use basis, only: total_basis_length
        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: walker_length
        use parallel
        use sort, only: qsort
        use system, only: sys_t

        integer, intent(in) :: determ_type
        integer, intent(in) :: determ_target_size
        type(sys_t), intent(in) :: sys
        type(semi_stoch_t), intent(inout) :: determ

        integer :: i, ierr

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

        ! Code to be added here to generate deterministic space.
        
        ! Let each process hold the number of deterministic states on each process.
        call mpi_allgather(determ%sizes(iproc), 1, mpi_integer, determ%sizes, nprocs, mpi_integer, MPI_COMM_WORLD)
        determ%tot_size = sum(determ%sizes)

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
        allocate(determ%dets(total_basis_length, determ%tot_size), stat=ierr)
        call check_allocate('determ%dets', size(determ%dets), ierr)
        call mpi_allgatherv(determ%temp_dets(:,1:determ%sizes(iproc)), determ%sizes(iproc), mpi_det_integer, &
                            determ%dets, determ%tot_size, determ%displs, mpi_det_integer, MPI_COMM_WORLD)

        call create_determ_hash_table(determ)

        call create_determ_hamil(sys, determ)

        call add_determ_dets_to_walker_dets(sys, determ)

        ! We don't need this temporary space anymore.
        deallocate(determ%temp_dets, stat=ierr)
        call check_deallocate('determ%temp_dets', ierr)

    end subroutine init_semi_stochastic

    subroutine create_determ_hash_table(determ)

        use basis, only: total_basis_length
        use checking, only: check_allocate, check_deallocate
        use hashing, only: murmurhash_bit_string

        type(semi_stoch_t), intent(inout) :: determ

        integer :: i, iz, hash, ierr
        integer, allocatable :: nclash(:)

        associate(ht => determ%hash_table)

            ht%seed = 37
            ! For now just let there be as many hash values as deterministic states.
            ht%nhash = determ%tot_size

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
                ht%hash_ptr(i) = ht%hash_ptr(i-1) + nclash(i)
            end do
            ht%hash_ptr(ht%nhash+1) = determ%tot_size + 1

            ! Now loop over all states again and this time fill in the hash table.
            nclash = 0
            do i = 1, determ%tot_size
                hash = modulo(murmurhash_bit_string(determ%dets(:,i), total_basis_length, ht%seed), ht%nhash) + 1
                iz = ht%hash_ptr(i) + nclash(hash)
                ht%ind(iz) = i
                nclash(hash) = nclash(hash) + 1
            end do

            deallocate(nclash, stat=ierr)
            call check_deallocate('nclash', ierr)

        end associate

    end subroutine create_determ_hash_table

    subroutine create_determ_hamil(sys, determ)

        use checking, only: check_allocate
        use csr, only: init_csrp
        use hamiltonian, only: get_hmatel
        use parallel, only: iproc
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(semi_stoch_t), intent(inout) :: determ

        integer :: i, j, nnz, imode, ierr
        real(p) :: hmatel

        associate(hamil => determ%hamil)

            ! For imode = 1 count the number of non-zero elements, for imode = 2 store them.
            do imode = 1, 2
                ! Over all deterministic states on all processes (all rows).
                do i = 1, determ%tot_size
                    ! Over all deterministic states on this process (all columns).
                    do j = 1, determ%sizes(iproc)
                        hmatel = get_hmatel(sys, determ%dets(:,i), determ%temp_dets(:,j))
                        if (abs(hmatel) > depsilon) then
                            nnz = nnz + 1
                            if (imode == 2) then
                                hamil%mat(nnz) = hmatel
                                hamil%col_ind(nnz) = j
                                if (hamil%row_ptr(i) == 0) hamil%row_ptr(i) = nnz
                            end if
                        end if
                    end do
                end do

                ! Allocate the CSR type components.
                if (imode == 1) then
                    call init_csrp(hamil, determ%tot_size, nnz)
                    hamil%row_ptr = 0
                end if

            end do

        end associate

    end subroutine create_determ_hamil

    subroutine add_determ_dets_to_walker_dets(sys, determ)

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

        determ%flags = 0

        istart = 1
        iend = tot_walkers
        do i = 1, determ%sizes(iproc)
            call binary_search(walker_dets, determ%temp_dets(:,i), istart, iend, hit, pos)
            if (hit) then
                ! This deterministic state is already in walker_dets. We simply
                ! need to set the deterministic flag.
                determ%flags(pos) = 1
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

        use basis, only: total_basis_length
        use hashing, only: murmurhash_bit_string

        type(determ_hash_t), intent(in) :: ht
        integer(i0), intent(in) :: dets(:,:)
        integer(i0), intent(in) :: f(total_basis_length)
        integer :: iz, hash
        logical :: is_determ

        is_determ = .false.

        hash = modulo(murmurhash_bit_string(f, total_basis_length, ht%seed), ht%nhash) + 1
        do iz = ht%hash_ptr(hash), ht%hash_ptr(hash+1)-1
            if (all(f == dets(:,ht%ind(iz)))) then
                is_determ = .true.
                return
            end if
        end do

    end function check_if_determ

end module semi_stoch
