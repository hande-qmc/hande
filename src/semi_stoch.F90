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
end type semi_stoch_t

type(semi_stoch_t) :: qmc_determ

contains

    subroutine init_semi_stochastic(determ_type, determ_target_size, determ)

        ! Create a semi_stoch_t object which holds all of the necessary
        ! information to perform a semi-stochastic calculation. The type of
        ! deterministic space is determined by determ_type.

        use basis, only: total_basis_length
        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: walker_length
        use parallel

        integer, intent(in) :: determ_type
        integer, intent(in) :: determ_target_size
        type(semi_stoch_t), intent(inout) :: determ

        integer :: i, ierr

        ! Create the temporary space for enumerating the deterinistic space and
        ! also the arrays to hold the space sizes and displacements.
        allocate(determ%temp_dets(total_basis_length, walker_length), stat=ierr)
        call check_allocate('determ%temp_dets', size(determ%temp_dets), ierr)
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

        ! Array to hold all deterministic states from all processes.
        allocate(determ%dets(total_basis_length, determ%tot_size), stat=ierr)
        call check_allocate('determ%dets', size(determ%dets), ierr)
        call mpi_allgatherv(determ%temp_dets(:,1:determ%sizes(iproc)), determ%sizes(iproc), mpi_det_integer, &
                            determ%dets, determ%tot_size, determ%displs, mpi_det_integer, MPI_COMM_WORLD)

        call create_determ_hash_table(determ)

        ! We don't need this temporary space anymore.
        deallocate(determ%temp_dets, stat=ierr)
        call check_deallocate('determ%temp_dets', ierr)

    end subroutine init_semi_stochastic

    subroutine create_determ_hash_table(determ)

        use basis, only: total_basis_length
        use checking, only: check_allocate, check_deallocate
        use hashing, only: murmurhash_bit_string

        type(semi_stoch_t) :: determ

        integer :: i, iz, hash, ierr
        integer, allocatable :: nclash(:)

        associate(ht => determ%hash_table)

            ht%seed = 37
            ! For now just let there be as many hash values as deterministic states.
            ht%nhash = determ%tot_size

            allocate(ht%ind(determ%tot_size), stat=ierr) 
            call check_allocate('determ%hash_table%ind', determ%tot_size, ierr)
            allocate(ht%hash_ptr(ht%nhash), stat=ierr) 
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

end module semi_stoch
