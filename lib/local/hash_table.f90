module hash_table

    ! This follows the hash table storage scheme described in "Linear-scaling
    ! and parallelizable algorithms for stochastic quantum chemistry", GH Booth,
    ! SD Smart and A Alavi, http://arxiv.org/abs/1305.6981.

    ! The idea revolves around the following key points:

    ! * want the data stored to be as contiguous as possible, so iterating
    !   through the data stored does not require iterating through the entire
    !   data array and is also cache friendly (because we don't jump around the
    !   data arrays).  For QMC this means that propogation remains cheap.
    ! * want data to be stored in the same slot in the data array for the
    !   lifetime it is stored, so lookups are O(1).  For QMC this reduces the
    !   cost of annihilation from an O(N_d N_s) to an O(N_s) operation, where
    !   N_s and N_d and the number of spawned psips and N_d the number of
    !   occupied determinants.

    ! In order to achieve this we:

    ! * hash the data label being stored and take the mod of it to reduce it to
    !   within a (comparitively) small range.
    ! * store the index of the data label in the hash table, which is shaped to
    !   allow for hash collisions.
    ! * store the data label and associated information at that index.
    ! * in the case of hash collisions, we must iterate through the list of
    !   collisions and compare to the stored data labels in order to find the
    !   desired entry.
    ! * use a circular array to store a (small) number of free slots which stop
    !   interwoven with stored data (due to entries being deleted).  These free
    !   slots are re-used preferentially in order to keep the data stored
    !   as contiguously as possible.

    use const
    use, intrinsic :: iso_c_binding, only: c_int

    implicit none

    type hash_table_t

        ! Number of slots in the hash table.
        integer :: nslots
        ! Maximum number of entries (ie hash collisions) per slot in the hash
        ! table.
        integer :: nentries

        ! Number of elements in the array forming the data element to be hashed.
        integer :: data_len
        ! Number of elements in the payload/info array, where each payload array
        ! is associated with a given item of hashed data.
        integer :: payload_len

        ! The index for the location in data_label and paload for a given data
        ! item, d, is table(1:,mod(hash(d),nslots)), where hash is the hash function.
        ! In the event of a hash collision, one must search over the number of
        ! entries (given by the 0-th element) to find the desired item.
        integer, allocatable :: table(:,:) ! (0:entries, 0:nslots-1)

        ! data_label(:,i) is the i-th (unique) data item being stored.
        integer(int_s), pointer :: data_label(:,:) => null() ! (:,nentries*nslots)
        ! Using an internal array for data label (ie allocation) or pointing to
        ! an external array?
        logical, private :: extern_data_label = .false.

        ! payload(:,i) i$ the 'payload' or information associated with
        ! the i-th data item in data_label.
        integer(int_p), pointer :: payload(:,:) => null() ! (:,nentries*nslots)
        ! Using an internal array for payload (ie allocation) or pointing to
        ! an external array?
        logical, private :: extern_payload = .false.

        ! We usually want to iterate over the data items, so it is best if they
        ! are (as far as possible) contiguous.
        ! If not, we can insert an entry into the empty entry.  Typically (we
        ! hope!) we have more room than the number of items we're actually
        ! storing.  If there are no slots in the occupied section of data_label
        ! and payload, then we use the next free entry at the head of the data
        ! structure.

        ! Highest occupied entry in the data_label and payload arrays.
        integer :: head

        ! Maximum number of free slots we might wish to access before iterating
        ! over data_label again (at which point we can easily update the list of
        ! free entries).
        integer :: max_free_reqd
        ! The index of the element in free_entries which grants the next free
        ! entry in data_label and payload arrays.
        integer :: next_free_entry
        ! Circular array.  Non-zero elements correspond to free entries in the
        ! data_label and corresponding payload arrays.
        integer, allocatable :: free_entries(:) ! (0:max_free_reqd-1)

        ! Seed for the hash function.
        integer(c_int) :: seed

    end type hash_table_t

    type hash_table_pos_t
        ! 'slot' in the hash table (hash_table_t%table), ie the hash (mod number
        ! of slots in the hash table) of the data label.
        integer :: islot
        ! 'entry' in the hash table (hash_table_t%table), ie the index of the
        ! data label amongst all known data labels which share the same hash.
        integer :: ientry
        ! index of the data label in the hash_table_t%data_label and
        ! hash_table_t%payload arrays, ie the value stored by
        ! hash_table_t%table(ientry, islot).
        integer :: indx
    end type hash_table_pos_t

    ! Error code if hash table is full.
    integer, parameter :: HT_ERR_FULL = 1
    ! Error code if there have been more collisions than the hash table can
    ! handle (i.e. more than nentries need to share the same slot in the hash
    ! table).
    integer, parameter :: HT_ERR_COLLISIONS = 2

    contains

!--- Allocation/deallocation ---

        subroutine alloc_hash_table(nslots, nentries, data_len, payload_len, max_free_reqd, seed, ht, data_label, payload)

            ! Allocate hash table and initialise storage.

            ! In:
            !    nslots, nentries, data_len, payload_len, max_free_reqd, seed: see
            !        descriptions in hash_table_t.
            !    data_label (optional), payload (optional): external arrays to
            !        use for the data label and payload arrays (see descriptions
            !        in hash_table_t) rather than internal stores.
            ! Out:
            !    ht: hash table.  On output all array attributes are allocated
            !        and all scalar attributes are set to their corresponding
            !        input arguments.  The hash table (ht%table) is set to
            !        0 (see reset_hash_table), ready for use.

            use checking, only: check_allocate

            integer, intent(in) :: nslots, nentries, data_len, payload_len, max_free_reqd, seed
            integer(int_s), intent(in), optional, target :: data_label(:,:)
            integer(int_p), intent(in), optional, target :: payload(:,:)
            type(hash_table_t), intent(out) :: ht

            integer :: ierr

            ht%nslots = nslots
            ht%nentries = nentries
            ht%data_len = data_len
            ht%payload_len = payload_len
            ht%max_free_reqd = max_free_reqd
            ht%seed = seed

            allocate(ht%table(0:ht%nentries,ht%nslots), stat=ierr)
            call check_allocate('ht%table', size(ht%table), ierr)
            if (present(data_label)) then
                ht%data_label => data_label
                ht%extern_data_label = .true.
            else
                allocate(ht%data_label(ht%data_len,ht%nentries*ht%nslots))
                call check_allocate('ht%data_label', size(ht%data_label), ierr)
            end if
            if (present(payload)) then
                ht%payload => payload
                ht%extern_payload = .true.
            else
                allocate(ht%payload(ht%payload_len,ht%nentries*ht%nslots), stat=ierr)
                call check_allocate('ht%payload', size(ht%payload), ierr)
            end if
            allocate(ht%free_entries(0:ht%max_free_reqd-1), stat=ierr)
            call check_allocate('ht%free_entries', size(ht%free_entries), ierr)

            call reset_hash_table(ht)

        end subroutine alloc_hash_table

        subroutine free_hash_table(ht)

            ! Deallocate hash table and set scalar attributes to 0.

            ! In/Out:
            !    ht: hash table.  On output all array attributes are deallocated/
            !        and all scalar attributes are set to 0.

            use checking, only: check_deallocate

            type(hash_table_t), intent(inout) :: ht

            integer :: ierr

            if (allocated(ht%table)) then
                deallocate(ht%table)
                call check_deallocate('ht%table', ierr)
            end if
            if (associated(ht%data_label)) then
                if (ht%extern_data_label) then
                    ht%data_label => null()
                else
                    deallocate(ht%data_label)
                    call check_deallocate('ht%data_label', ierr)
                end if
            end if
            if (associated(ht%payload)) then
                if (ht%extern_payload) then
                    ht%payload => null()
                else
                    deallocate(ht%payload)
                    call check_deallocate('ht%payload', ierr)
                end if
            end if
            if (allocated(ht%free_entries)) then
                deallocate(ht%free_entries)
                call check_deallocate('ht%free_entries', ierr)
            end if

            ht%nslots = 0
            ht%nentries = 0
            ht%data_len = 0
            ht%payload_len = 0
            ht%head = 0
            ht%max_free_reqd = 0
            ht%seed = 0
            ht%next_free_entry = 0
            ht%extern_data_label = .false.
            ht%extern_payload = .false.

        end subroutine free_hash_table

        elemental subroutine reset_hash_table(ht)

            ! Reset hash table.

            ! In/Out:
            !    ht: hash table.  On output the table, head, next_free_entry
            !        and free_entries attributes are all set to 0.  This means
            !        that no data is stored in the hash table and new entries 
            !        are added by consecutively by incrementing head.

            type(hash_table_t), intent(inout) :: ht

            integer :: i

            ht%table(0,:) = 0
            ht%head = 0
            ht%next_free_entry = 0
            ht%free_entries = 0

        end subroutine reset_hash_table

!--- Remove/add known free entries in occupied part of data_load/payload arrays ---

        elemental subroutine take_hash_table_free_entry(ht, indx, err_code)

            ! Get the index of the next entry of the data arrays to be used.

            ! In/Out:
            !    ht: hash table.  On output, the entry given by indx is
            !        registered as being used (i.e. removed from the
            !        ht%free_entries or ht%head is incremented).
            ! Out:
            !    indx: index of free entry in ht%data_label and ht%payload to be
            !        used.  The index is taken from the free_entries_array if
            !        there is a index available, otherwise the index is the
            !        entry directly above the current head (top) entry of
            !        ht%data_label and ht%payload arrays.
            !    err_code: error code.  Non-zero if an error is encountered.
            !        See error codes defined at module-level.

            type(hash_table_t), intent(inout) :: ht
            integer, intent(out) :: indx, err_code

            err_code = 0
            indx = 0
            if (ht%max_free_reqd > 0) then
                indx = ht%free_entries(ht%next_free_entry)
                ht%free_entries(ht%next_free_entry) = 0
                ht%next_free_entry = modulo(ht%next_free_entry-1,ht%max_free_reqd)
            end if

            if (indx == 0) then
                ! Have run out of free slots that we know about in the
                ! used part of the storage.  Now just add to the head of the
                ! occupied space...
                ht%head = ht%head + 1
                indx = ht%head
                if (ht%head > ubound(ht%data_label,dim=2)) err_code = HT_ERR_FULL
            end if

        end subroutine take_hash_table_free_entry

        elemental subroutine register_hash_table_free_entry(ht, indx)

            ! Return an entry in the data arrays to the hash table for reuse.

            ! In:
            !    indx: index of ht%data_label and ht%payload which can be
            !        reused.
            ! In/Out:
            !    ht: hash table.  On output, ht%free_entries is updated with
            !        the specified index.

            integer, intent(in) :: indx
            type(hash_table_t), intent(inout) :: ht

            ht%next_free_entry = modulo(ht%next_free_entry+1, ht%max_free_reqd)
            ht%free_entries(ht%next_free_entry) = indx

        end subroutine register_hash_table_free_entry

        elemental subroutine delete_hash_table_entry(ht, pos)

            ! Remove an entry from the hash table.

            ! In:
            !    pos: object containing the position of the entry in both the
            !        ht%table and ht%data_label/ht%payload arrays to be
            !        removed.
            ! In/Out:
            !    ht: hash table.  On exit, the entry described by pos has been
            !        deleted (NOTE: only in the hash table lookup, the data
            !        label and payload have, for efficiency, not been
            !        overwritten) and the element used in the data label and
            !        payload arrays has been registered for re-use.

            type(hash_table_pos_t), intent(in) :: pos
            type(hash_table_t), intent(inout) :: ht

            ht%table(pos%ientry:pos%islot-1,pos%islot) = ht%table(pos%ientry+1:pos%islot,pos%islot)
            ht%table(0,pos%islot) = ht%table(0,pos%islot) - 1

            call register_hash_table_free_entry(ht, pos%indx)

        end subroutine delete_hash_table_entry

!--- Add to hash table ---

        elemental subroutine assign_hash_table_entry(ht, slot, pos, err_code)

            ! Assign a new entry to a given location in the hash table.

            ! In:
            !    slot: the slot in the hash table to which a new entry is
            !        to be assigned.
            ! In/Out:
            !    ht: hash table.  On exit, a new entry in the specified slot in
            !        the hash table is assigned an available element of the
            !        data_label and payload arrays.
            ! Out:
            !    pos: object containing the position of the entry in both the
            !        ht%table and ht%data_label/ht%payload arrays.
            !    err_code: error code.  Non-zero if an error is encountered.
            !        See error codes defined at module-level.

            use utils, only: int_fmt

            type(hash_table_t), intent(inout) :: ht
            integer, intent(in) :: slot
            type(hash_table_pos_t), intent(out) :: pos
            integer, intent(out) :: err_code

            pos%islot = slot
            if (ht%table(0,pos%islot) == ht%nentries) then
                err_code = HT_ERR_COLLISIONS
            else
                call take_hash_table_free_entry(ht, pos%indx, err_code)
                if (err_code == 0) then
                    pos%ientry = ht%table(0,pos%islot) + 1
                    ht%table(0,pos%islot) = pos%ientry
                    ht%table(pos%ientry,pos%islot) = pos%indx
                end if
            end if

        end subroutine assign_hash_table_entry

!--- Query hashed store ---

        subroutine lookup_hash_table_entry(ht, label, pos, hit)

            ! Find the position/location of a data item.

            ! In:
            !    ht: hash table.
            !    label: data label to find in the hash table.
            ! Out:
            !    pos: position of the label in ht%table and in the
            !        ht%table/ht%payload arrays.  WARNING: if hit is false then
            !        pos%ientry and pos%indx contain incorrect information but
            !        pos%slot is the 'slot' in the hash table in which the
            !        label should be placed.
            !    hit: true if label is found in the hash table.

            use hashing, only: murmurhash_bit_string

            type(hash_table_t), intent(in) :: ht
            integer(i0), intent(in) :: label(:)
            type(hash_table_pos_t), intent(out) :: pos
            logical, intent(out) :: hit

            integer :: i

            hit = .false.
            pos%islot = modulo(murmurhash_bit_string(label, size(label), ht%seed),ht%nslots) + 1

            ! Need to search over elements in the table with this hash%nslots
            ! value (i.e. search over hash collisions).
            do i = 1, ht%table(0,pos%islot)
                pos%indx = ht%table(i,pos%islot)
                ! If data_label points to an external array, then the first
                ! dimension might exceed ht%data_len.  We do, however, assume
                ! the user has only passed in an array of size ht%data_len...
                if (all(ht%data_label(:ht%data_len,pos%indx) == int(label,int_s))) then
                    pos%ientry = i
                    hit = .true.
                    exit
                end if
            end do

        end subroutine lookup_hash_table_entry

!--- Bulk operations ---

        subroutine accumulate_in_hash_table(ht, labels, payload, err_code)

            ! Sum data belonging to the same label using a hash table.

            ! In:
            !    labels: label of each item in payload.
            !    payload: set of data to be summed over.
            ! In/Out:
            !    ht: hash table.  Must be allocated on input but all existing
            !        information contained by ht is wiped.  On output, contains
            !        one entry per data label and the payload for each entry is
            !        the sum of the items in payload corresponding to that data
            !        label.
            ! Out:
            !    err_code: error code.  Non-zero if an error is encountered.
            !        See error codes defined at module-level.

            type(hash_table_t), intent(inout) :: ht
            integer(i0), intent(in) :: labels(:,:)
            integer(int_p), intent(in) :: payload(:,:)
            integer, intent(out) :: err_code

            type(hash_table_pos_t) :: pos
            logical :: hit
            integer :: i

            call reset_hash_table(ht)
            err_code = 0

            do i = 1, ubound(labels, dim=2)
                call lookup_hash_table_entry(ht, labels(:,i), pos, hit)
                if (hit) then
                    ht%payload(:,pos%indx) = ht%payload(:,pos%indx) + payload(:,i)
                else
                    call assign_hash_table_entry(ht, pos%islot, pos, err_code)
                    ht%payload(:,pos%indx) = payload(:,i)
                    if (err_code /= 0) exit
                end if
            end do

        end subroutine accumulate_in_hash_table

!--- Hash table manipulation utilities ---

! TODO: expand/resize ht%table

end module hash_table
