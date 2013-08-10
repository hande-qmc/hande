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
        integer(i0), allocatable :: data_label(:,:) ! (:,nentries*nslots)

        ! payload(:,i) is the 'payload' or information associated with
        ! the i-th data item in data_label.
        integer, allocatable :: payload(:,:) ! (:,nentries*nslots)

        ! We usually want to iterate over the data items, so it is best if they
        ! are (as far as possible) contiguous.
        ! If not, we can insert elements into the empty element.  Typically (we
        ! hope!) we have more room than the number of items we're actually
        ! storing.  If there are no slots in the occupied section of data_label
        ! and payload, then we use the next free element at the head of the data
        ! structure.

        ! Highest occupied element in the data_label and payload arrays.
        integer :: head

        ! Maximum number of free slots we might wish to access before iterating
        ! over data_label again (at which point we can easily update the list of
        ! free entries).
        integer :: max_free_reqd
        ! The index of the element in free_elements which grants the next free
        ! entry in data_label and payload arrays.
        integer :: next_free_element
        ! Circular array.  Non-zero elements correspond to free elements in the
        ! data_label and corresponding payload arrays.
        integer, allocatable :: free_elements(:) ! (0:max_free_reqd-1)

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

    contains

!--- Allocation/deallocation ---

        subroutine alloc_hash_table(nslots, nentries, data_len, payload_len, max_free_reqd, ht)

            ! Allocate hash table and initialise storage.

            ! In:
            !    nslots, nentries, data_len, payload_len, max_free_reqd: see
            !        descriptions in hash_table_t.
            ! Out:
            !    ht: hash table.  On output all array attributes are allocated
            !        and all scalar attributes are set to their corresponding
            !        input arguments.  The hash table (ht%table) is set to
            !        0 (see reset_hash_table), ready for use.

            use checking, only: check_allocate

            integer, intent(in) :: nslots, nentries, data_len, payload_len, max_free_reqd
            type(hash_table_t), intent(out) :: ht

            integer :: ierr

            ht%nslots = nslots
            ht%nentries = nentries
            ht%data_len = data_len
            ht%payload_len = payload_len
            ht%max_free_reqd = max_free_reqd

            allocate(ht%table(ht%nentries,ht%nslots), stat=ierr)
            call check_allocate('ht%table', size(ht%table), ierr)
            allocate(ht%data_label(ht%data_len,ht%nentries*ht%nslots))
            call check_allocate('ht%data_label', size(ht%data_label), ierr)
            allocate(ht%payload(ht%payload_len,ht%nentries*ht%nslots), stat=ierr)
            call check_allocate('ht%payload', size(ht%payload), ierr)
            allocate(ht%free_elements(ht%max_free_reqd), stat=ierr)
            call check_allocate('ht%free_elements', size(ht%free_elements), ierr)

            call reset_hash_table(ht)

        end subroutine alloc_hash_table

        subroutine free_hash_table(ht)

            ! Deallocate hash table and set scalar attributes to 0.

            ! In/Out:
            !    ht: hash table.  On output all array attributes are deallocated
            !        and all scalar attributes are set to 0.

            use checking, only: check_deallocate

            type(hash_table_t), intent(inout) :: ht

            integer :: ierr

            if (allocated(ht%table)) then
                deallocate(ht%table)
                call check_deallocate('ht%table', ierr)
            end if
            if (allocated(ht%data_label)) then
                deallocate(ht%data_label)
                call check_deallocate('ht%data_label', ierr)
            end if
            if (allocated(ht%payload)) then
                deallocate(ht%payload)
                call check_deallocate('ht%payload', ierr)
            end if
            if (allocated(ht%free_elements)) then
                deallocate(ht%free_elements)
                call check_deallocate('ht%free_elements', ierr)
            end if

            ht%nslots = 0
            ht%nentries = 0
            ht%data_len = 0
            ht%payload_len = 0
            ht%head = 0
            ht%max_free_reqd = 0
            ht%next_free_element = 0

        end subroutine free_hash_table

        elemental subroutine reset_hash_table(ht)

            ! Reset hash table.

            ! In/Out:
            !    ht: hash table.  On output the table, head, next_free_element
            !        and free_elements attributes are all set to 0.  This means
            !        that no data is stored in the hash table and new elements
            !        are added by consecutively by incrementing head.

            type(hash_table_t), intent(inout) :: ht

            integer :: i

            ht%head = 0
            ht%next_free_element = 0
            ht%free_elements = 0
            ht%table = 0

        end subroutine reset_hash_table

!--- Remove/add known free slots in occupied part of data_load/payload arrays ---

        elemental subroutine take_hash_table_free_slot(ht, indx)

            ! Get the next index of the data arrays to be used.

            ! In/Out:
            !    ht: hash table.  On output, the entry given by indx is
            !        registered as being used (i.e. removed from the
            !        ht%free_elements or ht%head is incremented).
            ! Out:
            !    indx: index of free entry in ht%data_label and ht%payload to be
            !        used.  The index is taken from the free_elements_array if
            !        there is a index available, otherwise the index is the
            !        entry directly above the current head (top) entry of
            !        ht%data_label and ht%payload arrays.

            ! WARNING: this does not check that the returned index is inside the
            ! available space (ie less than ht%head)...

            type(hash_table_t), intent(inout) :: ht
            integer, intent(out) :: indx

            indx = ht%free_elements(ht%next_free_element)
            ht%free_elements(ht%next_free_element) = 0
            ht%next_free_element = modulo(ht%next_free_element-1,ht%max_free_reqd)

            if (indx == 0) then
                ! Have run out of free slots that we know about in the
                ! used part of the storage.  Now just add to the head of the
                ! occupied space...
                ht%head = ht%head + 1
                indx = ht%head
            end if

        end subroutine take_hash_table_free_slot

        elemental subroutine register_hash_table_free_slot(ht, indx)

            ! Return an index of the data arrays for reuse.

            ! In:
            !    indx: index of ht%data_label and ht%payload which can be
            !        reused.
            ! In/Out:
            !    ht: hash table.  On output, ht%free_elements is updated with
            !        the specified index.

            integer, intent(in) :: indx
            type(hash_table_t), intent(inout) :: ht

            ht%next_free_element = modulo(ht%next_free_element+1, ht%max_free_reqd)
            ht%free_elements(ht%next_free_element) = indx

        end subroutine register_hash_table_free_slot

        pure subroutine delete_hash_table_entry(ht, pos)

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

            integer :: i, indx

            do i = pos%ientry+1, ht%table(0,pos%islot)
                ht%table(i-1,pos%islot) = ht%table(i,pos%islot)
            end do
            ht%table(0,pos%islot) = ht%table(0,pos%islot) - 1

            call register_hash_table_free_slot(ht, pos%indx)

        end subroutine delete_hash_table_entry

!--- Query hashed store ---

        subroutine lookup_hash_table_element(ht, label, register, pos, err_flag)

            ! Find the position/location of a data item.

            ! In:
            !    label: data label to find in the hash table.
            !    register: if true and the data label cannot be found in the
            !        hash table, then return the position where label should be
            !        stored, in the process updating the hash table to remove
            !        the given position from being used by another data item.
            ! In/Out:
            !    ht: hash table.  Modified only if register is true.
            ! Out:
            !    pos: position of the label in ht%table and in the
            !        ht%table/ht%payload arrays.  If err_flag is 1 and register
            !        is not true, then pos%indx is *not* set.
            !    err_flag: 0 indicates no error.  1 indicates label could not be
            !        found.  If register is true, 2 indicates that the data
            !        storage (ht%data_label and ht%payload) is completely full
            !        and 3 indicates the label could not be found and ht%table
            !        has run out of entries for the hash of label (ie we have
            !        had ht%nentries+1 hash collisions).

            use hashing, only: murmurhash_bit_string

            type(hash_table_t), intent(inout) :: ht
            integer(i0), intent(in) :: label(:)
            logical, intent(in) :: register
            type(hash_table_pos_t), intent(out) :: pos
            integer, intent(out) :: err_flag

            integer :: i

            pos%islot = modulo(murmurhash_bit_string(label, size(label)),ht%nslots)

            ! Need to search over elements in with this hash%nslots value.
            do i = 1, ht%table(0,pos%islot)
                pos%indx = ht%table(i,pos%islot)
                if (all(ht%data_label(:,pos%indx) == label)) exit
            end do
            pos%ientry = i
            if (i /= ht%table(0,pos%islot)+1) then
                err_flag = 0
            else if (register) then
                ! Get a free entry and say we're using it!
                call take_hash_table_free_slot(ht, pos%indx)
                if (pos%indx > size(ht%data_label,dim=2)) then
                    err_flag = 2
                else if (i > ht%nentries) then
                    err_flag = 3
                else
                    err_flag = 1
                    ht%table(0,pos%islot) = pos%ientry
                    ht%table(pos%ientry,pos%islot) = pos%indx
                end if
            else
                err_flag = 1
            end if

        end subroutine lookup_hash_table_element

!--- Hash table manipulation utilities ---

! TODO: expand/resize ht%table.

end module hash_table
