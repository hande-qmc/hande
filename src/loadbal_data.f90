module loadbal_data

implicit none

integer, allocatable :: proc_map(:)
integer :: p_map_size

contains

    subroutine initialise_proc_map(load_balancing_slots, proc_map_d, p_map_size)

        ! Determinants are assigned to processors using proc_map(0:load_balancing_slots*nprocs-1) which
        ! is initialised so that its entries cyclically contain processors 0,..,nprocs-1

        ! In:
        !    load_balancing_slots: number of "slots" we subdivide interval by on each processor
        ! Out:
        !    proc_map_d: will become proc_map upon allocation.
        !    p_map_size: length of proc_map

        use parallel, only : nprocs

        integer, intent(in) :: load_balancing_slots
        integer, intent(out) :: p_map_size
        integer , allocatable, intent(out) :: proc_map_d(:)
        integer :: i

        p_map_size = load_balancing_slots*nprocs

        allocate(proc_map(0:p_map_size-1))

        do i = 0, p_map_size - 1
            proc_map_d(i) = modulo(i, nprocs)
        end do

    end subroutine initialise_proc_map

end module loadbal_data

