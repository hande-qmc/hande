module loadbal_data

! module containing useful global data which also initialises proc_map
    use const
    use fciqmc_data
    use parallel, only : nprocs, iproc

    implicit none 

    integer :: num_slots=20
    integer, allocatable :: proc_map(:)
    integer :: p_map_size

contains 
       
    subroutine initialise_proc_map(proc_map_d, p_map_size, num_slots)
        
        ! Determinants are assigned to processors using proc_map(num_slots*nprocs) which 
        ! is initialised so that its entries cyclically contain processors 0,..,nprocs-1
        ! In/out: 
        ! proc_map: array which maps determinants to processors
        !       proc_map(modulo(hash(d),num_slots*nprocs)=processor
        ! num_slots: number of "slots" we subdivide interval by on each processor

        use parallel, only : nprocs

        integer :: i
        integer, intent(inout) :: p_map_size
        integer, intent(in) :: num_slots
        integer , allocatable, intent(out) :: proc_map_d(:)
        
        p_map_size=num_slots*nprocs
        
        allocate(proc_map(p_map_size))
        
        do i=0,p_map_size-1
            proc_map_d(i+1)=modulo(i, nprocs)
        end do

    end subroutine initialise_proc_map 

end module loadbal_data

