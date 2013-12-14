module loadbal_data

! Module containing useful global data which also initialises proc_map

use const
use fciqmc_data, only: load_balancing_slots
use parallel, only : nprocs, iproc

implicit none 

integer, allocatable :: proc_map(:)
integer :: p_map_size

contains 
       
    subroutine initialise_proc_map(proc_map_d, p_map_size, load_balancing_slots)
        
        ! Determinants are assigned to processors using proc_map(0:load_balancing_slots*nprocs-1) which 
        ! is initialised so that its entries cyclically contain processors 0,..,nprocs-1

        ! In/out: 
        ! proc_map_d: Will become proc_map upon allocation
        ! load_balancing_slots: number of "slots" we subdivide interval by on each processor
        ! p_map_size: length of proc_map

        use parallel, only : nprocs

        integer, intent(out) :: p_map_size
        integer, intent(in) :: load_balancing_slots
        integer , allocatable, intent(out) :: proc_map_d(:)
        integer :: i

        p_map_size=load_balancing_slots*nprocs
        
        allocate(proc_map(0:p_map_size-1))
        
        do i=0,p_map_size-1
            proc_map_d(i)=modulo(i, nprocs)
        end do

    end subroutine initialise_proc_map 

end module loadbal_data

