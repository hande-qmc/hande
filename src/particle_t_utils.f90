module particle_t_utils

! Allocation and deallocation of particle_t objects.  (Helpful to have this in its own
! module to avoid lots of code in qmc_data and to avoid circular dependencies.)

implicit none

contains

    subroutine alloc_particle_t(max_nstates, tensor_label_len, pl)

        ! In:
        !    max_nstates: maximum number of states that can be held in pl.  Sets the
        !       second dimension of the states, pops and dat arrays.
        !    tensor_label_len: number of integers which make up the tensor label bit string.
        ! In/Out:
        !    pl: particle_t object.  On input pl%nspaces and pl%info_size must be set.  On
        !       output all allocatable components are appropriately allocated.

        use qmc_data, only: particle_t
        use checking, only: check_allocate
        use parallel, only: nprocs

        integer, intent(in) :: max_nstates, tensor_label_len
        type(particle_t), intent(inout) :: pl
        integer :: ierr

        allocate(pl%nparticles(pl%nspaces), stat=ierr)
        call check_allocate('pl%nparticles', pl%nspaces, ierr)

        allocate(pl%tot_nparticles(pl%nspaces), stat=ierr)
        call check_allocate('pl%tot_nparticles', pl%nspaces, ierr)

        allocate(pl%states(tensor_label_len,max_nstates), stat=ierr)
        call check_allocate('pl%states', tensor_label_len*max_nstates, ierr)

        allocate(pl%pops(pl%nspaces,max_nstates), stat=ierr)
        call check_allocate('pl%pops', pl%nspaces*max_nstates, ierr)

        allocate(pl%dat(pl%nspaces+pl%info_size,max_nstates), stat=ierr)
        call check_allocate('pl%dat', size(pl%dat), ierr)

        allocate(pl%nparticles_proc(pl%nspaces, nprocs), stat=ierr)
        call check_allocate('pl%nparticles_proc', nprocs*pl%nspaces, ierr)

    end subroutine alloc_particle_t

    subroutine dealloc_particle_t(pl)

        ! In/Out:
        !    pl: particle_t object.  On output all allocatable components are deallocated.

        use qmc_data, only: particle_t
        use checking, only: check_deallocate

        type(particle_t), intent(inout) :: pl
        integer :: ierr

        deallocate(pl%nparticles, stat=ierr)
        call check_deallocate('pl%nparticles', ierr)

        deallocate(pl%tot_nparticles, stat=ierr)
        call check_deallocate('pl%tot_nparticles', ierr)

        deallocate(pl%states, stat=ierr)
        call check_deallocate('pl%states', ierr)

        deallocate(pl%pops, stat=ierr)
        call check_deallocate('pl%pops', ierr)

        deallocate(pl%dat, stat=ierr)
        call check_deallocate('pl%dat', ierr)

        deallocate(pl%nparticles_proc, stat=ierr)
        call check_deallocate('pl%nparticles_proc', ierr)

    end subroutine dealloc_particle_t

end module particle_t_utils
