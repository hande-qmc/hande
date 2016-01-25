module particle_t_utils

! Allocation and deallocation of particle_t objects.  (Helpful to have this in its own
! module to avoid lots of code in qmc_data and to avoid circular dependencies.)

implicit none

contains

    subroutine init_particle_t(max_nstates, nwalker_int_extra, tensor_label_len, real_amplitudes, real32, pl, verbose_output)

        ! In:
        !    max_nstates: maximum number of states that can be held in pl.  Sets the
        !       second dimension of the states, pops and dat arrays.
        !    nwalker_int_extra: the number of additional (32-bit) integers used per state
        !       outside of particle_t, e.g. for semi_stoch_t%determ.
        !    tensor_label_len: number of integers which make up the tensor label bit string.
        !    real_amplitudes: if true, use real ampltiudes (via fixed precision) rather
        !       than integer amplitudes.
        !    real32: if true, use the fractional precision used when POP_SIZE=32 even if
        !       POP_SIZE=64.
        !    verbose_output (optional): if true (default: false) print out information about
        !       the memory allocation of pl.
        ! In/Out:
        !    pl: particle_t object.  On input pl%nspaces and pl%info_size must be set.  On
        !       output all allocatable components are appropriately allocated.

        use const, only: int_p, i0_length, int_p_length, p, sp
        use errors, only: warning
        use parallel, only: nprocs, parent
        use utils, only: int_fmt

        use qmc_data, only: particle_t
        use checking, only: check_allocate

        integer, intent(in) :: max_nstates, tensor_label_len
        integer, intent(in) :: nwalker_int_extra ! (e.g.) for semi_stoch_t%determ
        logical, intent(in) :: real_amplitudes, real32
        type(particle_t), intent(inout) :: pl
        logical, intent(in), optional :: verbose_output
        integer :: ierr, nwalker_int_p, nwalker_real, size_main_walker, pop_bit_shift, max_nstates_elements
        logical :: verbose

        verbose = .true.
        if (present(verbose_output)) verbose = verbose_output

        ! Set the real encoding shift, depending on whether 32 or 64-bit integers
        ! are being used.
        if (real_amplitudes) then
            if (bit_size(0_int_p) == 64) then
                if (real32) then
                    ! Use same space for fractional part as 32-bit integers.
                    pop_bit_shift = 11
                else
                    ! Allow a maximum population of 2^32, and a minimum fractional
                    ! part of 2^-31.
                    pop_bit_shift = 31
                end if
            else if (bit_size(0_int_p) == 32) then
                ! Allow a maximum population of 2^20, and a minimum fractional
                ! part of 2^-11.
                pop_bit_shift = 11
                if (parent) then
                    call warning('init_particle_t', &
                        'You are using 32-bit walker populations with real amplitudes.'//new_line('')// &
                        ' The maximum population size on a given determinant is 2^20=1048576.&
                        & Errors will occur if this is exceeded.'//new_line('')//&
                        ' Compile HANDE with the CPPFLAG -DPOP_SIZE=64 to use 64-bit populations.', 2)
                end if
            end if
        else
            ! Allow no fractional part for walker populations.
            pop_bit_shift = 0
        end if
        ! Store 2**pop_bit_shift for ease.
        pl%pop_real_factor = 2_int_p**(int(pop_bit_shift, int_p))

        ! The the number of bits occupied by each determinant in the main
        ! walker list is given by string_len*i0_length+nwalker_int_extra*32+
        ! nwalker_int_p*int_p_length+nwalker_real*32 (*64 if double precision).
        ! The number of bytes is simply 1/8 this.
        nwalker_int_p = pl%nspaces ! for populations
        nwalker_real = pl%nspaces + pl%info_size ! for <D_i|O|D_i> and info storage.
        if (p == sp) then
            ! SINGLE_PRECISION
            size_main_walker = tensor_label_len*i0_length/8 + nwalker_int_p*int_p_length/8 + &
                               nwalker_int_extra*4 + nwalker_real*4
        else
            size_main_walker = tensor_label_len*i0_length/8 + nwalker_int_p*int_p_length/8 + &
                               nwalker_int_extra*4 + nwalker_real*8
        end if
        max_nstates_elements = max_nstates
        if (max_nstates_elements < 0) then
            ! Given in MB.  Convert.  Note: important to avoid overflow in the
            ! conversion!
            max_nstates_elements = int((-real(max_nstates_elements,p)*10**6)/size_main_walker)
        end if

        if (parent .and. verbose) then
            write (6,'(1X,a53,f7.2)') 'Memory allocated per core for main walker list (MB): ', &
                                      size_main_walker*real(max_nstates_elements,p)/10**6
            write (6,'(1X,a48,'//int_fmt(max_nstates_elements,1)//')') &
                    'Number of elements per core in main walker list:', max_nstates_elements
        end if

        allocate(pl%nparticles(pl%nspaces), stat=ierr)
        call check_allocate('pl%nparticles', pl%nspaces, ierr)

        allocate(pl%tot_nparticles(pl%nspaces), stat=ierr)
        call check_allocate('pl%tot_nparticles', pl%nspaces, ierr)

        allocate(pl%states(tensor_label_len,max_nstates_elements), stat=ierr)
        call check_allocate('pl%states', tensor_label_len*max_nstates_elements, ierr)

        allocate(pl%pops(pl%nspaces,max_nstates_elements), stat=ierr)
        call check_allocate('pl%pops', pl%nspaces*max_nstates_elements, ierr)

        allocate(pl%dat(pl%nspaces+pl%info_size,max_nstates_elements), stat=ierr)
        call check_allocate('pl%dat', (pl%nspaces+pl%info_size)*max_nstates_elements, ierr)

        allocate(pl%nparticles_proc(pl%nspaces, nprocs), stat=ierr)
        call check_allocate('pl%nparticles_proc', nprocs*pl%nspaces, ierr)

    end subroutine init_particle_t

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

    subroutine move_particle_t(old_particles, new_particles)

        ! Move the particle list from old_particles to new_particles. Using move_alloc avoids
        ! allocating another (potentially very large) main particle list.

        ! In/Out:
        !   old_particles: contains particle list on entry.  Deallocated on exit.
        ! Out:
        !   new_particles: on exit contains original contents of old_particles.

        use qmc_data, only: particle_t

        type(particle_t), intent(inout) :: old_particles
        type(particle_t), intent(out) :: new_particles

        new_particles%nstates = old_particles%nstates
        call move_alloc(old_particles%nparticles, new_particles%nparticles)
        call move_alloc(old_particles%tot_nparticles, new_particles%tot_nparticles)
        call move_alloc(old_particles%nparticles_proc, new_particles%nparticles_proc)
        new_particles%nspaces = old_particles%nspaces
        new_particles%info_size = old_particles%info_size
        new_particles%pop_real_factor = old_particles%pop_real_factor
        call move_alloc(old_particles%states, new_particles%states)
        call move_alloc(old_particles%pops, new_particles%pops)
        call move_alloc(old_particles%dat, new_particles%dat)

    end subroutine move_particle_t

end module particle_t_utils
