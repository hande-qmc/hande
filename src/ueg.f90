module ueg_system

! Core routines for the uniform electron gas system.

use const

implicit none

contains

!-------
! Initialisation, utilities and finalisation

    subroutine init_ueg_proc_pointers(ndim, sys_ueg)

        ! Initialise UEG procedure pointers

        ! In:
        !    ndim: dimensionality of the UEG.
        ! In/Out:
        !    sys_ueg: sys_ueg_t object.  On output, the procedure pointers are set.

        use system, only: sys_ueg_t
        use errors, only: stop_all

        integer, intent(in) :: ndim
        type(sys_ueg_t), intent(inout) :: sys_ueg

        ! Set pointers to integral routines
        select case(ndim)
        case(2)
            sys_ueg%coulomb_int => coulomb_int_ueg_2d
        case(3)
            sys_ueg%coulomb_int => coulomb_int_ueg_3d
        case default
            call stop_all('init_ueg_proc_pointers', 'Can only do 2D and 3D UEG.')
        end select

        ! For now, we don't treat exchange integrals differently.
        sys_ueg%exchange_int => sys_ueg%coulomb_int

    end subroutine init_ueg_proc_pointers

    subroutine init_ueg_indexing(sys)

        ! Create arrays and data for index mapping needed for UEG.

        ! In/Out:
        !    sys: UEG system to be studied.  On output, the basis indexing components are set.

        use system, only: sys_t

        use checking, only: check_allocate
        use utils, only: tri_ind

        type(sys_t), intent(inout) :: sys

        integer :: ierr, i, j, a, ind, N_kx, k_min(sys%lattice%ndim), bit_pos, bit_el, k(3)
        integer :: k1, k2, k3, ktest(sys%lattice%ndim), kija(sys%lattice%ndim)

        sys%ueg%basis%kmax = ceiling(sqrt(2*sys%ueg%ecutoff))

        N_kx = 2*sys%ueg%basis%kmax+1

        allocate(sys%ueg%basis%offset_inds(sys%lattice%ndim), stat=ierr)
        call check_allocate('sys%ueg%basis%offset_inds', sys%lattice%ndim, ierr)
        forall (i=1:sys%lattice%ndim) sys%ueg%basis%offset_inds(i) = N_kx**(i-1)

        ! Wish the indexing array to be 1-indexed.
        k_min = -sys%ueg%basis%kmax ! Bottom corner of grid.
        sys%ueg%basis%offset = -dot_product(sys%ueg%basis%offset_inds, k_min) + 1

        allocate(sys%ueg%basis%lookup(N_kx**sys%lattice%ndim), stat=ierr)
        call check_allocate('sys%ueg%basis%lookup', N_kx**sys%lattice%ndim, ierr)

        ! sys%ueg%basis%lookup should be -1 for any wavevector that is in the
        ! square/cubic grid defined by sys%ueg%basis%kmax but not in the actual basis
        ! set described by ecutoff.
        sys%ueg%basis%lookup = -1

        ! Now fill in the values for the alpha orbitals which are in the basis.
        associate(bfns=>sys%basis%basis_fns, ueg_basis=>sys%ueg%basis)
            forall (i=1:sys%basis%nbasis:2)
                ueg_basis%lookup(dot_product(bfns(i)%l, ueg_basis%offset_inds) + ueg_basis%offset) = i
            end forall
        end associate

        ! Now fill in the values for permitted k_a in an excitation
        ! k_i,k_j->k_a,k_b, given a choice of k_i ad k_j and requiring k_b is in
        ! the basis.

        k = 0
        k(1:sys%lattice%ndim) = 2*sys%ueg%basis%kmax
        allocate(sys%ueg%ternary_conserve(0:sys%basis%string_len, -k(1):k(1), -k(2):k(2), -k(3):k(3)), stat=ierr)
        call check_allocate('sys%ueg%ternary_conserve', size(sys%ueg%ternary_conserve), ierr)
        sys%ueg%ternary_conserve = 0_i0
        !$omp parallel do default(none) shared(k,sys) &
        !$omp private(k1,k2,k3,a,kija,ktest,bit_pos,bit_el)
        do k3 = -k(3), k(3)
            if (sys%lattice%ndim == 3) ktest(3) = k3
            do k2 = -k(2), k(2)
               if (sys%lattice%ndim >= 2) ktest(2) = k2
               do k1 = -k(1), k(1)
                   ktest(1) = k1
                   ! If this is still slow, we can improve matters by
                   ! restricting the range of a such that there must be at least one b
                   ! (i.e. all components of k_i+k_j-k_a must lie within +/- ! sys%ueg%basis%kmax,
                   ! thus providing lower and upper bounds for a).
                   do a = 1, sys%basis%nbasis-1, 2 ! only bother with alpha orbitals
                       kija = ktest - sys%basis%basis_fns(a)%l
                       if (real(dot_product(kija,kija),p)/2 - sys%ueg%ecutoff < 1.e-8) then
                           ! There exists an allowed b in the basis!
                           sys%ueg%ternary_conserve(0,k1,k2,k3) = sys%ueg%ternary_conserve(0,k1,k2,k3) + 1
                           bit_pos = sys%basis%bit_lookup(1, a)
                           bit_el = sys%basis%bit_lookup(2, a)
                           sys%ueg%ternary_conserve(bit_el,k1,k2,k3) = ibset(sys%ueg%ternary_conserve(bit_el,k1,k2,k3), bit_pos)
                       end if
                   end do
               end do
            end do
        end do
        !$omp end parallel do

    end subroutine init_ueg_indexing

    pure function ueg_basis_index(ueg_basis, k, spin) result(indx)

        ! In:
        !    ueg_basis: UEG basis-lookup info.
        !    k: wavevector in units of 2\pi/L.
        !    spin: 1 for alpha orbital, -1 for beta orbital
        ! Returns:
        !    Index of basis function in the (energy-ordered) sys%basis%basis_fns array.
        !    Set to < 0 if the spin-orbital described by k and spin is not in the
        !    basis set.

        use ueg_types, only: ueg_basis_t

        integer :: indx
        type(ueg_basis_t), intent(in) :: ueg_basis
        integer, intent(in) :: k(:), spin

        if (minval(k) < -ueg_basis%kmax .or. maxval(k) > ueg_basis%kmax) then
            indx = -1
        else
            ! ueg_basis%lookup contains the mapping between a wavevector in
            ! a given basis and its entry in the energy-ordered list of basis
            ! functions.
            indx = ueg_basis%lookup(dot_product(k,ueg_basis%offset_inds) + ueg_basis%offset)
            ! ueg_basis%lookup only contains entries for the alpha spin-orbital.
            ! The corresponding beta orbital is the next entry in the basis_fns
            ! array.
            if (spin < 0) indx = indx + 1
        end if

    end function ueg_basis_index

    subroutine end_ueg_indexing(sys_ueg)

        ! Clean up UEG index arrays.

        ! In/Out:
        !    sys_ueg_t: On output, all allocatable components are deallocated.

        use checking, only: check_deallocate
        use system, only: sys_ueg_t

        type(sys_ueg_t), intent(inout) :: sys_ueg

        integer :: ierr

        if (allocated(sys_ueg%basis%lookup)) then
            deallocate(sys_ueg%basis%lookup, stat=ierr)
            call check_deallocate('sys_ueg%basis%lookup', ierr)
        end if
        if (allocated(sys_ueg%basis%offset_inds)) then
            deallocate(sys_ueg%basis%offset_inds, stat=ierr)
            call check_deallocate('sys_ueg%basis%offset_inds', ierr)
        end if
        if (allocated(sys_ueg%ternary_conserve)) then
            deallocate(sys_ueg%ternary_conserve, stat=ierr)
            call check_deallocate('sys_ueg%ternary_conserve', ierr)
        end if

    end subroutine end_ueg_indexing

!-------
! Integrals

    pure function get_two_e_int_ueg(sys, i, j, a, b) result(intgrl)

        ! In:
        !    sys: system being studied.
        !    i,j:  index of the spin-orbital from which an electron is excited in
        !          the reference determinant.
        !    a,b:  index of the spin-orbital into which an electron is excited in
        !          the excited determinant.
        !
        ! Returns:
        !   The anti-symmetrized integral < ij || ab >.

        ! Warning: assume i,j /= a,b (ie not asking for < ij || ij > or < ij || ji >).

        use system, only: sys_t

        real(p) :: intgrl
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b

        intgrl = 0.0_p

        ! Crystal momentum conserved?
        associate(bg=>sys%basis)
            if (all(bg%basis_fns(i)%l + bg%basis_fns(j)%l - bg%basis_fns(a)%l - bg%basis_fns(b)%l == 0)) then

                ! Spin conserved?

                ! Coulomb
                if (bg%basis_fns(i)%ms == bg%basis_fns(a)%ms .and.  bg%basis_fns(j)%ms == bg%basis_fns(b)%ms) &
                    intgrl = intgrl + sys%ueg%coulomb_int(sys%lattice%box_length(1), bg, i, a)

                ! Exchange
                if (bg%basis_fns(i)%ms == bg%basis_fns(b)%ms .and.  bg%basis_fns(j)%ms == bg%basis_fns(a)%ms) &
                    intgrl = intgrl - sys%ueg%coulomb_int(sys%lattice%box_length(1), bg, i, b)

            end if
        end associate

    end function get_two_e_int_ueg

    pure function coulomb_int_ueg_2d(cell_param, basis_set, i, a) result(intgrl)

        ! In:
        !    cell_param: the length of one side of the simulation cell.
        !    basis_set: basis set of system being studied.
        !    i: index of spin-orbital basis function.
        !    a: index of spin-orbital basis function.
        !
        ! Returns:
        !    The Coulumb integral < i j | a b > = 2\pi/\Omega|k_i - k_a| for the 2D
        !    UEG.  Note that we assume i, j, a and b are such that spin and
        !    crystal momentum is conserved and hence the integral is not zero by
        !    symmetry.  We also assume that i,j /= a,b (ie the integral is not
        !    a Hartree integral).

        use basis, only: basis_t

        real(p) :: intgrl
        real(p), intent(in) :: cell_param
        type(basis_t), intent(in) :: basis_set
        integer, intent(in) :: i, a
        integer :: q(2)

        ! Wavevectors are stored in units of 2\pi/L, where L is the length of
        ! the cell.  As we only deal with cubic simulation cells (i.e. \Omega = L^2),
        ! the integral hence becomes 1/(L|q|), where q = k_i - k_a.

        q = basis_set%basis_fns(i)%l - basis_set%basis_fns(a)%l
        intgrl = 1.0_p/(cell_param*sqrt(real(dot_product(q,q),p)))

    end function coulomb_int_ueg_2d

    pure function coulomb_int_ueg_3d(cell_param, basis_set, i, a) result(intgrl)

        ! In:
        !    cell_param: the length of one side of the simulation cell.
        !    basis_set: basis set of system being studied.
        !    i: index of spin-orbital basis function.
        !    a: index of spin-orbital basis function.
        !
        ! Returns:
        !    The Coulumb integral < i j | a b > = 1/|k_i - k_a|^2 for the 3D
        !    UEG.  Note that we assume i, j, a and b are such that spin and
        !    crystal momentum is conserved and hence the integral is not zero by
        !    symmetry.  We also assume that i,j /= a,b (ie the integral is not
        !    a Hartree integral).

        use basis, only: basis_t

        real(p) :: intgrl
        real(p), intent(in) :: cell_param
        type(basis_t), intent(in) :: basis_set
        integer, intent(in) :: i, a
        integer :: q(3)

        ! Wavevectors are stored in units of 2\pi/L, where L is the length of
        ! the cell.  As we only deal with cubic simulation cells (i.e. \Omega = L^3),
        ! the integral hence becomes 1/(\pi.L.q^2), where q = k_i - k_a.

        q = basis_set%basis_fns(i)%l - basis_set%basis_fns(a)%l
        intgrl = 1.0_p/(pi*cell_param*dot_product(q,q))

    end function coulomb_int_ueg_3d

    subroutine set_derived_ueg_properties(sys)

        ! Set derived UEG properties e.g. the Fermi energy, and wavevector.

        ! In/Out:
        !    sys: sys_t object, on output all derived UEG quantities are
        !       set.

        use system, only: sys_t

        type(sys_t), intent(inout) :: sys

        integer :: pol_factor

        ! Polarisation factor = 2 for polarised system, 1 for unpolarised.
        pol_factor = 1 + abs((sys%nalpha-sys%nbeta)/sys%nel)

        ! Only deal with fully spin (un)polarised system, so that kf^{up} =
        ! kf^{down}.
        ! Fermi wavevector.
        sys%ueg%kf = (9.0_dp*pol_factor*pi/(4.0_dp*sys%ueg%r_s**3))**(1.0_dp/3.0_dp)
        ! Fermi Energy.
        sys%ueg%ef = 0.5 * sys%ueg%kf**2

        write (6, '(1X,a14, f13.10)') "Fermi Energy: ", sys%ueg%ef

    end subroutine set_derived_ueg_properties

end module ueg_system
