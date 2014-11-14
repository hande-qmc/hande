module ueg_system

! Core routines for the uniform electron gas system.

use const
use ueg_types, only: ueg_basis_t

implicit none

! UEG-specific basis lookup tables, etc.
type(ueg_basis_t) :: ueg_basis

! When creating an arbitrary excitation, k_i,k_j->k_a,k_b, we must conserve
! crystal momentum, k_i+k_j-k_a-k_b=0.  Hence once we've chosen k_i, k_j and
! k_a, k_b is uniquely defined.  Further, once we've chosen k_i and k_j and if
! we require k_b to exist in the basis, then only certain values of k_a are
! permitted.  ueg_ternary_conserve(0,k1,k2,k3) gives how many k_a are permitted
! for k_i+k_j = (k1,k2,k3) and ueg_ternary_conserve(1:,k1,k2,k3) gives a bit
! string with only bytes set corresponding to permitted k_a values.  Note only
! basis functions corresponding to *alpha* orbitals are set.
! For systems with dimensionality lower than 3, the higher ki values are set to
! 0, i.e. dimensions:
! (0:string_len,-N:N,0,0) (1D)
! (0:string_len,-N:N,-N:N,0) (2D)
! (0:string_len,-N:N,-N:N,-N:N) (3D)
! NOTE: this contains values of k_i+k_j which cannot be formed by the basis with
! the energy cutoff.  Memory can be saved by not using a cubic array for
! k_i+k_j...
integer(i0), allocatable :: ueg_ternary_conserve(:,:,:,:)

abstract interface
    ! UEG-specific integral procedure pointers.
    ! The integral routines are different for 2D and UEG.  Abstract them using
    ! procedure pointers.
    pure function i_int_ueg(cell_param, b, i, a) result(intgrl)
        use basis, only: basis_t
        use const, only: p
        real(p) :: intgrl
        real(p), intent(in) :: cell_param
        type(basis_t), intent(in) :: b
        integer, intent(in) :: i, a
    end function i_int_ueg
end interface

procedure(i_int_ueg), pointer :: coulomb_int_ueg => null()
procedure(i_int_ueg), pointer :: exchange_int_ueg => null()

contains

!-------
! Initialisation, utilities and finalisation

    subroutine init_ueg_proc_pointers(ndim)

        ! Initialise UEG procedure pointers

        ! In:
        !    ndim: dimensionality of the UEG.

        use system
        use errors, only: stop_all

        integer, intent(in) :: ndim

        ! Set pointers to integral routines
        select case(ndim)
        case(2)
            coulomb_int_ueg => coulomb_int_ueg_2d
        case(3)
            coulomb_int_ueg => coulomb_int_ueg_3d
        case default
            call stop_all('init_ueg_proc_pointers', 'Can only do 2D and 3D UEG.')
        end select

        ! For now, we don't treat exchange integrals differently.
        exchange_int_ueg => coulomb_int_ueg

    end subroutine init_ueg_proc_pointers

    subroutine init_ueg_indexing(sys)

        ! Create arrays and data for index mapping needed for UEG.

        ! In:
        !    sys: UEG system to be studied.

        use system, only: sys_t

        use checking, only: check_allocate
        use utils, only: tri_ind

        type(sys_t), intent(in) :: sys

        integer :: ierr, i, j, a, ind, N_kx, k_min(sys%lattice%ndim), bit_pos, bit_el, k(3)
        integer :: k1, k2, k3, ktest(sys%lattice%ndim), kija(sys%lattice%ndim)

        ueg_basis%kmax = ceiling(sqrt(2*sys%ueg%ecutoff))

        N_kx = 2*ueg_basis%kmax+1

        allocate(ueg_basis%offset_inds(sys%lattice%ndim), stat=ierr)
        call check_allocate('ueg_basis%offset_inds', sys%lattice%ndim, ierr)
        forall (i=1:sys%lattice%ndim) ueg_basis%offset_inds(i) = N_kx**(i-1)

        ! Wish the indexing array to be 1-indexed.
        k_min = -ueg_basis%kmax ! Bottom corner of grid.
        ueg_basis%offset = -dot_product(ueg_basis%offset_inds, k_min) + 1

        allocate(ueg_basis%lookup(N_kx**sys%lattice%ndim), stat=ierr)
        call check_allocate('ueg_basis%lookup', N_kx**sys%lattice%ndim, ierr)

        ! ueg_basis%lookup should be -1 for any wavevector that is in the
        ! square/cubic grid defined by ueg_basis%kmax but not in the actual basis
        ! set described by ecutoff.
        ueg_basis%lookup = -1

        ! Now fill in the values for the alpha orbitals which are in the basis.
        forall (i=1:sys%basis%nbasis:2) 
            ueg_basis%lookup(dot_product(sys%basis%basis_fns(i)%l, ueg_basis%offset_inds) + ueg_basis%offset) = i
        end forall

        ! Now fill in the values for permitted k_a in an excitation
        ! k_i,k_j->k_a,k_b, given a choice of k_i ad k_j and requiring k_b is in
        ! the basis.

        k = 0
        k(1:sys%lattice%ndim) = 2*ueg_basis%kmax
        allocate(ueg_ternary_conserve(0:sys%basis%string_len, -k(1):k(1), -k(2):k(2), -k(3):k(3)), stat=ierr)
        call check_allocate('ueg_ternary_conserve', size(ueg_ternary_conserve), ierr)
        ueg_ternary_conserve = 0_i0
        !$omp parallel do default(none) shared(k,sys,ueg_ternary_conserve) &
        !$omp private(k1,k2,k3,a,kija,ktest,bit_pos,bit_el)
        do k3 = -k(3), k(3)
            if (sys%lattice%ndim == 3) ktest(3) = k3
            do k2 = -k(2), k(2)
               if (sys%lattice%ndim >= 2) ktest(2) = k2
               do k1 = -k(1), k(1)
                   ktest(1) = k1
                   ! If this is still slow, we can improve matters by
                   ! restricting the range of a such that there must be at least one b
                   ! (i.e. all components of k_i+k_j-k_a must lie within +/- ! ueg_basis%kmax,
                   ! thus providing lower and upper bounds for a).
                   do a = 1, sys%basis%nbasis-1, 2 ! only bother with alpha orbitals
                       kija = ktest - sys%basis%basis_fns(a)%l
                       if (real(dot_product(kija,kija),p)/2 - sys%ueg%ecutoff < 1.e-8) then
                           ! There exists an allowed b in the basis!
                           ueg_ternary_conserve(0,k1,k2,k3) = ueg_ternary_conserve(0,k1,k2,k3) + 1
                           bit_pos = sys%basis%bit_lookup(1, a)
                           bit_el = sys%basis%bit_lookup(2, a)
                           ueg_ternary_conserve(bit_el,k1,k2,k3) = ibset(ueg_ternary_conserve(bit_el,k1,k2,k3), bit_pos)
                       end if
                   end do
               end do
            end do
        end do
        !$omp end parallel do

    end subroutine init_ueg_indexing

    pure function ueg_basis_index(k, spin) result(indx)

        ! In:
        !    k: wavevector in units of 2\pi/L.
        !    spin: 1 for alpha orbital, -1 for beta orbital
        ! Returns:
        !    Index of basis function in the (energy-ordered) sys%basis%basis_fns array.
        !    Set to < 0 if the spin-orbital described by k and spin is not in the
        !    basis set.

        use system

        integer :: indx
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

    subroutine end_ueg_indexing()

        ! Clean up UEG index arrays.

        use checking, only: check_deallocate

        integer :: ierr

        if (allocated(ueg_basis%lookup)) then
            deallocate(ueg_basis%lookup, stat=ierr)
            call check_deallocate('ueg_basis%lookup', ierr)
        end if
        if (allocated(ueg_basis%offset_inds)) then
            deallocate(ueg_basis%offset_inds, stat=ierr)
            call check_deallocate('ueg_basis%offset_inds', ierr)
        end if
        if (allocated(ueg_ternary_conserve)) then
            deallocate(ueg_ternary_conserve, stat=ierr)
            call check_deallocate('ueg_ternary_conserve', ierr)
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
                    intgrl = intgrl + coulomb_int_ueg(sys%lattice%box_length(1), bg, i, a)

                ! Exchange
                if (bg%basis_fns(i)%ms == bg%basis_fns(b)%ms .and.  bg%basis_fns(j)%ms == bg%basis_fns(a)%ms) &
                    intgrl = intgrl - coulomb_int_ueg(sys%lattice%box_length(1), bg, i, b)

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

end module ueg_system
