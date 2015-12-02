module real_lattice

! Real space formulation of the Hubbard model.

use const

implicit none

contains

    subroutine init_real_space(sys)

        ! Initialise real space Hubbard model and Heisenberg model: find and store
        ! the matrix elements < i | T | j > where i and j are real space basis functions.

        ! In/Out:
        !    sys: system to be studied.  On output the symmetry components are set.

        use basis, only: set_orb
        use calc, only: doing_dmqmc_calc, dmqmc_energy_squared
        use determinants, only: decode_det
        use system
        use bit_utils
        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all

        type(sys_t), intent(inout) :: sys
        integer :: i, j, k, ierr, pos, ind, ivec, v, isystem
        integer :: r(sys%lattice%ndim), bit_element, bit_pos, site_index
        logical :: diag_connection

        integer, allocatable :: lvecs(:,:)
        integer :: lattice_size(3)

        sys%nsym = 1
        sys%sym0 = 1
        sys%sym_max = 1
        sys%nsym_tot = 1
        sys%sym0_tot = 1
        sys%sym_max_tot = 1

        associate(sl=>sys%lattice, sr=>sys%real_lattice)

            sr%t_self_images = any(abs(sl%box_length-1.0_p) < depsilon)

            allocate(sr%tmat(sys%basis%string_len,sys%basis%nbasis), stat=ierr)
            call check_allocate('sr%tmat',sys%basis%string_len*sys%basis%nbasis,ierr)
            allocate(sr%connected_orbs(sys%basis%string_len,sys%basis%nbasis), stat=ierr)
            call check_allocate('sr%connected_orbs',sys%basis%string_len*sys%basis%nbasis,ierr)
            allocate(lvecs(sl%ndim,3**sl%ndim), stat=ierr)
            call check_allocate('lvecs', size(lvecs), ierr)
            if (sl%triangular_lattice) then
                allocate(sr%connected_sites(0:3*sl%ndim,sys%basis%nbasis), stat=ierr)
                call check_allocate('sr%connected_sites', size(sr%connected_sites), ierr)
            else
                allocate(sr%connected_sites(0:2*sl%ndim,sys%basis%nbasis), stat=ierr)
                call check_allocate('sr%connected_sites', size(sr%connected_sites), ierr)
            end if
            if (doing_dmqmc_calc(dmqmc_energy_squared)) then
                allocate(sr%next_nearest_orbs(sys%basis%nbasis,sys%basis%nbasis), stat=ierr)
                call check_allocate('sr%next_nearest_orbs',sys%basis%nbasis*sys%basis%nbasis,ierr)
            end if

            sr%tmat = 0_i0
            sr%connected_orbs = 0_i0

            ! For the Hubbard model, each orbital can have spin up or down, so
            ! basis_fns(i) refers to alternating alpha and beta orbitals.
            ! In the do loop we therefore loop over every *second* orbital (because
            ! spin must be the same for orbitals to be connected in this case).
            ! For Heisenberg and Chung--Landau models, we just want to loop over
            ! every component of basis_fns, so we set isystem = 1
            select case(sys%system)
            case(heisenberg, chung_landau)
                isystem = 1
            case default
                isystem = 2
            end select

            call enumerate_lattice_vectors(sl, lvecs)

            ! Construct how the sl%lattice is connected.
            diag_connection = .false. ! For sl%ndim /= 2.
            do i = 1, sys%basis%nbasis-(isystem-1), isystem
                do j = i, sys%basis%nbasis-(isystem-1), isystem
                    ! Loop only over one spin: the other spin is identical so can be
                    ! filled in automatically.
                    ! All matrix elements between different spins are zero
                    ! Allow j=i in case i is its own periodic image.
                    r = sys%basis%basis_fns(i)%l - sys%basis%basis_fns(j)%l
                    do ivec = 1, 3**sl%ndim
                        ! For the triangular sl%lattice, there are extra diagonal bonds between pairs
                        ! of sites which obey this condition.
                        if (sl%ndim == 2) then
                            diag_connection = all((r-lvecs(:,ivec)) == (/1,1/)) .or. &
                                              all((r-lvecs(:,ivec)) == (/-1,-1/))
                        end if
                        if (sum(abs(r-lvecs(:,ivec))) == 1 .or. &
                            (sl%triangular_lattice .and. diag_connection)) then
                            ! i and j are on sites which are nearest neighbours
                            if (all(lvecs(:,ivec) == 0)) then
                                ! Nearest neighbours within unit cell.
                                call set_orb(sys%basis%bit_lookup,sr%tmat(:,i),j)
                                if (isystem == 2) call set_orb(sys%basis%bit_lookup,sr%tmat(:,i+1),j+1)
                            else if (.not. sr%finite_cluster) then ! if we want inf. sl%lattice
                                ! Nearest neighbours due to periodic boundaries.
                                call set_orb(sys%basis%bit_lookup,sr%tmat(:,j),i)
                                if (isystem == 2) call set_orb(sys%basis%bit_lookup,sr%tmat(:,j+1),i+1)
                                ! else we just want connections to other cells to
                                ! stay as 0
                            end if

                            ! If we only want a discrete molecule and the sl%lattice
                            ! vector connecting the 2 sites is the 0-vector then the
                            ! 2 sites are connected in a unit cell and thus are
                            ! actually connected. (If they "connect" across cell
                            ! boundaries then they are not connected for a single
                            ! molecule).
                            if ( (sr%finite_cluster .and. all(lvecs(:,ivec) == 0)) .or. &
                                 .not. sr%finite_cluster) then
                                if (i /= j) then
                                    ! sr%connected_orbs does not contain self-connections
                                    ! due to the periodic boundary conditions.
                                    call set_orb(sys%basis%bit_lookup,sr%connected_orbs(:,i),j)
                                    if (isystem == 2) call set_orb(sys%basis%bit_lookup,sr%connected_orbs(:,i+1),j+1)
                                    call set_orb(sys%basis%bit_lookup,sr%connected_orbs(:,j),i)
                                    if (isystem == 2) call set_orb(sys%basis%bit_lookup,sr%connected_orbs(:,j+1),i+1)
                                end if
                            end if
                        end if

                    end do
                end do
            end do

            if (allocated(sr%next_nearest_orbs)) call create_next_nearest_orbs(sys%basis, sr)

            ! Decode sr%connected_orbs to store list of connections.
            sr%connected_sites = 0
            do i = 1, sys%basis%nbasis
                v = 0
                do ind = 1, sys%basis%string_len
                    do pos = 0, i0_end
                        if (btest(sr%connected_orbs(ind,i), pos)) then
                            v = v + 1
                            sr%connected_sites(v, i) = sys%basis%basis_lookup(pos, ind)
                        end if
                    end do
                end do
                sr%connected_sites(0,i) = v
            end do

        end associate

        select case(sys%system)
        case (heisenberg)
            ! This is the number of bonds for an arbitrary lattice
            ! with periodic or fixed end boundary conditions.
            sys%heisenberg%nbonds = sum(sys%real_lattice%connected_sites(0,:))/2
            ! Find lattice_mask for a gerenal bipartite lattice.
            if (sys%lattice%bipartite_lattice) then
                allocate (sys%heisenberg%lattice_mask(sys%basis%string_len), stat=ierr)
                associate(lattice_mask=>sys%heisenberg%lattice_mask)
                    call check_allocate('lattice_mask',sys%basis%string_len,ierr)
                    ! lattice_size is such that any loops over higher dimensions
                    ! than that of the model are single iterations but allows us to not have
                    ! to handle each dimension separately.
                    lattice_size = 1
                    lattice_size(1) = ceiling(sys%lattice%box_length(1))
                    if (sys%lattice%ndim > 1) lattice_size(2) = ceiling(sys%lattice%box_length(2))
                    if (sys%lattice%ndim > 2) lattice_size(3) = ceiling(sys%lattice%box_length(3))

                    lattice_mask = 0_i0
                    do k = 1, lattice_size(3)
                        do j = 1, lattice_size(2)
                            do i = 1, lattice_size(1),2
                                site_index = (lattice_size(2)*lattice_size(1))*(k-1) + &
                                              lattice_size(1)*(j-1) + mod(j+k,2) + i
                                bit_pos = sys%basis%bit_lookup(1, site_index)
                                bit_element = sys%basis%bit_lookup(2, site_index)
                                lattice_mask(bit_element) = ibset(lattice_mask(bit_element), bit_pos)
                            end do
                        end do
                    end do
                end associate
            end if
        end select

        deallocate(lvecs, stat=ierr)
        call check_deallocate('lvecs', ierr)

    end subroutine init_real_space

    subroutine end_real_space(sh)

        ! Clean up real_lattice specific allocations.

        ! In/Out:
        !    sh: Heisenberg system object to be deallocated.

        use checking, only: check_deallocate
        use system, only: sys_heisenberg_t

        type(sys_heisenberg_t), intent(inout) :: sh

        integer :: ierr

        if (allocated(sh%lattice_mask)) then
            deallocate(sh%lattice_mask, stat=ierr)
            call check_deallocate('sh%lattice_mask', ierr)
        end if

    end subroutine end_real_space

    subroutine enumerate_lattice_vectors(sl, lvecs)

        ! Enumerate combinations of lattice vectors which define the Wigner--Seitz cell.

        ! In:
        !    sl: sys_lattice_t object describing the real-space lattice.
        ! Out:
        !    lvecs: all possible primitive combinations of the above lattice vectors,
        !         where the amplitude for each lattice vector can be either -1, 0 or +1.
        !         lvec(:,i) stores the i'th such combination.  lvecs must have dimension
        !         (at least) (ndim, 3^ndim).

        use system, only: sys_lattice_t

        type(sys_lattice_t), intent(in) :: sl
        integer, intent(out) :: lvecs(:,:)

        integer :: i, j, k

        select case(sl%ndim)
        case(1)
            do i = -1, 1
                lvecs(:,i+2) = i*sl%lattice(:,1)
            end do
        case(2)
            do i = -1, 1
                do j = -1, 1
                    lvecs(:,j+2+3*(i+1)) = i*sl%lattice(:,1) + j*sl%lattice(:,2)
                end do
            end do
        case(3)
            do i = -1, 1
                do j = -1, 1
                    do k = -1, 1
                        lvecs(:,k+2+3*(j+1)+9*(i+1)) = i*sl%lattice(:,1) + j*sl%lattice(:,2) + k*sl%lattice(:,3)
                    end do
                end do
            end do
        end select

    end subroutine enumerate_lattice_vectors

    elemental function get_one_e_int_real(sys, i, j) result(one_e_int)

        ! In:
        !    sys: system being studied.
        !    i: index of a real-space basis function.
        !    j: index of a real-space basis function.
        ! Returns:
        !    <phi1 | T | phi2> where T is the kinetic energy operator.

        use system, only: sys_t

        real(p) :: one_e_int
        type(sys_t), intent(in) :: sys
        Integer, intent(in) ::  i,j
        integer :: ind, pos

        one_e_int = 0.0_p

        ! Need to check if i and j are on sites which are nearest neighbours
        ! either directly or due to periodic boundary conditions.
        pos = sys%basis%bit_lookup(1,j)
        ind = sys%basis%bit_lookup(2,j)
        ! Test if i <-> j.  If so there's a kinetic interaction.
        if (btest(sys%real_lattice%tmat(ind,i),pos)) one_e_int = one_e_int - sys%hubbard%t
        pos = sys%basis%bit_lookup(1,i)
        ind = sys%basis%bit_lookup(2,i)
        ! Test if i <-> j.  If so there's a kinetic interaction.
        if (btest(sys%real_lattice%tmat(ind,j),pos)) one_e_int = one_e_int - sys%hubbard%t

    end function get_one_e_int_real

    pure function get_coulomb_matel_real(sys, f) result(umatel)

        ! In:
        !    sys: system being studied.
        !    f(string_len): bit string representation of the Slater
        !        determinant, D.
        ! Returns:
        !    The matrix element < D | U | D >
        !    Note < D1 | U | D2 > = 0 if D1/=D2 within the real space
        !    formulation of the Hubbard model.

        use basis
        use system, only: sys_t
        use bit_utils, only: count_set_bits

        real(p) :: umatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer :: i
        integer(i0) :: b

        ! < D | U | D > = U*number of doubly occupied sites.
        associate(basis=>sys%basis)
            ! 1. Find the bit string representing the occupied beta orbitals.
            ! 2. Right shift it by one place.  The beta orbitals now line up with
            !    alpha orbitals.
            ! 3. AND the shifted beta bit string with the original bit string
            !    representing the list of occupied orbitals in the determinant.
            ! 4. The non-zero bits represent a sites which have both alpha and beta
            !    orbitals occupied.
            ! 5. Hence < D | U | D >.
            umatel = 0.0_p
            do i = 1, basis%string_len
                b = iand(f(i), basis%beta_mask)
                umatel = umatel + count_set_bits(iand(f(i), ishft(b,-1)))
            end do
            umatel = sys%hubbard%u*umatel
        end associate

    end function get_coulomb_matel_real

    subroutine create_next_nearest_orbs(basis, sr)

        ! Create the list of next nearest orbitals for each orbital.

        ! In:
        !    basis: basis set info.
        ! In/Out:
        !    sr: sys_real_lattice_t.  On input sr%connected_orbs must be set and
        !        sr%next_nearest_orbs must be allocated.  On output
        !        sr%next_nearest_orbs is filled in.

        use basis_types, only: basis_t
        use parallel
        use system, only: sys_real_lattice_t

        type(basis_t), intent(in) :: basis
        type(sys_real_lattice_t), intent(inout) :: sr
        integer :: ibasis, jbasis, kbasis
        integer :: bit_position, bit_element

        sr%next_nearest_orbs = 0_i0

        do ibasis = 1, basis%nbasis
            do jbasis = 1, basis%nbasis
                bit_position = basis%bit_lookup(1,jbasis)
                bit_element = basis%bit_lookup(2,jbasis)
                if (btest(sr%connected_orbs(bit_element,ibasis),bit_position)) then
                    do kbasis = 1, basis%nbasis
                        bit_position = basis%bit_lookup(1,kbasis)
                        bit_element = basis%bit_lookup(2,kbasis)
                        if (btest(sr%connected_orbs(bit_element,jbasis),bit_position)) then
                            sr%next_nearest_orbs(ibasis,kbasis) = sr%next_nearest_orbs(ibasis,kbasis)+1
                        end if
                    end do
                end if
            end do
            sr%next_nearest_orbs(ibasis,ibasis) = 0_i0
        end do

    end subroutine create_next_nearest_orbs

    subroutine find_translational_symmetry_vecs(sys, sym_vecs, nsym)

        ! This routine will find all symmetry vectors for the lattice provided
        ! and return them in sym_vecs.

        ! In:
        !     sys: system being studied.
        ! In/Out:
        !     sym_vecs: An array which on output will hold all translational
        !         symmetry vectors. Should be deallocated on input.
        ! Out:
        !     nsym: The total number of symmetry vectors.

        use checking, only: check_allocate
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        real(p), allocatable, intent(inout) :: sym_vecs(:,:)
        integer, intent(out) :: nsym
        integer :: i, j, k, ierr
        integer :: nvecs(3)
        real(p) :: v(sys%lattice%ndim), test_vec(sys%lattice%ndim)
        integer :: scale_fac

        ! The maximum number of translational symmetry vectors is nsites (for
        ! the case of a non-tilted lattice), so allocate this much storage.
        allocate(sym_vecs(sys%lattice%ndim,sys%lattice%nsites),stat=ierr)
        call check_allocate('sym_vecs',sys%lattice%ndim*sys%lattice%nsites,ierr)
        sym_vecs = 0

        ! The number of symmetry vectors in each direction.
        nvecs = 0
        ! The total number of symmetry vectors.
        nsym = 0

        do i = 1, sys%lattice%ndim
            scale_fac = maxval(abs(sys%lattice%lattice(:,i)))
            v = real(sys%lattice%lattice(:,i),p)/real(scale_fac,p)

            do j = 1, scale_fac-1
                test_vec = v*j
                ! If test_vec has all integer components.
                if (all(.not. (abs(test_vec-real(nint(test_vec),p)) > 0.0_p) )) then
                    ! If this condition is obeyed then this is a symmetry vector,
                    ! so store it.
                    nvecs(i) = nvecs(i) + 1
                    nsym = nsym + 1
                    sym_vecs(:,nsym) = test_vec
                end if
            end do
        end do

        ! Next, add all combinations of the above generated vectors to form a closed group.

        ! Add all pairs of the above vectors.
        do i = 1, nvecs(1)
            do j = nvecs(1)+1, sum(nvecs)
                nsym = nsym + 1
                sym_vecs(:,nsym) = sym_vecs(:,i)+sym_vecs(:,j)
            end do
        end do
        do i = nvecs(1)+1, nvecs(1)+nvecs(2)
            do j = nvecs(1)+nvecs(2)+1, sum(nvecs)
                nsym = nsym + 1
                sym_vecs(:,nsym) = sym_vecs(:,i)+sym_vecs(:,j)
            end do
        end do

        ! Add all triples of the above vectors.
        do i = 1, nvecs(1)
            do j = nvecs(1)+1, nvecs(1)+nvecs(2)
                do k = nvecs(1)+nvecs(2)+1, sum(nvecs)
                    nsym = nsym + 1
                    sym_vecs(:,nsym) = sym_vecs(:,i)+sym_vecs(:,j)+sym_vecs(:,k)
                end do
            end do
        end do

        ! Include the identity transformation vector in the first slot.
        sym_vecs(:,2:nsym+1) = sym_vecs(:,1:nsym)
        sym_vecs(:,1) = 0
        nsym = nsym + 1

    end subroutine find_translational_symmetry_vecs

    subroutine map_vec_to_cell(nbasis, basis_fns, ndim, lvecs, r)

        ! Map a vector, r, outside from outside to inside the simulation cell.
        ! This subroutine assumes that the site specified by r is outside the cell
        ! by no more than one lattice vector, along each lattice vector.

        ! In:
        !    nbasis: number of basis functions
        !    basis_fns: set of one-particle (spin) basis functions.
        !    ndim: dimensionality of the lattice.
        !    lvecs: all 3**ndim possible lattice vectors in the nearest 'shell'
        !       (ie all integer combinations from -1 to 1 for each lattice vector).
        ! In/Out:
        !    r: On output the site specified by r (in units of the lattice
        !       sites) is mapped into the equivalent site inside the simulation
        !       cell.

        use basis_types, only: basis_fn_t

        integer, intent(in) :: nbasis, ndim, lvecs(ndim, 3**ndim)
        type(basis_fn_t), intent(in) :: basis_fns(:)
        integer, intent(inout) :: r(ndim) 
        integer :: v(ndim)
        integer :: i, j

        do i = 1, 3**ndim
            ! Add all combinations of lattice vectors (stored in lvecs).
            v = r + lvecs(:,i)
            do j = 1, nbasis
                ! Loop over all basis functions and check if the shifted vector is
                ! now the same as any of these vectors. If so, it is in the cell,
                ! so keep it and return.
                if (all(v == basis_fns(j)%l)) then
                    r = v
                    return
                end if
            end do
        end do

    end subroutine map_vec_to_cell

end module real_lattice
