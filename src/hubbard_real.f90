module hubbard_real

! Real space formulation of the Hubbard model.

use const

implicit none

! The kinetic term is constant in the real space formulation:
! only the connectivity of the lattice matters.
! tmat(:,i) is a bit string.  The j-th bit corresponding to a basis function
! (as given by bit_lookup) is set if i and j are connected.
! We need to distinguish between connections within the cell and those due to
! periodic boundaries.  We do this by the following strategy:
!   a) j>i.
!          If the j-th bit is set then i and j are connected within the crystal
!          cell.
!   b) j<=i.
!          If the i-th bit of tmat(:,j) is set, then i and j are connected due
!          to periodic boundary conditions.
! This may seem like a somewhat arbitrary choice, but it enables for the
! correct evaluation of the kinetic energy using bit operations.
! Further it enables us to pick up cases such as the 2x2 (non-tilted) system,
! where a site is connected to a different site and that site's periodic image.
integer(i0), allocatable :: tmat(:,:) ! (basis_length, nbasis)

! Orbitals i and j are connected if the j-th bit of connected_orbs(:,i) is
! set.  This is a bit like tmat but without a bit set for a site being its own
! periodic image.  This is useful in FCIQMC for generating random
! excitations.
integer(i0), allocatable :: connected_orbs(:,:) ! (basis_length, nbasis)

! connected_sites(:,i) contains the list of sites connected to site i (ie is the
! decoded/non-bit list form of connected_orbs).
! If connected_orbs(j,i) is 0 then it means there are fewer than 2ndim sites
! that are connected to i that are not a periodic image of i.
! For the triangular lattice, there are 3ndim bonds, and ndim must equal 2,
! so each site is connected to 6.
integer, allocatable :: connected_sites(:,:) ! (2ndim, nbasis) or (6, nbasis)

! True if any site is its own periodic image.
! This is the case if one dimension (or more) has only one site per crystalisystem
! cell.  If so then the an orbital can incur a kinetic interaction with itself.
! This is the only way that the integral < i | T | i >, where i is a basis
! function centred on a lattice site, can be non-zero.
logical :: t_self_images

! True if we are actually only modelling a finite system (e.g. a H_2 molecule)
! False if we are modelling an infinite lattice
! The code is set up to model inifinite lattices by default, however in order
! to model only a finite "cluster" of sites, all one need do is set the 
! connection matrix elements corresponding to connections accross cell 
! boundaries (i.e. periodic boundary conditions) to 0
logical :: finite_cluster = .false. ! default to infinite crystals


contains

    subroutine init_real_space()

        ! Initialise real space Hubbard model and Heisenberg model: find and store 
        ! the matrix elements < i | T | j > where i and j are real space basis functions.

        use basis, only: nbasis, bit_lookup, basis_lookup, basis_length, basis_fns, set_orb
        use determinants, only: decode_det
        use system, only: lattice, ndim, box_length, system_type
        use system, only: heisenberg, triangular_lattice
        use bit_utils
        use checking, only: check_allocate
        use errors, only: stop_all
        use parallel, only: parent

        integer :: i, j, k, ierr, pos, ind, ivec, v, isystem
        integer :: basis_find, row_1, row_2
        integer :: r(ndim)

        integer :: lvecs(ndim, 3**ndim)
        integer :: difference_vec(2), shifted_vec(2), unshifted_vec(2)
        
        logical :: connected = .false.

        t_self_images = any(abs(box_length-1.0_p) < depsilon)

        allocate(tmat(basis_length,nbasis), stat=ierr)
        call check_allocate('tmat',basis_length*nbasis,ierr)
        allocate(connected_orbs(basis_length,nbasis), stat=ierr)
        call check_allocate('connected_orbs',basis_length*nbasis,ierr)
        if (triangular_lattice) then
            allocate(connected_sites(6,nbasis), stat=ierr)
            call check_allocate('connected_sites',basis_length*6,ierr)
        else
            allocate(connected_sites(2*ndim,nbasis), stat=ierr)
            call check_allocate('connected_sites',basis_length*2*ndim,ierr)
        end if

        tmat = 0
        connected_orbs = 0
        
        ! For the Hubbard model, each orbital can have spin up or down, so 
        ! basis_fns(i) refers to alternating alpha and beta orbitals.
        ! In the do loop we therefore loop over every *second* orbital (because  
        ! spin must be the same for orbitals to be connected in this case).
        ! For Heisenberg, we just want to loop over every component of
        ! basis_fns, so we set isystem = 1
        if (system_type == heisenberg) then
            isystem = 1
        else
            isystem = 2
        endif
            

        ! Form all lattice vectors
        select case(ndim)
        case(1)
            do i = -1, 1
                lvecs(:,i+2) = i*lattice(:,1)
            end do
        case(2)
            do i = -1, 1
                do j = -1, 1
                    lvecs(:,j+2+3*(i+1)) = i*lattice(:,1) + j*lattice(:,2)
                end do
            end do
        case(3)
            do i = -1, 1
                do j = -1, 1
                    do k = -1, 1
                        lvecs(:,k+2+3*(j+1)+9*(i+1)) = i*lattice(:,1) + j*lattice(:,2) + k*lattice(:,3)
                    end do
                end do
            end do
        end select

        ! Construct how the lattice is connected.
        do i = 1, nbasis-(isystem-1), isystem
            do j = i, nbasis-(isystem-1), isystem
                ! Loop only over one spin: the other spin is identical so can be
                ! filled in automatically.
                ! All matrix elements between different spins are zero
                ! Allow j=i in case i is its own periodic image.
                r = basis_fns(i)%l - basis_fns(j)%l
                do ivec = 1, 3**ndim
                    if (sum(abs(r-lvecs(:,ivec))) == 1) then
                        ! i and j are on sites which are nearest neighbours 
                        if (all(lvecs(:,ivec) == 0)) then
                            ! Nearest neighbours within unit cell.
                            call set_orb(tmat(:,i),j)
                            if (isystem == 2) call set_orb(tmat(:,i+1),j+1)
                        else if (.not. finite_cluster) then ! if we want inf. lattice
                            ! Nearest neighbours due to periodic boundaries.
                            call set_orb(tmat(:,j),i)
                            if (isystem == 2) call set_orb(tmat(:,j+1),i+1) 
                            ! else we just want connections to other cells to
                            ! stay as 0 
                        end if        
                       
                        ! If we only want a discrete molecule and the lattice
                        ! vector connecting the 2 sites is the 0-vector then the
                        ! 2 sites are connected in a unit cell and thus are
                        ! actually connected. (If they "connect" accross cell
                        ! boundaries then they are not connected for a single
                        ! molecule).
                        if ( (finite_cluster .and. all(lvecs(:,ivec) == 0)) .or. &
                             .not. finite_cluster) then
                            if (i /= j) then
                                ! connected_orbs does not contain self-connections 
                                ! due to the periodic boundary conditions.
                                call set_orb(connected_orbs(:,i),j)
                                if (isystem == 2) call set_orb(connected_orbs(:,i+1),j+1)                      
                                call set_orb(connected_orbs(:,j),i)
                                if (isystem == 2) call set_orb(connected_orbs(:,j+1),i+1)
                            end if
                        end if
                    end if
                    
                    if (triangular_lattice) then
                        
                        ! We want to find all the connected sites for the triangular lattice.
                        ! We already have some connections - all connections from the
                        ! square lattice are still present in the triangular lattice,
                        ! just with some potential extra ones. We find these below.
                        
                        ! We can treat a triangular lattice as a 2d rectangular lattice
                        ! by taking a triangular lattice (with n rows of m columns) and
                        ! shifting every other row across, to get a n by m rectangular
                        ! lattice. The only difference to a rectangular lattice are the
                        ! extra connections. So long as we keep these, we still have a
                        ! have a triangular lattice for all purposes.
                        
                        ! The picture below shows (very roughly!) the various connections
                        ! for a (4 by 4) triangular lattice. The extra conectons are the 
                        ! zig-zag ones.
                        !    _ _ _
                        !   |/|/|/|
                        !   |\|\|\|
                        !   |/|/|/|
                        !   |\|\|\|

                        ! Once all the appropriate vectors have been taken away,
                        ! we will want their difference to be the vector below. In this case,
                        ! the two sites will be connected.
                        difference_vec = (/1,0/)
                        ! Depending on whether we are on an odd or even numbered row (if we
                        ! number the rows 1,2,3...) the extra connected sites will be in different
                        ! positions relative to the site - either both to the right of it
                        ! or both to the left of it (see picture above).
                        ! It is important to distinguish between these two cases, by taking
                        ! mod(row_1,2) and mod(row_2,2) to see which of the two cases each row is.
                        row_1 = basis_fns(i)%l(1)-basis_fns(1)%l(1)
                        row_2 = basis_fns(j)%l(1)-basis_fns(1)%l(1)
                        
                        ! If the two sites are on these different types of rows, the corresponding
                        ! sites may be connected.
                        ! If on same types, they won't be (atleast not via the new connections
                        ! which we are adding)...
                        if (mod(row_1,2) == 0 .and. mod(row_2,2) == 0) then
                            connected = .false.
                        ! If the first sites is on the correct row type, we shift it, and also
                        ! take the lvec away to allow for periodic boundaries.
                        else if (mod(row_1,2) == 0) then
                            shifted_vec = basis_fns(i)%l
                            shifted_vec(2) = shifted_vec(2) - 1
                            shifted_vec = shifted_vec - lvecs(:,ivec)
                            unshifted_vec = basis_fns(j)%l
                            connected = .true.
                        ! If the second site is on the correct row type, shift this.
                        else if (mod(row_2,2) == 0) then
                            shifted_vec = basis_fns(j)%l
                            shifted_vec(2) = shifted_vec(2) - 1
                            shifted_vec = shifted_vec - lvecs(:,ivec)
                            unshifted_vec = basis_fns(i)%l
                            connected = .true.
                        end if

                        ! If the two are on the correct rows to have a possible extra
                        ! connection on the triangular lattice then...
                        if (connected) then
                            ! If connected not through boundary conditions:
                            if (all(lvecs(:,ivec) == 0)) then
                                ! If, after the shifting, the sites are vertically above each
                                ! other by one position, then there is an extra connection.       
                                if (sum(abs(shifted_vec-unshifted_vec)-difference_vec) == 0) then
                                    call set_orb(connected_orbs(:,i),j)
                                    if (isystem == 2) call set_orb(connected_orbs(:,i+1),j+1)                      
                                    call set_orb(connected_orbs(:,j),i)
                                    if (isystem == 2) call set_orb(connected_orbs(:,j+1),i+1)
                                    end if
                            ! If connected through boundary conditions, and boundary conditions
                            ! are turned on:
                            else if (.not.finite_cluster) then
                                if (sum(abs(shifted_vec-unshifted_vec)-difference_vec) == 0) then
                                    call set_orb(connected_orbs(:,i),j)
                                    if (isystem == 2) call set_orb(connected_orbs(:,i+1),j+1)                      
                                    call set_orb(connected_orbs(:,j),i)
                                    if (isystem == 2) call set_orb(connected_orbs(:,j+1),i+1)
                                end if
                            end if 
                        end if
                        ! For triangular lattice, just set tmat and connected orbs to be the same...
                        tmat = connected_orbs
                    end if
                    
                end do
            end do
        end do

        ! Decode connected_orbs to store list of connections.
        connected_sites = 0
        do i = 1, nbasis
            v = 0
            do ind = 1, basis_length
                do pos = 0, i0_end
                    if (btest(connected_orbs(ind,i), pos)) then
                        v = v + 1
                        connected_sites(v, i) = basis_lookup(pos, ind)
                    end if
                end do
            end do
        end do
        
    end subroutine init_real_space

    subroutine end_real_space()

        ! Clean up hubbard_real specific allocations.

        use checking, only: check_deallocate

        integer :: ierr

        if (allocated(tmat)) then
            deallocate(tmat, stat=ierr)
            call check_deallocate('tmat',ierr)
        end if
        if (allocated(connected_orbs)) then
            deallocate(connected_orbs, stat=ierr)
            call check_deallocate('connected_orbs',ierr)
        end if
        if (allocated(connected_sites)) then
            deallocate(connected_sites, stat=ierr)
            call check_deallocate('connected_sites',ierr)
        end if

    end subroutine end_real_space

    elemental function get_one_e_int_real(i, j) result(one_e_int)

        ! In:
        !    i: index of a real-space basis function.
        !    j: index of a real-space basis function.
        ! Returns:
        !    <phi1 | T | phi2> where T is the kinetic energy operator.

        use basis, only: basis_fns, bit_lookup
        use system, only: hubt

        real(p) :: one_e_int
        Integer, intent(in) ::  i,j
        integer :: ind, pos

        one_e_int = 0.0_p

        ! Need to check if i and j are on sites which are nearest neighbours
        ! either directly or due to periodic boundary conditions.
        pos = bit_lookup(1,j)
        ind = bit_lookup(2,j)
        ! Test if i <-> j.  If so there's a kinetic interaction.
        if (btest(tmat(ind,i),pos)) one_e_int = one_e_int - hubt
        pos = bit_lookup(1,i)
        ind = bit_lookup(2,i)
        ! Test if i <-> j.  If so there's a kinetic interaction.
        if (btest(tmat(ind,j),pos)) one_e_int = one_e_int - hubt

    end function get_one_e_int_real

    pure function get_coulomb_matel_real(f) result(umatel)

        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant, D.
        ! Returns:
        !    The matrix element < D | U | D >
        !    Note < D1 | U | D2 > = 0 if D1/=D2 within the real space
        !    formulation of the Hubbard model.

        use basis
        use system, only: hubu
        use bit_utils, only: count_set_bits
        use determinants, only: beta_mask, separate_strings

        real(p) :: umatel
        integer(i0), intent(in) :: f(basis_length)
        integer :: i
        integer(i0) :: b

        ! < D | U | D > = U*number of doubly occupied sites.
        if (separate_strings) then
            ! Just need to AND the alpha string with the beta string.
            umatel = sum(count_set_bits(iand(f(:basis_length/2),f(basis_length/2+1:))))
        else
            ! 1. Find the bit string representing the occupied beta orbitals.
            ! 2. Right shift it by one place.  The beta orbitals now line up with
            !    alpha orbitals.
            ! 3. AND the shifted beta bit string with the original bit string
            !    representing the list of occupied orbitals in the determinant.
            ! 4. The non-zero bits represent a sites which have both alpha and beta
            !    orbitals occupied.
            ! 5. Hence < D | U | D >.
            umatel = 0.0_p
            do i = 1, basis_length
                b = iand(f(i), beta_mask)
                umatel = umatel + count_set_bits(iand(f(i), ishft(b,-1)))
            end do
        end if
        umatel = hubu*umatel

    end function get_coulomb_matel_real

end module hubbard_real
