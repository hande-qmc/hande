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

! True if any site is its own periodic image.
! This is the case if one dimension (or more) has only one site per crystal
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

    subroutine init_real_space_hub()

        ! Initialise real space Hubbard model: find and store the matrix
        ! elements < i | T | j > where i and j are real space basis functions.

        use basis, only: nbasis, bit_lookup, basis_length, basis_fns, set_orb
        use system, only: lattice, ndim, box_length
        use bit_utils
        use checking, only: check_allocate
        use errors, only: stop_all

        integer :: i, j, k, ierr, pos, ind, ivec
        integer :: r(ndim)

        integer :: lvecs(ndim, 3**ndim)

        t_self_images = any(abs(box_length-1.0_p) < depsilon)

        allocate(tmat(basis_length,nbasis), stat=ierr)
        call check_allocate('tmat',basis_length*nbasis,ierr)
        allocate(connected_orbs(basis_length,nbasis), stat=ierr)
        call check_allocate('connected_orbs',basis_length*nbasis,ierr)

        tmat = 0
        connected_orbs = 0

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
        do i = 1, nbasis-1, 2
            do j = i, nbasis-1, 2
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
                            call set_orb(tmat(:,i+1),j+1)
                            ! Nearest neighbours due to periodic boundaries.
                        else if(.not. finite_cluster) then ! if we want inf. lattice
                            call set_orb(tmat(:,j),i)
                            call set_orb(tmat(:,j+1),i+1) 
                            ! else we just want connections to other cells to
                            ! stay as 0 
                        end if        
                       
                        ! if we only want a discrete molecule and the lattice
                        ! vector connecting the 2 sites is the 0-vector then the
                        ! 2 sites are connected in a unit cell and thus are
                        ! actually connected. (If they "connect" accross cell
                        ! boundaries then they are not connected for a single
                        ! molecule).
                        if( (finite_cluster .and. all(lvecs(:,ivec) == 0)) .or. &
                             .not. finite_cluster) then
                            call set_orb(connected_orbs(:,i),j)
                            call set_orb(connected_orbs(:,i+1),j+1)                      
                            call set_orb(connected_orbs(:,j),i)
                            call set_orb(connected_orbs(:,j+1),i+1)
                        end if
 
                    end if
                end do
            end do
        end do

        ! Information for the random excitation generators.
        do i = 1, nbasis
            pos = bit_lookup(1,i)
            ind = bit_lookup(2,i)
            connected_orbs(ind,i) = ibclr(connected_orbs(ind,i),pos)
        end do

    end subroutine init_real_space_hub

    subroutine end_real_space_hub()

        ! Clean up hubbard_real specific allocations.

        integer :: ierr

        if (allocated(tmat)) deallocate(tmat, stat=ierr)

    end subroutine end_real_space_hub

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
        use determinants, only: beta_mask

        real(p) :: umatel
        integer(i0), intent(in) :: f(basis_length)
        integer :: i
        integer(i0) :: b

        ! < D | U | D > = U*number of doubly occupied sites.
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
        umatel = hubu*umatel

    end function get_coulomb_matel_real

end module hubbard_real
