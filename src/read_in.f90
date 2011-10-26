module read_in_system

! Module for reading in and manipulating integral files which entirely define
! the Hamiltonian corresponding to a 'real' system.

use const

implicit none

contains

    subroutine read_in_fcidump

        ! Read in a FCIDUMP file, which contains the kinetic and coulomb
        ! integrals, one-particle eigenvalues and other system information.
        ! File format (partially) defined in Comp. Phys. Commun. 54 (1989) 75.
        ! See also notes below.

        use basis, only: nbasis, basis_fns, init_basis_fn
        use basis, only: write_basis_fn
        use molecular_integrals, only: init_molecular_integrals, store_one_e_int_mol, &
                                       store_two_e_int_mol
        use point_group_symmetry, only: init_pg_symmetry
        use system, only: fcidump, uhf

        use utils, only: get_free_unit
        use checking, only: check_allocate
        use errors, only: stop_all

        use, intrinsic :: iso_fortran_env

        ! System data
        ! We don't know how many orbitals we have until we read in the FCI
        ! namelist, so have to hardcode the array sizes.
        ! It's reasonably safe to assume that we'll never use more than 1000
        ! orbitals!
        integer :: norb, nelec, ms2, orbsym(1000), isym, syml(1000), symlz(1000)

        ! Integrals
        integer :: i, j, a, b
        real(dp) :: x

        ! reading in...
        integer :: ir, ios, ierr
        logical :: t_exists

        namelist /FCI/ norb, nelec, ms2, orbsym, uhf, isym, syml, symlz

        ! avoid annoying compiler warnings over unused variables in FCI namelist
        ! that are present for NECI compatibility.
        isym = 0
        syml = 0
        symlz = 0

        ! FCIDUMP file format is as follows:

        ! &FCI             ! FCI namelist.  See below.
        ! /                ! / terminates a namelist.  Most compilers also
        !                  ! implement the extension where &END is used to
        !                  ! terminate the namelist instead.
        ! x i a j b        ! x is a float, i, j, a and b are integers.

        ! &FCI namelist:
        !  * NORB: number of orbitals in the basis.  See note on basis indices below.
        !  * NELEC: number of electrons in system.
        !  * MS2: spin polarisation.
        !  * ORBSYM: array containing symmetry label of each orbital.  See
        !    symmetry notes below and in pg_symmetry.
        !  * UHF: true if FCIDUMP file was produced from an unrestricted
        !    Hartree-Fock calculation.  See note on basis indices below.
        !  * ISYM: currently unused.  Defined solely for compatibility with NECI
        !    FCIDUMP files.  Gives the symmetry of the wavefunction formed by
        !    occupied the NELEC lowest energy spin-orbitals.
        !  * SYML: currently unused.  Defined solely for compatibility with NECI
        !    FCIDUMP files.  Array containing L (angular momentum) for each orbital.
        !    Set to -1 if L is not a good quantum number.
        !  * SYMLZ: currently unused.  Defined solely for compatibility with NECI
        !    FCIDUMP files.  Array containing Lz (angular momentum along the
        !    z-axis) for each orbital.
        !    For example d_xz would have L=2 and Lz=1, and dyz L=2, Lz=-1.

        ! Integrals:
        !  * if i = j = a = b = 0, E_core = x , where E_core contains the
        !    nuclear-nuclear and other non-electron contributions to the
        !    Hamiltonian.
        !  * if a = j = b = 0, \epsilon_i = x, the single-particle eigenvalue
        !    of the i-th orbital.
        !  * if j = b = 0, < i | h | a > = x, the one-body Hamiltonian matrix element
        !    between the i-th and a-th orbitals, where h = T+V_ext.
        !  * otherwise < i j | 1/r_12 | a b > = x, the Coulomb integral between
        !    the i-a co-density and the j-b codensity.  Note the Coulomb
        !    integrals are given in Chemists' notation, a rare instance that
        !    Chemists are wrong and Physicists correct.

        ! Basis indices:
        ! RHF: All indices are in terms of spatial orbitals.  NORB is the
        ! number of spatial orbitals.
        ! UHF: All indices are in terms of spin orbitals.  NORB is the
        ! number of spin orbitals.
        ! Basis functions (as stored by basis_fns) are always stored as spin
        ! orbitals (the memory saving involved in storing only spatial orbitals
        ! is not worth the additional overhead/headache, as FCIQMC involves
        ! working in spin orbitals).  Integrals are expensive to store, so we
        ! store them in as compressed format as possible.

        ! Symmetry:
        ! Molecular orbitals are defined by the D2h point group (or a subgroup
        ! thereof)by the quantum chemistry packages (QChem, MOLPRO) used to
        ! produce FCIDUMP files , so we need only concern ourselves with Abelian
        ! symmetries.
        ! If ORBSYM(i) = 0, then the symmetry of the i-th orbital is not
        ! well-defined.  In this case, we can only resort to turning off all
        ! symmetry (i.e. set all orbitals to be totally symmetric).  Note that
        ! this has memory implications for the integral storage.
        ! ORBSYM(i) = S+1, where S is the symmetry label defining the
        ! irreducible representation spanned by the i-th orbital.
        ! See notes in pg_symmetry about the symmetry label for Abelian point
        ! groups.

        ir = get_free_unit()
        inquire(file=fcidump, exist=t_exists)
        if (.not.t_exists) call stop_all('read_in_fcidump', 'FCIDUMP does not &
                                                           &exist:'//trim(fcidump))
        open (ir, file=fcidump, status='old', form='formatted')

        ! read system data
        read (ir, FCI)

        if (uhf) then
            nbasis = norb
        else
            nbasis = 2*norb
        end if

        allocate(basis_fns(nbasis), stat=ierr)
        call check_allocate('basis_fns', nbasis, ierr)

        ! Set up basis functions.
        do i = 1, norb
            if (uhf) then
                if (mod(i,2) == 0) then
                    call init_basis_fn(basis_fns(i), sym=orbsym(i)-1, ms=-1)
                else
                    call init_basis_fn(basis_fns(i), sym=orbsym(i)-1, ms=1)
                end if
                ! Assume orbitals are ordered appropriately in FCIDUMP...
                basis_fns(i)%spatial_index = (i+1)/2
            else
                ! Need to initialise both up- and down-spin basis functions.
                call init_basis_fn(basis_fns(2*i-1), sym=orbsym(i)-1, ms=-1)
                call init_basis_fn(basis_fns(2*i), sym=orbsym(i)-1, ms=-1)
                basis_fns(2*i-1)%spatial_index = i
                basis_fns(2*i)%spatial_index = i
            end if
        end do

        ! Was a symmetry found for all basis functions?  If not, then we must
        ! turn symmetry off.
        if (minval(orbsym(:norb)) < 0) then
            write (6,'(1X,a62)') 'Unconverged symmetry found.  Turning point group symmetry off.'
            forall (i=1:nbasis) basis_fns(i)%sym = 0
        end if

        ! Set up symmetry information.
        call init_pg_symmetry()

        ! Initialise integral stores.
        call init_molecular_integrals()

        ! read integrals and eigenvalues
        ios = 0
        do
            ! loop over lines.
            read (ir,*, iostat=ios) x, i, a, j, b
            if (ios == iostat_end) exit ! reached end of file
            if (ios /= 0) call stop_all('read_input','Problem reading input.')
            if (.not.uhf) then
                ! Working in spin-orbitals but FCIDUMP is in spatial orbitals.
                ! Need to only store integrals in one spin-channel in RHF, so
                ! can just with (e.g.) beta orbitals.
                i = i*2
                j = j*2
                a = a*2
                b = b*2
            end if
            if (ios < 0) exit ! end of file
            if (i == 0 .and. j == 0 .and. a == 0 .and. b == 0) then
                ! Ecore
                Ecore = x
            else if (j == 0 .and. a == 0 .and. b == 0) then
                ! \epsilon_i
                if (uhf) then
                    basis_fns(i)%sp_eigv = x
                else
                    basis_fns(i-1)%sp_eigv = x
                    basis_fns(i)%sp_eigv = x
                end if
            else if (j == 0 .and. b == 0) then
                ! < i | h | a >
                call store_one_e_int_mol(i, a, x)
            else
                ! < i j | 1/r_12 | a b >
                call store_two_e_int_mol(i, j, a, b, x)
            end if
        end do

        close(ir, status='keep')

        do i = 1, nbasis
            call write_basis_fn(basis_fns(i), new_line=.true.)
        end do

    end subroutine read_in_fcidump

end module read_in_system
