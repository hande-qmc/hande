module read_in

use const

implicit none

contains

    subroutine read_in_fcidump

        ! Read in a FCIDUMP file, which contains the kinetic and coulomb
        ! integrals, one-particle eigenvalues and other system information.
        ! File format (partially) defined in Comp. Phys. Commun. 54 (1989) 75.
        ! See also notes below.

        use system, only: fcidump

        use utils, only: get_free_unit
        use errors, only: stop_all

        ! System data
        integer :: norb, nelec, ms2, orbsym(1000), isym
        logical :: uhf = .false.

        ! Integrals
        integer :: i, j, a, b
        real(dp) :: x

        ! reading in...
        integer :: ir, ios
        logical :: t_exists

        namelist /FCI/ norb, nelec, ms2, orbsym, uhf

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
        !    symmetry notes.
        !  * UHF: true if FCIDUMP file was produced from an unrestricted
        !    Hartree-Fock calculation.  See note on basis indices below.
        !  * ISYM: unused.  Defined solely for compatibility with NECI FCIDUMP files.

        ! Integrals:
        !  * if i = j = a = b = 0, E_core = x , where E_core contains the
        !    nuclear-nuclear and other non-electron contributions to the
        !    Hamiltonian.
        !  * if a = j = b = 0, \epsilon_i = x, the single-particle eigenvalue
        !    of the i-th orbital.
        !  * if j = b = 0, < i | T | a > = x, the one-body Hamiltonian matrix element
        !    between the i-th and a-th orbitals.
        !  * otherwise < i j | 1/r_12 | a b > = x, the Coulomb integral between
        !    the i-a co-density and the j-b codensity.  Note the Coulomb
        !    integrals are given in Chemist's notation, a rare instance that
        !    Chemists are wrong and Physicists correct.

        ! Basis indices:
        ! RHF: All indices are in terms of spatial orbitals.  NORBS is the
        ! number of spatial orbitals.
        ! UHF: All indices are in terms of spin orbitals.  NORBS is the
        ! number of spin orbitals.
        ! Basis functions (as stored by basis_fns) are always stored as spin
        ! orbitals (the memory saving involved in storing only spatial orbitals
        ! is not worth the additional overhead/headache, as FCIQMC involves
        ! working in spin orbitals).  Integrals are expensive to store, so we
        ! store them in as compressed format as possible.

        ir = get_free_unit()
        inquire(file=fcidump, exist=t_exists)
        if (.not.t_exists) call stop_all('read_in_fcidump', 'FCIDUMP does not &
                                                           &exist:'//trim(fcidump))
        open (ir, file=fcidump, status='old', form='formatted', iostat=ios)

        ! read system data
        read (ir, FCI)

        if (uhf) then
            nbasis = norbs
        else
            nbasis = 2*norbs
        end if

        allocate(basis_fns(nbasis), stat=ierr)
        call check_allocate('basis_fns', nbasis, ierr)

        ! avoid annoying compiler warnings over unused variables in FCI namelist
        ! that are present for NECI compatibility.
        isym = 0

        ! read integrals and eigenvalues
        do 
            ! loop over lines.
            read (ir,*) x, i, a, j, b
            if (ios < 0) exit ! end of file
        end do

        if (ios > 0) call stop_all('read_input','Problem reading input.')
        close(ir, status='keep')

    end subroutine read_in_fcidump

end module read_in
