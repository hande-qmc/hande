module read_in_system

! Module for reading in and manipulating integral files which entirely define
! the Hamiltonian corresponding to a 'real' system.

use const

implicit none

contains

    subroutine read_in_fcidump(store_info, cas_info)

        ! Read in a FCIDUMP file, which contains the kinetic and coulomb
        ! integrals, one-particle eigenvalues and other system information.
        ! File format (partially) defined in Comp. Phys. Commun. 54 (1989) 75.
        ! See also notes below.
        !
        ! In:
        !    store_info (optional): if true (default) then store the data read
        !    in.  Otherwise the basis defined by the FCIDUMP file is simply
        !    printed out.
        !    cas_info (optional): if present, then defines the complete active
        !    space.  Defailt: (N,M), where N is the total number of electrons
        !    and M is the number of active *spatial* orbitals.  The default thus
        !    uses the entire space available.

        use basis, only: basis_fn, nbasis, basis_fns, end_basis_fns, write_basis_fn, &
                         write_basis_fn_header
        use molecular_integrals, only: init_molecular_integrals, store_one_body_int_mol, &
                                       store_two_body_int_mol, zero_one_body_int_store, &
                                       get_one_body_int_mol, one_e_h_integrals, coulomb_integrals
        use point_group_symmetry, only: init_pg_symmetry
        use system, only: fcidump, uhf, ecore, nel

        use utils, only: get_free_unit, tri_ind_reorder, int_fmt
        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all

        use, intrinsic :: iso_fortran_env

        logical, intent(in), optional :: store_info
        integer, optional :: cas_info(2)
        logical :: t_store
        integer :: cas(2)

        ! System data
        ! We don't know how many orbitals we have until we read in the FCI
        ! namelist, so have to hardcode the array sizes.
        ! It's reasonably safe to assume that we'll never use more than 1000
        ! orbitals!
        integer :: norb, nelec, ms2, orbsym(1000), isym, syml(1000), symlz(1000)

        ! all basis functions, including inactive ones.
        type(basis_fn), allocatable :: all_basis_fns(:)

        ! Integrals
        integer :: i, j, a, b, ii, jj, aa, bb
        real(dp) :: x

        ! reading in...
        integer :: ir, ios, ierr
        logical :: t_exists
        integer :: active_basis_offset, rhf_fac
        integer, allocatable :: seen_ijij(:), seen_iaib(:,:)
        logical, allocatable :: seen_iha(:)

        namelist /FCI/ norb, nelec, ms2, orbsym, uhf, isym, syml, symlz

        ! avoid annoying compiler warnings over unused variables in FCI namelist
        ! that are present for NECI compatibility.
        isym = 0
        syml = 0
        symlz = 0

        if (present(store_info)) then
            t_store = store_info
        else
            t_store = .true.
        end if

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
        !    NOTE:
        !         We assume that in UHF calculations the number of spin-up basis
        !         functions is equal to the number of spin-down basis functions.
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
            rhf_fac = 1
        else
            nbasis = 2*norb
            rhf_fac = 2  ! need to double count some integrals.
        end if

        allocate(all_basis_fns(nbasis), stat=ierr)
        call check_allocate('all_basis_fns', nbasis, ierr)

        ! Set up basis functions, including those which are subsequently frozen.
        call init_basis_fns_read_in(norb, uhf, orbsym(1:norb), all_basis_fns) 

        if (present(cas_info)) then
            if (any(cas_info < 1)) then
                cas = (/ nel, nbasis/2 /)
            else
                cas = cas_info
            end if
        else
            cas = (/ nel, nbasis/2 /)
        end if

        ! From CAS work out the start of the active basis functions, the number
        ! of active basis functions and the number of active electrons.
        active_basis_offset = nel-cas(1) ! number of core *spin* orbitals
        ! Note that nbasis is spin-orbitals whereas cas(2)=M is in spatial orbitals
        ! (as we use the conventional CAS definition).
        nbasis = min(nbasis, 2*cas(2))
        nel = nel - ( nel - cas(1) )
        if (uhf) then
            norb =  nbasis
        else
            norb = nbasis/2
        end if

        allocate(basis_fns(nbasis), stat=ierr)
        call check_allocate('basis_fns', nbasis, ierr)

        ! Set up basis functions.
        call init_basis_fns_read_in(norb, uhf, orbsym(active_basis_offset/rhf_fac+1:norb+active_basis_offset/rhf_fac), basis_fns) 

        ! Was a symmetry found for all basis functions?  If not, then we must
        ! turn symmetry off.
        if (minval(orbsym(:norb)) < 0) then
            write (6,'(1X,a62)') 'Unconverged symmetry found.  Turning point group symmetry off.'
            forall (i=1:nbasis) basis_fns(i)%sym = 0
        end if

        ! Set up symmetry information.
        if (t_store) call init_pg_symmetry()

        ! Initialise integral stores.
        if (t_store) call init_molecular_integrals()

        ! Freezing core orbitals amounts to changing the Hamiltonian.  In
        ! particular:

        ! E_core = E_nuc + \sum_{i=1}^{N_fr} <i|h|i> + \sum_{i<j}^{N_fr} <ij|ij> - <ij|ji>
        ! <a|h'|b> = <a|h|b> + \sum_{i=1}^{N_fr} <ia|ib> - <ia|bi>

        ! where we use i,j to refer to frozen core orbitals and a,b to refer to
        ! active orbitals.

        ! See J. Chem. Phys. 62, 4764 (1975), An Introduction to
        ! Configuration Interaction Theory by C. David Sherrill 
        ! (http://vergil.chemistry.gatech.edu/notes/ci/ci.html) and
        ! Alex Thom's thesis.

        ! We don't wish to incur the overhead of ever storing integrals
        ! involving frozen orbitals so instead we simply accumulate the
        ! quantities above.  Hence we must zero them:

        Ecore = 0.0_p
        if (t_store) call zero_one_body_int_store(one_e_h_integrals)

        ! Now, there is no guarantee that FCIDUMP files will include all
        ! permutation symmetry and so we must avoid double-counting when
        ! accumulating Ecore and <a|h'|b>.

        ! * <i|h|i> can only occur once (and requires a factor of 2 in RHF
        !   calculations), hence there is no risk of double counting.
        ! * <ij|ij> and <ij|ji>: store whether seen in a triangular array.
        ! * <ia|ib> and <ia|bi>: store whether seen in a triangular array for
        !   each i.

        ! In the 'seen' arrays, +1 indicates the Coulomb integral has already
        ! been seen and +2 indicates the exchange integral has already been
        ! seen.  We only vaguely attempt to save memory in these arrays (ie one
        ! can do better if needed---see integral stores in molecular_integrals
        ! for inspiration).

        ! We similarly need to remember if we've seen <i|h|a> or <a|h|i> as we
        ! only store one of the pair.

        allocate(seen_iha((nbasis*(nbasis+1))/2), stat=ierr)
        call check_allocate('seen_iha', size(seen_ijij), ierr)
        allocate(seen_ijij((active_basis_offset*(active_basis_offset+1))/2), stat=ierr)
        call check_allocate('seen_ijij', size(seen_ijij), ierr)
        allocate(seen_iaib(active_basis_offset,(nbasis*(nbasis+1))/2), stat=ierr)
        call check_allocate('seen_iaib', size(seen_iaib), ierr)
        seen_iha = .false.
        seen_ijij = 0
        seen_iaib = 0

        ! Freezing virtual orbitals amounts to simply not exciting into them,
        ! which can easily be enforced by removing them from the basis set.

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

            if (i > 0 .and. j == 0 .and. a == 0 .and. b == 0) then
                ! \epsilon_i --- temporarily store for all basis functions,
                ! including inactive (frozen) orbitals.
                if (uhf) then
                    all_basis_fns(i)%sp_eigv = x
                else
                    all_basis_fns(i-1)%sp_eigv = x
                    all_basis_fns(i)%sp_eigv = x
                end if
            end if

            ! Adjust indices to take into account frozen core orbitals.
            ii = i - active_basis_offset
            jj = j - active_basis_offset
            aa = a - active_basis_offset
            bb = b - active_basis_offset

            if (max(ii,jj,aa,bb) <= nbasis) then

                ! Have integrals involving only core or active orbitals.
                if (i == 0 .and. j == 0 .and. a == 0 .and. b == 0) then

                    ! Nuclear energy.
                    Ecore = Ecore + x

                else if (ii > 0 .and. j == 0 .and. a == 0 .and. b == 0) then

                    ! \epsilon_i
                    if (uhf) then
                        basis_fns(ii)%sp_eigv = x
                    else
                        basis_fns(ii-1)%sp_eigv = x
                        basis_fns(ii)%sp_eigv = x
                    end if

                else if (j == 0 .and. b == 0) then

                    ! < i | h | a >
                    if (t_store) then
                        if (ii < 1 .and. ii == aa) then
                            ! Have <i|h|i> from a core orbital.  Add
                            ! contribution to Ecore.
                            Ecore = Ecore + x*rhf_fac
                        else if (all( (/ ii, aa /) > 0)) then
                            if (.not.seen_iha(tri_ind_reorder(ii,aa))) then
                                x = x + get_one_body_int_mol(one_e_h_integrals, ii, aa)
                                call store_one_body_int_mol(ii, aa, x, one_e_h_integrals)
                                seen_iha(tri_ind_reorder(ii,aa)) = .true.
                            end if
                        end if
                    end if

                else

                    ! < i j | 1/r_12 | a b >
                    if (t_store) then
                        if (all( (/ ii,jj,aa,bb /) <= 0)) then
                            ! Have <ij|ab> involving only core orbitals.
                            if (ii == aa .and. jj == bb .and. ii == jj) then
                                ! RHF calculations: need to include <i,up i,down|i,up i,down>
                                if (.not.uhf .and. mod(seen_ijij(tri_ind_reorder(i,j)),2) == 0) then
                                    Ecore = Ecore + x
                                    seen_ijij(tri_ind_reorder(i,j)) = seen_ijij(tri_ind_reorder(i,j)) + 1
                                end if
                            else if (ii == aa .and. jj == bb .and. ii < jj) then
                                if (mod(seen_ijij(tri_ind_reorder(i,j)),2) == 0) then
                                    Ecore = Ecore + x*rhf_fac
                                    seen_ijij(tri_ind_reorder(i,j)) = seen_ijij(tri_ind_reorder(i,j)) + 1
                                end if
                            else if (ii == bb .and. jj == aa .and. ii < jj) then
                                if (seen_ijij(tri_ind_reorder(i,j)) < 2) then
                                    Ecore = Ecore - x
                                    seen_ijij(tri_ind_reorder(i,j)) = seen_ijij(tri_ind_reorder(i,j)) + 2
                                end if
                            end if
                        else if (min(ii,jj,aa,bb) <= 0) then
                            ! Have an integral involving some core and some
                            ! active orbitals.
                            if (ii == aa .and. min(jj,bb) >= 1) then
                                ! Update <j|h|b> with contribution <ij|ib>.
                                if (mod(seen_iaib(i, tri_ind_reorder(jj,bb)),2) == 0) then
                                    x = get_one_body_int_mol(one_e_h_integrals, jj, bb) + x*rhf_fac
                                    call store_one_body_int_mol(jj, bb, x, one_e_h_integrals)
                                    seen_iaib(i, tri_ind_reorder(jj,bb)) = seen_iaib(i, tri_ind_reorder(jj,bb)) + 1
                                end if
                            else if (jj == bb .and. min(ii,aa) >= 1) then
                                ! Update <i|h|a> with contribution <ij|aj>
                                if (mod(seen_iaib(j, tri_ind_reorder(ii,aa)),2) == 0) then
                                    x = get_one_body_int_mol(one_e_h_integrals, ii, aa) + x*rhf_fac
                                    call store_one_body_int_mol(ii, aa, x, one_e_h_integrals)
                                    seen_iaib(j, tri_ind_reorder(ii,aa)) = seen_iaib(j, tri_ind_reorder(ii,aa)) + 1
                                end if
                            else if (ii == bb .and. min(jj,aa) >= 1) then
                                ! Update <j|h|a> with contribution <ij|ai>.
                                if (seen_iaib(i, tri_ind_reorder(jj,aa)) < 2) then
                                    x = get_one_body_int_mol(one_e_h_integrals, jj, aa) - x
                                    call store_one_body_int_mol(jj, aa, x, one_e_h_integrals)
                                    seen_iaib(i, tri_ind_reorder(jj,aa)) = seen_iaib(j, tri_ind_reorder(jj,aa)) + 2
                                end if
                            else if (jj == aa .and. min(ii,bb) >= 1) then
                                ! Update <i|h|b> with contribution <ij|jb>
                                if (seen_iaib(j, tri_ind_reorder(ii,bb)) < 2) then
                                    x = get_one_body_int_mol(one_e_h_integrals, ii, bb) - x
                                    call store_one_body_int_mol(ii, bb, x, one_e_h_integrals)
                                    seen_iaib(j, tri_ind_reorder(ii,bb)) = seen_iaib(j, tri_ind_reorder(ii,bb)) + 2
                                end if
                            end if
                        else if (all( (/ ii,jj,aa,bb /) > 0)) then
                            ! Have <ij|ab> involving active orbitals.
                            call store_two_body_int_mol(ii, jj, aa, bb, x, coulomb_integrals)
                        end if
                    end if

                end if

            end if

        end do

        deallocate(seen_iha, stat=ierr)
        call check_deallocate('seen_iha', ierr)
        deallocate(seen_ijij, stat=ierr)
        call check_deallocate('seen_ijij', ierr)
        deallocate(seen_iaib, stat=ierr)
        call check_deallocate('seen_iaib', ierr)

        close(ir, status='keep')

        if (size(basis_fns) /= size(all_basis_fns)) then
            ! We froze some orbitals...
            ! Print out entire original basis.
            call write_basis_fn_header()
            do i = 1, size(all_basis_fns)
                call write_basis_fn(all_basis_fns(i), ind=i, new_line=.true.)
            end do
            write (6,'(/,1X,"Freezing...",/,1X,"Using complete active space: (",' &
                      //int_fmt(cas(1),0)//',",",'//int_fmt(cas(2),0)//',")",/)') cas
        end if

        do i = 1, size(all_basis_fns)
            deallocate(all_basis_fns(i)%l, stat=ierr)
            call check_deallocate('all_basis_fns_element',ierr)
        end do
        deallocate(all_basis_fns, stat=ierr)
        call check_deallocate('all_basis_fns',ierr)

        call write_basis_fn_header()
        do i = 1, nbasis
            call write_basis_fn(basis_fns(i), ind=i, new_line=.true.)
        end do
        write (6,'(/,1X,a8,f18.12)') 'E_core =', Ecore

        if (.not.t_store) then
            ! Should tidy up and deallocate everything we allocated.
            call end_basis_fns()
        end if

    end subroutine read_in_fcidump

    subroutine init_basis_fns_read_in(norb, uhf, orbsym, basis_arr)

        ! Initialise list of basis functions based upon the FCI namelist`data
        ! read in from an FCIDUMP file.

        ! In:
        !    norb: number of orbitals/functions in the FCIDUMP file.
        !    uhf: true if FCIDUMP file was produced from an unrestricted (UHF)
        !         calculation.  norb refers to spatial orbitals in a RHF-based
        !         FCIDUMP file and spin orbitals in a UHF-based FCIDUMP file.
        !    orbsym: list of symmetries of each orbital/function.
        ! In/Out:
        !    basis_arr: array of basis functions.  On entry it is allocated to
        !    (at least) the size of the basis.  On exit the individual basis_fn
        !    elements have been initialised and symmetry and spatial_index
        !    information assigned.

        use basis, only: basis_fn, init_basis_fn

        integer, intent(in) :: norb
        logical, intent(in) :: uhf
        integer, intent(in) :: orbsym(:)
        type(basis_fn), intent(inout) :: basis_arr(:)

        integer :: i

        do i = 1, norb
            if (uhf) then
                if (mod(i,2) == 0) then
                    call init_basis_fn(basis_arr(i), sym=orbsym(i)-1, ms=-1)
                else
                    call init_basis_fn(basis_arr(i), sym=orbsym(i)-1, ms=1)
                end if
                ! Assume orbitals are ordered appropriately in FCIDUMP...
                basis_arr(i)%spatial_index = (i+1)/2
            else
                ! Need to initialise both up- and down-spin basis functions.
                call init_basis_fn(basis_arr(2*i-1), sym=orbsym(i)-1, ms=-1)
                call init_basis_fn(basis_arr(2*i), sym=orbsym(i)-1, ms=1)
                basis_arr(2*i-1)%spatial_index = i
                basis_arr(2*i)%spatial_index = i
            end if
        end do

    end subroutine init_basis_fns_read_in

end module read_in_system
