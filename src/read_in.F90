module read_in_system

! Module for reading in and manipulating integral files which entirely define
! the Hamiltonian corresponding to a 'real' system.

use const

implicit none

contains

    subroutine read_in_integrals(store_info, cas_info)

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
        use molecular_integrals
        use point_group_symmetry, only: init_pg_symmetry
        use system, only: fcidump, uhf, ecore, nel, nvirt, dipole_int_file, dipole_frozen_core

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use parallel
        use ranking, only: insertion_rank_rp
        use utils, only: get_free_unit, tri_ind_reorder, int_fmt

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
        integer :: i, j, a, b, ii, jj, aa, bb, orbs(4), active(2), core(2), ia, ic, iorb
        real(dp) :: x

        ! reading in...
        integer :: ir, ios, ierr
        logical :: t_exists
        integer :: active_basis_offset, rhf_fac
        integer, allocatable :: seen_ijij(:), seen_iaib(:,:), sp_eigv_rank(:), sp_fcidump_rank(:)
        logical, allocatable :: seen_iha(:)
        real(p), allocatable :: sp_eigv(:)
        logical :: not_found_sp_eigv

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

        ! Only do i/o on root processor.
        if (parent) then
            ir = get_free_unit()
            inquire(file=fcidump, exist=t_exists)
            if (.not.t_exists) call stop_all('read_in_integrals', 'FCIDUMP does not &
                                                               &exist:'//trim(fcidump))
            open (ir, file=fcidump, status='old', form='formatted')

            ! read system data
            read (ir, FCI)
        end if

#ifdef PARALLEL
        ! Distribute FCI namelist.
        call MPI_BCast(norb, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_BCast(nelec, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_BCast(ms2, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_BCast(orbsym, 1000, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_BCast(isym, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_BCast(uhf, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
        call MPI_BCast(syml, 1000, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        call MPI_BCast(symlz, 1000, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
#endif

        ! NOTE: nbasis is currently the number of spin-orbitals in the FCIDUMP file.
        ! This is changed to the number of spin-orbitals active the calculation
        ! later on.

        if (uhf) then
            nbasis = norb
            rhf_fac = 1
        else
            nbasis = 2*norb
            rhf_fac = 2  ! need to double count some integrals.
        end if

        if (present(cas_info)) then
            if (any(cas_info < 1)) then
                cas = (/ nel, nbasis/2 /)
            else
                cas = cas_info
            end if
        else
            cas = (/ nel, nbasis/2 /)
        end if

        ! Read in FCIDUMP file to get single-particle eigenvalues.
        not_found_sp_eigv = .true.
        allocate(sp_eigv(norb), stat=ierr)
        call check_allocate('sp_eigv', norb, ierr)
        ios = 0
        if (parent) then
            do
                ! loop over lines.
                read (ir,*, iostat=ios) x, i, a, j, b
                if (ios == iostat_end) exit ! reached end of file
                if (ios /= 0) call stop_all('read_in_integrals','Problem reading integrals file: '//trim(FCIDUMP))
                if (i > 0 .and. j == 0 .and. a == 0 .and. b == 0) then
                    ! \epsilon_i --- temporarily store for all basis functions,
                    ! including inactive (frozen) orbitals.
                    not_found_sp_eigv = .false.
                    sp_eigv(i) = x
                end if
            end do

            if (not_found_sp_eigv) &
                call stop_all('read_in_integrals',fcidump//' file does not contain &
                              &single-particle eigenvalues.  Please implement &
                              &calculating them from the integrals.')
        end if

#ifdef PARALLEL
        call MPI_BCast(sp_eigv, norb, mpi_preal, root, MPI_COMM_WORLD, ierr)
#endif

        ! Rank basis functions by single-particle energy.
        ! Note that we use a *stable* ranking algorithm.
        allocate(sp_eigv_rank(0:norb), stat=ierr)
        call check_allocate('sp_eigv_rank', norb+1, ierr)
        allocate(sp_fcidump_rank(0:norb), stat=ierr)
        call check_allocate('sp_fcidump_rank', norb+1, ierr)

        if (uhf) then
            ! Cannot simply naively sort all basis functions by energy as that
            ! causes basis functions to be labelled with the incorrect spin.
            ! Further we assume that basis functions have alternating alpha, beta
            ! spins.  To overcome this we sort by energy within each spin and then
            ! merge the two ranking orders.
            ! Yes, this the calls to insertion_rank_rp do result in array copies,
            ! but the arrays are small and the elegance and brevity of the code is
            ! more than aadequate compensation.
            call insertion_rank_rp(sp_eigv(1::2), sp_eigv_rank(1::2), tolerance=depsilon)
            call insertion_rank_rp(sp_eigv(2::2), sp_eigv_rank(2::2), tolerance=depsilon)
            ! Interweave so the basis functions remain with alternating spins.
            forall (i = 1:nbasis-1:2)
                sp_eigv_rank(i) = 2*sp_eigv_rank(i) - 1
                sp_eigv_rank(i+1) = 2*sp_eigv_rank(i+1)
            end forall
        else
            ! RHF case is much easier, as spin channels are degenerate and have
            ! the same spatial components.
            call insertion_rank_rp(sp_eigv(1:), sp_eigv_rank(1:), tolerance=depsilon)
        end if

        ! the 0 index is used as a null entry in the FCIDUMP file, so need
        ! to handle that...
        sp_eigv_rank(0) = 0

        ! sp_eigv_rank(i) = i' takes us from the {i'} basis (as ordered in
        ! the FCIDUMP file) to the {i} basis (ordered by energy).  We also
        ! need the inverse.
        do i = 0, norb
            do j = 0, norb
                if (sp_eigv_rank(j) == i) then
                    sp_fcidump_rank(i) = j
                    exit
                end if
            end do
        end do

        ! Set up basis functions, including those which are subsequently frozen.
        allocate(all_basis_fns(nbasis), stat=ierr)
        call check_allocate('all_basis_fns', nbasis, ierr)
        call init_basis_fns_read_in(norb, uhf, orbsym, sp_eigv, sp_eigv_rank(1:), all_basis_fns)

        ! From CAS work out the start of the active basis functions, the number
        ! of active basis functions and the number of active electrons.

        ! NOTE: this sets nbasis to be the number of spin orbitals in the active
        ! basis.

        active_basis_offset = nel-cas(1) ! number of core *spin* orbitals
        ! Note that nbasis is spin-orbitals whereas cas(2)=M is in spatial orbitals
        ! (as we use the conventional CAS definition).
        nbasis = min(nbasis, 2*cas(2))
        nel = nel - ( nel - cas(1) )
        nvirt = nbasis - nel
        if (uhf) then
            norb =  nbasis
        else
            norb = nbasis/2
        end if

        ! Set up basis functions used in calculation.
        allocate(basis_fns(nbasis), stat=ierr)
        call check_allocate('basis_fns', nbasis, ierr)
        call init_basis_fns_read_in(norb, uhf, orbsym, sp_eigv, sp_eigv_rank(1+active_basis_offset/rhf_fac:), basis_fns)

        deallocate(sp_eigv, stat=ierr)
        call check_deallocate('sp_eigv', ierr)

        ! Was a symmetry found for all basis functions?  If not, then we must
        ! turn symmetry off.
        if (minval(basis_fns(:)%sym) < 0) then
            if (parent) write (6,'(1X,a62)') 'Unconverged symmetry found.  Turning point group symmetry off.'
            forall (i=1:nbasis) basis_fns(i)%sym = 0
        end if

        ! Set up symmetry information.
        if (t_store) call init_pg_symmetry()

        ! Initialise integral stores.
        if (t_store) call init_molecular_integrals()

        if (parent) then
            ! Now, read in FCIDUMP again to get the integrals.
            rewind(ir)
            ! Have to re-read the FCI namelist.
            ! ***WARNING***
            ! This will, e.g. overwrite norb, but we don't need those variables any
            ! more.
            read(ir,FCI)
        end if

        ! Freezing core orbitals amounts to changing the Hamiltonian.  In
        ! particular:

        ! E_core = E_nuc + \sum_{i=1}^{N_fr} <i|h|i> + \sum_{i<j}^{N_fr} <ij|ij> - <ij|ji>  [1]
        ! <a|h'|b> = <a|h|b> + \sum_{i=1}^{N_fr} <ia|ib> - <ia|bi>                          [2]

        ! where we use i,j to refer to frozen core spin-orbitals and a,b to refer to
        ! active spin-orbitals.

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

        if (parent) then

            allocate(seen_iha((nbasis*(nbasis+1))/2), stat=ierr)
            call check_allocate('seen_iha', size(seen_ijij), ierr)
            allocate(seen_ijij((active_basis_offset*(active_basis_offset+1))/2), stat=ierr)
            call check_allocate('seen_ijij', size(seen_ijij), ierr)
            allocate(seen_iaib(-active_basis_offset+1:0,(nbasis*(nbasis+1))/2), stat=ierr)
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
                if (ios /= 0) call stop_all('read_in_integrals','Problem reading integrals file: '//trim(FCIDUMP))

                ! Working in spin orbitals but FCIDUMP is in spatial orbitals in RHF
                ! calculations and spin orbitals in UHF calculations, and te basis
                ! is ! not necessarily ordered by energy.  We wish to work in spin
                ! orbitals ordered by energy.
                ! Need to only store integrals in one spin-channel in RHF, so
                ! can just get away with referring (e.g.) beta orbitals, hence the
                ! factor of 2 in RHF.  The integral store routines convert the
                ! indices further to compress the integral stores as much as
                ! possible.
                i = rhf_fac*sp_fcidump_rank(i)
                j = rhf_fac*sp_fcidump_rank(j)
                a = rhf_fac*sp_fcidump_rank(a)
                b = rhf_fac*sp_fcidump_rank(b)

                ! Adjust indices to take into account frozen core orbitals and to
                ! convert to an energy ordering.
                ii = i - active_basis_offset
                jj = j - active_basis_offset
                aa = a - active_basis_offset
                bb = b - active_basis_offset

                if (max(ii,jj,aa,bb) <= nbasis) then

                    ! Have integrals involving only core or active orbitals.
                    if (i == 0 .and. j == 0 .and. a == 0 .and. b == 0) then

                        ! Nuclear energy.
                        Ecore = Ecore + x

                    else if (i > 0 .and. j == 0 .and. a == 0 .and. b == 0) then

                        ! \epsilon_i
                        ! Already dealt with in previous pass over FCIDUMP.

                    else if (j == 0 .and. b == 0) then

                        ! < i | h | a >
                        if (t_store) then
                            if (ii < 1 .and. ii == aa) then
                                ! Have <i|h|i> from a core orbital.  Add
                                ! contribution to Ecore.
                                ! If RHF need to include <i,up|h|i,up> and
                                ! <i,down|h|i,down>.
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
                            orbs = (/ ii, jj, aa, bb /)
                            select case(count(orbs > 0))
                            case(0)
                                ! Have <ij|ab> involving only core orbitals.
                                ! We should only run over i<j but there's no
                                ! guarantee that the FCIDUMP file contains the
                                ! permutations we expect (ie it might contain
                                ! <ji|ji> but not <ij|ij>, where i<j), so
                                ! instead we let the seen_ijij array make sure
                                ! we only use a unique integral once, no matter
                                ! which permutation(s) occur in the FCIDUMP
                                ! file.
                                if (ii == aa .and. jj == bb .and. ii == jj) then
                                    if (.not.uhf .and. mod(seen_ijij(tri_ind_reorder(i,j)),2) == 0) then
                                        ! RHF calculations: need to include <i,up i,down|i,up i,down>.
                                        Ecore = Ecore + x
                                        seen_ijij(tri_ind_reorder(i,j)) = seen_ijij(tri_ind_reorder(i,j)) + 1
                                    end if
                                else if (ii == aa .and. jj == bb .and. ii /= jj) then
                                    ! <ij|ij>, i/=j
                                    if (mod(seen_ijij(tri_ind_reorder(i,j)),2) == 0) then
                                        ! If RHF, then need to include:
                                        !   <i,up j,up|i,up, j,up>
                                        !   <i,up j,down|i,up, j,down>
                                        !   <i,down j,up|i,down, j,up>
                                        !   <i,down j,down|i,down, j,down>
                                        Ecore = Ecore + x*rhf_fac**2
                                        seen_ijij(tri_ind_reorder(i,j)) = seen_ijij(tri_ind_reorder(i,j)) + 1
                                    end if
                                else if (ii == bb .and. jj == aa .and. ii /= jj) then
                                    ! <ij|ji>, i/=j
                                    if (seen_ijij(tri_ind_reorder(i,j)) < 2) then
                                        ! If RHF, then need to include:
                                        !   <i,up j,up|j,up, i,up>
                                        !   <i,down j,down|j,down, i,down>
                                        Ecore = Ecore - rhf_fac*x
                                        seen_ijij(tri_ind_reorder(i,j)) = seen_ijij(tri_ind_reorder(i,j)) + 2
                                    end if
                                end if
                            case(2)
                                ! Have an integral involving two core and two active orbitals.
                                ic = 1
                                ia = 1
                                do iorb = 1, 4
                                    if (orbs(iorb) > 0) then
                                        active(ia) = orbs(iorb)
                                        ia = ia +1
                                    else
                                        core(ic) = orbs(iorb)
                                        ic = ic + 1
                                    end if
                                end do
                                if (core(1) == core(2)) then
                                    ! Have integral of type < i a | i b > or < i a | b i >,
                                    ! where i is a core orbital and a and b are
                                    ! active orbitals.
                                    if ((ii == core(1) .and. aa == core(1)) .or. (jj == core(1) .and. bb == core(1))) then
                                        ! < i a | i b >
                                        if (mod(seen_iaib(core(1), tri_ind_reorder(active(1),active(2))),2) == 0) then
                                            ! Update <a|h|b> with contribution <ia|ib>.
                                            x = get_one_body_int_mol(one_e_h_integrals, active(1), active(2)) + x*rhf_fac
                                            call store_one_body_int_mol(active(1), active(2), x, one_e_h_integrals)
                                            seen_iaib(core(1), tri_ind_reorder(active(1),active(2))) = &
                                                seen_iaib(core(1), tri_ind_reorder(active(1),active(2))) + 1
                                        end if
                                    else
                                        ! < i a | b i > (or a permutation thereof)
                                        if (seen_iaib(core(1), tri_ind_reorder(active(1),active(2))) < 2) then
                                            ! Update <j|h|a> with contribution <ij|ai>.
                                            x = get_one_body_int_mol(one_e_h_integrals, active(1), active(2)) - x
                                            call store_one_body_int_mol(active(1), active(2), x, one_e_h_integrals)
                                            seen_iaib(core(1), tri_ind_reorder(active(1),active(2))) = &
                                                seen_iaib(core(1), tri_ind_reorder(active(1),active(2))) + 2
                                        end if
                                    end if
                                end if
                            case(4)
                                ! Have <ij|ab> involving active orbitals.
                                call store_two_body_int_mol(ii, jj, aa, bb, x, coulomb_integrals)
                            end select
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

        end if

#ifdef PARALLEL
        call MPI_BCast(ECore, 1, mpi_preal, root, MPI_COMM_WORLD, ierr)
#endif
        call broadcast_one_body_int(one_e_h_integrals, root)
        call broadcast_two_body_int(coulomb_integrals, root)

        if (size(basis_fns) /= size(all_basis_fns) .and. parent) then
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

        if (parent) then
            call write_basis_fn_header()
            do i = 1, nbasis
                call write_basis_fn(basis_fns(i), ind=i, new_line=.true.)
            end do
            write (6,'(/,1X,a8,f18.12)') 'E_core =', Ecore
        end if

        if (t_store .and. dipole_int_file /= '') then
            call read_in_one_body(dipole_int_file, sp_fcidump_rank(1:), active_basis_offset, &
                                  one_body_op_integrals, dipole_frozen_core)
        end if

        deallocate(sp_eigv_rank, stat=ierr)
        call check_deallocate('sp_eigv_rank', ierr)
        deallocate(sp_fcidump_rank, stat=ierr)
        call check_deallocate('sp_fcidump_rank', ierr)

        if (.not.t_store) then
            ! Should tidy up and deallocate everything we allocated.
            call end_basis_fns()
        end if

    end subroutine read_in_integrals

    subroutine init_basis_fns_read_in(norb, uhf, orbsym, sp_eigv, sp_eigv_rank, basis_arr)

        ! Initialise list of basis functions based upon the FCI namelist`data
        ! read in from an FCIDUMP file.

        ! In:
        !    norb: number of orbitals/functions to initialise.
        !    uhf: true if FCIDUMP file was produced from an unrestricted (UHF)
        !         calculation.  norb refers to spatial orbitals in a RHF-based
        !         FCIDUMP file and spin orbitals in a UHF-based FCIDUMP file.
        !    orbsym: list of symmetries of each orbital/function.
        !    sp_eigv: single-particle eigenvalues of each orbital/function.
        !    sp_eigv_rank: rank of each orbital/function by the single-particle
        !         eigenvalue.
        !
        ! Note: input arrays *must* have at least norb set elements.  Only the
        ! first norb elements are used.  In RHF systems the up and down spin
        ! orbitals have identical eigenvalues and symmetries (due to sharing the
        ! same spatial orbitals).
        !
        ! In/Out:
        !    basis_arr: array of (spin) basis functions.  On entry it is
        !    allocated to (at least) the size of the basis.  On exit the
        !    individual basis_fn elements have been initialised and symmetry and
        !    spatial_index information assigned.

        use basis, only: basis_fn, init_basis_fn

        integer, intent(in) :: norb
        logical, intent(in) :: uhf
        integer, intent(in) :: orbsym(:)
        real(p), intent(in) :: sp_eigv(:)
        integer, intent(in) :: sp_eigv_rank(:)
        type(basis_fn), intent(inout) :: basis_arr(:)

        integer :: i, rank

        do i = 1, norb
            rank = sp_eigv_rank(i)
            if (uhf) then
                if (mod(i,2) == 0) then
                    call init_basis_fn(basis_arr(i), sym=orbsym(rank)-1, ms=-1)
                else
                    call init_basis_fn(basis_arr(i), sym=orbsym(rank)-1, ms=1)
                end if
                ! Assume orbitals are ordered appropriately in FCIDUMP...
                basis_arr(i)%spatial_index = (i+1)/2
                basis_arr(i)%sp_eigv = sp_eigv(rank)
            else
                ! Need to initialise both up- and down-spin basis functions.
                call init_basis_fn(basis_arr(2*i), sym=orbsym(rank)-1, ms=-1)
                call init_basis_fn(basis_arr(2*i-1), sym=orbsym(rank)-1, ms=1)
                basis_arr(2*i-1)%spatial_index = i
                basis_arr(2*i)%spatial_index = i
                basis_arr(2*i-1)%sp_eigv = sp_eigv(rank)
                basis_arr(2*i)%sp_eigv = sp_eigv(rank)
            end if
        end do

    end subroutine init_basis_fns_read_in

    subroutine read_in_one_body(integral_file, sp_fcidump_rank, active_basis_offset, store, frozen_core_expectation)

        ! Read in an integral file containing the integrals <i|O|a>, where O is
        ! a one-body operator which is not part of the Hamiltonian.

        ! The integral file consists soley of lines with the format:
        !    x i a
        ! where <i|O|a> = x.

        ! We assume the orbitals have the same indexing as used in the FCIDUMP
        ! file and i,a are in spatial (spin) indices if produced by a RHF (UHF)
        ! calculation.

        ! In:
        !    integral_file: file containing integrals.
        !    sp_fcidump_rank: ranking array which converts index in the
        !         integrals file(s) to the (energy-ordered) index used in HANDE,
        !         i.e. sp_fcidump_rank(a) = i, where a is the a-th
        !         orbital according to the ordering used in the integral files
        !         and i is the i-th orbital by energy ordering.
        !    active_basis_offset: number of frozen core orbitals.
        ! Out:
        !    store: one-body store of integrals given in integral_file..
        !    frozen_core_expectation: contribution to <\Psi|O|\Psi> from the
        !         frozen core orbitals.

        use basis, only: nbasis
        use system, only: uhf
        use point_group_symmetry, only: cross_product_pg_basis
        use molecular_integrals, only: one_body, init_one_body_int_store,              &
                                       end_one_body_int_store, store_one_body_int_mol, &
                                       broadcast_one_body_int

        use errors, only: stop_all
        use parallel, only: parent, root
        use utils, only: get_free_unit

        use, intrinsic :: iso_fortran_env, only: iostat_end

        character(*), intent(in) :: integral_file
        integer, intent(in) :: sp_fcidump_rank(:), active_basis_offset
        type(one_body), intent(out) :: store
        real(p), intent(out) :: frozen_core_expectation

        integer :: ir, op_sym, ios, i, a, ii, aa, rhf_fac
        logical :: t_exists

        real(p) :: x

        if (uhf) then
            rhf_fac = 1
        else
            rhf_fac = 2  ! need to double count some integrals.
        end if

        ! Only do i/o on root processor.
        if (parent) then
            ! We don't know the symmetry of the operator.
            ! However, we do know that a non-zero integral must have a totally
            ! symmetric integrand *and* we know the symmetries of all the orbitals.
            ir = get_free_unit()
            inquire(file=integral_file, exist=t_exists)
            if (.not.t_exists) call stop_all('read_in_one_body', 'Integral file does not &
                                                               &exist:'//trim(integral_file))
            open (ir, file=integral_file, status='old', form='formatted')

            do
                read (ir,*, iostat=ios) x, i, a
                if (ios == iostat_end) exit ! reached end of file
                if (ios /= 0) call stop_all('read_in_one_body','Problem reading integrals file: '//trim(integral_file))
                ii = rhf_fac*sp_fcidump_rank(i) - active_basis_offset
                aa = rhf_fac*sp_fcidump_rank(a) - active_basis_offset
                if (abs(x) > depsilon .and. ii > 0 .and. aa > 0) exit ! Found the first non-zero integral in active space.
            end do
            ! We only use Abelian symmetries so all representations are their own
            ! inverse.
            op_sym = cross_product_pg_basis(ii,aa)
        else
            ! We'll broadcast the symmetry and the integrals to all other
            ! processors later.
            op_sym = -1
        end if

        ! Allocate integral store on *all* processors.
        if (allocated(store%integrals)) call end_one_body_int_store(store)
        call init_one_body_int_store(op_sym, store)

        ! In addition to reading in the integrals, we must also calculate the
        ! contribution from the core (frozen) orbitals.
        !
        ! Consider |\Psi> = c_i |D_i>, where {|D_i>} includes the frozen core
        ! electrons and Einstein summation is used throughout. Then:
        !
        !     <\Psi | O | \Psi> = <D_i|O|D_j> c_i^* c_j
        !
        ! where, for spin-orbitals {|a>}:
        !
        !     <D_i|O|D_j> = | <a|O|a> if |D_i>=|D_j>
        !                   | <a|O|b> if |D_i> and |D_j> are related by the excitation a->b
        !                   | 0       otherwise
        !
        !     <\Psi |O | \Psi > = <a_c|O|a_c> c_i^* c_i + <D_i'|O|D_j'> c_i^* c_j
        !
        ! where {|D_i'>} now only incude active occupied orbitals and {|a_c>} are
        ! the frozen core electrons.  As {|a_c>} is indentical for all
        ! determinants and |\Psi> is normalised, the summation over i in the first
        ! term is unity and hence:
        !
        !     <\Psi |O | \Psi > = <a_c|O|a_c> + <D_i'|O|D_j'> c_i^* c_j
        !
        ! The contribution from the frozen core electrons is hence just the sum
        ! over the diagonal integrals.

        ! And back to root...
        frozen_core_expectation = 0.0_p
        if (parent) then
            rewind(ir)
            do
                read (ir,*, iostat=ios) x, i, a
                if (ios == iostat_end) exit ! reached end of file
                if (ios /= 0) call stop_all('read_in_one_body','Problem reading integrals file: '//trim(integral_file))
                ii = rhf_fac*sp_fcidump_rank(i) - active_basis_offset
                aa = rhf_fac*sp_fcidump_rank(a) - active_basis_offset
                if (ii < 1 .and. ii == aa) then
                    frozen_core_expectation = frozen_core_expectation + rhf_fac*x
                else if (min(ii,aa) >= 1 .and. max(ii,aa) <= nbasis) then
                    call store_one_body_int_mol(ii, aa, x, store)
                end if
            end do
        end if

        ! And now send info everywhere...
#ifdef PARALLEL
        call MPI_BCast(frozen_core_expectation, 1, mpi_preal, root, MPI_COMM_WORLD, ierr)
#endif
        call broadcast_one_body_int(store, root)

    end subroutine read_in_one_body

end module read_in_system
