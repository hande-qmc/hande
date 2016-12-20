module qmc_io

! QMC Output procedures.

use const
use spawn_data, only: spawn_t

implicit none

private
public :: write_qmc_report_header, write_qmc_report, write_dmqmc_report_header, write_dmqmc_report
public :: write_qmc_var, write_column_title

interface write_qmc_var
    module procedure :: write_qmc_var_int_32
    module procedure :: write_qmc_var_int_64
    module procedure :: write_qmc_var_real_sp
    module procedure :: write_qmc_var_real_dp
end interface

contains

    subroutine write_qmc_report_header(ntypes, cmplx_est, rdm_energy, nattempts, io_unit)

        ! In:
        !    ntypes: number of particle types being sampled.
        ! In (optional):
        !    cmplx_est: Print out information of cmplx_estlex estimators.
        !    rdm_energy: Print out energy calculated from rdm.
        !    nattempts: Print out number of attempts made each iteration
        !       (for CCMC).
        !    io_unit: io unit to write to.

        use calc, only: doing_calc, hfs_fciqmc_calc

        integer, intent(in) :: ntypes
        logical, optional, intent(in) :: cmplx_est, rdm_energy, nattempts
        integer, intent(in), optional :: io_unit

        logical :: cmplx_est_set
        integer :: i, nreplicas, iunit
        character(20) :: column_title

        iunit = 6
        if (present(io_unit)) iunit = io_unit

        cmplx_est_set = .false.
        if (present(cmplx_est)) cmplx_est_set = cmplx_est

        ! Data table info.
        write (iunit,'(1X,"Information printed out every QMC report loop:",/)')
        write (iunit,'(1X,"Shift: the energy offset calculated at the end of the report loop.")')
        write (iunit,'(1X,"H_0j: <D_0|H|D_j>, Hamiltonian matrix element.")')
        write (iunit,'(1X,"N_j: population of Hamiltonian particles on determinant D_j.")')
        if (doing_calc(hfs_fciqmc_calc)) then
            write (iunit,'(1X,"O_0j: <D_0|O|D_j>, operator matrix element.")')
            write (iunit,'(1X,a67)') "N'_j: population of Hellmann--Feynman particles on determinant D_j."
            write (iunit,'(1X,"# HF psips: current total population of Hellmann--Feynman particles.")')
        end if

        write (iunit,'(1X,"# H psips: current total population of Hamiltonian particles.")')
        write (iunit,'(1X,"# states: number of many-particle states occupied.")')
        write (iunit,'(1X,"spawn_events: number of successful spawning events across all processors.")')
        write (iunit,'(1X,"R_spawn: average rate of spawning across all processors.")')
        write (iunit,'(1X,"time: average time per Monte Carlo cycle.",/)')
        write (iunit,'(1X,"Note that all particle populations are averaged over the report loop.",/)')

        if (cmplx_est_set) then
            nreplicas = ntypes/2
        else
            nreplicas = ntypes
        end if

        if (nreplicas > 1) then
            ! Label replicas
            write (iunit,'(1X,"#",1X)', advance='no')
            call write_column_title(iunit, '', int_val=.true.)
            do i = 1, nreplicas
                write (column_title, '("Replica ",i0)') i
                call write_column_title(iunit, trim(column_title))
                call write_column_title(iunit, '')
                call write_column_title(iunit, '')
                call write_column_title(iunit, '')
                if (cmplx_est_set) then
                    call write_column_title(iunit, '')
                    call write_column_title(iunit, '')
                end if
            end do
            write (iunit,'()')
        end if

        write (iunit,'(1X,"#",1X)', advance='no')
        call write_column_title(iunit, 'iterations', int_val=.true., justify=1)
        ! NOTE: HFS and complex are not currently compatible.
        if (cmplx_est_set) then
            do i = 1, ntypes, 2
                call write_column_title(iunit, 'Shift')
                call write_column_title(iunit, 'Re{\sum H_0j N_j}')
                call write_column_title(iunit, 'Im{\sum H_0j N_j}')
                call write_column_title(iunit, 'Re{N_0}')
                call write_column_title(iunit, 'Im{N_0}')
                call write_column_title(iunit, '# H psips')
            end do
        else
            do i = 1, ntypes
                call write_column_title(iunit, 'Shift')
                call write_column_title(iunit, '\sum H_0j N_j')
                call write_column_title(iunit, 'N_0')
                if (doing_calc(hfs_fciqmc_calc)) then
                    call write_column_title(iunit, 'HF shift')
                    call write_column_title(iunit, '\sum O_0j N_j')
                    call write_column_title(iunit, "\sum H_0j N'_j")
                    call write_column_title(iunit, "N'_0")
                end if
                call write_column_title(iunit, '# H psips')
                if (doing_calc(hfs_fciqmc_calc)) call write_column_title(iunit, '# HF psips')
            end do
        end if
        if (present(rdm_energy)) then
            if (rdm_energy) call write_column_title(iunit, 'RDM Energy')
        end if
        call write_column_title(iunit, '# states', int_val=.true., justify=1)
        call write_column_title(iunit, '# spawn_events', int_val=.true., justify=1)
        if (present(nattempts)) then
            if (nattempts) call write_column_title(iunit, '# attempts', int_val=.true., justify=1)
        end if
        call write_column_title(iunit, 'R_spawn', low_prec_val=.true.)
        call write_column_title(iunit, '  time', low_prec_val=.true., justify=2)

        write (iunit,'()')

    end subroutine write_qmc_report_header

    subroutine write_dmqmc_report_header(ntypes, dmqmc_in, max_excit, estimates)

        ! Write header for DMQMC specific information.

        ! In:
        !    ntypes: number of particle types being sampled.
        !    dmqmc_in: input options for dmqmc calculations.
        !    max_excit: maximum number of excitations in system.

        use calc, only: doing_calc, hfs_fciqmc_calc, dmqmc_calc, doing_dmqmc_calc
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t
        use calc, only: dmqmc_energy, dmqmc_energy_squared, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_correlation, dmqmc_full_r2, dmqmc_rdm_r2, dmqmc_kinetic_energy
        use calc, only: dmqmc_H0_energy, dmqmc_potential_energy, dmqmc_HI_energy
        use utils, only: int_fmt

        integer, intent(in) :: ntypes
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer, intent(in) :: max_excit
        type(dmqmc_estimates_t), intent(in) :: estimates

        integer :: i, j, iunit
        character(16) :: excit_header
        character(10) :: header_iidx, header_jidx

        iunit = 6

        write (iunit,'(1X,"Information printed out every QMC report loop:",/)')
        write (iunit,'(1X,"Shift: the energy offset calculated at the end of the report loop.")')
        if (doing_dmqmc_calc(dmqmc_full_r2)) then
            write (iunit, '(1X, "Trace: The current total population on the diagonal elements of the &
                                 &first replica of the density matrix.")')
            write (iunit, '(1X, "Trace 2: The current total population on the diagonal elements of the &
                                 &second replica of the density matrix.")')
        else
            write (iunit, '(1X, "Trace: The current total population on the diagonal elements of the &
                                 &density matrix.")')
        end if
        if (doing_dmqmc_calc(dmqmc_full_r2)) then
            write (iunit, '(1X, "Full S2: The numerator of the estimator for the Renyi entropy of the &
                                  &full system.")')
        end if
        if (doing_dmqmc_calc(dmqmc_energy)) then
            write (iunit, '(1X, "\sum\rho_{ij}H_{ji}: The numerator of the estimator for the expectation &
                                 &value of the energy.")')
        end if
        if (doing_dmqmc_calc(dmqmc_energy_squared)) then
            write (iunit, '(1X, "\sum\rho_{ij}H2{ji}: The numerator of the estimator for the expectation &
                                 &value of the energy squared.")')
        end if
        if (doing_dmqmc_calc(dmqmc_correlation)) then
            write (iunit, '(1X, "\sum\rho_{ij}S_{ji}: The numerator of the estimator for the expectation &
                                 &value of the spin correlation function.")')
        end if
        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
            write (iunit, '(1X, "\sum\rho_{ij}M2{ji}: The numerator of the estimator for the expectation &
                                 &value of the staggered magnetisation.")')
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            write (iunit, '(1x, "RDM(n) S2: The numerator of the estimator for the Renyi entropy of RDM n.")')
        end if
        if (dmqmc_in%rdm%calc_inst_rdm) then
            write (iunit, '(1x, "RDM(n) trace m: The current total population on the diagonal of replica m &
                                  &of RDM n.")')
        end if
        if (dmqmc_in%calc_excit_dist) then
            write (iunit, '(1x, "Excit. level n: The fraction of particles on excitation level n of the &
                             &density matrix.")')
        end if
        if (dmqmc_in%calc_mom_dist) write (iunit, '(1x, "n_k: The numerator of the estimator for the &
                                                     &momentum distribution at momentum k")')

        write (iunit,'(1X,"# particles: current total population of Hamiltonian particles.")')
        write (iunit,'(1X,"# states: number of many-particle states occupied.")')
        write (iunit,'(1X,"# spawn_events: number of successful spawning events across all processors.")')
        write (iunit,'(1X,"R_spawn: average rate of spawning across all processors.")')
        write (iunit,'(1X,"time: average time per Monte Carlo cycle.",/)')
        write (iunit,'(1X,"Note that all particle populations are averaged over the report loop.",/)')

        ! Header of data table.
        write (iunit,'(1X,"#",1X)', advance='no')
        call write_column_title(iunit, 'iterations', int_val=.true., justify=1)
        call write_column_title(iunit, 'Instant shift')
        call write_column_title(iunit, 'Trace')
        if (doing_dmqmc_calc(dmqmc_full_r2)) then
            call write_column_title(iunit, 'Trace 2')
            call write_column_title(iunit, 'Full S2')
        end if
        if (doing_dmqmc_calc(dmqmc_energy)) then
            call write_column_title(iunit, '\sum\rho_{ij}H_{ji}')
        end if
        if (doing_dmqmc_calc(dmqmc_energy_squared)) then
            call write_column_title(iunit, '\sum\rho_{ij}H2{ji}')
        end if
        if (doing_dmqmc_calc(dmqmc_correlation)) then
            call write_column_title(iunit, '\sum\rho_{ij}S_{ji}')
        end if
        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
            call write_column_title(iunit, '\sum\rho_{ij}M2{ji}')
        end if
        if (doing_dmqmc_calc(dmqmc_kinetic_energy)) then
            call write_column_title(iunit, '\sum\rho_{ij}T_{ji}')
        end if
        if (doing_dmqmc_calc(dmqmc_H0_energy)) then
            call write_column_title(iunit, '\sum\rho_{ij}H0{ji}')
        end if
        if (doing_dmqmc_calc(dmqmc_HI_energy)) then
            call write_column_title(iunit, '\sum\rho_{ij}HI{ji}')
        end if
        if (doing_dmqmc_calc(dmqmc_potential_energy)) then
            call write_column_title(iunit, '\sum\rho_{ij}U_{ji}')
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            do i = 1, dmqmc_in%rdm%nrdms
                write(header_iidx, '('//int_fmt(i,0)//')') i
                call write_column_title(iunit, 'RDM'//trim(header_iidx)//' S2')
            end do
        end if
        if (dmqmc_in%rdm%calc_inst_rdm) then
            do i = 1, dmqmc_in%rdm%nrdms
                do j = 1, ntypes
                    write(header_iidx, '('//int_fmt(i,0)//')') i
                    write(header_jidx, '('//int_fmt(j,0)//')') j
                    call write_column_title(iunit, 'RDM'//trim(header_iidx)//' trace '//trim(header_jidx)//'')
                end do
            end do
        end if
        if (dmqmc_in%calc_excit_dist) then
            do i = 0, max_excit
                write (excit_header, '("Excit. level",1X,'//int_fmt(i,0)//')') i
                call write_column_title(iunit, excit_header)
            end do
        end if
        if (dmqmc_in%calc_mom_dist) then
            call write_momentum_array_header('   n_', .true., estimates%mom_dist%kpoints)
        end if

        call write_column_title(iunit, '# particles')
        call write_column_title(iunit, '# states', int_val=.true., justify=1)
        call write_column_title(iunit, '# spawn_events', int_val=.true., justify=1)
        call write_column_title(iunit, 'R_spawn', low_prec_val=.true.)
        call write_column_title(iunit, '  time', low_prec_val=.true., justify=2)

        write (iunit,'()')

    end subroutine write_dmqmc_report_header

    subroutine write_column_title(io, key, int_val, low_prec_val, justify, sep)

        ! Write column headers using the same column width as write_qmc_var.

        ! In:
        !    io: unit to write to.
        !    int_val: integer value column.
        !    low_prec_val: low-precision real column.
        !    justify: 1: right justify, 2: no justification, otherwise: left justify.

        ! Note: initial space and newline is left to caller.

        integer, intent(in) :: io
        character(*), intent(in) :: key
        logical, intent(in), optional :: int_val, low_prec_val
        integer, intent(in), optional :: justify
        character(1), intent(in), optional :: sep
        logical :: int_val_loc, low_prec_val_loc
        character(20) :: key_str
        integer :: justify_loc, str_len
        character(1) :: sep_loc

        int_val_loc = .false.
        low_prec_val_loc = .false.
        justify_loc = 0
        sep_loc = ' '
        if (present(int_val)) int_val_loc = int_val
        if (present(low_prec_val)) low_prec_val_loc = low_prec_val
        if (present(justify)) justify_loc = justify
        if (present(sep)) sep_loc = sep

        key_str = key
        str_len = len(key_str)
        if (int_val_loc) then
            str_len = 14
        end if
        if (low_prec_val_loc) str_len = 8
        select case(justify_loc)
        case(1)
            ! Right justify.
            write (io,'("'//adjustr(key_str(:str_len))//'", a1,1X)', advance='no') sep_loc
        case(2)
            ! No justification.
            write (io,'("'//key_str(:str_len)//'", a1,1X)', advance='no') sep_loc
        case default
            ! Left justify.
            ! an extra space so the column lines up with the first non-sign entry in the value.
            write (io,'(1X,"'//adjustl(key_str(:str_len-1))//'", a1, 1X)', advance='no') sep_loc
        end select

    end subroutine write_column_title

    subroutine write_qmc_report(qmc_in, qs, ireport, ntot_particles, elapsed_time, comment, non_blocking_comm, io_unit, cmplx_est, &
                                rdm_energy, nattempts)

        ! Write the report line at the end of a report loop.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    qs: QMC state (containing shift and various estimators).
        !    ireport: index of the report loop.
        !    ntot_particles: total number of particles in main walker list.
        !    elapsed_time: time taken for the report loop.
        !    comment: if true, then prefix the line with a #.
        !    non_blocking_comm: true if using non-blocking communications
        ! In (optional):
        !    io_unit: io unit to write report to.
        !    cmplx_est: if true, doing calculation with real and imaginary walkers
        !       so need to print extra parameters.
        !    rdm_energy: Print energy calculated from RDM.
        !    nattempts: Print number of selection attempts made. Only used in CCMC.

        use calc, only: doing_calc, hfs_fciqmc_calc
        use qmc_data, only: qmc_in_t, qmc_state_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        integer, intent(in) :: ireport
        real(dp), intent(in) :: ntot_particles(:)
        real, intent(in) :: elapsed_time
        logical, intent(in) :: comment, non_blocking_comm
        logical, intent(in), optional :: cmplx_est, rdm_energy, nattempts
        integer, intent(in), optional :: io_unit

        logical :: cmplx_est_set, nattempts_set
        integer :: mc_cycles, ntypes, i, iunit

        ntypes = size(ntot_particles)
        iunit = 6
        if (present(io_unit)) iunit = io_unit
        ! For non-blocking communications we print out the nth report loop
        ! after the (n+1)st iteration. Adjust mc_cycles accordingly
        if (.not. non_blocking_comm) then
            mc_cycles = ireport*qmc_in%ncycles
        else
            mc_cycles = (ireport-1)*qmc_in%ncycles
        end if

        cmplx_est_set = .false.
        if (present(cmplx_est)) cmplx_est_set = cmplx_est
        nattempts_set = .false.
        if (present(nattempts)) nattempts_set = nattempts

        if (comment) then
            write (iunit,'(1X,"#",1X)', advance='no')
        else
            write (iunit,'(3X)', advance='no')
        end if

        call write_qmc_var(iunit, qs%mc_cycles_done+mc_cycles)

        ! NOTE: HFS and complex are not currently compatible.
        if (cmplx_est_set) then
            do i = 1, ntypes, 2
                call write_qmc_var(iunit, qs%shift(i))
                call write_qmc_var(iunit, real(qs%estimators(i)%proj_energy_comp, p))
                call write_qmc_var(iunit, aimag(qs%estimators(i)%proj_energy_comp))
                call write_qmc_var(iunit, real(qs%estimators(i)%D0_population_comp, p))
                call write_qmc_var(iunit, aimag(qs%estimators(i)%D0_population_comp))
                call write_qmc_var(iunit, ntot_particles(i)+ntot_particles(i+1))
            end do
        else
            do i = 1, ntypes
                call write_qmc_var(iunit, qs%shift(i))
                call write_qmc_var(iunit, qs%estimators(i)%proj_energy)
                call write_qmc_var(iunit, qs%estimators(i)%D0_population)
                call write_qmc_var(iunit, ntot_particles(i))
            end do
        end if

        if (present(rdm_energy)) then
            if (rdm_energy) call write_qmc_var(iunit, qs%estimators(1)%rdm_energy/qs%estimators(1)%rdm_trace)
        end if

        call write_qmc_var(iunit, qs%estimators(1)%tot_nstates)
        call write_qmc_var(iunit, qs%estimators(1)%tot_nspawn_events)
        if (nattempts_set) then
            call write_qmc_var(iunit, qs%estimators(1)%nattempts)
        end if

        call write_qmc_var(iunit, qs%spawn_store%rspawn, low_prec=.true.)
        call write_qmc_var(iunit, elapsed_time/qmc_in%ncycles, low_prec=.true.)

        write (iunit,'()')

    end subroutine write_qmc_report

    subroutine write_dmqmc_report(qmc_in, qs, ireport, ntot_particles, elapsed_time, comment, &
                                  dmqmc_in, dmqmc_estimates)

        ! Write the report line at the end of a report loop.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    qs: QMC state (containing shift and various estimators).
        !    ireport: index of the report loop.
        !    ntot_particles: total number of particles in main walker list.
        !    elapsed_time: time taken for the report loop.
        !    comment: if true, then prefix the line with a #.
        !    dmqmc_in: input options relating to DMQMC.
        !    dmqmc_estimates: type containing all DMQMC estimates to be printed.

        use calc, only: doing_calc, dmqmc_calc, doing_dmqmc_calc
        use calc, only: dmqmc_energy, dmqmc_energy_squared, dmqmc_full_r2, dmqmc_rdm_r2
        use calc, only: dmqmc_correlation, dmqmc_staggered_magnetisation, dmqmc_kinetic_energy
        use calc, only: dmqmc_H0_energy, dmqmc_potential_energy, dmqmc_HI_energy
        use qmc_data, only: qmc_in_t, qmc_state_t
        use dmqmc_data

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        integer, intent(in) :: ireport
        real(dp), intent(in) :: ntot_particles(:)
        real, intent(in) :: elapsed_time
        logical, intent(in) :: comment
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(dmqmc_estimates_t), intent(in) :: dmqmc_estimates

        integer :: mc_cycles, i, j, ntypes, iunit

        iunit = 6

        ntypes = size(ntot_particles)

        mc_cycles = ireport*qmc_in%ncycles

        if (comment) then
            write (iunit,'(1X,"#",1X)', advance='no')
        else
            write (iunit,'(3X)', advance='no')
        end if

        call write_qmc_var(iunit, qs%mc_cycles_done+mc_cycles-qmc_in%ncycles)
        call write_qmc_var(iunit, qs%shift(1))
        call write_qmc_var(iunit, dmqmc_estimates%trace(1))
        ! The trace on the second replica.
        if (doing_dmqmc_calc(dmqmc_full_r2)) then
            call write_qmc_var(iunit, dmqmc_estimates%trace(2))
        end if

        ! Renyi-2 entropy for the full density matrix.
        if (doing_dmqmc_calc(dmqmc_full_r2)) then
            call write_qmc_var(iunit, dmqmc_estimates%numerators(full_r2_ind))
        end if

        ! Energy.
        if (doing_dmqmc_calc(dmqmc_energy)) then
            call write_qmc_var(iunit, dmqmc_estimates%numerators(energy_ind))
        end if

        ! Energy squared.
        if (doing_dmqmc_calc(dmqmc_energy_squared)) then
            call write_qmc_var(iunit, dmqmc_estimates%numerators(energy_squared_ind))
        end if

        ! Correlation function.
        if (doing_dmqmc_calc(dmqmc_correlation)) then
            call write_qmc_var(iunit, dmqmc_estimates%numerators(correlation_fn_ind))
        end if

        ! Staggered magnetisation.
        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
            call write_qmc_var(iunit, dmqmc_estimates%numerators(staggered_mag_ind))
        end if

        ! Kinetic energy
        if (doing_dmqmc_calc(dmqmc_kinetic_energy)) then
            call write_qmc_var(iunit, dmqmc_estimates%numerators(kinetic_ind))
        end if

        ! H^0 energy, where H = H^0 + V.
        if (doing_dmqmc_calc(dmqmc_H0_energy)) then
            call write_qmc_var(iunit, dmqmc_estimates%numerators(H0_ind))
        end if

        ! H^I energy, where H^I = exp(-(beta-tau)/2 H^0) H exp(-(beta-tau)/2. H^0).
        if (doing_dmqmc_calc(dmqmc_HI_energy)) then
            call write_qmc_var(iunit, dmqmc_estimates%numerators(HI_ind))
        end if

        ! Potential energy.
        if (doing_dmqmc_calc(dmqmc_potential_energy)) then
            call write_qmc_var(iunit, dmqmc_estimates%numerators(potential_ind))
        end if

        ! Renyi-2 entropy for all RDMs being sampled.
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            do i = 1, dmqmc_in%rdm%nrdms
                call write_qmc_var(iunit, dmqmc_estimates%inst_rdm%renyi_2(i))
            end do
        end if

        ! Traces for instantaneous RDM estimates.
        if (dmqmc_in%rdm%calc_inst_rdm) then
            do i = 1, dmqmc_in%rdm%nrdms
                do j = 1, ntypes
                    call write_qmc_var(iunit, dmqmc_estimates%inst_rdm%traces(j,i))
                end do
            end do
        end if

        ! The distribution of walkers on different excitation levels of the
        ! density matrix.
        if (dmqmc_in%calc_excit_dist) then
            do i = 0, ubound(dmqmc_estimates%excit_dist,1)
                call write_qmc_var(iunit, dmqmc_estimates%excit_dist(i)/ntot_particles(1))
            end do
        end if

        if (dmqmc_in%calc_mom_dist) then
            call write_momentum_array(dmqmc_estimates%mom_dist%n_k, dmqmc_estimates%mom_dist%kpoints, .true.)
        end if

        call write_qmc_var(iunit, ntot_particles(1))
        call write_qmc_var(iunit, qs%estimators(1)%tot_nstates)
        call write_qmc_var(iunit, qs%estimators(1)%tot_nspawn_events)
        call write_qmc_var(iunit, qs%spawn_store%rspawn, low_prec=.true.)
        call write_qmc_var(iunit, elapsed_time/qmc_in%ncycles, low_prec=.true.)

        write (iunit,'()')

    end subroutine write_dmqmc_report

    ! -- write_qmc_var --
    ! Write a QMC variable to a fixed-width column (header created in write_qmc_report_header).
    ! In:
    !    io: unit to write to.
    !    val: value to print out.
    !    low_prec (optional, real values only, default false): only write 4 decimal places instead of 10.

    subroutine write_qmc_var_int_32(io, val, sep)

        use const, only: int_32

        integer, intent(in) :: io
        integer(int_32), intent(in) :: val
        character(1), intent(in), optional :: sep
        character(1) :: sep_loc

        sep_loc = ' '
        if (present(sep)) sep_loc = sep
        write (io, '(i14,a1,1X)', advance='no') val, sep_loc

    end subroutine write_qmc_var_int_32

    subroutine write_qmc_var_int_64(io, val, sep)

        use const, only: int_64

        integer, intent(in) :: io
        integer(int_64), intent(in) :: val
        character(1), intent(in), optional :: sep
        character(1) :: sep_loc

        sep_loc = ' '
        if (present(sep)) sep_loc = sep
        write (io, '(i14,a1,1X)', advance='no') val, sep_loc

    end subroutine write_qmc_var_int_64

    subroutine write_qmc_var_real_sp(io, val, low_prec, sep)

        use const, only: sp, dp

        integer, intent(in) :: io
        real(sp), intent(in) :: val
        logical, intent(in), optional :: low_prec
        character(1), intent(in), optional :: sep
        character(1) :: sep_loc

        sep_loc = ' '
        if (present(sep)) sep_loc = sep
        ! Just forward to dp version to keep identical formatting.
        call write_qmc_var(io, real(val, dp), low_prec, sep_loc)

    end subroutine write_qmc_var_real_sp

    subroutine write_qmc_var_real_dp(io, val, low_prec, sep)

        use const, only: dp

        integer, intent(in) :: io
        real(dp), intent(in) :: val
        logical, intent(in), optional :: low_prec
        logical :: low_prec_loc
        character(1), intent(in), optional :: sep
        character(1) :: sep_loc

        sep_loc = ' '
        if (present(sep)) sep_loc = sep

        low_prec_loc = .false.
        if (present(low_prec)) low_prec_loc = low_prec

        if (low_prec_loc) then
            write (io, '(f8.4,a1,1X)', advance='no') val, sep_loc
        else
            write (io, '(es17.10,a1,4X)', advance='no') val, sep_loc
        end if

    end subroutine write_qmc_var_real_dp

    subroutine write_momentum_array_header(header, compress_k, kpoints)

        ! Write header for momentum space array estimator.

        ! In:
        !    header: base column header
        !    compress_k: if true only unique kpoints (in terms of the single-particle
        !        eigenvalues) will be printed with with results averaged over degeneracies.
        !    kvals: unique kpoints at which the correlation function is evaluated at.
        !        To be used with compress_k is true. Otherwise the orbital index
        !        will be printed.

        use utils, only: int_fmt

        character(5), intent(in) :: header
        logical, intent(in) :: compress_k
        real(p), intent(in) :: kpoints(:)

        integer :: iorb, iunit
        real(p) :: k
        character(14) :: klabel

        iunit = 6

        if (compress_k) then
            k = kpoints(1)
            do iorb = 2, size(kpoints)
                if (kpoints(iorb) > k + depsilon) then
                    write(klabel, '(es14.8)') k
                    call write_column_title(iunit, trim(header)//klabel)
                    k = kpoints(iorb)
                end if
            end do
            write(klabel, '(es14.8)') k
            call write_column_title(iunit, trim(header)//klabel)
        else
            do iorb = 1, size(kpoints)
                write (klabel, '('//int_fmt(iorb)//')') iorb
                call write_column_title(iunit, trim(header)//klabel)
            end do
        end if

    end subroutine write_momentum_array_header

    subroutine write_momentum_array(f_k, grid_points, compress_k)

        ! Write header for momentum space array estimator.

        ! In:
        !    f_k: momentum space correlation function to be printed.
        !    grid_points: unique kpoints at which the correlation function is evaluated at.
        !        To be used with compress_k is true.
        !    compress_k: if true only unique kpoints (in magnitude) will
        !        be printed with with results averaged over degeneracies.

        real(p), intent(in) :: f_k(:), grid_points(:)
        logical, intent(in) :: compress_k

        integer :: i, deg, iunit
        real(p) :: kp, compressed

        iunit = 6

        kp = grid_points(1)
        compressed = f_k(1)
        deg = 1
        if (compress_k) then
            do i = 2, size(f_k)
                if (grid_points(i) > kp + depsilon) then
                    call write_qmc_var(iunit, compressed/deg)
                    kp = grid_points(i)
                    compressed = f_k(i)
                    deg = 1
                else
                    compressed = compressed + f_k(i)
                    deg = deg + 1
                end if
            end do
            call write_qmc_var(iunit, compressed/deg) ! last point
        else
            do i = 1, size(f_k)
                call write_qmc_var(iunit, compressed/deg)
            end do
        end if

    end subroutine write_momentum_array

end module qmc_io
