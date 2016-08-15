module qmc_io

! QMC Output procedures.

use const
use spawn_data, only: spawn_t

implicit none

! [todo] - fix compilation so we can make module private by default -- too many modules implicitly reuse the module level imports here.
!private
!public :: write_qmc_report_header, write_qmc_report

interface write_qmc_var
    module procedure :: write_qmc_var_int
    module procedure :: write_qmc_var_int_64
    module procedure :: write_qmc_var_real_sp
    module procedure :: write_qmc_var_real_dp
end interface

contains

    subroutine write_qmc_report_header(ntypes, cmplx_est)

        ! In:
        !    ntypes: number of particle types being sampled.
        ! In (optional):
        !    cmplx_est: Print out information of cmplx_estlex estimators.

        use calc, only: doing_calc, hfs_fciqmc_calc

        integer, intent(in) :: ntypes
        logical, optional, intent(in) :: cmplx_est
        logical :: cmplx_est_set

        cmplx_est_set = .false.
        if (present(cmplx_est)) cmplx_est_set = cmplx_est

        ! Data table info.
        write (6,'(1X,"Information printed out every QMC report loop:",/)')
        write (6,'(1X,"Shift: the energy offset calculated at the end of the report loop.")')
        write (6,'(1X,"H_0j: <D_0|H|D_j>, Hamiltonian matrix element.")')
        write (6,'(1X,"N_j: population of Hamiltonian particles on determinant D_j.")')
        if (doing_calc(hfs_fciqmc_calc)) then
            write (6,'(1X,"O_0j: <D_0|O|D_j>, operator matrix element.")')
            write (6,'(1X,a67)') "N'_j: population of Hellmann--Feynman particles on determinant D_j."
            write (6,'(1X,"# HF psips: current total population of Hellmann--Feynman particles.")')
        end if

        write (6,'(1X,"# H psips: current total population of Hamiltonian particles.")')
        write (6,'(1X,"# states: number of many-particle states occupied.")')
        write (6,'(1X,"spawn_events: number of successful spawning events across all processors.")')
        write (6,'(1X,"R_spawn: average rate of spawning across all processors.")')
        write (6,'(1X,"time: average time per Monte Carlo cycle.",/)')
        write (6,'(1X,"Note that all particle populations are averaged over the report loop.",/)')

        write (6,'(1X,"#",1X)', advance='no')
        call write_column_title(6, 'iterations', int_val=.true., justify=1)
        call write_column_title(6, 'Shift')
        ! NOTE: HFS and complex are not currently compatible.
        if (cmplx_est_set) then
            call write_column_title(6, 'Re{\sum H_0j N_j}')
            call write_column_title(6, 'Im{\sum H_0j N_j}')
            call write_column_title(6, 'Re{N_0}')
            call write_column_title(6, 'Im{N_0}')
        else
            call write_column_title(6, '\sum H_0j N_j')
            call write_column_title(6, 'N_0')
            if (doing_calc(hfs_fciqmc_calc)) then
                call write_column_title(6, 'HF shift')
                call write_column_title(6, '\sum O_0j N_j')
                call write_column_title(6, "\sum H_0j N'_j")
                call write_column_title(6, "N'_0")
                call write_column_title(6, '# H psips')
            end if
            call write_column_title(6, '# H psips')
            if (doing_calc(hfs_fciqmc_calc)) call write_column_title(6, '# HF psips')
        end if
        call write_column_title(6, '# states', int_val=.true., justify=1)
        call write_column_title(6, '# spawn_events', int_val=.true., justify=1)
        call write_column_title(6, 'R_spawn', low_prec_val=.true.)
        call write_column_title(6, '  time', low_prec_val=.true., justify=2)

        write (6,'()')

        contains

            subroutine write_column_title(io, key, int_val, low_prec_val, justify)

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
                logical :: int_val_loc, low_prec_val_loc
                character(17) :: key_str
                integer :: justify_loc, str_len

                int_val_loc = .false.
                low_prec_val_loc = .false.
                justify_loc = 0
                if (present(int_val)) int_val_loc = int_val
                if (present(low_prec_val)) low_prec_val_loc = low_prec_val
                if (present(justify)) justify_loc = justify

                key_str = key
                str_len = len(key_str)
                if (int_val_loc) then
                    str_len = 14
                end if
                if (low_prec_val_loc) str_len = 8
                select case(justify_loc)
                case(1)
                    ! Right justify.
                    write (io,'("'//adjustr(key_str(:str_len))//'", 2X)', advance='no')
                case(2)
                    ! No justification.
                    write (io,'("'//key_str(:str_len)//'", 2X)', advance='no')
                case default
                    ! Left justify.
                    ! an extra space so the column lines up with the first non-sign entry in the value.
                    write (io,'(1X,"'//adjustl(key_str(:str_len-1))//'", 2X)', advance='no')
                end select

            end subroutine write_column_title

    end subroutine write_qmc_report_header

    subroutine write_qmc_report(qmc_in, qs, ireport, ntot_particles, elapsed_time, comment, non_blocking_comm, cmplx_est)

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
        !    cmplx_est: if true, doing calculation with real and imaginary walkers
        !       so need to print extra parameters.

        use calc, only: doing_calc, hfs_fciqmc_calc
        use qmc_data, only: qmc_in_t, qmc_state_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        integer, intent(in) :: ireport
        real(dp), intent(in) :: ntot_particles(:)
        real, intent(in) :: elapsed_time
        logical, intent(in) :: comment, non_blocking_comm
        logical, intent(in), optional :: cmplx_est

        logical :: cmplx_est_set
        integer :: mc_cycles, i, j, ntypes

        ntypes = size(ntot_particles)

        ! For non-blocking communications we print out the nth report loop
        ! after the (n+1)st iteration. Adjust mc_cycles accordingly
        if (.not. non_blocking_comm) then
            mc_cycles = ireport*qmc_in%ncycles
        else
            mc_cycles = (ireport-1)*qmc_in%ncycles
        end if

        cmplx_est_set = .false.
        if (present(cmplx_est)) cmplx_est_set = cmplx_est

        if (comment) then
            write (6,'(1X,"#",1X)', advance='no')
        else
            write (6,'(3X)', advance='no')
        end if

        call write_qmc_var(6, qs%mc_cycles_done+mc_cycles)
        call write_qmc_var(6, qs%shift(1))

        ! NOTE: HFS and complex are not currently compatible.
        if (cmplx_est_set) then
            call write_qmc_var(6, real(qs%estimators%proj_energy_comp, p))
            call write_qmc_var(6, aimag(qs%estimators%proj_energy_comp))
            call write_qmc_var(6, real(qs%estimators%D0_population_comp, p))
            call write_qmc_var(6, aimag(qs%estimators%D0_population_comp))
            call write_qmc_var(6, ntot_particles(1)+ntot_particles(2))
        else
            call write_qmc_var(6, qs%estimators%proj_energy)
            call write_qmc_var(6, qs%estimators%D0_population)
            if (doing_calc(hfs_fciqmc_calc)) then
            end if
            call write_qmc_var(6, ntot_particles(1))
            if (doing_calc(hfs_fciqmc_calc)) call write_qmc_var(6, ntot_particles(2))
        end if

        call write_qmc_var(6, qs%estimators%tot_nstates)
        call write_qmc_var(6, qs%estimators%tot_nspawn_events)

        call write_qmc_var(6, qs%spawn_store%rspawn, low_prec=.true.)
        call write_qmc_var(6, elapsed_time/qmc_in%ncycles, low_prec=.true.)

        write (6,'()')

    end subroutine write_qmc_report

    ! -- write_qmc_var --
    ! Write a QMC variable to a fixed-width column (header created in write_qmc_report_header).
    ! In:
    !    io: unit to write to.
    !    val: value to print out.
    !    low_prec (optional, real values only, default false): only write 4 decimal places instead of 10.

    subroutine write_qmc_var_int(io, val)

        integer, intent(in) :: io, val
        write (io, '(i14,2X)', advance='no') val

    end subroutine write_qmc_var_int

    subroutine write_qmc_var_int_64(io, val)

        use const, only: int_64
        integer(int_64), intent(in) :: io, val
        write (io, '(i14,2X)', advance='no') val

    end subroutine write_qmc_var_int_64

    subroutine write_qmc_var_real_sp(io, val, low_prec)

        use const, only: sp, dp

        integer, intent(in) :: io
        real(sp), intent(in) :: val
        logical, intent(in), optional :: low_prec
        logical :: low_prec_loc

        ! Just forward to dp version to keep identical formatting.
        call write_qmc_var(io, real(val, dp), low_prec)

    end subroutine write_qmc_var_real_sp

    subroutine write_qmc_var_real_dp(io, val, low_prec)

        use const, only: dp

        integer, intent(in) :: io
        real(dp), intent(in) :: val
        logical, intent(in), optional :: low_prec
        logical :: low_prec_loc

        low_prec_loc = .false.
        if (present(low_prec)) low_prec_loc = low_prec

        if (low_prec_loc) then
            write (io, '(f8.4,2X)', advance='no') val
        else
            write (io, '(es17.10,2X)', advance='no') val
        end if

    end subroutine write_qmc_var_real_dp

end module qmc_io
