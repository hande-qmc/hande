module sg_correction_energies

    ! Module for calculating the finite size correction to the energy via S(G)

    use sg_correction, only: sg_correction_t
    use qmc_data, only: qmc_state_t
    use system, only: sys_t
    use const, only: p

    implicit none
    
    contains

        subroutine calculate_uncorrected_sg_energy(qs, sys)

                ! Average the S(G) values at the end of each report loop and calculate 
                ! the (uncorrected) correlation energy.

                ! In:
                !   sys: System (sys_t) object. Contains lattice box length information. 
                ! In/Out:
                !   qs: QMC state (qmc_state_t) object. Contains various information about the QMC state. 
            
                type(qmc_state_t), intent(inout) :: qs 
                type(sys_t), intent(in) :: sys

                real(p) :: mag2
                integer :: i
                 
                ! The normalisation by N0 and averages over the number of cycles per repoort loop 
                ! that were previously done here are now in energy_evaluation.
                
                do i = 1, size(qs%sg_correction%sg)
                    
                    ! Prefactor converts from dimensionless to real units
                    mag2 = ((2 * acos(-1.0_p)) / sys%lattice%box_length(1))**2 * real(qs%sg_correction%q_vectors(1, i)**2 & 
                            + qs%sg_correction%q_vectors(2, i)**2 &
                            + qs%sg_correction%q_vectors(3, i)**2, p)

                    ! Skip |G| = 0
                    if (mag2 /= 0.0_p) then 
                        qs%estimators%energy_sg = qs%estimators%energy_sg & 
                            + (4 * acos(-1.0_p) / mag2) * qs%sg_correction%sg(i) 
                    end if
                end do 

        end subroutine calculate_uncorrected_sg_energy

        subroutine calculate_corrected_sg_energy(qs, sys)

                type(qmc_state_t), intent(inout) :: qs 
                type(sys_t), intent(in) :: sys
                
                real(p) :: mag2, window, Gc, h
                integer :: i, Nk, Mk
                real(p) :: contribution, term1, term2, term3, term4
                real(p) :: real_unit_conversion
                
                ! Calculate short range part with non-interpolated q_vectors and S(G)
                qs%estimators%corrected_energy_sg = 0.0_p

                do i = 1, size(qs%sg_correction%q_vectors, 2) 

                    mag2 = ((2 * acos(-1.0)) / sys%lattice%box_length(1))**2 * real(qs%sg_correction%q_vectors(1, i)**2 & 
                            + qs%sg_correction%q_vectors(2, i)**2 &
                            + qs%sg_correction%q_vectors(3, i)**2, p)
                    
                    Gc = qs%sg_correction%unique_q_mags(2) ! Largest non-zero magnitude
                    window = window_function(sqrt(mag2), Gc)
                    
                    ! Skip |G| = 0. Short range contriubtion to corrected energy
                    if (mag2 /= 0.0_p) then 
                        qs%estimators%corrected_energy_sg = qs%estimators%corrected_energy_sg & 
                            + (4 * acos(-1.0_p) / mag2) * (1 - window) * qs%sg_correction%sg(i) 
                    end if

                end do

                ! Long range contriubtion (correction) to energy 
                Nk = 1 ! Number of k-points per cell. Need to find where this is given 
                Mk = Nk**3
                
                ! Simpson quadrature. See thesis for derivation of integrand. 
                real_unit_conversion = ((2 * acos(-1.0_p)) / sys%lattice%box_length(1))
                h = (qs%sg_correction%q_spline_vals(qs%sg_correction%n_spline_points) &
                  - qs%sg_correction%q_spline_vals(1)) / (qs%sg_correction%n_spline_points - 1)
                term1 = qs%sg_correction%sg_spline_vals(1) * window_function(qs%sg_correction%q_spline_vals(1), Gc)
                term4 = qs%sg_correction%sg_spline_vals(qs%sg_correction%n_spline_points) & 
                      * window_function(qs%sg_correction%q_spline_vals(qs%sg_correction%n_spline_points), Gc)
                term2 = 0.0_p
                term3 = 0.0_p

                do i = 2, qs%sg_correction%n_spline_points - 1 
                    
                    contribution = qs%sg_correction%sg_spline_vals(i) &
                                 * window_function(qs%sg_correction%q_spline_vals(i), Gc)

                    if (mod(i, 2) == 0) then 
                        term2 = term2 + contribution  
                    else if (mod(i, 2) /= 0) then 
                        term3 = term3 + contribution
                    end if 
                end do

                qs%estimators%corrected_energy_sg = qs%estimators%corrected_energy_sg &
                    + real_unit_conversion * Mk * sys%lattice%box_length(1)**3 &
                    * (2 / acos(-1.0_p)) * (h / 3.0_p) * (term1 + 4 * term2 + 2 * term3 + term4)

            end subroutine calculate_corrected_sg_energy

            pure function window_function(magnitude, Gc) result (window)

                real(p), intent(in) :: magnitude 
                real(p), intent(in) :: Gc 

                real(p) :: window 

                ! Heaviside function
                if ((Gc - magnitude) >= 0.0_p) then 
                    window = 0.5_p * (cos(acos(-1.0_p) * (magnitude / Gc)) + 1.0_p) 
                else if ((Gc - magnitude) < 0.0_p) then 
                    window = 0.0_p
                end if
                
            end function window_function

end module sg_correction_energies

