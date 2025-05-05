module sg_correction_energies

    ! Module for calculating the finite size correction to the energy via S(G)

    use sg_correction, only: sg_correction_t
    use qmc_data, only: qmc_state_t
    use system, only: sys_t
    use const, only: p

    implicit none
    
    contains

        subroutine calculate_uncorrected_sg_energy(sg, box_length, q_vectors, energy_sg)

                ! Calculate the uncorrected energy via S(G) for comparison with the projected energy to validate the method.    
                
                ! In:
                !   sg: The now full array of S(G) values.
                !   box_length: Side length of the simulation cell.
                !   q_vectors: Physical momentum transfer vectors of the system.
                ! In/Out:
                !   energy_sg: Correlation energy of the system calculated via S(G)
                
                real(p), intent(in) :: sg(:) 
                real(p), intent(in) :: box_length 
                integer, intent(in) :: q_vectors(:, :)
                real(p), intent(inout) :: energy_sg(:)

                real(p) :: mag2
                integer :: i
                 
                ! The normalisation by N0 and averages over the number of cycles per repoort loop are in energy_evaluation.
                
                do i = 1, size(sg)
                    
                    ! Prefactor converts from dimensionless to real units
                    mag2 = ((2 * acos(-1.0_p)) / box_length)**2 * real(q_vectors(1, i)**2 & 
                                               + q_vectors(2, i)**2 + q_vectors(3, i)**2, p)

                    ! Skip |G| = 0
                    if (mag2 /= 0.0_p) then 
                        energy_sg = energy_sg + (4 * acos(-1.0_p) / mag2) * sg(i) 
                    end if
                end do 

        end subroutine calculate_uncorrected_sg_energy

        subroutine calculate_corrected_sg_energy(corrected_energy_sg, q_vectors, box_length, unique_q_mags, sg, q_spline_vals, &
                                                 sg_spline_vals, n_spline_points)
                
                ! Calculate the long range correction to the energy via S(G)
                
                ! In: 
                !   q_vectors: List of unique momentum transfers within the system. 
                !   box_length: Side length of the simulation cell. 
                !   unique_q_mags: Unique magnitudes of the q vectors. 
                !   sg: Array of S(G) values with index corresponding to the q vector which will give that value.
                !   q_spline_vals: Uniform distribution of momentum transfer magnitudes in the range q_min to q_max.  
                !   sg_spline_vals: Values of S(G) evaluated at the q_spline_vals points using the cubic spline. 
                !   n_spline_points: Total number of points which were evaluated on the spline. 
                ! In/Out: 
                !   corrected_energy_sg: Correlation energy of the system with the long range correction added. 
                
                real(p), intent(inout) :: corrected_energy_sg(:) 
                integer, intent(in) :: q_vectors(:, :) 
                real(p), intent(in) :: box_length 
                real(p), intent(in) :: unique_q_mags(:)
                real(p), intent(in) :: sg(:)
                real(p), intent(in) :: q_spline_vals(:), sg_spline_vals(:)
                integer, intent(in) :: n_spline_points

                real(p) :: mag2, window, Gc, h
                integer :: i, Nk, Mk
                real(p) :: contribution, term1, term2, term3, term4
                real(p) :: real_unit_conversion
                
                ! Calculate short range part with non-interpolated q_vectors and S(G)
                corrected_energy_sg = 0.0_p
                do i = 1, size(q_vectors, 2) 
                    mag2 = ((2 * acos(-1.0)) / box_length)**2 * real(q_vectors(1, i)**2 &
                                             + q_vectors(2, i)**2 + q_vectors(3, i)**2, p)
                    
                    Gc = unique_q_mags(2) ! Largest non-zero magnitude
                    window = window_function(sqrt(mag2), Gc)
                    
                    ! Skip |G| = 0. Short range contriubtion to corrected energy
                    if (mag2 /= 0.0_p) then 
                        corrected_energy_sg = corrected_energy_sg + (4 * acos(-1.0_p) / mag2) * (1 - window) * sg(i) 
                    end if

                end do

                ! Long range contriubtion (correction) to energy 
                Nk = 1 ! Number of k-points per cell. Hard coded as 1 for UEG. Will need to change to be more general.
                Mk = Nk**3
                
                ! Composite Simpson quadrature to evaluate long range correction 
                real_unit_conversion = ((2 * acos(-1.0_p)) / box_length)
                h = (q_spline_vals(n_spline_points) - q_spline_vals(1)) / (n_spline_points - 1)
                term1 = sg_spline_vals(1) * window_function(q_spline_vals(1), Gc)
                term4 = sg_spline_vals(n_spline_points) * window_function(q_spline_vals(n_spline_points), Gc)
                term2 = 0.0_p
                term3 = 0.0_p

                do i = 2, (n_spline_points - 1) 
                    
                    contribution = sg_spline_vals(i) * window_function(q_spline_vals(i), Gc)

                    if (mod(i, 2) == 0) then 
                        term2 = term2 + contribution  
                    else if (mod(i, 2) /= 0) then 
                        term3 = term3 + contribution
                    end if 
                end do

                corrected_energy_sg = corrected_energy_sg + real_unit_conversion * Mk * box_length**3 &
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

