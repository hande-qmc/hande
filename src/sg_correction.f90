module sg_correction

    ! Module for calculating the finite size correction to the correlation energy

    use system, only: sys_t
    use qmc_data, only: qmc_state_t, qmc_in_t
    use excitations, only: excit_t
    use checking, only: check_allocate
    use sort, only: insert_sort_unique_vector
    use const, only: p

    implicit none

    contains

        subroutine get_q_vectors(sys, qs)

            ! Find all unique pairs of momentum transfer vectors (q_vectors) between basis functions 
            ! and allocate sg_array which will contain the corresponding S(G) contribution of each
            ! q_vector.

            ! In:
            !   sys: System (sys_t) object. Contains basis set information.
            ! In/Out:
            !   qs: QMC state (qmc_state_t) object. Contains various information about the QMC state. 
            
            type(sys_t), intent(in) :: sys 
            type(qmc_state_t), intent(inout) :: qs 

            integer :: ispace, ierr, ki, ka, q(3)

            do ispace = 1, qs%psip_list%nspaces 

                ! Get all unique q vectors
                do ki = 1, sys%basis%nbasis
                    do ka = 1, sys%basis%nbasis 
                        q = sys%basis%basis_fns(ki)%l - sys%basis%basis_fns(ka)%l
                        call insert_sort_unique_vector(qs%estimators(ispace)%q_vectors, q)
                    end do 
                end do 

                ! Allocate corresponding S(G) array 
                allocate(qs%estimators(ispace)%sg(size(qs%estimators(ispace)%q_vectors, 2)), stat=ierr)
                call check_allocate('sg', size(qs%estimators(ispace)%q_vectors, 2), ierr)

                qs%estimators(ispace)%sg = 0.0_p

            end do

        end subroutine get_q_vectors

        subroutine generate_unit_fibonacci_sphere(n_points, points)

            ! Generate a fibonacci unit sphere with n_points points. 

            ! In:
            !   n_points: Number of points on the fibonacci sphere. 
            ! Out:
            !   points: Array of 3D vectors containing the fibonacci sphere points.
                
            integer, intent(in) :: n_points
            real(p), intent(out) :: points(3, n_points)

            integer :: n, k, idx
            real(p) :: phi, lat, lon

            phi = (1.0_p + sqrt(5.0_p)) / 2.0_p 
            n = (n_points - 1) / 2

            do k = -n, n 
               
                idx = k + (n + 1) ! Maps k values to correct index

                lat = asin(2.0_p * k / (2.0_p * n + 1.0_p))
                lon = 2.0_p * acos(-1.0_p) * k / phi    

                points(1, idx) = cos(lat) * cos(lon)
                points(2, idx) = cos(lat) * sin(lon)
                points(3, idx) = sin(lat)

            end do
            
        end subroutine generate_unit_fibonacci_sphere

        subroutine generate_unique_mag_fibonacci_sphere(qs, tolerance, fibonacci_points_unit, n_points)

            ! Calculates the magnitudes of the stored q_vectors and creates a list 
            ! containing only unique magnitudes. 

            ! In:
            !   tolerance: Real number which dictates to what extent 2 magnitudes close in value 
            !              are the same.
            !   n_points: Number of points on each fibonacci sphere.
            !   fibonacci_points_unit: Fibonacci sphere points on a unit sphere. 
            ! In/Out:
            !   qs: QMC state (qmc_state_t) object. Contains various information about the QMC state.

            type(qmc_state_t), intent(inout) :: qs
            real(p), intent(in) :: tolerance 
            real(p), intent(in) :: fibonacci_points_unit(:, :)
            integer, intent(in) :: n_points

            integer :: i 
            integer :: n_unique

            ! Get magnitudes of q_vectors
            allocate(qs%estimators(1)%q_mags(size(qs%estimators(1)%q_vectors)))
            do i = 1, size(qs%estimators(1)%q_vectors, 2)
                qs%estimators(1)%q_mags(i) = sqrt(real(qs%estimators(1)%q_vectors(1, i)**2 + &
                                                       qs%estimators(1)%q_vectors(2, i)**2 + &
                                                       qs%estimators(1)%q_vectors(3, i)**2, p))
            end do

            ! Find unique magnitudes only 
            allocate(qs%estimators(1)%unique_q_mags(size(qs%estimators(1)%q_vectors)))
            n_unique = 0 
            do i = 1, size(qs%estimators(1)%q_vectors, 2)
                if (i == 1 .or. abs(qs%estimators(1)%q_mags(i) - qs%estimators(1)%q_mags(i - 1)) > tolerance) then 
                    n_unique = n_unique + 1 
                    qs%estimators(1)%unique_q_mags(n_unique) = qs%estimators(1)%q_mags(i)
                end if
            end do 
            
            ! Remove excess size
            qs%estimators(1)%unique_q_mags = qs%estimators(1)%unique_q_mags(1:n_unique)

            ! Generate 1 fibonacci sphere per unique magnitude
            allocate(qs%estimators(1)%fibonacci_spheres(3, n_points, n_unique))
            do i = 1, n_unique 
                qs%estimators(1)%fibonacci_spheres(:, :, i) = qs%estimators(1)%unique_q_mags(i) * fibonacci_points_unit
            end do

        end subroutine generate_unique_mag_fibonacci_sphere

        subroutine calculate_sg_contribution(sys, connection, qs, real_population, ispace)

            ! Calculate the contribution to S(G) of the current determinant or cluster.

            ! In:
            !   sys: System (sys_t) object. Contains basis set information.
            !   connection: Excitation (excit_t) object. Contains information about the excitation between 
            !               the current determinant or cluster and their children. 
            !   ispace: Integer identifying the current walker space. This is always 1 for S(G) correction (at present).
            !   real_population :: Real number which gives the population of the current determinant or cluster.
            ! In/Out:
            !   qs: QMC state (qmc_state_t) object. Contains various information about the QMC state.

            type(sys_t), intent(in) :: sys
            type(excit_t), intent(in) :: connection
            integer, intent(in) :: ispace
            real(p), intent(in) :: real_population
            type(qmc_state_t), intent(inout) :: qs

            real(p) :: N0, omega
            integer :: km(3), kn(3), kp(3), kq(3), G(3)
            integer :: sm, sn, sp, sq
            integer :: idx
            integer :: phase

            omega = sys%lattice%box_length(1)**3
            phase = merge(-1, 1, connection%perm)
            N0 = real(qs%psip_list%pops(ispace, 1), p) / qs%psip_list%pop_real_factor
    
            ! to_orb contains occupied orbitals, from_orb contains virtual orbitals
            km = sys%basis%basis_fns(connection%to_orb(1))%l 
            kn = sys%basis%basis_fns(connection%to_orb(2))%l
            kp = sys%basis%basis_fns(connection%from_orb(1))%l 
            kq = sys%basis%basis_fns(connection%from_orb(2))%l
            sm = sys%basis%basis_fns(connection%to_orb(1))%ms
            sn = sys%basis%basis_fns(connection%to_orb(2))%ms 
            sp = sys%basis%basis_fns(connection%from_orb(1))%ms 
            sq = sys%basis%basis_fns(connection%from_orb(2))%ms

            ! Term 1
            G = km - kp
            if (sm == sp .and. sn == sq .and. all(kn - kq == -G)) then 
                idx = find_q_index(qs%estimators(ispace)%q_vectors, G)
                if (idx /= 0) then 
                    qs%estimators(ispace)%sg(idx) = qs%estimators(ispace)%sg(idx) &
                        + (phase * real_population) / (N0 * omega)
                end if 
            end if

            ! Term 2 
            G = km - kq 
            if (sm == sq .and. sn == sp .and. all(kn - kp == - G)) then
                idx = find_q_index(qs%estimators(ispace)%q_vectors, G)
                if (idx /= 0) then
                    qs%estimators(ispace)%sg(idx) = qs%estimators(ispace)%sg(idx) &
                        - (phase * real_population) / (N0 * omega)
                end if 
            end if

        end subroutine calculate_sg_contribution

        subroutine calculate_uncorrected_sg_energy(qs, qmc_in, sys)

            ! Average the S(G) values at the end of each report loop and calculate 
            ! the (uncorrected) correlation energy.

            ! In:
            !   qmc_in: Input options which realte to the QMC method, (qmc_in_t) type.
            !   sys: System (sys_t) object. Contains lattice box length information. 
            ! In/Out:
            !   qs: QMC state (qmc_state_t) object. Contains various information about the QMC state. 
        
            type(qmc_state_t), intent(inout) :: qs 
            type(qmc_in_t), intent(in) :: qmc_in
            type(sys_t), intent(in) :: sys

            real(p) :: mag2
            integer :: i
            
            do i = 1, size(qs%estimators(1)%sg)
                qs%estimators(1)%sg(i) = qs%estimators(1)%sg(i) / qmc_in%ncycles
                
                ! Prefactor converts from dimensionless to real units
                mag2 = ((2 * acos(-1.0)) / sys%lattice%box_length(1))**2 * real(qs%estimators(1)%q_vectors(1, i)**2 & 
                        + qs%estimators(1)%q_vectors(2, i)**2 &
                        + qs%estimators(1)%q_vectors(3, i)**2)

                ! Skip |G| = 0
                if (mag2 /= 0.0_p) then 
                    qs%estimators(1)%energy_sg = qs%estimators(1)%energy_sg & 
                        + (4 * acos(-1.0) / mag2) * qs%estimators(1)%sg(i)
                end if
            end do 

        end subroutine calculate_uncorrected_sg_energy

        subroutine interpolate_sg_fibonacci(qs)

            ! Interpolate the calculated S(G) values of each unique magnitude onto the fibonacci
            ! sphere points using inverse distance weighting interpolation. 
            ! See: https://en.wikipedia.org/wiki/Inverse_distance_weighting.

            ! In/Out:
            !   qs: QMC state (qmc_state_t) object. Contains various information about the QMC state. 
        
            type(qmc_state_t), intent(inout) :: qs 
            
            real(p) :: locality = 6.0_p
            real(p) :: d, w, sum_w, sum_wS
            real(p) :: tolerance = 1e-8_p
            integer :: i, j, k

            ! Interpolate for each fibonacci sphere 
            do i = 1, size(qs%estimators(1)%unique_q_mags)
                
                ! Skip interpolation for zero magnitude sphere and just set S(G) at all points to zero
                if (abs(qs%estimators(1)%unique_q_mags(i)) < tolerance) then 
                    qs%estimators(1)%sg_fibonacci(:, i) = 0.0_p
                    cycle 
                end if 
                
                ! Loop over fibonacci points 
                do j = 1, size(qs%estimators(1)%fibonacci_spheres, 2)

                    sum_w = 0.0_p 
                    sum_wS = 0.0_p

                    ! Sum over S(G) values / q_vectors
                    do k = 1, size(qs%estimators(1)%q_vectors, 2)

                        d = sqrt((qs%estimators(1)%fibonacci_spheres(1, j, i) - qs%estimators(1)%q_vectors(1, k))**2 &
                               + (qs%estimators(1)%fibonacci_spheres(2, j, i) - qs%estimators(1)%q_vectors(2, k))**2 &
                               + (qs%estimators(1)%fibonacci_spheres(3, j, i) - qs%estimators(1)%q_vectors(3, k))**2)

                        if (d < tolerance) then
                            sum_wS = qs%estimators(1)%sg(k)
                            sum_w = 1.0_p
                            exit 
                        end if 

                        w = 1.0_p / d**locality 
                        sum_w = sum_w + w 
                        sum_wS = sum_wS + (w * qs%estimators(1)%sg(k))

                    end do 

                    qs%estimators(1)%sg_fibonacci(j, i) = sum_wS / sum_w 

                end do 
            end do

        end subroutine interpolate_sg_fibonacci

        subroutine spherical_average_sg(qs)

            type (qmc_state_t), intent(inout) :: qs

            real(p) :: sum_sg
            integer :: i, j

            do i = 1, size(qs%estimators(1)%sg_fibonacci, 2)
                
                sum_sg = 0.0_p

                do j = 1, size(qs%estimators(1)%sg_fibonacci, 1)
                    sum_sg = sum_sg + qs%estimators(1)%sg_fibonacci(j, i)
                end do 

                qs%estimators(1)%spherically_averaged_sg(i) = sum_sg / size(qs%estimators(1)%sg_fibonacci, 1)
            end do 
       
        end subroutine spherical_average_sg

        subroutine cubic_spline_sg(qs)

            type(qmc_state_t), intent(inout) :: qs

            integer :: n_points, i, w
            real(p), allocatable :: q_points(:), sg_points(:)
            real(p), allocatable :: h(:)
            real(p), allocatable :: lower_diag(:), main_diag(:), upper_diag(:) ! Tridiagonal matrix
            real(p), allocatable :: m_vector(:), right_vector(:)

            n_points = size(qs%estimators(1)%spherically_averaged_sg)
            allocate(q_points(n_points), sg_points(n_points))
            sg_points = qs%estimators(1)%spherically_averaged_sg
            q_points = qs%estimators(1)%unique_q_mags

            allocate(h(n_points - 1), lower_diag(n_points), main_diag(n_points), upper_diag(n_points), &
                     right_vector(n_points), m_vector(n_points))
            
            do i = 1, n_points - 1 
                h(i) = q_points(i + 1) - q_points(i)
            end do 
            
            main_diag = 0.0_p
            lower_diag = 0.0_p
            upper_diag = 0.0_p
            right_vector = 0.0_p

            ! Natural spline boundary conditions. S''(G) = 0.  
            main_diag(1) = 1.0_p
            main_diag(n_points) = 1.0_p
            right_vector(1) = 0.0_p
            right_vector(n_points) = 0.0_p
            
            ! Build matrix
            do i = 2, n_points - 1 
                lower_diag(i) = h(i - 1) ! lower_diag(1) = 0 as there is no matrix element (1, 0)
                main_diag(i) = 2.0_p * (h(i-1) + h(i)) ! Boundary conditions give main_diag(1) and main_diag(n_points) = 1
                upper_diag(i) = h(i) ! upper_diag(1) = 0 as there is no matrix element (0, 1)
                right_vector(i) = 6.0_p * ((sg_points(i + 1) - sg_points(i)) / h(i) &
                                       - (sg_points(i) - sg_points(i - 1)) / h(i - 1)) ! Boundary conditions give right_vector(1) 
                                                                                       ! and right_vector(n_points) = 0
            end do

            ! Thomas algorithm
            do i = 2, n_points 
                w = lower_diag(i) / main_diag(i - 1)
                main_diag(i) = main_diag(i) - w * upper_diag(i - 1)
                right_vector(i) = right_vector(i) - w * right_vector(i - 1)
            end do 

            ! Back substitution 
            m_vector(n_points) = right_vector(n_points) / main_diag(n_points)
            do i = n_points - 1, 1, -1 ! Reverse loop
                m_vector(i) = (right_vector(i) - upper_diag(i) * m_vector(i + 1)) / main_diag(i)
            end do 

            ! Compute coefficients 
            do i = 1, n_points - 1
                qs%estimators(1)%spline_coeffs_a(i) = (m_vector(i + 1) - m_vector(i)) / (6.0_p * h(i))
                qs%estimators(1)%spline_coeffs_b(i) = m_vector(i) / 2.0_p
                qs%estimators(1)%spline_coeffs_c(i) = (sg_points(i + 1) - sg_points(i)) / h(i) &
                                                    - (h(i) / 6.0_p) * (2 * m_vector(i) + m_vector(i + 1))
                qs%estimators(1)%spline_coeffs_d(i) = sg_points(i)
                
            end do 

        end subroutine cubic_spline_sg

        pure function evaluate_spline(q, n_knots, qs) result(sg)

            real(p), intent(in) :: q 
            integer, intent(in) :: n_knots 
            type(qmc_state_t), intent(in) :: qs 

            integer :: i 
            real(p) :: delta
            real(p) :: sg

            ! Find which interval contains q 
            do i = 1, n_knots - 1 
                if (q < qs%estimators(1)%unique_q_mags(i + 1)) exit 
            end do 

            ! Evaluate spline in this interval 
            delta = q - qs%estimators(1)%unique_q_mags(i)
            sg = qs%estimators(1)%spline_coeffs_a(i) * delta**3 + qs%estimators(1)%spline_coeffs_b(i) * delta**2 &
                                                                + qs%estimators(1)%spline_coeffs_c(i) * delta &
                                                                + qs%estimators(1)%spline_coeffs_d(i)

        end function evaluate_spline

        subroutine calculate_corrected_sg_energy(qs, sys, q_spline_vals, sg_spline_vals, n_spline_points)

            type(qmc_state_t), intent(inout) :: qs 
            type(sys_t), intent(in) :: sys
            integer, intent(in) :: n_spline_points
            real(p), intent(in) :: q_spline_vals(n_spline_points), sg_spline_vals(n_spline_points) 
            
            real(p) :: mag2, window, Gc, h
            integer :: i, Nk, Mk
            real(p) :: contribution, term1, term2, term3, term4
            real(p) :: real_unit_conversion
            
            ! Calculate short range part with non-interpolated q_vectors and S(G)
            qs%estimators(1)%corrected_energy_sg = 0.0_p

            do i = 1, size(qs%estimators(1)%q_vectors, 2) 

                mag2 = ((2 * acos(-1.0)) / sys%lattice%box_length(1))**2 * real(qs%estimators(1)%q_vectors(1, i)**2 & 
                        + qs%estimators(1)%q_vectors(2, i)**2 &
                        + qs%estimators(1)%q_vectors(3, i)**2)
                
                Gc = qs%estimators(1)%unique_q_mags(2) ! Largest non-zero magnitude
                window = window_function(sqrt(mag2), Gc)
                
                ! Skip |G| = 0. Short range contriubtion to corrected energy
                if (mag2 /= 0.0_p) then 
                    qs%estimators(1)%corrected_energy_sg = qs%estimators(1)%corrected_energy_sg & 
                        + (4 * acos(-1.0) / mag2) * (1 - window) * qs%estimators(1)%sg(i)
                end if

            end do

            ! Long range contriubtion (correction) to energy 
            Nk = 1 ! Number of k-points per cell. Need to find where this is given 
            Mk = Nk**3
            
            ! Simpson quadrature. See thesis for derivation of integrand. 
            real_unit_conversion = ((2 * acos(-1.0)) / sys%lattice%box_length(1))
            h = (q_spline_vals(n_spline_points) - q_spline_vals(1)) / (n_spline_points - 1)
            term1 = sg_spline_vals(1) * window_function(q_spline_vals(1), Gc)
            term4 = sg_spline_vals(n_spline_points) * window_function(q_spline_vals(n_spline_points), Gc)
            term2 = 0.0_p
            term3 = 0.0_p

            do i = 2, n_spline_points - 1 
                
                contribution = sg_spline_vals(i) * window_function(q_spline_vals(i), Gc)

                if (mod(i, 2) == 0) then 
                    term2 = term2 + contribution  
                else if (mod(i, 2) /= 0) then 
                    term3 = term3 + contribution
                end if 
            end do

            qs%estimators(1)%corrected_energy_sg = qs%estimators(1)%corrected_energy_sg &
                + real_unit_conversion * Mk * sys%lattice%box_length(1)**3 * (2 / acos(-1.0)) * (h / 3.0_p) &
                * (term1 + 4 * term2 + 2 * term3 + term4) 

        end subroutine calculate_corrected_sg_energy

        subroutine dump_sg_data(qs)

            ! Dump information about q_vectors and their contribution S(G) to file SG.

            ! In:
            !   qs: QMC state (qmc_state_t) object. Contains various information about the QMC state.

            type(qmc_state_t), intent(in) :: qs 

            integer :: sg_unit, ierr 
            integer :: i

            open(newunit=sg_unit, file="SG", status="replace", action="write", iostat=ierr)

            write(sg_unit, *) 'G vectors and their S(G) contribution: '
            do i = 1, size(qs%estimators(1)%sg)
                write(sg_unit, *) qs%estimators(1)%q_vectors(1, i), &
                                  qs%estimators(1)%q_vectors(2, i), &
                                  qs%estimators(1)%q_vectors(3, i), &
                                  qs%estimators(1)%sg(i) 
            end do

            write(sg_unit,*) 'n_unique =', size(qs%estimators(1)%unique_q_mags)
            write(sg_unit,*) 'n_avg    =', size(qs%estimators(1)%spherically_averaged_sg)
            
            write(sg_unit, *) '|G| and their spherically averaged S(G) values: '
            do i = 1, size(qs%estimators(1)%spherically_averaged_sg)
                write(sg_unit, *) qs%estimators(1)%unique_q_mags(i), &
                                  qs%estimators(1)%spherically_averaged_sg(i)
            end do 

            close(sg_unit)

        end subroutine dump_sg_data
        
        pure function find_q_index(q_vectors, G) result(idx)

            ! Compare a vector G to list of vectors q_vectors and find the 
            ! index of q_vectors where q_vectors(idx) == G. If not found 
            ! return 0.

            ! In:
            !   q_vectors: List of 3D vectors.
            !   G: One 3D vector.
            ! Out:
            !   idx: Index for which q_vectors(idx) == G. 0 if not Found. 

            integer, intent(in) :: q_vectors(:, :)
            integer, intent(in) :: G(3)
            integer :: idx, n, i
            
            ! Total number of q_vectors in list
            n = size(q_vectors, 2)
            idx = 0 ! If no match is found idx is returned as is

            do i = 1, n 
                if (all(q_vectors(:, i) == G)) then
                    idx = i 
                    return 
                end if
            end do 

        end function find_q_index

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

end module sg_correction
