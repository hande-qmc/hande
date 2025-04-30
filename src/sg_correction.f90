module sg_correction

    ! Module for calculating quantities relating to S(G) required to calculate finite size correction to the energy.

    use system, only: sys_t
    use excitations, only: excit_t
    use checking, only: check_allocate
    use sort, only: insert_sort_unique_vector
    use const, only: p

    implicit none

    type sg_correction_t

        integer, allocatable :: q_vectors(:, :)

        integer, allocatable :: q_map(:, :, :)

        real(p), allocatable :: q_mags(:)
        real(p), allocatable :: unique_q_mags(:)

        real(p), allocatable :: fibonacci_points_unit(:, :)
        real(p), allocatable :: fibonacci_spheres(:, :, :)

        real(p), allocatable :: sg(:)

        real(p), allocatable :: sg_fibonacci(:, :)

        real(p), allocatable :: spherically_averaged_sg(:)

        real(p), allocatable :: spline_coeffs_a(:)
        real(p), allocatable :: spline_coeffs_b(:)
        real(p), allocatable :: spline_coeffs_c(:)
        real(p), allocatable :: spline_coeffs_d(:)

        real(p) :: tolerance = 1e-8_p 
        integer :: n_fibonacci_points = 199, n_spline_points = 100

        real(p) :: q_min, q_max 
        real(p), allocatable :: q_spline_vals(:), sg_spline_vals(:)

    end type sg_correction_t

    contains

        subroutine get_q_vectors(sys, sg_correction)

            ! Find all unique pairs of momentum transfer vectors (q_vectors) between basis functions 
            ! and allocate sg_array which will contain the corresponding S(G) contribution of each
            ! q_vector.

            ! In:
            !   sys: System (sys_t) object. Contains basis set information.
            ! In/Out:
            !   sg_correction: (sg_correction_t) object. Contains information required to calculated energies via S(G). 
            
            type(sys_t), intent(in) :: sys 
            type(sg_correction_t), intent(inout) :: sg_correction

            integer :: i, ierr, ki, ka, q(3)
            integer :: qx_min, qx_max, qy_min, qy_max, qz_min, qz_max

            ! Get all unique q vectors
            do ki = 1, sys%basis%nbasis
                do ka = 1, sys%basis%nbasis 
                    q = sys%basis%basis_fns(ki)%l - sys%basis%basis_fns(ka)%l
                    call insert_sort_unique_vector(sg_correction%q_vectors, q)
                end do 
            end do 

            ! Allocate corresponding S(G) array 
            allocate(sg_correction%sg(size(sg_correction%q_vectors, 2)), stat=ierr)
            call check_allocate('sg', size(sg_correction%q_vectors, 2), ierr)

            sg_correction%sg = 0.0_p

            ! Construct lookup map for the index of the q_vectors 
            qx_min = minval(sg_correction%q_vectors(1, :))
            qx_max = maxval(sg_correction%q_vectors(1, :))
            qy_min = minval(sg_correction%q_vectors(2, :))
            qy_max = maxval(sg_correction%q_vectors(2, :))
            qz_min = minval(sg_correction%q_vectors(3, :))
            qz_max = maxval(sg_correction%q_vectors(3, :))
            allocate(sg_correction%q_map(qx_min:qx_max, qy_min:qy_max, qz_min:qz_max))
            sg_correction%q_map = 0.0_p ! If q vector not in the list of q_vectors the map should return 0 

            do i = 1, size(sg_correction%q_vectors, 2)
                q = sg_correction%q_vectors(:, i)
                sg_correction%q_map(q(1), q(2), q(3)) = i
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

        subroutine generate_unique_mag_fibonacci_sphere(sg_correction, fibonacci_points_unit, n_points)

            ! Calculates the magnitudes of the stored q_vectors and creates a list 
            ! containing only unique magnitudes. 

            ! In:
            !   n_points: Number of points on each fibonacci sphere.
            !   fibonacci_points_unit: Fibonacci sphere points on a unit sphere. 
            ! In/Out:
            !   sg_correction: (sg_correction_t) object. Contains information required to calculated energies via S(G). 

            type(sg_correction_t), intent(inout) :: sg_correction
            real(p), intent(in) :: fibonacci_points_unit(:, :)
            integer, intent(in) :: n_points

            integer :: i 
            integer :: n_unique

            ! Get magnitudes of q_vectors
            allocate(sg_correction%q_mags(size(sg_correction%q_vectors)))
            do i = 1, size(sg_correction%q_vectors, 2)
                sg_correction%q_mags(i) = sqrt(real(sg_correction%q_vectors(1, i)**2 + &
                                                       sg_correction%q_vectors(2, i)**2 + &
                                                       sg_correction%q_vectors(3, i)**2, p))
            end do

            ! Find unique magnitudes only 
            allocate(sg_correction%unique_q_mags(size(sg_correction%q_vectors)))
            n_unique = 0 
            do i = 1, size(sg_correction%q_vectors, 2)
                if (i == 1 .or. abs(sg_correction%q_mags(i) - sg_correction%q_mags(i - 1)) > sg_correction%tolerance) then 
                    n_unique = n_unique + 1 
                    sg_correction%unique_q_mags(n_unique) = sg_correction%q_mags(i)
                end if
            end do 
            
            ! Remove excess size
            sg_correction%unique_q_mags = sg_correction%unique_q_mags(1:n_unique)

            ! Generate 1 fibonacci sphere per unique magnitude
            allocate(sg_correction%fibonacci_spheres(3, n_points, n_unique))
            do i = 1, n_unique 
                sg_correction%fibonacci_spheres(:, :, i) = sg_correction%unique_q_mags(i) * fibonacci_points_unit
            end do

        end subroutine generate_unique_mag_fibonacci_sphere

        subroutine calculate_sg_contribution(sys, connection, sg_correction, real_population)

            ! Calculate the contribution to S(G) of the current determinant or cluster.

            ! In:
            !   sys: System (sys_t) object. Contains basis set information.
            !   connection: Excitation (excit_t) object. Contains information about the excitation between 
            !               the current determinant or cluster and their children. 
            !   real_population :: Real number which gives the population of the current determinant or cluster.
            ! In/Out:
            !   sg_correction: (sg_correction_t) object. Contains information required to calculated energies via S(G). 


            type(sys_t), intent(in) :: sys
            type(excit_t), intent(in) :: connection
            real(p), intent(in) :: real_population
            type(sg_correction_t), intent(inout) :: sg_correction

            real(p) :: omega
            integer :: km(3), kn(3), kp(3), kq(3), G(3)
            integer :: sm, sn, sp, sq
            integer :: idx
            integer :: phase

            omega = sys%lattice%box_length(1)**3
            phase = merge(-1, 1, connection%perm)
    
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
                idx = sg_correction%q_map(G(1), G(2), G(3))
                if (idx /= 0) then 
                    sg_correction%sg(idx) = sg_correction%sg(idx) &
                        + (phase * real_population) / omega
                end if 
            end if

            ! Term 2 
            G = km - kq 
            if (sm == sq .and. sn == sp .and. all(kn - kp == - G)) then
                idx = sg_correction%q_map(G(1), G(2), G(3))
                if (idx /= 0) then
                    sg_correction%sg(idx) = sg_correction%sg(idx) &
                        - (phase * real_population) / omega
                end if 
            end if

        end subroutine calculate_sg_contribution
 
        subroutine interpolate_sg_fibonacci(sg_correction)

            ! Interpolate the calculated S(G) values of each unique magnitude onto the fibonacci
            ! sphere points using inverse distance weighting interpolation. 
            ! See: https://en.wikipedia.org/wiki/Inverse_distance_weighting.

            ! In/Out:
            !   sg_correction: (sg_correction_t) object. Contains information required to calculated energies via S(G). 
 
            type(sg_correction_t), intent(inout) :: sg_correction 
            
            real(p) :: dx, dy, dz, d2, w, sum_w, sum_wS
            real(p) :: tolerance = 1e-8_p
            integer :: i, j, k

            ! Interpolate for each fibonacci sphere 
            do i = 1, size(sg_correction%unique_q_mags)
                
                ! Skip interpolation for zero magnitude sphere and just set S(G) at all points to zero
                if (abs(sg_correction%unique_q_mags(i)) < tolerance) then 
                    sg_correction%sg_fibonacci(:, i) = 0.0_p
                    cycle 
                end if 
                
                ! Loop over fibonacci points 
                do j = 1, size(sg_correction%fibonacci_spheres, 2)

                    sum_w = 0.0_p 
                    sum_wS = 0.0_p

                    ! Sum over S(G) values / q_vectors
                    do k = 1, size(sg_correction%q_vectors, 2)

                        dx = sg_correction%fibonacci_spheres(1, j, i) - sg_correction%q_vectors(1, k)
                        dy = sg_correction%fibonacci_spheres(2, j, i) - sg_correction%q_vectors(2, k)
                        dz = sg_correction%fibonacci_spheres(3, j, i) - sg_correction%q_vectors(3, k)
                        d2 = dx*dx + dy*dy * dz*dz

                        if (d2 < tolerance) then
                            sum_wS = sg_correction%sg(k)
                            sum_w = 1.0_p
                            exit 
                        end if 
                        
                        ! Exponent of weight w gives measure of locality, i.e., how much the interpolated value is 
                        ! influenced by values closer or further away. Higher locality is more influenced by nearer points.
                        ! If locality == 6 then (1 / d)**6 = (1 / d2)**3.  
                        w = 1.0_p / (d2 * d2 * d2) 
                        sum_w = sum_w + w 
                        sum_wS = sum_wS + (w * sg_correction%sg(k))

                    end do 

                    sg_correction%sg_fibonacci(j, i) = sum_wS / sum_w 

                end do 
            end do

        end subroutine interpolate_sg_fibonacci

        subroutine spherical_average_sg(sg_correction)

            ! In/Out:
            !   sg_correction: (sg_correction_t) object. Contains information required to calculated energies via S(G).

            type (sg_correction_t), intent(inout) :: sg_correction

            real(p) :: sum_sg
            integer :: i, j

            do i = 1, size(sg_correction%sg_fibonacci, 2)
                
                sum_sg = 0.0_p

                do j = 1, size(sg_correction%sg_fibonacci, 1)
                    sum_sg = sum_sg + sg_correction%sg_fibonacci(j, i)
                end do 

                sg_correction%spherically_averaged_sg(i) = sum_sg / size(sg_correction%sg_fibonacci, 1)
            end do 
       
        end subroutine spherical_average_sg

        subroutine cubic_spline_sg(sg_correction)

            ! In/Out:
            !   sg_correction: (sg_correction_t) object. Contains information required to calculated energies via S(G).

            type(sg_correction_t), intent(inout) :: sg_correction

            integer :: n_points, i, w
            real(p), allocatable :: q_points(:), sg_points(:)
            real(p), allocatable :: h(:)
            real(p), allocatable :: lower_diag(:), main_diag(:), upper_diag(:) ! Tridiagonal matrix
            real(p), allocatable :: m_vector(:), right_vector(:)

            n_points = size(sg_correction%spherically_averaged_sg)
            allocate(q_points(n_points), sg_points(n_points))
            sg_points = sg_correction%spherically_averaged_sg
            q_points = sg_correction%unique_q_mags

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
                sg_correction%spline_coeffs_a(i) = (m_vector(i + 1) - m_vector(i)) / (6.0_p * h(i))
                sg_correction%spline_coeffs_b(i) = m_vector(i) / 2.0_p
                sg_correction%spline_coeffs_c(i) = (sg_points(i + 1) - sg_points(i)) / h(i) &
                                                    - (h(i) / 6.0_p) * (2 * m_vector(i) + m_vector(i + 1))
                sg_correction%spline_coeffs_d(i) = sg_points(i) 
            end do 

        end subroutine cubic_spline_sg

        pure function evaluate_spline(q, n_knots, sg_correction) result(sg)

            real(p), intent(in) :: q 
            integer, intent(in) :: n_knots 
            type(sg_correction_t), intent(in) :: sg_correction

            integer :: i 
            real(p) :: delta
            real(p) :: sg

            ! Find which interval contains q 
            do i = 1, n_knots - 1 
                if (q < sg_correction%unique_q_mags(i + 1)) exit 
            end do 

            ! Evaluate spline in this interval 
            delta = q - sg_correction%unique_q_mags(i)
            sg = sg_correction%spline_coeffs_a(i) * delta**3 + sg_correction%spline_coeffs_b(i) * delta**2 &
                                                                + sg_correction%spline_coeffs_c(i) * delta &
                                                                + sg_correction%spline_coeffs_d(i)

        end function evaluate_spline

end module sg_correction
