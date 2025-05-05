module sg_correction

    ! Module for calculating quantities relating to S(G) required to calculate finite size correction to the energy.

    use system, only: sys_t
    use basis, only: basis_fn_t
    use excitations, only: excit_t
    use checking, only: check_allocate
    use sort, only: insert_sort_unique_vector
    use const, only: p

    implicit none

    type sg_correction_t

        integer, allocatable :: q_vectors(:, :)
        
        integer :: qx_min, qx_max, qy_min, qy_max, qz_min, qz_max
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

        subroutine get_q_vectors(nbasis, basis_fns, q_vectors, q_map, sg, qx_min, qx_max, qy_min, qy_max, qz_min, qz_max)

            ! Find all unique pairs of momentum transfer vectors (q_vectors) between basis functions 
            ! and allocate sg_array which will contain the corresponding S(G) contribution of each
            ! q_vector.

            ! In:
            !   nbasis: Number of basis functions in the system. 
            !   basis_fns: Object of type (basis_fn_t), gives information about the basis functions. 
            ! Out:
            !   q_vectors: List of unique momentum transfers within the system. 
            !   q_map: Lookup map for easy access to the indicies corresponding to a certain q_vector.
            !   sg: Array of S(G) values with index corresponding to the q vector which will give that value.
            !   qx/y/x_min/max: Min and max values of the q vectors in each direction.
            
            integer, intent(in) :: nbasis 
            type(basis_fn_t), intent(in) :: basis_fns(:)
            integer, allocatable, intent(out) :: q_vectors(:, :)
            integer, allocatable, intent(out) :: q_map(:, :, :)
            real(p), allocatable, intent(out) :: sg(:)
            integer, intent(out) :: qx_min, qx_max, qy_min, qy_max, qz_min, qz_max

            integer :: i, ierr, ki, ka, q(3), n_q
            

            ! Get all unique q vectors
            do ki = 1, nbasis
                do ka = 1, nbasis 
                    q = basis_fns(ki)%l - basis_fns(ka)%l
                    call insert_sort_unique_vector(q_vectors, q)
                end do 
            end do 

            ! Allocate corresponding S(G) array 
            allocate(sg(size(q_vectors, 2)), stat=ierr)
            call check_allocate('sg', size(q_vectors, 2), ierr)
            sg = 0.0_p

            ! Construct lookup map for the index of the q_vectors 
            qx_min = minval(q_vectors(1, :))
            qx_max = maxval(q_vectors(1, :))
            qy_min = minval(q_vectors(2, :))
            qy_max = maxval(q_vectors(2, :))
            qz_min = minval(q_vectors(3, :))
            qz_max = maxval(q_vectors(3, :))
            allocate(q_map(qx_min:qx_max, qy_min:qy_max, qz_min:qz_max))
            q_map = 0 ! If q vector not in the list of q_vectors the map should return 0 

            do i = 1, size(q_vectors, 2)
                q = q_vectors(:, i)
                q_map(q(1), q(2), q(3)) = i
            end do 

        end subroutine get_q_vectors

        subroutine generate_unit_fibonacci_sphere(n_points, points)

            ! Generate a Fibonacci sphere with n_points points of unit magnitude. 

            ! In:
            !   n_points: Number of points on the fibonacci sphere. 
            ! Out:
            !   points: Array of 3D vectors containing the fibonacci sphere points.
                
            integer, intent(in) :: n_points
            real(p), allocatable, intent(out) :: points(:, :) 

            integer :: n, k, idx
            real(p) :: phi, lat, lon

            allocate(points(3, n_points))

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

        subroutine generate_unique_mag_fibonacci_sphere(q_mags, q_vectors, unique_q_mags, tolerance, n_fibonacci_points, &
                                                        fibonacci_spheres, fibonacci_points_unit)

            ! Calculates the magnitudes of the stored q_vectors and creates a list 
            ! containing only unique magnitudes. 

            ! In:
            !   q_vectors: List of unique momentum transfers within the system. 
            !   tolerance: Determines to what extent a magnitude of q is or isnt 'unique'.
            !   n_fibonacci_points: Number of points on each Fibonacci sphere .
            !   fibonacci_points_unit: Location of points on the unit magnitude Fibonacci spheres.
            !   
            ! Out:
            !   q_mags: Magnitudes of the physical momentum transfer vectors (q_vectors).
            !   unique_q_mags: Unique magnitudes of the physical (i.e., calculated from the basis functions) momentum tranfer 
            !                  vectors.
            !   fibonacci_spheres: Fibonacci spheres scaled to each unique magnitude of q.

            integer, intent(in) :: q_vectors(:, :)
            real(p), intent(in) :: tolerance 
            integer, intent(in) :: n_fibonacci_points 
            real(p), intent(in) :: fibonacci_points_unit(:, :)
            real(p), allocatable, intent(out) :: q_mags(:)
            real(p), allocatable, intent(out) :: unique_q_mags(:)
            real(p), allocatable, intent(out) :: fibonacci_spheres(:, :, :)
   
            integer :: i 
            integer :: n_unique

            ! Get magnitudes of q_vectors
            allocate(q_mags(size(q_vectors, 2)))
            do i = 1, size(q_vectors, 2)
                q_mags(i) = sqrt(real(q_vectors(1, i)**2 + &
                                      q_vectors(2, i)**2 + &
                                      q_vectors(3, i)**2, p))
            end do

            ! Find unique magnitudes only 
            allocate(unique_q_mags(size(q_vectors, 2)))
            n_unique = 0 
            do i = 1, size(q_vectors, 2)
                if (i == 1 .or. abs(q_mags(i) - q_mags(i - 1)) > tolerance) then 
                    n_unique = n_unique + 1 
                    unique_q_mags(n_unique) = q_mags(i)
                end if
            end do 
            
            ! Remove excess size
            unique_q_mags = unique_q_mags(1:n_unique)

            ! Generate 1 fibonacci sphere per unique magnitude
            allocate(fibonacci_spheres(3, n_fibonacci_points, n_unique))
            do i = 1, n_unique 
                fibonacci_spheres(:, :, i) = unique_q_mags(i) * fibonacci_points_unit
            end do

        end subroutine generate_unique_mag_fibonacci_sphere

        subroutine calculate_sg_contribution(box_length, connection, basis_fns, q_map, sg, real_population, q_vectors, qx_min, &
                                             qx_max, qy_min, qy_max, qz_min, qz_max)

            ! Calculate the contribution to S(G) of the current determinant or cluster.

            ! In:
            !   box_length: Side length of the simulation cell.
            !   connection: Excitation (excit_t) object. Contains information about the excitation between 
            !               the current determinant or cluster and their children.
            !   basis_fns: Object of type (basis_fn_t), gives information about the basis functions.
            !   q_map: Lookup map of the q_vectors for quick access to their indicies. 
            !   real_population : Population of the current determinant or cluster.
            !   q_vectors: Physical momentum transfer vectors of the system.
            !   qx/y/x_min/max: Min and max values of the q vectors in each direction.
            ! In/Out:
            !   sg: Array of S(G) values to which the current determinant or cluster may or may not contribute.
            
            real(p), intent(in) :: box_length 
            type(excit_t), intent(in) :: connection 
            type(basis_fn_t), intent(in) :: basis_fns(:)
            integer, intent(in) :: q_vectors(:, :)
            real(p), intent(inout) :: sg(:)
            real(p), intent(in) :: real_population
            integer, intent(in) :: qx_min, qx_max, qy_min, qy_max, qz_min, qz_max
            integer, intent(in) :: q_map(qx_min:qx_max, qy_min:qy_max, qz_min:qz_max)

            real(p) :: omega
            integer :: km(3), kn(3), kp(3), kq(3), G(3)
            integer :: sm, sn, sp, sq
            integer :: idx
            integer :: phase

            omega = box_length**3
            phase = merge(-1, 1, connection%perm)
    
            ! to_orb contains occupied orbitals, from_orb contains virtual orbitals
            km = basis_fns(connection%to_orb(1))%l 
            kn = basis_fns(connection%to_orb(2))%l
            kp = basis_fns(connection%from_orb(1))%l 
            kq = basis_fns(connection%from_orb(2))%l
            sm = basis_fns(connection%to_orb(1))%ms
            sn = basis_fns(connection%to_orb(2))%ms 
            sp = basis_fns(connection%from_orb(1))%ms 
            sq = basis_fns(connection%from_orb(2))%ms
 
            ! Term 1
            G = km - kp
            if (sm == sp .and. sn == sq .and. all(kn - kq == -G)) then
                idx = q_map(G(1), G(2), G(3))
                if (idx /= 0) then
                    sg(idx) = sg(idx) + (phase * real_population) / omega
                end if 
            end if

            ! Term 2 
            G = km - kq 
            if (sm == sq .and. sn == sp .and. all(kn - kp == - G)) then
                idx = q_map(G(1), G(2), G(3))
                if (idx /= 0) then
                    sg(idx) = sg(idx) - (phase * real_population) / omega
                end if 
            end if

        end subroutine calculate_sg_contribution
 
        subroutine interpolate_sg_fibonacci(n_fibonacci_points, fibonacci_spheres, unique_q_mags, sg_fibonacci, q_vectors, sg)

            ! Interpolate the calculated S(G) values of each unique magnitude onto the fibonacci
            ! sphere points using inverse distance weighting interpolation.

            ! In:
            !   n_fibonacci_points: Number of points on each Fibonacci spehere.
            !   fibonacci_spheres: Points on each Fibonacci sphere of unique magnitude.
            !   unique_q_mags: Unique magnitudes of the physical (i.e., calculated from the basis functions) momentum tranfer 
            !                  vectors.
            !   q_vectors: The physical momentum tranfer vectors.
            !   sg: Corresponding values of S(G)
            ! Out:
            !   sg_fibonacci: S(G) values interpolated from the physical vector points to those on the Fibonacci spheres.
            
            integer, intent(in) :: n_fibonacci_points
            real(p), intent(in) :: fibonacci_spheres(:, :, :)
            real(p), intent(in) :: unique_q_mags(:) 
            integer, intent(in) :: q_vectors(:, :)
            real(p), intent(in) :: sg(:)
            real(p), allocatable, intent(out) :: sg_fibonacci(:, :)
            
            real(p) :: dx, dy, dz, d2, w, sum_w, sum_wS
            real(p) :: tolerance = 1e-8_p
            integer :: i, j, k

            allocate(sg_fibonacci(n_fibonacci_points, size(unique_q_mags)))

            ! Interpolate for each fibonacci sphere 
            do i = 1, size(unique_q_mags)
                
                ! Skip interpolation for zero magnitude sphere and just set S(G) at all points to zero
                if (abs(unique_q_mags(i)) < tolerance) then 
                    sg_fibonacci(:, i) = 0.0_p
                    cycle 
                end if 
                
                ! Loop over fibonacci points 
                do j = 1, size(fibonacci_spheres, 2)

                    sum_w = 0.0_p 
                    sum_wS = 0.0_p

                    ! Sum over S(G) values / q_vectors
                    do k = 1, size(q_vectors, 2)

                        dx = fibonacci_spheres(1, j, i) - q_vectors(1, k)
                        dy = fibonacci_spheres(2, j, i) - q_vectors(2, k)
                        dz = fibonacci_spheres(3, j, i) - q_vectors(3, k)
                        d2 = dx*dx + dy*dy + dz*dz

                        if (d2 < tolerance) then
                            sum_wS = sg(k)
                            sum_w = 1.0_p
                            exit 
                        end if 
                        
                        ! Exponent of weight w gives measure of locality, i.e., how much the interpolated value is 
                        ! influenced by values closer or further away. Higher locality is more influenced by nearer points.
                        ! If locality == 6 then (1 / d)**6 = (1 / d2)**3.  
                        w = 1.0_p / (d2 * d2 * d2) 
                        sum_w = sum_w + w 
                        sum_wS = sum_wS + (w * sg(k))

                    end do 

                    sg_fibonacci(j, i) = sum_wS / sum_w 

                end do 
            end do

        end subroutine interpolate_sg_fibonacci

        subroutine spherical_average_sg(sg_fibonacci, spherically_averaged_sg, unique_q_mags)

            ! Takes a simple spherical average of the interpolated S(G) values across each Fibonacci sphere of unique magnitude.

            ! In:
            !   sg_fibonacci: S(G) values interpolated from the physical vector points to those on the Fibonacci spheres. 
            !   unique_q_mags: Unique magnitudes of the physical (i.e., calculated from the basis functions) momentum tranfer 
            !                  vectors. 

            ! Out:
            !   spherically_averaged_sg: Spherical average of the interpolated S(G) valaues.
            
            real(p), intent(in) :: sg_fibonacci(:, :)
            real(p), allocatable, intent(out) :: spherically_averaged_sg(:)
            real(p), intent(in) :: unique_q_mags(:)

            real(p) :: sum_sg
            integer :: i, j

            allocate(spherically_averaged_sg(size(unique_q_mags)))

            do i = 1, size(sg_fibonacci, 2)
                
                sum_sg = 0.0_p

                do j = 1, size(sg_fibonacci, 1)
                    sum_sg = sum_sg + sg_fibonacci(j, i)
                end do 

                spherically_averaged_sg(i) = sum_sg / size(sg_fibonacci, 1)
            end do 
       
        end subroutine spherical_average_sg

        subroutine cubic_spline_sg(spherically_averaged_sg, unique_q_mags, spline_coeffs_a, spline_coeffs_b, spline_coeffs_c, &
                                   spline_coeffs_d)

            ! Constructs a continous cubic spline through all avaliable data points, those being the unique q magnitudes 
            ! and the corresponding spherically averaged values of S(G). Spline currently uses natural boundary conditions 
            ! but this is potentially open to experimentation. The Thomas algorithm is used to find the coefficients of the spline 
            ! in each region between knots.

            ! In:
            !   spherically_averaged_sg: Spherical average of the interpolated S(G) valaues.
            !   unique_q_mags: Unique magnitudes of the physical (i.e., calculated from the basis functions) momentum tranfer 
            !                  vectors. 
            ! Out:
            !   spline_coeffs_a/b/c/d: Vectors of coefficients for each region (section between points) of the spline. A refers to 
            !                          the cubic term, B the quadratic, etc.
            
            real(p), intent(in) :: spherically_averaged_sg(:)
            real(p), intent(in) :: unique_q_mags(:)

            real(p), allocatable, intent(out) :: spline_coeffs_a(:)
            real(p), allocatable, intent(out) :: spline_coeffs_b(:)
            real(p), allocatable, intent(out) :: spline_coeffs_c(:)
            real(p), allocatable, intent(out) :: spline_coeffs_d(:)

            integer :: n_spline_points, i, w
            real(p), allocatable :: h(:)
            real(p), allocatable :: lower_diag(:), main_diag(:), upper_diag(:) 
            real(p), allocatable :: m_vector(:), right_vector(:)

            allocate(spline_coeffs_a(size(unique_q_mags) - 1))
            allocate(spline_coeffs_b(size(unique_q_mags) - 1))
            allocate(spline_coeffs_c(size(unique_q_mags) - 1))
            allocate(spline_coeffs_d(size(unique_q_mags) - 1))

            n_spline_points = size(spherically_averaged_sg)
            allocate(h(n_spline_points - 1), lower_diag(n_spline_points), main_diag(n_spline_points), &
                     upper_diag(n_spline_points), right_vector(n_spline_points), m_vector(n_spline_points))
            
            do i = 1, n_spline_points - 1 
                h(i) = unique_q_mags(i + 1) - unique_q_mags(i)
            end do 
            
            main_diag = 0.0_p
            lower_diag = 0.0_p
            upper_diag = 0.0_p
            right_vector = 0.0_p

            ! Natural spline boundary conditions. S''(G) = 0.  
            main_diag(1) = 1.0_p
            main_diag(n_spline_points) = 1.0_p
            right_vector(1) = 0.0_p
            right_vector(n_spline_points) = 0.0_p
            
            ! Build matrix
            do i = 2, n_spline_points - 1 
                lower_diag(i) = h(i - 1) ! lower_diag(1) = 0 as there is no matrix element (1, 0)
                main_diag(i) = 2.0_p * (h(i-1) + h(i)) ! Boundary conditions give main_diag(1) and main_diag(n_spline_points) = 1
                upper_diag(i) = h(i) ! upper_diag(1) = 0 as there is no matrix element (0, 1)
                right_vector(i) = 6.0_p * ((spherically_averaged_sg(i + 1) - spherically_averaged_sg(i)) / h(i) &
                                       - (spherically_averaged_sg(i) - spherically_averaged_sg(i - 1)) / h(i - 1)) ! Boundary conditions give right_vector(1) 
                                                                                                        ! and right_vector(n_spline_points) = 0
            end do

            ! Thomas algorithm
            do i = 2, n_spline_points 
                w = lower_diag(i) / main_diag(i - 1)
                main_diag(i) = main_diag(i) - w * upper_diag(i - 1)
                right_vector(i) = right_vector(i) - w * right_vector(i - 1)
            end do 
            m_vector(n_spline_points) = right_vector(n_spline_points) / main_diag(n_spline_points)
            do i = n_spline_points - 1, 1, -1 ! Reverse loop
                m_vector(i) = (right_vector(i) - upper_diag(i) * m_vector(i + 1)) / main_diag(i)
            end do 

            ! Compute coefficients 
            do i = 1, n_spline_points - 1
                spline_coeffs_a(i) = (m_vector(i + 1) - m_vector(i)) / (6.0_p * h(i))
                spline_coeffs_b(i) = m_vector(i) / 2.0_p
                spline_coeffs_c(i) = (spherically_averaged_sg(i + 1) - spherically_averaged_sg(i)) / h(i) &
                                                    - (h(i) / 6.0_p) * (2 * m_vector(i) + m_vector(i + 1))
                spline_coeffs_d(i) = spherically_averaged_sg(i) 
            end do 

        end subroutine cubic_spline_sg

        pure function evaluate_spline(q, n_knots, sg_correction) result(sg)

            ! Takes in a q magnitude and returns the value of S(G) at this point according to the cubic spline 

            ! In:
            !   sg_correction: (sg_correction_t) object. Contains information required to calculated energies via S(G).
 
            real(p), intent(in) :: q 
            integer, intent(in) :: n_knots 
            type(sg_correction_t), intent(in) :: sg_correction

            integer :: i 
            real(p) :: delta
            real(p) :: sg

            ! Find which interval/region between knots contains the given q magnitude 
            do i = 1, n_knots - 1 
                if (q <= sg_correction%unique_q_mags(i + 1)) exit 
            end do 

            ! Evaluate value of S(G) for the given q using the spline in this interval 
            delta = q - sg_correction%unique_q_mags(i)
            sg = sg_correction%spline_coeffs_a(i) * delta**3 + sg_correction%spline_coeffs_b(i) * delta**2 &
                                                                + sg_correction%spline_coeffs_c(i) * delta &
                                                                + sg_correction%spline_coeffs_d(i)

        end function evaluate_spline

end module sg_correction
