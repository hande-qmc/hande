module sg_correction     

! Module performing the finite size error correction method of Liao et al (doi.org/10.1063/1.4964307). The correlation energy may be
! written in terms of the Coulomb kernel in reciprocal space and the transition structure factor S(G):
!   \sum_{\textbf{G}}v(\textbf{G})S(\textbf{G})
! where S(G) is:
!   S(\textbf{G}) = \sum_{i < j}\sum_{a < b}C_{ij}^{ab}[C_i^a(\textbf{G})C_b^{j*}(\textbf{G}) -
!   C_i^b(\textbf{G})C_a^{j*}(\textbf{G})].
! The C_i^a(\textbf{G}) are the Fourier coefficients of the orbital overlap density, and C_{ij}^{ab} is the amplitude of
! determinant |D_{ij}^{ab}\rangle.  
! Note that this is a different quantity from the static structure factor S(q) implemented for the UEG with DMQMC.
! The set of \textbf{G} vectors are the orbital momentum transfer vectors of the system, and finite size errors in the 
! correlation energy arise from having the minimum value of \textbf{G} != 0. This error can be written as the missing part of the
! above correlation energy sum:
! \Delta E_{FS} = \sum_{\textbf{G} < G_{min}}v(\textbf{G})S(\textbf{G}) \propto int_0^{G_{min}} (v(\textbf{G})S(\textbf{G})).
! We calculate S(G) for all avaliable momentum transfer vectors G in the system, interpolate these values across a number of
! Fibonacci spheres using inverse distance weighting to obtain an isotropic distribution of S(G) in momentum space. The inversion
! distance weighting procedure is taken from: doi.org/10.1145/800186.810616. The spherical average 
! of each sphere is taken to yield S_{avg}(|G|). A cubic spline is constructed with S_{avg}(|G|) as the knots, code for this spline
! and the associated 'Thomas algorithm' are adapted from Numerical Recipes 3rd Edition (9780521880688). The resulting spline 
! can then be integrated in the region [0, Gmin] using composite Simpson's quadrature to obtain a correction for the finite size error.
! This method is implemented for both CCMC and FCIQMC, the only difference being that in CCMC calculate_sg_contribution is called
! within the do_ccmc_accumulation subroutine, i.e., the contribution is totalled for every excitor, whereas in FCIQMC it is simply
! called once per determinant.

use system, only: sys_t, sys_lattice_t
use basis, only: basis_t, basis_fn_t
use excitations, only: excit_t
use checking, only: check_allocate, check_deallocate
use sort, only: insert_sort_unique_vector
use const, only: p

implicit none

! Derived type used for calculation of quantities relating to S(G) and the finite size error correction
type sg_correction_t

! List of unique pairs of momentum transfer vectors between basis_fns(i) and basis_fns(a). 2D array with shape (3, N), where N
! is the number of q_vectors. Note that these are labelled by q for consistency with the codebase, whereas Liao et al refer to them
! by G. 
integer, allocatable :: q_vectors(:, :)

! Stores the lowest and highest values of each entry x, y, z from all of q_vectors. These are used to assign the size of q_map
! correctly in get_q_vectors.
integer :: qx_min, qx_max, qy_min, qy_max, qz_min, qz_max

! Lookup table for quick access to the index of a momentum transfer vector q in q_vectors. 3D array with extent in each dimension
! determined at run time based on qx_min, qx_max, qy_min, qy_max, qz_min, qz_max. These indexes are needed to determine at what
! position in S(G) a contribution from G should go. Quick lookup is required as the calculate_sg_contribution subroutine is called
! in the loops over determinants/clusters.  
integer, allocatable :: q_map(:, :, :)

! Magnitudes of the q_vectors. 1D array with length N. 
real(p), allocatable :: q_mags(:)

! Defines to what extent two magnitudes are the same. If magnitudes differ by less than tolerance they are considered to be
! the same magnitude. 
real(p) :: tolerance = 1e-8_p

! Contains only the unique entries of q_mags. The length of this list (the number of unique magnitudes) determines the number of
! Fibonacci spheres that will be used and the entries determine the magnitudes of each sphere. 
real(p), allocatable :: unique_q_mags(:)

! The number of Fibonacci points generated for each Fibonacci sphere which are used in the interpolation via inverse distance
! weighting subroutine. This number has been decided upon via convergence testing, but more experimentation should be carried out
! here. 
integer :: n_fibonacci_points = 199

! 2D array containing a list of positions on a Fibonacci sphere of unit magnitude. Shape is (3, n_fibonacci_points). This array is
! populated in generate_unit_fibonacci_sphere.
real(p), allocatable :: fibonacci_points_unit(:, :)

! 3D array containing list of positions on multiple Fibonacci spheres. The magnitudes of each sphere are the entries in
! unique_q_mags. Array shape is (3, n_fibonacci_points, n_unique) where n_unique is number of unique magnitudes. This array is
! populated in generate_unique_mag_fibonacci_sphere.
real(p), allocatable :: fibonacci_spheres(:, :, :)

! 1D array which stores the values of the transition structure factor S(G) of length size(q_vectors). As mentioned the contributions
! to this list are calculated in calculate_sg_contribution.
real(p), allocatable :: sg(:)

! Values of the transition structure factor S(G) interpolated onto the Fibonacci sphere positions. 2D array with dimension
! (n_fibonacci_points, n_unique) as we have S(G) values for Fibonacci spheres of all magnitudes, array filled in
! interpolate_sg_fibonacci subroutine.
real(p), allocatable :: sg_fibonacci(:, :)

! List of S(G) values obtained by averaging the Fibonacci spheres. Length is size(sg_fibonacci, 1), that is, the number of Fibonacci
! spheres or equivalently the number of unique_q_mags. These values are calculated in spherical_average_sg.
real(p), allocatable :: spherically_averaged_sg(:)

! Number of points at which to evaluate the cubic spline built from spherically_averaged_sg. This value is again the result of
! convergence testing but further experimentation should be done with this.
integer :: n_spline_points = 100

! Coefficients ai, bi, ci, di for interval i of the cubic spline. The spline interpolates in the intervals between the values of
! spherically_averaged_sg, therefore the length of each list is size(spherically_averaged_sg - 1), or equivalently again the number
! of unique_q_mags - 1. Entries computed in cubic_spline_sg.
real(p), allocatable :: spline_coeffs_a(:)
real(p), allocatable :: spline_coeffs_b(:)
real(p), allocatable :: spline_coeffs_c(:)
real(p), allocatable :: spline_coeffs_d(:)

! Maximum and minimum magnitudes of unique_q_mags. q_min should be zero but this causes division by 0 when calculating the
! correction term. These values are assigned at the end of the report loop in the main FCIQMC and CCMC drivers. This should be moved
! to be in the S(G) initialisation area before any loop as these do not change. 
real(p) :: q_min, q_max 

! List containing uniformly distributed values of q between q_min and q_max. The number of these is set by n_spline_points. Again
! this should be moved out of the report loop.
real(p), allocatable :: q_spline_vals(:)

! List of values of S(G) evaluated at the above q_spline_vals using the spline. The length is n_spline_points.
real (p), allocatable :: sg_spline_vals(:)

end type sg_correction_t

    contains
        
        subroutine initialise_sg_correction(sys, sg_correction)

            ! Wrapper to collect and setup any data needed for performing the finite size correction via S(G). Collects the system momentum
            ! transfers, generates the points on a Fibonacci sphere of unit magnitude, and generates the additional Fibonacci
            ! spheres according to the unique momentum transfer magnitudes.
            
            ! In:
            !   sys: sys_t type which holds information about the system basis set.
            ! In/Out:
            !   sg_correction: sg_correction_t derived type which stores relevant data for the correction scheme.

            type(sys_t), intent(in) :: sys
            type(sg_correction_t), intent(inout) :: sg_correction
 
            call get_q_vectors(sys%basis%nbasis, sys%basis%basis_fns, sg_correction%q_vectors, sg_correction%q_map, &
                               sg_correction%sg, sg_correction%qx_min, sg_correction%qx_max, sg_correction%qy_min, &
                               sg_correction%qy_max, sg_correction%qz_min, sg_correction%qz_max)
            call generate_unit_fibonacci_sphere(sg_correction%n_fibonacci_points, sg_correction%fibonacci_points_unit)
            call generate_unique_mag_fibonacci_sphere(sg_correction%q_mags, sg_correction%q_vectors, &
                                                      sg_correction%unique_q_mags, sg_correction%tolerance, &
                                                      sg_correction%n_fibonacci_points, sg_correction%fibonacci_spheres, &
                                                      sg_correction%fibonacci_points_unit)

        end subroutine initialise_sg_correction

        subroutine get_q_vectors(nbasis, basis_fns, q_vectors, q_map, sg, qx_min, qx_max, qy_min, qy_max, qz_min, qz_max)

            ! Find all unique pairs of momentum transfer vectors (q_vectors) between basis functions 
            ! and allocate sg_array which will contain the corresponding S(G) contribution of each
            ! q_vector. The vectors obtained here are integer vectors with units of 2pi/L. When using these in calculation of
            ! energies they must be converted to physical units.

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
            allocate(q_map(qx_min:qx_max, qy_min:qy_max, qz_min:qz_max), stat=ierr)
            call check_allocate('q_map', size(q_map), ierr)
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
   
            integer :: i, ierr 
            integer :: n_unique

            ! Get magnitudes of q_vectors
            allocate(q_mags(size(q_vectors, 2)), stat=ierr)
            call check_allocate('q_mags', size(q_vectors, 2), ierr)
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
        
        subroutine calculate_sg_contribution(lattice, connection, basis_fns, sg_correction, real_population)

            ! Calculate the contribution to S(G) of the current determinant or cluster.

            ! In:
            !   lattice: sys_lattice_t object which contains information about the system lattice, in particular the side length of
            !   the simulation cell.
            !   connection: Excitation (excit_t) object. Contains information about the excitation between 
            !               the current determinant or cluster and their children.
            !   basis_fns: Object of type (basis_fn_t), gives information about the basis functions.
            !   real_population: Population of the current determinant or cluster.
            ! In/Out:
            !   sg_correction: sg_correction_t derived type which stores relevant data for the correction scheme.
           
            type(sys_lattice_t), intent(in) :: lattice
            type(excit_t), intent(in) :: connection 
            type(basis_fn_t), intent(in) :: basis_fns(:)
            type(sg_correction_t), intent(inout) :: sg_correction
            real(p), intent(in) :: real_population

            real(p) :: omega
            integer :: G(3)
            integer :: idx
            integer :: phase

            omega = lattice%box_length(1) * lattice%box_length(2) * lattice%box_length(3) 
            phase = 1
            if (connection%perm) phase = -phase
            
            ! The transition structure factor is S(\mathbf{G}) = C_{ij}^{ab} [C_i^a(\mathbf{G})C_b^{j*}(\mathbf{G}) 
            ! - C_i^b(\mathbf{G})C_a^{j*}(\mathbf{G})] for a given set of orbitals i, j, k, l.
            ! The orbitals in the UEG are plane waves, and as such the Fourier coefficients of the orbital overlap 
            ! density C_i^a(\textbf{G}) become delta conditions on both the spin and k-vectors of the orbitals.
            ! This gives rise to a determinants contribution to S(G) being calculated as seen below.
           
            !(\delta_{\sigma_m, \sigma_p}\delta_{\mathbf{k}_m -\mathbf{k}_p, \mathbf{G}})
            ! * (\delta_{\sigma_q, \sigma_n}\delta_{\mathbf{k}_n -\mathbf{k}_q, -\mathbf{G}})
            G = basis_fns(connection%to_orb(1))%l - basis_fns(connection%from_orb(1))%l
            if (basis_fns(connection%to_orb(1))%ms == basis_fns(connection%from_orb(1))%ms .and. &
                basis_fns(connection%to_orb(2))%ms == basis_fns(connection%from_orb(2))%ms .and. &
                all(basis_fns(connection%to_orb(2))%l - basis_fns(connection%from_orb(2))%l == -G)) then
                idx = sg_correction%q_map(G(1), G(2), G(3))
                if (idx /= 0) then
                    sg_correction%sg(idx) = sg_correction%sg(idx) + (phase * real_population) / omega
                end if 
            end if
            
            !- (\delta_{\sigma_m, \sigma_q}\delta_{\mathbf{k}_m -\mathbf{k}_q, \mathbf{G}})
            ! * (\delta_{\sigma_p, \sigma_n}\delta_{\mathbf{k}_n -\mathbf{k}_p, -\mathbf{G}})
            G = basis_fns(connection%to_orb(1))%l  - basis_fns(connection%from_orb(2))%l 
            if (basis_fns(connection%to_orb(1))%ms == basis_fns(connection%from_orb(2))%ms .and. &
                basis_fns(connection%to_orb(2))%ms == basis_fns(connection%from_orb(1))%ms .and. &
                all(basis_fns(connection%to_orb(2))%l - basis_fns(connection%from_orb(1))%l == -G)) then
                idx = sg_correction%q_map(G(1), G(2), G(3))
                if (idx /= 0) then
                    sg_correction%sg(idx) = sg_correction%sg(idx) - (phase * real_population) / omega
                end if 
            end if

        end subroutine calculate_sg_contribution
        
        subroutine calculate_sg_correction(sg_correction, lattice, energy_sg, corrected_energy_sg, ireport, nreport)

            ! Wrapper to process the gathered S(G) data at the end of each report loop. Calls the interpolation and spline
            ! subroutines, as well as the routines in sg_correction_energies.
            
            ! In:
            !   lattice: sys_lattice_t object which contains information about the system lattice, in particular the side length of
            !   the simulation cell.
            !   ireport: Current report loop number. 
            !   nreport: Total number of report loops.
            ! In/Out:
            !   sg_correction: sg_correction_t derived type which stores relevant data for the correction scheme.
            !   energy_sg: Estimator of the uncorrected energy. Elements of array are equivalent to different walker 
            !   spaces such as with complex systems or when using replica tricks. Since UEG populations are always real 
            !   however, nspaces should always be 1 when using sg_correction here. 
            !   corrected_energy_sg: Estimator of the corrected energy. Array elements same as above.

            type(sg_correction_t), intent(inout) :: sg_correction 
            type(sys_lattice_t), intent(in) :: lattice
            real(p), intent(inout) :: energy_sg(:), corrected_energy_sg(:)
            integer :: ireport, nreport, unit_sg_out
            character(len = 6) :: filename

            integer :: i

            corrected_energy_sg(:) = 0.0_p

            call interpolate_sg_fibonacci(sg_correction%n_fibonacci_points, sg_correction%fibonacci_spheres, &
                                          sg_correction%unique_q_mags, sg_correction%sg_fibonacci, & 
                                          sg_correction%q_vectors, sg_correction%sg)
            call spherical_average_sg(sg_correction%sg_fibonacci, sg_correction%spherically_averaged_sg, &
                                      sg_correction%unique_q_mags)
            call cubic_spline_sg(sg_correction%spherically_averaged_sg, sg_correction%unique_q_mags, &
                                 sg_correction%spline_coeffs_a, sg_correction%spline_coeffs_b, &
                                 sg_correction%spline_coeffs_c, sg_correction%spline_coeffs_d)

            sg_correction%q_min = sg_correction%unique_q_mags(1) + 1.0e-6_p ! Non-zero q_min avoids divergence in long range integral
            sg_correction%q_max = sg_correction%unique_q_mags(size(sg_correction%unique_q_mags))
            if (.not. allocated(sg_correction%q_spline_vals)) then 
                allocate(sg_correction%q_spline_vals(sg_correction%n_spline_points))
                allocate(sg_correction%sg_spline_vals(sg_correction%n_spline_points))
            end if

            do i = 1, sg_correction%n_spline_points
                ! Uniform distribution of q magnitudes
                sg_correction%q_spline_vals(i) = sg_correction%q_min + (i - 1) &
                                                  * (sg_correction%q_max - sg_correction%q_min) &
                                                  / (sg_correction%n_spline_points - 1) 
                sg_correction%sg_spline_vals(i) = evaluate_spline(sg_correction%q_spline_vals(i), & 
                                                                     size(sg_correction%unique_q_mags), & 
                                                                     sg_correction)
            end do

            call calculate_corrected_sg_energy(corrected_energy_sg, sg_correction%q_vectors, &
                                               lattice%box_length, sg_correction%unique_q_mags, &
                                               sg_correction%sg, sg_correction%q_spline_vals, &
                                               sg_correction%sg_spline_vals, sg_correction%n_spline_points)

            sg_correction%sg = 0.0_p

        end subroutine calculate_sg_correction
        
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
            
            real(p) :: dx, dy, dz, d2, d2_inverse, w, sum_w, sum_wS
            real(p) :: fx, fy, fz
            real(p), allocatable :: qx(:), qy(:), qz(:)
            real(p) :: tolerance = 1e-8_p
            integer :: i, j, k, ierr, n_unqiue_q_mags, n_q_vectors
            logical, allocatable :: is_zero_mag_sphere(:)

            n_unqiue_q_mags = size(unique_q_mags)
            n_q_vectors = size(q_vectors, 2)

            allocate(sg_fibonacci(n_fibonacci_points, n_unqiue_q_mags), stat=ierr)
            call check_allocate('sg_fibonacci', n_unqiue_q_mags, ierr)

            allocate(qx(n_q_vectors), stat=ierr)
            call check_allocate('qx', n_q_vectors, ierr)
            allocate(qy(n_q_vectors), stat=ierr)
            call check_allocate('qy', n_q_vectors, ierr)
            allocate(qz(n_q_vectors), stat=ierr)
            call check_allocate('qz', n_q_vectors, ierr)
            allocate(is_zero_mag_sphere(n_unqiue_q_mags), stat=ierr)
            call check_allocate('is_zero_mag_sphere', n_unqiue_q_mags, ierr)
            
            qx = q_vectors(1, :)
            qy = q_vectors(2, :)
            qz = q_vectors(3, :)
            is_zero_mag_sphere = abs(unique_q_mags) < tolerance

            ! Interpolate for each fibonacci sphere 
            do i = 1, n_unqiue_q_mags
                
                ! Skip interpolation for zero magnitude sphere and just set S(G) at all points to zero
                if (is_zero_mag_sphere(i)) then 
                    sg_fibonacci(:, i) = 0.0_p
                    cycle 
                end if 
                
                ! Loop over fibonacci points 
                do j = 1, n_fibonacci_points

                    sum_w = 0.0_p 
                    sum_wS = 0.0_p
                    
                    fx = fibonacci_spheres(1, j, i)
                    fy = fibonacci_spheres(2, j, i)
                    fz = fibonacci_spheres(3, j, i) 

                    ! Sum over S(G) values / q_vectors
                    do k = 1, n_q_vectors

                        dx = fx - qx(k)
                        dy = fy - qy(k)
                        dz = fz - qz(k)
                        d2 = dx*dx + dy*dy + dz*dz

                        !if (d2 < tolerance) then
                        !    sum_wS = sg(k)
                        !    sum_w = 1.0_p
                        !    exit 
                        !end if 
                        d2 = d2 + tolerance
                        
                        ! Exponent of weight w gives measure of locality, i.e., how much the interpolated value is 
                        ! influenced by values closer or further away. Higher locality is more influenced by nearer points.
                        ! If locality == 6 then (1 / d)**6 = (1 / d2)**3. 
                        d2_inverse = 1.0_p / d2 
                        w = d2_inverse * d2_inverse * d2_inverse
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
            ! but this is may need to be modified for larger systems. The Thomas algorithm  is used to find the coefficients 
            ! of the spline in each region between knots. Code here is adapted from Numerical Recipes 3rd Edition (9780521880688).

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

            integer :: n_spline_points, i, ierr
            real(p) :: w
            real(p), allocatable :: h(:)
            real(p), allocatable :: lower_diag(:), main_diag(:), upper_diag(:) 
            real(p), allocatable :: m_vector(:), right_vector(:)

            allocate(spline_coeffs_a(size(unique_q_mags) - 1), stat=ierr)
            call check_allocate('spline_coeffs_a', size(spline_coeffs_a), ierr)
            allocate(spline_coeffs_b(size(unique_q_mags) - 1), stat=ierr)
            call check_allocate('spline_coeffs_b', size(spline_coeffs_b), ierr)
            allocate(spline_coeffs_c(size(unique_q_mags) - 1), stat=ierr)
            call check_allocate('spline_coeffs_c', size(spline_coeffs_c), ierr)
            allocate(spline_coeffs_d(size(unique_q_mags) - 1), stat=ierr)
            call check_allocate('spline_coeffs_d', size(spline_coeffs_d), ierr)

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
                main_diag(i) = 2.0_p * (h(i-1) + h(i)) 
                upper_diag(i) = h(i) ! upper_diag(1) = 0 as there is no matrix element (0, 1)
                right_vector(i) = 6.0_p * ((spherically_averaged_sg(i + 1) - spherically_averaged_sg(i)) / h(i) &
                                       - (spherically_averaged_sg(i) - spherically_averaged_sg(i - 1)) / h(i - 1))  
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
        
        subroutine calculate_corrected_sg_energy(corrected_energy_sg, q_vectors, box_lengths, unique_q_mags, sg, q_spline_vals, &
                                                 sg_spline_vals, n_spline_points)
                
                ! Calculate the long range correction to the energy via S(G)
                
                ! In: 
                !   q_vectors: List of unique momentum transfers within the system. 
                !   box_lengths: Side lengths of the simulation cell. 
                !   unique_q_mags: Unique magnitudes of the q vectors. 
                !   sg: Array of S(G) values with index corresponding to the q vector which will give that value.
                !   q_spline_vals: Uniform distribution of momentum transfer magnitudes in the range q_min to q_max.  
                !   sg_spline_vals: Values of S(G) evaluated at the q_spline_vals points using the cubic spline. 
                !   n_spline_points: Total number of points which were evaluated on the spline. 
                ! In/Out: 
                !   corrected_energy_sg: Correlation energy of the system with the long range correction added. 
                
                real(p), intent(inout) :: corrected_energy_sg(:) 
                integer, intent(in) :: q_vectors(:, :) 
                real(p), intent(in) :: box_lengths(3) 
                real(p), intent(in) :: unique_q_mags(:)
                real(p), intent(in) :: sg(:)
                real(p), intent(in) :: q_spline_vals(:), sg_spline_vals(:)
                integer, intent(in) :: n_spline_points

                real(p) :: mag2, window, Gc, h
                real(p) :: Gx, Gy, Gz, omega, prefactor, L, real_unit_conversion
                integer :: i
                real(p) :: contribution, term1, term2, term3, term4

                Gc = unique_q_mags(2) ! Largest non-zero magnitude
                omega = box_lengths(1) * box_lengths(2) * box_lengths(3)
                prefactor = (2.0_p * omega) / (acos(-1.0_p))
                L = omega**(1.0_p / 3.0_p)
                real_unit_conversion = (2 * acos(-1.0_p)) / L 
                
                ! Calculate short range part with non-interpolated q_vectors and S(G)
                corrected_energy_sg = 0.0_p
                do i = 1, size(q_vectors, 2)

                    ! Momentum transfers are in units of 2\pi/L so we must convert to real units for calculation of the energy,
                    ! hence the prefactor ((2 * acos(-1.0_p)) / L)**2.
                    Gx = ((2 * acos(-1.0_p)) / box_lengths(1)) * q_vectors(1, i)
                    Gy = ((2 * acos(-1.0_p)) / box_lengths(2)) * q_vectors(2, i)
                    Gz = ((2 * acos(-1.0_p)) / box_lengths(3)) * q_vectors(3, i)
                    mag2 = Gx * Gx + Gy * Gy + Gz * Gz
                    window = window_function(sqrt(mag2), Gc)

                    ! Skip |G| = 0. Short range contriubtion to corrected energy
                    if (mag2 /= 0.0_p) then 
                        corrected_energy_sg = corrected_energy_sg + (4 * acos(-1.0_p) / mag2) * (1 - window) * sg(i) 
                    end if

                end do

                ! Long range contriubtion (correction) to energy 
                ! Composite Simpson quadrature to evaluate long range correction 
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

                corrected_energy_sg = corrected_energy_sg + prefactor * (h / 3.0_p) & 
                                    * real_unit_conversion * (term1 + 4 * term2 + 2 * term3 + term4)

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

    subroutine finalise_sg_correction(sg)

        ! Deallocate all S(G) correction arrays 

        ! In/Out: 
        !   sg: sg_correction_t type which stores all S(G) correction associated arrays and variables

        type(sg_correction_t), intent(inout) :: sg

        integer :: ierr
        
        if (allocated(sg%q_vectors)) then
            deallocate(sg%q_vectors, stat=ierr)
            call check_deallocate('sg%q_vectors', ierr)
        end if
        if (allocated(sg%q_map)) then
            deallocate(sg%q_map, stat=ierr)
            call check_deallocate('sg%q_map', ierr)
        end if
        if (allocated(sg%q_mags)) then
            deallocate(sg%q_mags, stat=ierr)
            call check_deallocate('sg%q_mags', ierr)
        end if
        if (allocated(sg%unique_q_mags)) then
            deallocate(sg%unique_q_mags, stat=ierr)
            call check_deallocate('sg%unique_q_mags', ierr)
        end if
        if (allocated(sg%fibonacci_points_unit)) then
            deallocate(sg%fibonacci_points_unit, stat=ierr)
            call check_deallocate('sg%fibonacci_points_unit', ierr)
        end if
        if (allocated(sg%fibonacci_spheres)) then
            deallocate(sg%fibonacci_spheres, stat=ierr)
            call check_deallocate('sg%fibonacci_spheres', ierr)
        end if
        if (allocated(sg%sg)) then
            deallocate(sg%sg, stat=ierr)
            call check_deallocate('sg%sg', ierr)
        end if
        if (allocated(sg%sg_fibonacci)) then
            deallocate(sg%sg_fibonacci, stat=ierr)
            call check_deallocate('sg%sg_fibonacci', ierr)
        end if
        if (allocated(sg%spherically_averaged_sg)) then
            deallocate(sg%spherically_averaged_sg, stat=ierr)
            call check_deallocate('sg%spherically_averaged_sg', ierr)
        end if
        if (allocated(sg%spline_coeffs_a)) then
            deallocate(sg%spline_coeffs_a, stat=ierr)
            call check_deallocate('sg%spline_coeffs_a', ierr)
        end if
        if (allocated(sg%spline_coeffs_b)) then
            deallocate(sg%spline_coeffs_b, stat=ierr)
            call check_deallocate('sg%spline_coeffs_b', ierr)
        end if
        if (allocated(sg%spline_coeffs_c)) then
            deallocate(sg%spline_coeffs_c, stat=ierr)
            call check_deallocate('sg%spline_coeffs_c', ierr)
        end if
        if (allocated(sg%spline_coeffs_d)) then
            deallocate(sg%spline_coeffs_d, stat=ierr)
            call check_deallocate('sg%spline_coeffs_d', ierr)
        end if
        if (allocated(sg%q_spline_vals)) then
            deallocate(sg%q_spline_vals, stat=ierr) 
            call check_deallocate('sg%q_spline_vals', ierr)
        end if
        if (allocated(sg%sg_spline_vals)) then
            deallocate(sg%sg_spline_vals, stat=ierr) 
            call check_deallocate('sg%sg_spline_vals', ierr)
        end if

    end subroutine finalise_sg_correction


end module sg_correction
