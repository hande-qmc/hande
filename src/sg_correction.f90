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

        subroutine generate_fibonacci_sphere(n_points, points)

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
            
        end subroutine generate_fibonacci_sphere

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
            real(p), intent(in) :: fibonacci_points_unit(:, :)
            real(p), intent(in) :: tolerance
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
            do i = 1, size(qs%estimators(1)%q_vectors)
                if (qs%estimators(1)%q_mags(i) > tolerance) then ! Do not add zero magnitude to list 
                    if (i == 1 .or. abs(qs%estimators(1)%q_mags(i) - qs%estimators(1)%q_mags(i - 1)) > tolerance) then 
                        n_unique = n_unique + 1 
                        qs%estimators(1)%unique_q_mags(n_unique) = qs%estimators(1)%q_mags(i)
                    end if
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

        subroutine calculate_sg_energy(qs, qmc_in, sys)

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

                mag2 = ((2 * acos(-1.0)) / sys%lattice%box_length(1))**2 * real(qs%estimators(1)%q_vectors(1, i)**2 & 
                        + qs%estimators(1)%q_vectors(2, i)**2 &
                        + qs%estimators(1)%q_vectors(3, i)**2)

                ! Skip |G| = 0
                if (mag2 /= 0.0_p) then 
                    qs%estimators(1)%energy_sg = qs%estimators(1)%energy_sg & 
                        + (4 * acos(-1.0) / mag2) * qs%estimators(1)%sg(i)
                end if
            end do 

        end subroutine calculate_sg_energy

        subroutine dump_sg_data(qs)

            ! Dump information about q_vectors and their contribution S(G) to file SG.

            ! In:
            !   qs: QMC state (qmc_state_t) object. Contains various information about the QMC state.

            type(qmc_state_t), intent(in) :: qs 

            integer :: sg_unit, ierr 
            integer :: i

            open(newunit=sg_unit, file="SG", status="replace", action="write", iostat=ierr)

            do i = 1, size(qs%estimators(1)%sg, 1)
                write(sg_unit, *) &
                    qs%estimators(1)%q_vectors(1, i), &
                    qs%estimators(1)%q_vectors(2, i), &
                    qs%estimators(1)%q_vectors(3, i), &
                    qs%estimators(1)%sg(i) 
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

end module sg_correction
