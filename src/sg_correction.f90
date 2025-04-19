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

            ! Find all unique pairs of momentum transfer vectors between basis functions 
            ! and allocate the corresponding array to store S(G) values

            ! In:
            !   sys: Information about system, specifically the number of basis functions and 
            !        the basis functions themselves
            ! In/Out:
            !   qs: Information about state, nspaces, q_vectors are access through here, and 
            !       the sg arrat is allocated and zerod
            
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

        subroutine calculate_sg_contribution(sys, connection, qs, real_population, ispace)

            ! Calculate the contribution to S(G) of the current determinant

            ! In:
            !   sys: Information about system, the cell length and basis functions are used 
            !   connection: Used to identify the occupied and virtual orbitals as well as phase 
            !               change required due to permutation
            !   ispace: Identifies current walker space, this should always be 1 (i.e., the real part)
            !   real_population :: Current determinant population 
            ! In/Out:
            !   qs: State information, used to access population data and q_vectors. Updates to sg array 
            !       go through here

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

            ! Average the calculated S(G) values at the end of the report loop and calculate 
            ! the uncorrected correlation energy via S(G)

            ! In:
            !   qmc_in: Read in qmc options, gives the number of cycles per report loop 
            !   sys: System information, provides cell length 
            ! In/Out: 
            !   qs: State information, gives S(G) values which are averaged, as well as 
            !       access to the q_vectors 
        
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

            ! Dump information about each momentum transfer vector, its contribution to S(G)
            ! and the current uncorrected correlation energy 

            ! In:
            !   qs: State information, used to access sg array and q_vectors

            type(qmc_state_t), intent(in) :: qs 

            integer :: sg_unit, ierr 
            integer :: i

            open(newunit=sg_unit, file="SG", status="replace", action="write", iostat=ierr)
            write(sg_unit, *) 'energy by S(G): ', qs%estimators(1)%energy_sg

            do i = 1, size(qs%estimators(1)%sg, 1)
                write(sg_unit, *) &
                    qs%estimators(1)%q_vectors(1, i), &
                    qs%estimators(1)%q_vectors(2, i), &
                    qs%estimators(1)%q_vectors(3, i), &
                    qs%estimators(1)%sg(i) 
            end do

            close(sg_unit)

        end subroutine dump_sg_data
        
        ! This can probably be a lot faster
        pure function find_q_index(q_vectors, G) result(idx)

            ! Given a momentum tranfser vector G, compare it to a full list
            ! of vectors q_vectors and find the index of the q vector which
            ! matches G, if none are found return 0

            ! In:
            !   q_vectors: Full list of 3D vectors to compare to
            !   G: A single 3D vector to check
            ! Out:
            !   idx: Index of the q vector which matches G 

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
