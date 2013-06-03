program analyse_rdm

    use dSFMT_interface

    implicit none

    integer :: seed = 328576
    integer :: rdm_unit, stat, lwork, ierr, info
    real(dp), allocatable :: work(:)
    real(dp) :: rand_var
    integer :: rdm_size
    integer :: i, j, k
    logical :: does_exist
    character(11) :: filename = 'average_rdm'

    real(dp) :: next_gaussian_rv
    logical :: gaussian_rv_left

    real(dp), allocatable :: rdm_elements(:,:), rdm_with_noise(:,:), rdm_errors(:,:)
    real(dp), allocatable :: eigv(:)
    real(dp) :: temp_store(2)
    real(dp) :: vn_entropy, renyi_2

    inquire(file=filename, exist=does_exist)

    if (.not. does_exist) then
        write(6,'(a21)') "No RDM file detected."
        stop 999
    end if

    call dSFMT_init(seed)

    rdm_unit = 10
    open(rdm_unit, file=filename, status='old')
    read(rdm_unit, *, iostat=stat) rdm_size

    allocate(rdm_elements(rdm_size, rdm_size))
    allocate(rdm_with_noise(rdm_size, rdm_size))
    allocate(rdm_errors(rdm_size, rdm_size))
    allocate(eigv(rdm_size))
    rdm_elements = 0.0_dp
    rdm_errors = 0.0_dp
    eigv = 0.0_dp

    ! Read in the RDM elements and errors to the above arrays.
    do i = 1, rdm_size
        do j = i, rdm_size
            read(rdm_unit, *, iostat=stat) temp_store
            rdm_elements(i,j) = temp_store(1)
            rdm_elements(j,i) = temp_store(1)
            rdm_errors(i,j) = temp_store(2)
        end do
    end do

    ! Temporary space, as dsyev will destory the upper half of the input matrix.
    rdm_with_noise = rdm_elements

    ! Find the optimal size of the workspace.
    allocate(work(1), stat=ierr)
    call dsyev('N', 'U', rdm_size, rdm_with_noise, rdm_size, eigv, work, -1, info)
    lwork = nint(work(1))
    deallocate(work)
    allocate(work(lwork), stat=ierr)

    ! Now perform the diagonalisation.
    call dsyev('N', 'U', rdm_size, rdm_with_noise, rdm_size, eigv, work, lwork, info)

    ! Find the mean estimate.
    vn_entropy = calculate_von_neumann_entropy(eigv)
    write(6,'(a35,f12.7)') "# von Neumann entropy mean estimate =", vn_entropy

    ! Now, estimate the error on the estimate of the above VN entropy by adding Gaussian noise
    ! to the mean RDM elements.
    ! For testing, just do this 10000 times.
    do k = 1, 10000

        ! Loop over all elements and add the appropriate noise.
        do i = 1, rdm_size
            do j = i, rdm_size
                rand_var = generate_gaussian_rv(rdm_errors(i,j))
                rdm_with_noise(i,j) = rdm_elements(i,j) + rand_var
                rdm_with_noise(j,i) = rdm_with_noise(i,j)
            end do
        end do

        call dsyev('N', 'U', rdm_size, rdm_with_noise, rdm_size, eigv, work, lwork, info)

        ! Calculate and write out the estmate from the noisy RDM.
        vn_entropy = calculate_von_neumann_entropy(eigv)
        write(6,'(f12.7)') vn_entropy

    end do

    ! Next calculate Renyi 2 entropy.

    ! First, the mean estimate from the RDM without noise added.
    rdm_with_noise = rdm_elements
    call dsyev('N', 'U', rdm_size, rdm_with_noise, rdm_size, eigv, work, lwork, info)

    renyi_2 = calculate_renyi_2_entropy(eigv)
    write(6,'(a33,f12.7)') "# Renyi 2 entropy mean estimate =", renyi_2

    ! Then, once again, add noise 10000 times to estimate the error.
    do k = 1, 10000

        do i = 1, rdm_size
            do j = i, rdm_size
                rand_var = generate_gaussian_rv(rdm_errors(i,j))
                rdm_with_noise(i,j) = rdm_elements(i,j) + rand_var
                rdm_with_noise(j,i) = rdm_with_noise(i,j)
            end do
        end do

        call dsyev('N', 'U', rdm_size, rdm_with_noise, rdm_size, eigv, work, lwork, info)

        renyi_2 = calculate_renyi_2_entropy(eigv)
        write(6,'(f12.7)') renyi_2

    end do

contains

        function calculate_r2_error(rdm_size, rdm, rdm_error) result(error)

            ! This function calculate the error on Renyi 2 directly using the error
            ! propagation formula, rather than by adding noise. This is for checking.

            integer, intent(in) :: rdm_size
            real(dp), intent(in) :: rdm(rdm_size, rdm_size), rdm_error(rdm_size, rdm_size)
            real(dp) :: error, elem_sum, derivative
            integer :: i, j

            elem_sum = 0.0_dp
            do i = 1, rdm_size
                do j = 1, rdm_size
                    elem_sum = elem_sum + rdm(i,j)**2
                end do
            end do
            elem_sum = 2.0_dp/elem_sum

            error = 0.0_dp
            do i = 1, rdm_size
                do j = 1, rdm_size
                    derivative = elem_sum*rdm(i,j)
                    derivative = derivative/log(2.0_dp)
                    error = error + (derivative**2)*(rdm_error(i,j)**2)
                end do
            end do
            error = sqrt(error)

        end function calculate_r2_error

        function calculate_r2_direct(rdm_size, rdm) result(renyi_2)

            ! Function to calculate Renyi 2 from matrix elements rather than from
            ! the eigenvalues. This is for checking.

            integer, intent(in) :: rdm_size
            real(dp), intent(in) :: rdm(rdm_size, rdm_size)
            real(dp) :: renyi_2, elem_sum
            integer :: i, j

            renyi_2 = 0.0_dp
            do i = 1, rdm_size
                do j = 1, rdm_size
                    renyi_2 = renyi_2 + rdm(i,j)**2
                end do
            end do

            renyi_2 = -log(renyi_2)/log(2.0_dp)

        end function calculate_r2_direct

        function calculate_von_neumann_entropy(eigv) result(vn_entropy)

            ! Calculate VN entropy from a set of eigenvalues.

            real(dp), intent(out) :: eigv(:)
            real(dp) :: vn_entropy
            integer :: i

            vn_entropy = 0.0_dp
            do i = 1, ubound(eigv,1)
                ! Throw away all eigenvalues which are zero or negative.
                if (.not. (eigv(i) > 0.0_dp)) then
                    cycle
                end if
                vn_entropy = vn_entropy - eigv(i)*(log(eigv(i))/log(2.0_p))
            end do

        end function calculate_von_neumann_entropy

        function calculate_renyi_2_entropy(eigv) result(renyi_2)

            ! Calculate Renyi entropy from a set of eigenvalues.

            real(dp), intent(out) :: eigv(:)
            real(dp) :: renyi_2
            integer :: i

            renyi_2 = 0.0_dp
            do i = 1, ubound(eigv,1)
                renyi_2 = renyi_2 + eigv(i)**2
            end do
            renyi_2 = -log(renyi_2)/log(2.0_dp)

        end function calculate_renyi_2_entropy

        function generate_gaussian_rv(sd) result(gaussian_rv)

            ! This function returns a random variable from a Gaussian distribution with
            ! standard deviation sd (and mean of 0) as input by the user. It uses the
            ! polar form of the Box-Muller transform. This transform generates two
            ! Gaussian random variables each time it is performed, and so the second unused
            ! random variable is stored and used next time this function is called to save
            ! time.

            real(dp), intent(in) :: sd
            real(dp) :: gaussian_rv
            real(dp) :: u, v, s

            ! If there is still a Gaussian random variable left over from
            ! last time, simply return with this.
            if (gaussian_rv_left) then
                gaussian_rv = sd*next_gaussian_rv
                gaussian_rv_left = .false.
                return
            end if

            do
                ! Uniform random variable between 0 and 1.
                u =  genrand_real2()
                ! Uniform random variable between -1 and 1.
                u = 2.0_dp*u - 1
                v =  genrand_real2()
                v = 2.0_dp*v - 1
                s = u*u + v*v
                ! Only exit the loop if s < 1.0 (if u and v are within a circle with unit radius).
                ! If not, generate another pair.
                if (s < 1.0_dp) exit
            end do

            ! This produces a random variable from a Gaussian distribution with mean 0 and standard
            ! deviation 1.
            gaussian_rv = u*sqrt((-2.0_dp*log(s))/s)
            ! And this transforms the random variable to be from a Gaussian with mean
            ! 0 and standard deviation sd.
            gaussian_rv = sd*gaussian_rv

            ! Store the second produced random variable for next time.
            next_gaussian_rv = v*sqrt((-2.0_dp*log(s))/s)
            gaussian_rv_left = .true.

        end function generate_gaussian_rv

end program analyse_rdm
