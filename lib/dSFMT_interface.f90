module dSFMT_interface

! This contains a handy wrapper around a subset of of the functionality in
! dSFMT RNG (double precision SIMD-oriented Fast Mersenne Twister random number
! generator).

! See also the functions defined in dSFMT_wrapper.cpp.

use const
implicit none

! It is much faster to generate random numbers in blocks.
! genrand_real2 is a wrapper around accessing the random_store,
! filling it up again as necessary.

! Testing indicates that 50000 is a very good size for the array storing the
! random numbers.  Testing was done standalone, so undoubtedly influenced by
! cache size and this might be different for real-world applications, but it's easy to
! change to allocatable later on.
integer, parameter :: random_store_size=5*10**4
real(dp), save :: random_store(random_store_size)

integer, save :: current_element=1

contains

    subroutine dSFMT_init(seed)

        ! Initialise the dSFMT RNG and fill random_store with
        ! a block of random numbers in interval [0,1).
        !
        ! In:
        !    seed: seed for the RNG.

        integer, intent(in) :: seed

        call init_gen_rand(seed)

        call fill_array_close_open(random_store, random_store_size)

    end subroutine dSFMT_init

    function genrand_real2() result(r)

        ! Return:
        !    random number in interval [0,1).  Name comes from the function
        !    defined in the original Mersenne Twist implementation.

        real(dp) :: r

        if (current_element == random_store_size+1) then
            ! Run out of random numbers: get more.
            current_element = 1
            call fill_array_close_open(random_store, random_store_size)
        end if

        r = random_store(current_element)
        current_element = current_element + 1 


    end function genrand_real2

    subroutine test_rand()

        integer :: i
        real(dp) :: r

        call dSFMT_init(7)

        do i = 1,10**7
            r = genrand_real2()
        end do
        write (6,*) r
        stop

    end subroutine test_rand
        

end module dSFMT_interface
