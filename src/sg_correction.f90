module sg_correction

    ! Module for calculating the finite size correction to the correlation energy

    use system, only: sys_t

    implicit none

    contains
        
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
