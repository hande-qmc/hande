module alias

! A module using the alias method to choose an integer out of a discrete probability distribution.
! The naming convention here is such that the discrete probability p(i) is given as 
! weights(i)/total weight.

use const, only: i0, p

implicit none

contains

    ! [review] - JSS: name isn't immediately obvious.  Prec?
    function select_weighted_value_prec(rng, N, aliasU, aliasK) result(ret)

        ! Select an element, i=1..N with probability from pre-generated alias method weights.
        ! [review] - JSS: O(N) setup cost, O(1) to select?  Repetition below the arguments comments.
        ! This uses the alias method, requiring O(N) storage, and O(N) time.
        
        ! In:
        !    N: the number of objects to select from
        ! [review] - JSS: what are alias reals and alias integers?
        !    aliasU: a length N array of precomputed alias reals.
        !    aliasK: a length N array of precomputed alias integers
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    ret: the index of the element chosen.

        ! The 'alias method' allows one to select from a discrete probability distribution of N objects 
        ! (with object i having probability p_i) in O(1) time. There's an O(N) storage and O(N) setup cost 
        ! - a list of N reals (U_i) and N integers (K_i) which requires O(N) setup.

        ! Pick a random real number x, 0<=x<N.
        ! Let i=floor(x) and V=x-i. (so i is an integer and V the remainder).
        ! The randomly selected object, ret, will be ret=i if V<U_i and ret=K_i otherwise.

        ! Here's Knuth's exercise:
        ! Vol2: 3.4.1 Ex 7 [20] (A. J. Walker)
        ! Suppose we have a bunch of cubes of k different colors, say n_j cubes of color C_j for 1<=j<=k, and we have k boxes
        ! {B_1,...,B_k} each of which can hold exactly n cubes.  Furthermore n_1+...+n_k=kn, so the cubes will just fit in the
        ! boxes.  Prove (constructively) that there is always a way to put the cubes into the boxes so that each box contains at
        ! most two different colors of cubes; in fact there is a way to do it so that, whenever box B_j contains two colors, one of
        ! those colors is C_j.  Show how to use this principle to compute the U and K tables given a probability distribution
        ! (p_1,...p_k).

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        type(dSFMT_t), intent(inout) :: rng
        
        integer, intent(in) :: N
        integer :: ret

        real(p) :: aliasU(N)
        integer :: aliasK(N)
        
        real(p) :: x
        integer :: K 

        x = get_rand_close_open(rng)*N
        K = floor(x)
        x = x-K
        K = K+1
        if (x < aliasU(K)) then
            ret = K
        else
            ret = aliasK(K)
        end if

    end function select_weighted_value_prec

    ! [review] - JSS: is the alias method useful compared to simple binary search of the probabilities for one-off selections?
    ! [review] - VAN: How would you do the binary search to select an integer from a discrete probability distribution?
    function select_weighted_value(rng, N, weights, totweight) result(ret)

        ! Select an element, i=1..N with probability weights(i)/totweight.
        ! This uses the alias method, requiring O(N) storage, and O(N) time
        
        ! In:
        !    N: the number of objects to select from
        !    weights: a length N array of reals containing the weights of each element
        !    totweight: sum(weights(1:N))
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    ret: the index of the element chosen.

        ! See notes in select_weighted_value_prec

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        type(dSFMT_t), intent(inout) :: rng
        
        integer, intent(in) :: N
        real(p), intent(in) :: totweight, weights(N) 
        integer :: ret

        real(p) :: aliasU(N)
        integer :: aliasK(N)
        
        call generate_alias_tables(N, weights, totweight, aliasU, aliasK)        
        ret = select_weighted_value_prec(rng, N, aliasU, aliasK)
    end function select_weighted_value

    subroutine generate_alias_tables(N, weights, totweight, aliasU, aliasK)

        ! Generate an alias table for the alias method.
        ! This requires O(N) time and O(2N) scratch space
        !
        ! In:
        !    N: number of objects to select from
        !    weights: a length N array of reals containing the weights of each element
        !    totweight: sum(weights(1:N))
        ! Out:
        !    aliasU: a length N array of reals for the U table
        !    aliasK: a length N array of integers for the K table.

        !
        ! The alias method (a la Wikipedia)
        ! The distribution may be padded with additional probabilities /p_i / = 0 to 
        ! increase n to a convenient value, such as a power of two.

        ! To generate the table, first initialize /U_i / = /np_i /. While doing this, 
        ! divide the table entries into three categories:

        !  * The "overfull" group, where /U_i / > 1,
        !  * The "underfull" group, where /U_i / < 1, and
        !  * The "exactly full" group, where /U_i / = 1.

        ! If /U_i / = 1, the corresponding value K_i will never be consulted and is unimportant, 
        ! but a value of /K_i / = /i/ is sensible.

        ! As long as not all table entries are exactly full, repeat the following steps:

        ! 1. Arbitrarily choose an overfull entry /U_i / > 1 and an underfull
        !    entry /U_j / < 1. (If one of these exists, the other must, as well.)
        ! 2. Allocate the unused space in entry j to outcome i, by setting /K_j /
        !    = /i/.
        ! 3. Remove the space from entry i that entry j got from entry i by changing 
        ! /U_i / = /U_i / - (1 - /U_j /) = /U_i / + /U_j / - 1.
        ! 4. Entry j is now exactly full.
        ! 5. Assign entry i to the appropriate category based on the new value of
        !    U_i .  

        integer, intent(in) :: N
        real(p), intent(in) :: totweight, weights(N) 
        real(p), intent(out) :: aliasU(N)
        integer, intent(out) :: aliasK(N)

        ! Working space:  We need a list of the underfull and the overfull, and a copy of the weights
        integer :: underfull(N)
        integer :: overfull(N)
        integer :: i, nunder, nover, ov, un

        aliasU = weights * (N / totweight)
        nunder = 0
        nover = 0
        do i = 1, N
            if (aliasU(i) <= 1.0_p) then
                nunder = nunder + 1
                underfull(nunder) = i
            else ! account for aliasU(i)=1 as underfull
                nover = nover +1
                overfull(nover) = i
            end if
            ! [review] - JSS: in case of what?!
            ! [review] - VAN: Sensible to set it in case there is a tiny numerical error
            ! [review] - VAN: when matching up the last underfull and overfull Us. 
            aliasK(i) = i ! This is a sensible safe choice in case aliasK(i) is not set below.
        end do
        do while (nover > 0 .and. nunder > 0)
            ! match the last nover with the last nunder
            ! [review] - JSS: how arbitrary is this choice of over and under and how does it impact efficiency?
            ov = overfull(nover)
            un = underfull(nunder)
            ! put ov as the alternate for un
            aliasK(un) = ov
            ! un is now full, so we remove it from the un list
            nunder = nunder - 1
            ! remove that much probability from the ov's amount
            aliasU(ov) = aliasU(ov) - (1 - aliasU(un))
            if (aliasU(ov) < 1.0_p) then
                ! Move it to the under list
                nunder = nunder + 1
                underfull(nunder) = overfull(nover)
                nover = nover - 1
            end if
        end do

    end subroutine generate_alias_tables

end module alias  
