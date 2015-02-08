module stoch_utils

implicit none

public
private :: stochastic_round_int_32, stochastic_round_int_64

interface stochastic_round
    module procedure stochastic_round_int_32
    module procedure stochastic_round_int_64
end interface stochastic_round

contains

    subroutine stochastic_round_int_32(rng, population, cutoff, ntypes)

        ! For any values in population less than cutoff, round up to cutoff or
        ! down to zero. This is done such that the expectation value of the
        ! resulting populations is equal to the input values.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    population: populations to be stochastically rounded.
        !    cutoff: the value to round up to.
        !    ntypes: the number of values in population to apply this op to.

        use const, only: int_32, p
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(dSFMT_t), intent(inout) :: rng
        integer(int_32), intent(inout) :: population(:)
        integer(int_32), intent(in) :: cutoff
        integer, intent(in) :: ntypes
        integer :: itype
        real(p) :: r

        do itype = 1, ntypes
            if (abs(population(itype)) < cutoff .and. population(itype) /= 0_int_32) then
                r = get_rand_close_open(rng)*cutoff
                if (abs(population(itype)) > r) then
                    population(itype) = sign(cutoff, population(itype))
                else
                    population(itype) = 0_int_32
                end if
            end if
        end do

    end subroutine stochastic_round_int_32

    subroutine stochastic_round_int_64(rng, population, cutoff, ntypes)

        ! For any values in population less than cutoff, round up to cutoff or
        ! down to zero. This is done such that the expectation value of the
        ! resulting populations is equal to the input values.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    population: populations to be stochastically rounded.
        !    cutoff: the value to round up to.
        !    ntypes: the number of values in population to apply this op to.

        use const, only: int_64, p
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(dSFMT_t), intent(inout) :: rng
        integer(int_64), intent(inout) :: population(:)
        integer(int_64), intent(in) :: cutoff
        integer, intent(in) :: ntypes
        integer :: itype
        real(p) :: r

        do itype = 1, ntypes
            if (abs(population(itype)) < cutoff .and. population(itype) /= 0_int_64) then
                r = get_rand_close_open(rng)*cutoff
                if (abs(population(itype)) > r) then
                    population(itype) = sign(cutoff, population(itype))
                else
                    population(itype) = 0_int_64
                end if
            end if
        end do

    end subroutine stochastic_round_int_64

end module stoch_utils
