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

    function stochastic_round_spawned_particle(spawn_cutoff, pspawn, rng) result(nspawn)

        ! Stochastically round spawned particle whose population is encoded in
        ! pspawn so that the desired population is achieved on average.

        ! In:
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    pspawn: Encoded spawning probability.
        ! In/Out:
        !    rng: random number generator.
        ! Returns:
        !    nspawn: number of successful spawning events from given pspawn.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use const, only: p, int_p

        integer(int_p), intent(in) :: spawn_cutoff
        real(p), intent(in) :: pspawn
        type(dSFMT_t), intent(inout) :: rng

        integer(int_p) :: nspawn
        real(p) :: padd

        if (pspawn < spawn_cutoff) then

            ! If the spawning amplitude is below the minimum spawning event
            ! allowed, stochastically round it either down to zero or up
            ! to the cutoff.
            if (pspawn > get_rand_close_open(rng)*spawn_cutoff) then
                nspawn = spawn_cutoff
            else
                nspawn = 0_int_p
            end if

        else

            ! Need to take into account the possibilty of a spawning attempt
            ! producing multiple offspring...
            ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and
            ! then spawn a particle with probability pspawn-floor(pspawn).
            nspawn = int(pspawn, int_p)
            padd = pspawn - nspawn

            if (padd > get_rand_close_open(rng)) nspawn = nspawn + 1_int_p

        end if

    end function stochastic_round_spawned_particle

end module stoch_utils
