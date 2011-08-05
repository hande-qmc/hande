module folded_spectrum_system_choice


use const, only: p
use proc_pointers, only: system_gen_excit_ptr, system_sc0_ptr, system_decoder_ptr, system_update_proj_energy_ptr

use system, only: hubt


implicit none

real(p) :: system_hub_matel

contains
    subroutine initialise_hubbard_real_space_system()
        use spawning, only: gen_excit_hub_real
        use hamiltonian, only: slater_condon0_hub_real
        use energy_evaluation, only: update_proj_energy_hub_real

        use determinants, only: decode_det_occ
        implicit none
        system_hub_matel = hubt

        system_gen_excit_ptr => gen_excit_hub_real
        system_sc0_ptr       => slater_condon0_hub_real
    
        system_decoder_ptr   => decode_det_occ
        system_update_proj_energy_ptr => update_proj_energy_hub_real
    end subroutine initialise_hubbard_real_space_system





end module folded_spectrum_system_choice
