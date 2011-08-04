module folded_spectrum_system_choice


use spawning, only: spawn_hub_real
use hamiltonian, only: slater_condon0_hub_real



implicit none

    !real space hubbard system
    system_gen_excit_ptr => spawn_hub_real
    system_sc0_ptr       => slater_condon0_hub_real
    
    system_decoder_ptr   => decode_det_occ
    system_hub_matel     => hubt
    





end module folded_spectrum_system_choice
