#! /bin/bash

~/hubbard_fciqmc/tools/join_rdms.sh
~/hubbard_fciqmc/tools/average_rdm.py rdm > average_rdm
~/hubbard_fciqmc/tools/analyse_rdm/analyse.x > entropy_estimates
~/hubbard_fciqmc/tools/analyse_rdm/stats.py entropy_estimates > final_estimates
