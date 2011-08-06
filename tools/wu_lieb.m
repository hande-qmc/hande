(* See, e.g. http://dx.doi.org/10.1016/S0378-4371(02)01785-5 *)
(* E0 gives the exact ground state energy for a 1D Hubbard model at half-filling for a given U and N sites in the crystal cell in the thermodynamic limit *)
E0[N_,U_] := -4 N NIntegrate[BesselJ[0,w]BesselJ[1,w]/(w(1+Exp[w U/2])), {w,0,Infinity}]
