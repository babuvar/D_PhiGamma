#! /bin/tcsh -f

foreach i (charged charm mixed uds)

foreach j (11 13 15 17 19 21 23 25 27 7 9)

h2root 4s_mc_${i}_s10_exp${j}_op.hbk 4s_mc_${i}_s10_exp${j}_op.root

end

foreach k (31 33 35 37 39 41 43 45 47 49 51 55 61 63 65)

h2root 4s_mc_${i}_s0_exp${k}_op.hbk 4s_mc_${i}_s0_exp${k}_op.root

end

end
