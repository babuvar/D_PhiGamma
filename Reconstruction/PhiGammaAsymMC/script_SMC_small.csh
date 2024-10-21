#! /bin/tcsh -f
                     
cd /gpfs/fs02/belle/users/varghese/PhiGamma_Analysis/PhiGammaAsymMC

source ~/Env/cshrc.pre
setenv BASF_USER_INIT geant_init

#foreach FILE (D0toPhiGamma_573 D0BtoPhiGamma_573)   #Y(5s) 
#foreach FILE (D0toPhiGamma_3550 D0BtoPhiGamma_3550) #Y(4s)
#foreach FILE (D0toPhiGamma_3550_2 D0BtoPhiGamma_3550_2) #Y(4s)
#foreach FILE (D0toPhiGamma_7100 D0BtoPhiGamma_7100) #Y(4s)
#foreach FILE (D0toPhiGamma_1146 D0BtoPhiGamma_1146) #Y(5s)
#foreach FILE (D0toPhiGamma_5K D0BtoPhiGamma_5K) #Y(4s)
#foreach FILE (D0toPhiGamma_1K_Y5S D0BtoPhiGamma_1K_Y5S) #Y(5s)

foreach FILE (D0toPhiGamma_30K D0BtoPhiGamma_30K) #Y(4s)
#foreach FILE (D0toPhiGamma_30K_2 D0BtoPhiGamma_30K_2) #Y(4s)
#foreach FILE (D0toPhiGamma_5K_2 D0BtoPhiGamma_5K_2) #Y(4s)
#foreach FILE (D0toPhiGamma_6K_Y5S D0BtoPhiGamma_6K_Y5S) #Y(5s)

#foreach FILE (D0toPhiGamma_104000 D0BtoPhiGamma_104000) #Y(4s)
#foreach FILE (D0toPhiGamma_55000_Y5s D0BtoPhiGamma_55000_Y5s) #Y(5s)
#foreach FILE (D0toPhiGamma_52900_Y5s D0BtoPhiGamma_52900_Y5s) #Y(5s)
#foreach FILE (D0toPhiPi0_105800 D0BtoPhiPi0_105800)  #Y(5s)

#setenv FILE D0toPhiGamma_3550
#setenv FILE D0BtoPhiGamma_3550_2
#setenv FILE trial2


basf <<EOF >! ./log/log_$FILE.out


path create main
path create Skim

module register fix_mdst d0rad
path add_module main fix_mdst
path add_module Skim d0rad

path add_condition main >:0:Skim
path add_condition main =<:0:KILL

initialize

histogram define ./hbook/$FILE.hbk  

process_dir     /gpfs/fs02/belle/users/varghese/PhiGamma_Analysis/mcproduzh/gsim/mdst/$FILE/

terminate 

EOF

end




