#! /bin/tcsh -f
                     
cd /gpfs/fs02/belle/users/varghese/PhiGammaAsymMC

source ~/Env/cshrc.pre
setenv BASF_USER_INIT geant_init



foreach FILE (D0toPhiGamma_5200_0 D0toPhiGamma_5200_1 D0toPhiGamma_5200_2 D0toPhiGamma_5200_3 D0toPhiGamma_5200_4 D0toPhiGamma_5200_5 D0BtoPhiGamma_5200_0 D0BtoPhiGamma_5200_1 D0BtoPhiGamma_5200_2 D0BtoPhiGamma_5200_3 D0BtoPhiGamma_5200_4 D0BtoPhiGamma_5200_5)  #Y(4s)

#foreach FILE (D0toPhiGamma_550_Y5s_0 D0toPhiGamma_550_Y5s_1 D0toPhiGamma_550_Y5s_2 D0toPhiGamma_550_Y5s_3 D0toPhiGamma_550_Y5s_4 D0toPhiGamma_550_Y5s_5 D0BtoPhiGamma_550_Y5s_0 D0BtoPhiGamma_550_Y5s_1 D0BtoPhiGamma_550_Y5s_2 D0BtoPhiGamma_550_Y5s_3 D0BtoPhiGamma_550_Y5s_4 D0BtoPhiGamma_550_Y5s_5)


basf <<EOF >! ./log/log_$FILE.out


path create main
path create Skim

module register fix_mdst d0rad
path add_module main fix_mdst
path add_module Skim d0rad

path add_condition main >:0:Skim
path add_condition main =<:0:KILL

initialize

histogram define ./hbook_x6/$FILE.hbk  

process_dir  /gpfs/fs02/belle/users/varghese/mcproduzh/gsim/mdst/$FILE/

terminate 

EOF

end




