#! /bin/tcsh -f
                     
cd /gpfs/fs02/belle/users/varghese/PhiGammaAsymMC

source ~/Env/cshrc.pre
setenv BASF_USER_INIT geant_init

setenv FILE D0toPhiGamma_60400

basf <<EOF >! ./log/log_$FILE.out

path create main
path create Skim

module register fix_mdst d0rad
path add_module main fix_mdst
path add_module Skim d0rad

path add_condition main >:0:Skim
path add_condition main =<:0:KILL

initialize

histogram define ./hbooktest/$FILE.hbk  

process_dir /gpfs/fs02/belle/users/varghese/mcproduzh/gsim/mdst/D0toPhiGamma_60400 

terminate 

EOF






