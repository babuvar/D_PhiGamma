#! /bin/tcsh -f
                       
cd /gpfs/fs02/belle/users/varghese/PhiPi0CutOptandAsymGMC

source ~/Env/cshrc.pre
setenv BASF_USER_INIT geant_init

setenv FILE 4s_mc_$1_s$2_exp$3_op

basf <<EOF >! ./log$4/log_$FILE.out

path create main
path create Skim

module register fix_mdst d0rad
path add_module main fix_mdst
path add_module Skim d0rad

path add_condition main >:0:Skim
path add_condition main =<:0:KILL

initialize

histogram define ./hbook$4/$FILE.hbk  

`/gpfs/home/belle/nishida6/public/fileloc/jwicht_script/skim-process_event-mc.sh d0rad on_resonance $1 $2 $3`
 
terminate 

EOF






