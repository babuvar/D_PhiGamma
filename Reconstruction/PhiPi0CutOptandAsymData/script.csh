#! /bin/tcsh -f
                        
cd /gpfs/fs02/belle/users/varghese/PhiPi0CutOptandAsymData

source ~/Env/cshrc.pre
setenv BASF_USER_INIT geant_init

setenv FILE 4s_data_exp$1_op
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

`/gpfs/home/belle/nishida6/public/fileloc/jwicht_script/skim-process_event-data.sh d0rad on_resonance $1`

terminate 

EOF






