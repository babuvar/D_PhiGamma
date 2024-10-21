#! /bin/tcsh -f
                        
cd /gpfs/fs02/belle/users/varghese/PhiGammaAsymData

source ~/Env/cshrc.pre
setenv BASF_USER_INIT geant_init

setenv FILE 3s_data_$1_$2_op

basf <<EOF >! ./log/log_$FILE.out

path create main
path create Skim

module register fix_mdst d0rad
path add_module main fix_mdst
path add_module Skim d0rad

path add_condition main >:0:Skim
path add_condition main =<:0:KILL


initialize

histogram define ./hbook_3s/$FILE.hbk  

process_url http://bweb3/mdst.php?ex=49&rs=$1&re=$2&skm=HadronBorJ&dt=3S_scan&bl=caseB

terminate 

EOF






