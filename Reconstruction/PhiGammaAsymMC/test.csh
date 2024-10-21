#! /bin/tcsh -f
                        
cd /gpfs/fs02/belle/users/varghese/PhiGammaAsymMC

source ~/Env/cshrc.pre
setenv BASF_USER_INIT geant_init

setenv FILE $1_$2


echo "./log/log_$FILE.out"

echo "./hbook/$FILE.hbk"  

echo "process_dir /gpfs/fs02/belle/users/varghese/mcproduzh/gsim/mdst/$1/dir$2/"

echo




