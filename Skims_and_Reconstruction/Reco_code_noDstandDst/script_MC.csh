#! /bin/tcsh -f
                        
setenv BELLE_LEVEL b20090127_0910
setenv BASF_USER_IF basfsh.so
setenv BASF_NPROCESS 0
setenv BELLE_DEBUG opt
setenv USE_GRAND_REPROCESS_DATA 1
setenv BASF_USER_INIT geant_init

source /sw/belle/local/etc/cshrc_general


setenv FILE exp$3_run$1_$2_MC_st$5_$4


basf <<EOF >! ./log_MC$5/log_$FILE.out



path create main
path create Skim

module register fix_mdst dtohpi
path add_module main fix_mdst
path add_module Skim dtohpi 


path add_condition main >:0:Skim
path add_condition main =<:0:KILL



initialize

histogram define ./hbook_MC$5/$FILE.hbk  


 
process_url  http://bweb3/montecarlo.php?ex=$3&rs=$1&re=$2&ty=evtgen-$4&dt=on_resonance&bl=caseB&st=$5 
 

terminate 

EOF






