#! /bin/tcsh -f
                        
setenv BELLE_LEVEL b20090127_0910
setenv BASF_USER_IF basfsh.so
setenv BASF_NPROCESS 0
setenv BELLE_DEBUG opt
setenv USE_GRAND_REPROCESS_DATA 1
setenv BASF_USER_INIT geant_init

source /sw/belle/local/etc/cshrc_general



#setenv FILE exp$3_run$1_$2_MC_st0_$4
setenv FILE exp$3_run$1_$2_Data




basf <<EOF >! ./log_data_4s/log_$FILE.out


path create main
path create Skim

module register fix_mdst dtohpi
path add_module main fix_mdst
path add_module Skim dtohpi 


path add_condition main >:0:Skim
path add_condition main =<:0:KILL



initialize

histogram define ./hbook_data_4s/$FILE.hbk

process_url http://bweb3/mdst.php?ex=$3&rs=$1&re=$2&skm=HadronBorJ&dt=on_resonance&bl=caseB 



terminate 

EOF






