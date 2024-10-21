#! /bin/tcsh -f
                        
setenv BELLE_LEVEL b20090127_0910
setenv BASF_USER_IF basfsh.so
setenv BASF_NPROCESS 0
setenv BELLE_DEBUG opt
setenv USE_GRAND_REPROCESS_DATA 1
setenv BASF_USER_INIT geant_init

source /sw/belle/local/etc/cshrc_general

setenv FILE HadronBJ_4s_100K

basf <<EOF >! ./log/log_$FILE.out

path create main
path create Skim

module register fix_mdst dtohpi
path add_module main fix_mdst
path add_module Skim dtohpi 


path add_condition main >:0:Skim
path add_condition main =<:0:KILL



initialize

histogram define ./hbook/$FILE.hbk  


process_event /group/belle/bdata_b/dstprod/dat/e000055/HadronBJ/0127/on_resonance/00/HadronBJ-e000055r000003-b20090127_0910.mdst 100002

terminate 

EOF






