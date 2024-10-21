#! /bin/tcsh -f
                        
cd /home/belle/varghese/Reco_KPiPi
source  /home/belle/varghese/.cshrc
setenv BASF_USER_INIT geant_init

setenv FILE HadronBJ_4s_500K

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


process_event /group/belle/bdata_b/dstprod/dat/e000055/HadronBJ/0127/on_resonance/00/HadronBJ-e000055r000003-b20090127_0910.mdst 10002

terminate 

EOF






