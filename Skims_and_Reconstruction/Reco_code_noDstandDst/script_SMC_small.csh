#! /bin/tcsh -f

setenv BELLE_LEVEL b20090127_0910
setenv BASF_USER_IF basfsh.so
setenv BASF_NPROCESS 0
setenv BELLE_DEBUG opt
setenv USE_GRAND_REPROCESS_DATA 1
setenv BASF_USER_INIT geant_init

source /sw/belle/local/etc/cshrc_general


foreach FILE (DMtoPiZPiM_50000_Y4s   DPtoPiZPiP_50000_Y4s) #Y(4s)



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

process_dir  /home/belle/varghese/mcproduzh/gsim/mdst/$FILE/

terminate 

EOF

end




