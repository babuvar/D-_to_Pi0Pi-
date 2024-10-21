#! /bin/tcsh -f

cd /home/belle/varghese/Reco_code_noDstandDst
source  /home/belle/varghese/.cshrc
setenv BASF_USER_INIT geant_init

         
setenv FILE $1_$2

basf <<EOF >! ./log_SMC/log_$FILE.out

path create main
path create Skim

module register fix_mdst dtohpi
path add_module main fix_mdst
path add_module Skim dtohpi

path add_condition main >:0:Skim
path add_condition main =<:0:KILL

initialize

histogram define ./hbook_SMC/$FILE.hbk  

process_dir /home/belle/varghese/mcproduzh/gsim/mdst/$1/dir$2/
terminate 

EOF






