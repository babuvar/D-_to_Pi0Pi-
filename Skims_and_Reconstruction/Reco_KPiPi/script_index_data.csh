#! /bin/tcsh -f

cd /home/belle/varghese/Reco_KPiPi/ 

source /home/belle/varghese/.cshrc
setenv BASF_USER_INIT geant_init


setenv FILE $1



basf <<EOF >! ./log_index_DataFull_Y5s/log_$FILE.out


path create main
path create Skim

module register fix_mdst dtohpi
path add_module main fix_mdst
path add_module Skim dtohpi

path add_condition main >:0:Skim
path add_condition main =<:0:KILL

initialize

histogram define ./hbook_index_DataFull_Y5s/$FILE.hbk  

process_event  /home/belle/varghese/Reco_KPiPi/index_Data_Y5S/$FILE 0

terminate 

EOF






