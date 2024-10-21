#! /bin/tcsh -f
                        
cd /home/belle/varghese/Skim_combined/
source  /home/belle/varghese/.cshrc
setenv BASF_USER_INIT geant_init



#setenv FILE exp$3_run$1_$2_MC_st0_$4
setenv FILE exp$3_run$1_$2_Data

setenv com1 "module put_parameter dtohpi_skim parm\$1"
setenv com2 "module put_parameter dtohpi_skim parm2\$2"
setenv com3 "module put_parameter dtohpi_skim parm3\$3"




basf <<EOF >! ./log_Data_Y5S/log_$FILE.out


path create main
path create Skim

module register fix_mdst dtohpi_skim
path add_module main fix_mdst
path add_module Skim dtohpi_skim 


path add_condition main >:0:Skim
path add_condition main =<:0:KILL

$com1
$com2 
$com3


initialize

histogram define ./hbook_Data_Y5S/$FILE.hbk

process_url http://bweb3/mdst.php?ex=$3&rs=$1&re=$2&skm=HadronBorJ&dt=on_resonance&bl=caseB 



terminate 

EOF






