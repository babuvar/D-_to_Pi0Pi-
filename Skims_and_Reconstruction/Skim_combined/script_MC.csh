#! /bin/tcsh -f
                        
cd /gpfs/fs02/belle/users/varghese/Charged_D_Analysis/Skim_combined/
source ~/Env/cshrc.pre
setenv BASF_USER_INIT geant_init



setenv FILE exp$3_run$1_$2_MC_st$5_$4


setenv com1 "module put_parameter dtohpi_skim parm\$1"
setenv com2 "module put_parameter dtohpi_skim parm2\$2"
setenv com3 "module put_parameter dtohpi_skim parm3\$3"
setenv com4 "module put_parameter dtohpi_skim parm4\$4"
setenv com5 "module put_parameter dtohpi_skim parm5\$5"

basf <<EOF >! ./log_MC$5/log_$FILE.out



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
$com4
$com5


initialize

histogram define ./hbook_MC$5/$FILE.hbk  


 
process_url  http://bweb3/montecarlo.php?ex=$3&rs=$1&re=$2&ty=evtgen-$4&dt=on_resonance&bl=caseB&st=$5 
 

terminate 

EOF






