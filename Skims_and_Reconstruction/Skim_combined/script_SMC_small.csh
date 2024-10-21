#! /bin/tcsh -f

cd /gpfs/fs02/belle/users/varghese/Charged_D_Analysis/Skim_combined

source ~/Env/cshrc.pre
setenv BASF_USER_INIT geant_init


foreach FILE (DMtoPiZPiM_50000_Y4s   DPtoPiZPiP_50000_Y4s) #Y(4s)
#foreach FILE (DPtoKPiPi_50000_Y4s DMtoKPiPi_50000_Y4s) #Y(4s)


basf <<EOF >! ./log/log_$FILE.out


path create main
path create Skim

module register fix_mdst dtohpi_skim
path add_module main fix_mdst
path add_module Skim dtohpi_skim

path add_condition main >:0:Skim
path add_condition main =<:0:KILL

initialize

histogram define ./hbook/$FILE.hbk  

process_dir  /gpfs/fs02/belle/users/varghese/Charged_D_Analysis/mcproduzh/gsim/mdst/$FILE

terminate 

EOF

end




