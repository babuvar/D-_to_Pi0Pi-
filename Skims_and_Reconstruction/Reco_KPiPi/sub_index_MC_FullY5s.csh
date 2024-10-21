#! /bin/tcsh -f

foreach name (`grep -i "dex" IndexFileNames_FullMC_Y5s.txt`)

bsub -q l ./script_index_MC.csh $name
#echo "bsub -q l ./script_index_MC.csh $name"

end


