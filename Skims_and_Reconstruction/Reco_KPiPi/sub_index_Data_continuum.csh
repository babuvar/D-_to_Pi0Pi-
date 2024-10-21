#! /bin/tcsh -f

foreach name (`grep -i "dex" IndexFileNames_FullData_continuum.txt`)

bsub -q l ./script_index_data_continuum.csh $name
#echo "bsub -q l ./script_index_data_continuum.csh $name"

end


