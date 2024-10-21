#! /bin/tcsh -f

foreach name (`grep -i "dex" IndexFileNames_FullData_Y5s.txt`)

bsub -q l ./script_index_data.csh $name
#echo "bsub -q l ./script_index_data.csh $name"

end


