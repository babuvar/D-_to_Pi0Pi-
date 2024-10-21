#! /bin/tcsh -f

#setenv Stream 11

foreach Stream (10 11)

foreach Type (charged mixed charm uds)

echo "bsub -q b_index ./script_MC.csh        6       872     7       $Type    $Stream"
#bsub -q b_index ./script_MC.csh         807     864     17      $Type    $Stream


end

end

