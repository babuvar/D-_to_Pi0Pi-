#! /bin/bash

mkdir Files_done
mkdir RootFiles_by_exp


for VAR in 7 9 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 55 61 63 65 
do


hadd KPiPi_exp${VAR}_GMC_4s.root exp${VAR}_*root
mv KPiPi_exp${VAR}_GMC_4s.root RootFiles_by_exp
mv exp${VAR}_*root Files_done


done
