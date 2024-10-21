#! /bin/bash

mkdir Files_done
mkdir RootFiles_by_exp


for VAR in     43 53 67 69 71

do


hadd KPiPi_exp${VAR}_data_4s.root exp${VAR}_*root
mv KPiPi_exp${VAR}_data_4s.root RootFiles_by_exp
mv exp${VAR}_*root Files_done


done
