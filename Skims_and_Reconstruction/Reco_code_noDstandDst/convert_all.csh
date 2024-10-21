#! /bin/tcsh -f


#foreach i (hbook_*)

#foreach i (hbook_MC1  hbook_MC11 hbook_MC12 hbook_MC13 hbook_MC14 hbook_MC15 hbook_MC2 hbook_MC3 hbook_MC4 hbook_MC5)

foreach i (hbook_MC_5s_0 hbook_MC_5s_1 hbook_MC_5s_2 hbook_MC_5s_3 hbook_MC_5s_4 hbook_MC_5s_5)



cd $i

./convert.csh

cd ..

end


