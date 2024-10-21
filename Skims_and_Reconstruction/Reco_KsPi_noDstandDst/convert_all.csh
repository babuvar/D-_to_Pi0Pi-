#! /bin/tcsh -f


#foreach i (hbook_*)
foreach i (hbook_MC0 hbook_MC10)

cd $i

./convert.csh

cd ..

end


