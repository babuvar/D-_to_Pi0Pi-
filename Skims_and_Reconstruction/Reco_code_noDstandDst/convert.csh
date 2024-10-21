#! /bin/tcsh -f

mkdir root

foreach i (*.hbk)

h2root $i

end

foreach j (*.root)

mv $j root/

end

