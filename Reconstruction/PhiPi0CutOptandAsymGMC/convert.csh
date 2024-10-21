#! /bin/tcsh -f


foreach i (*.hbk)

h2root $i

end

foreach j (*.root)

mv $j root/

end

cd root
hadd Merged_GMC.root *root

