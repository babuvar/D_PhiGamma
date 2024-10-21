#! /bin/tcsh -f

mkdir root

foreach i (1s*.hbk)

h2root $i

end

foreach j (1s*.root)

mv $j root/

end

cd root
hadd Merged.root 1s*root

