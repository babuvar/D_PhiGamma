#! /bin/tcsh -f


foreach j (D0toPhiPi0_100000 D0BtoPhiPi0_100000 )

bsub -q b_a ./script.csh $j > results.txt 
#./script.csh $j

end

