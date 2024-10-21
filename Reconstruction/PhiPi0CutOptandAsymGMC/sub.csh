#! /bin/tcsh -f

foreach stream (0 1 2 3 4 5)
#foreach stream (1 2 5)
#foreach stream (5)

foreach i (charged charm mixed uds)
#foreach i (charged)

foreach j (11 13 15 17 19 21 23 25 27 7 9)
#foreach j (21)

bsub -q b_index ./script.csh $i 1$stream $j $stream

end

foreach k (31 33 35 37 39 41 43 45 47 49 51 55 61 63 65)
#foreach k (39)

bsub -q b_index ./script.csh $i $stream $k $stream
#./testscript.csh $i $stream $k $stream

end

end

end

