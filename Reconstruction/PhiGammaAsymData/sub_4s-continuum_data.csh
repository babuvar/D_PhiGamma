#! /bin/tcsh -f

foreach i ( 7 11 13 15 17 19 23 25 27 31 33 35 37 39 41 43 45 47 49 51 55 61 63 65 67 69 71)

bsub -q b_index ./script_4s-continuum_data.csh $i 

end

