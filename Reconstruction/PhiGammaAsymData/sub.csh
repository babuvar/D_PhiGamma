#! /bin/tcsh -f

foreach i ( 7 9 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 55 61 63 65 )

foreach j ( on_resonance)

bsub -q b_index ./script.csh $i $j  

end
            
end



foreach i (43 53 67 69 71 )

foreach j ( 5S_onresonance)

bsub -q b_index ./script.csh $i $j

end

end



