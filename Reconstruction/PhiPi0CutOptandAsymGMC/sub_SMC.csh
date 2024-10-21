#! /bin/tcsh -f



foreach i (D0toPhiPi0_1M D0BtoPhiPi0_1M)               #10 directories
#foreach i (D0toPhiGamma_500K D0BtoPhiGamma_500K)        # 5 directories
#foreach i (D0toPhiGamma_71000 D0BtoPhiGamma_71000)        # 5 directories
#foreach i (D0toPhiGamma_114650 D0BtoPhiGamma_114650)        # 5 directories
#foreach i (D0toPhiPi0_161520 D0BtoPhiPi0_161520)        # 6 directories
#foreach i (D0toPhiGamma_80760 D0BtoPhiGamma_80760)        # 5 directories


foreach j (1 2 3 4 5 6 7 8 9 10)
#foreach j (1 2 3 4 5)
#foreach j (1 2 3 4 5 6)


bsub -q b_a ./script_SMC.csh $i $j


end

end




