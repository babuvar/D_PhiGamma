#! /bin/tcsh -f



#foreach i (D0toPhiPi0_1M D0BtoPhiPi0_1M)               #10 directories
#foreach i (D0toPhiGamma_500K D0BtoPhiGamma_500K)        # 5 directories
#foreach i (D0toPhiGamma_71000 D0BtoPhiGamma_71000)        # 5 directories
#foreach i (D0toPhiGamma_114650 D0BtoPhiGamma_114650)        # 5 directories
#foreach i (D0toPhiPi0_161520 D0BtoPhiPi0_161520)        # 6 directories
#foreach i (D0toPhiGamma_80760 D0BtoPhiGamma_80760)        # 5 directories
#foreach i (D0toPhiGamma_30K_2 D0BtoPhiGamma_30K_2 )        # 5 directories



#foreach i (DtoPhiGamma_100K_Y5s D0BtoPhiGamma_100K_Y5s)        # 4 directories
#foreach i (DtoPhiGamma_100K_Y4s D0BtoPhiGamma_100K_Y4s)        # 4 directories

#foreach i (D0toPhiGamma_500K_Y5s D0BtoPhiGamma_500K_Y5s)        # 11 directories
foreach i (D0toPhiPi0_1M_Y5s D0BtoPhiPi0_1M_Y5s)        # 15 directories

#foreach j (1 2 3 4 5 6 7 8 9 10)
#foreach j (1 2 3 4 5)
#foreach j (1 2 3 4 5 6)

#foreach j (1 2 3 4)

#foreach j (1 2 3 4 5 6 7 8 9 10 11)
foreach j (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)


bsub -q b_a ./script_SMC.csh $i $j


end

end




