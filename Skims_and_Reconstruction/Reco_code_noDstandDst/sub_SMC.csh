#! /bin/tcsh -f



#foreach i (D0toPhiPi0_1M D0BtoPhiPi0_1M)               #10 directories
#foreach i (D0toPhiGamma_500K D0BtoPhiGamma_500K)        # 5 directories
#foreach i (D0toPhiGamma_71000 D0BtoPhiGamma_71000)        # 5 directories
#foreach i (D0toPhiGamma_114650 D0BtoPhiGamma_114650)        # 5 directories
#foreach i (D0toPhiPi0_161520 D0BtoPhiPi0_161520)        # 6 directories
#foreach i (D0toPhiGamma_80760 D0BtoPhiGamma_80760)        # 5 directories
#foreach i (DMtoPiZPiM_50000_Y4s DPtoPiZPiP_50000_Y4s)        # 11 directories
#foreach i (DMtoPiZPiM_1M_Y4s DPtoPiZPiP_1M_Y4s)        # 12 directories
foreach i (DtoPiZPi_2M_Y4s)        # 23 directories



#foreach j (1 2 3 4 5 6 7 8 9 10)
#foreach j (1 2 3 4 5)
#foreach j (1 2 3 4 5 6)
#foreach j (1 2 3 4 5 6 7 8 9 10 11)
#foreach j (1 2 3 4 5 6 7 8 9 10 11 12)
foreach j (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23)


bsub -q s ./script_SMC.csh $i $j
#./script_SMC.csh $i $j


end

end




