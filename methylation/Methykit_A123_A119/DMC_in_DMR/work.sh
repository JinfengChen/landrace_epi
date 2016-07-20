ln -s ../A123_A119.Methykit/Before_dedup_normalized/default/A123_A119.Methykit.CG.DMR.P0.05.txt ./
ln -s ../A123_A119.Methykit/Before_dedup_normalized/default/A123_A119.Methykit.CG.DMR.Q0.05.txt ./
ln -s ../../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.txt ./
bedtools intersect -f 1 -a A123_A119.Methykit.CG.DMR.P0.05.txt -b MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.txt -wao > A123_A119.Methykit.CG.DMR.P0.05.overlapped &
bedtools intersect -f 1 -a A123_A119.Methykit.CG.DMR.Q0.05.txt -b MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.txt -wao > A123_A119.Methykit.CG.DMR.Q0.05.overlapped &
awk '$18 >= 4' A123_A119.Methykit.CG.DMR.P0.05.overlapped > A123_A119.Methykit.CG.DMR.P0.05.overlapped_DMC4.txt
awk '$18 >= 4' A123_A119.Methykit.CG.DMR.Q0.05.overlapped > A123_A119.Methykit.CG.DMR.Q0.05.overlapped_DMC4.txt
ln -s ../../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.P0.05.Diff4mC.txt ./

#methykit DMC4
#116467
wc -l A123_A119.Methykit.CG.DMR.P0.05.overlapped_DMC4.txt
#112927
wc -l A123_A119.Methykit.CG.DMR.Q0.05.overlapped_DMC4.txt
#112927
bedtools intersect -f 1 -a A123_A119.Methykit.CG.DMR.P0.05.overlapped_DMC4.txt -b A123_A119.Methykit.CG.DMR.Q0.05.overlapped_DMC4.txt | awk '{print $1"_"$2}' | wc -l

#Fisher DMC4
#3546
bedtools intersect -f 1 -a Fisher_A123_A119R1/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.Q0.05.Diff4mC.txt -b Fisher_A123_A119R2/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.Q0.05.Diff4mC.txt | awk '{print $1"_"$2}' | uniq | sort | uniq | wc -l
#31596
bedtools intersect -f 1 -a Fisher_A123_A119R1/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.P0.05.Diff4mC.txt -b Fisher_A123_A119R2/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.P0.05.Diff4mC.txt | awk '{print $1"_"$2}' | uniq | sort | uniq | wc -l


#compared to A123 vs. A119R1
#39801
wc -l MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.P0.05.Diff4mC.txt
#4426
wc -l MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.Q0.05.Diff4mC.txt
#36250
bedtools intersect -f 1 -a A123_A119.Methykit.CG.DMR.P0.05.overlapped_DMC4.txt -b MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.P0.05.Diff4mC.txt| awk '{print $1"_"$2}' | uniq | sort | uniq | wc -l
#4426
bedtools intersect -f 1 -a A123_A119.Methykit.CG.DMR.P0.05.overlapped_DMC4.txt -b MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.Q0.05.Diff4mC.txt| awk '{print $1"_"$2}' | uniq | sort | uniq | wc -l

#compared to A123 vs. A119R2
#51873
wc -l MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.P0.05.Diff4mC.txt
#6638
wc -l MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.Q0.05.Diff4mC.txt
#33949
bedtools intersect -f 1 -a A123_A119.Methykit.CG.DMR.P0.05.overlapped_DMC4.txt -b MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.P0.05.Diff4mC.txt| awk '{print $1"_"$2}' | uniq | sort | uniq | wc -l
#5479
bedtools intersect -f 1 -a A123_A119.Methykit.CG.DMR.P0.05.overlapped_DMC4.txt -b MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.Q0.05.Diff4mC.txt| awk '{print $1"_"$2}' | uniq | sort | uniq | wc -l

#Methykit DMC4 vs. diff40
ln -s ../A123_A119.Methykit/Before_dedup_normalized/A123_A119.Methykit.CG.DMR.Q0.05.diff40.txt ./
#116467
wc -l A123_A119.Methykit.CG.DMR.P0.05.overlapped_DMC4.txt
#65129
wc -l A123_A119.Methykit.CG.DMR.Q0.05.diff40.txt
#36302
bedtools intersect -f 1 -a A123_A119.Methykit.CG.DMR.Q0.05.diff40.txt -b A123_A119.Methykit.CG.DMR.P0.05.overlapped_DMC4.txt | awk '{print $1"_"$2}' | uniq | sort | uniq | wc -l

