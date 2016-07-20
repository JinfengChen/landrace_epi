echo "Comparision DMR"
cut -f2- A123_A119.Methykit/default/A123_A119.Methykit.CG.DMR.txt | awk '$6<0.05 && ($7>40 || $7<-40)' > A123_A119.Methykit/default/A123_A119.Methykit.CG.DMR.diff40.txt
#6700
bedtools intersect -a A123_A119.Methykit/default/A123_A119.Methykit.CG.DMR.diff40.txt -b ../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.DMR_window.Chr2.txt | wc -l
#11091
wc -l A123_A119.Methykit/default/A123_A119.Methykit.CG.DMR.diff40.txt
#28651
wc -l ../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.DMR_window.Chr2.txt

#1915
bedtools intersect -a A123_A119.Methykit/filtered_normalized/A123_A119.Methykit.CG.DMR.diff40.txt -b ../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.DMR_window.Chr2.txt| cut -f2| uniq | wc -l
#2380
wc -l A123_A119.Methykit/filtered_normalized/A123_A119.Methykit.CG.DMR.diff40.txt


echo "Comparsion DMC between methykit and fisher before deduplication"
#1028231
cut -f2- A123_A119.Methykit.CG.DMC.txt | awk '$5<=0.05' > A123_A119.Methykit.CG.DMC.P0.05.txt
wc -l A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.P0.05.txt
#308368
cut -f2- A123_A119.Methykit.CG.DMC.txt | awk '$6<=0.05' > A123_A119.Methykit.CG.DMC.Q0.05.txt
wc -l A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.Q0.05.txt

#in Fisher_A123_A119
#700646
awk '$21<=0.05' MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.txt > MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.P0.05.txt
#144099
awk '$22<=0.05' MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.txt > MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.Q0.05.txt

#methykit compared to Fisher_A123_A119
#639348
bedtools intersect -f 1 -a A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.P0.05.txt -b ../Fisher_A123_A119/MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.P0.05.txt | wc -l
#141222
bedtools intersect -f 1 -a A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.Q0.05.txt -b ../Fisher_A123_A119/MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.Q0.05.txt | wc -l

#in Lulu_A123_A119
#351812
cat A123_A119_new_old_CG.qvalue | awk '$13 <= 0.05 && $14 <= 0.05 && $15 >= 0.05' | awk '{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$13}' > A123_A119_new_old_CG.qvalue.P0.05
wc -l A123_A119_new_old_CG.qvalue.P0.05

#methykit compared to Lulu_A123_A119
#337181
bedtools intersect -a A123_A119.Methykit.CG.DMC.P0.05.txt -b ../Lulu_A123_A119/A123_A119_new_old_CG.qvalue.P0.05 | wc -l
#253381
bedtools intersect -a A123_A119.Methykit.CG.DMC.Q0.05.txt -b ../Lulu_A123_A119/A123_A119_new_old_CG.qvalue.P0.05 | wc -l

echo "Comparision DMR between methykit and fisher before deduplication"
#1727736
cut -f2- A123_A119.Methykit.CG.DMR.txt | awk '$5<=0.05' > A123_A119.Methykit.CG.DMR.P0.05.txt
#65129
cut -f2- A123_A119.Methykit.CG.DMR.txt | awk '$5<=0.05 && ($7>40 || $7<-40)' > A123_A119.Methykit.CG.DMR.P0.05.diff40.txt

#1539868
cut -f2- A123_A119.Methykit.CG.DMR.txt | awk '$6<=0.05' > A123_A119.Methykit.CG.DMR.Q0.05.txt
#65129
cut -f2- A123_A119.Methykit.CG.DMR.txt | awk '$6<=0.05 && ($7>40 || $7<-40)' > A123_A119.Methykit.CG.DMR.Q0.05.diff40.txt

#1539868, some overlapping region are three, almost close to each other
bedtools intersect -f 0.9 -a A123_A119.Methykit.CG.DMR.P0.05.txt -b A123_A119.Methykit.CG.DMR.Q0.05.txt | awk '{print $1":"$2}' | uniq | sort | uniq | wc -l


#in Fisher_A123_A119
#622992
awk '$9<=0.05' MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.txt > MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.P0.05.txt
#347508
awk '$9<=0.05 && ($5/($5+$6)-$7/($7+$8) >= 0.4 || $5/($5+$6)-$7/($7+$8) <= -0.4)' MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.txt > MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.P0.05.diff40.txt
#350726
awk '$10<=0.05' MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.txt > MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Q0.05.txt
#255667
awk '$10<=0.05 && ($5/($5+$6)-$7/($7+$8) >= 0.4 || $5/($5+$6)-$7/($7+$8) <= -0.4)' MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.txt > MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Q0.05.diff40.txt

#methykit compared to Fisher_A123_A119
#165945
bedtools intersect -f 0.9 -a A123_A119.Methykit.CG.DMR.P0.05.txt -b ../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.P0.05.txt| awk '{print $1":"$2}' | uniq | sort | uniq | wc -l
#70191
bedtools intersect -f 0.9 -a A123_A119.Methykit.CG.DMR.P0.05.txt -b ../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Q0.05.txt| awk '{print $1":"$2}' | uniq | sort | uniq | wc -l
#414478
bedtools intersect -f 0.5 -a A123_A119.Methykit.CG.DMR.P0.05.txt -b ../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.P0.05.txt| awk '{print $1":"$2}' | uniq | sort | uniq | wc -l
#183486
bedtools intersect -f 0.5 -a A123_A119.Methykit.CG.DMR.P0.05.txt -b ../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Q0.05.txt| awk '{print $1":"$2}' | uniq | sort | uniq | wc -l
#diff methylation level by 0.4
#27765
bedtools intersect -f 0.9 -a A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMR.P0.05.diff40.txt -b ../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.P0.05.diff40.txt | awk '{print $1":"$2}' | uniq | sort | uniq | wc -l
#22696
bedtools intersect -f 0.9 -a A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMR.Q0.05.diff40.txt -b ../Fisher_A123_A119/MSU7_w200_s50.fisher_test.CG.P_adjusted_BH.Filter_by_DMC.P0.05.Diff4mC.txt | awk '{print $1":"$2}' | uniq | sort | uniq | wc -l

#methykit normalized vs. not normalized, DMC
#Pvalue
#not normalized,701673 
wc -l A123_A119.Methykit/Before_dedup_normalized/default/A123_A119.Methykit.CG.DMC.P0.05.txt
#normalized, 1028231, more sensitive
wc -l A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.P0.05.txt
#compare, 639502 
bedtools intersect -f 1 -a A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.P0.05.txt -b A123_A119.Methykit/Before_dedup_normalized/default/A123_A119.Methykit.CG.DMC.P0.05.txt | wc -l
#Qvalue
#not normalized, 179535
wc -l A123_A119.Methykit/Before_dedup_normalized/default/A123_A119.Methykit.CG.DMC.Q0.05.txt
#normalized, 308368
wc -l A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.Q0.05.txt
#compare, 175518
bedtools intersect -f 1 -a A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.Q0.05.txt -b A123_A119.Methykit/Before_dedup_normalized/default/A123_A119.Methykit.CG.DMC.Q0.05.txt | wc -l

#Fisher_A123_A119
#700646
wc -l ../Fisher_A123_A119/MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.P0.05.txt
#144099
wc -l ../Fisher_A123_A119/MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.Q0.05.txt
#compare
#639348
bedtools intersect -f 1 -a A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.P0.05.txt -b ../Fisher_A123_A119/MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.P0.05.txt | wc -l
#141222
bedtools intersect -f 1 -a A123_A119.Methykit/Before_dedup_normalized/Filer_normalized/A123_A119.Methykit.CG.DMC.Q0.05.txt -b ../Fisher_A123_A119/MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.Q0.05.txt | wc -l


echo "Methykit use mC reads and C reads in all the sample to do test, not number of mC and number of C"
echo "Either use if fisher or methykit, that's two different strategies"
echo "Fisher is a conserved and need to consider difference between replicate"

