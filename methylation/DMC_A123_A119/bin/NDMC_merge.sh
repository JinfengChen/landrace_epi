#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

module load bedtools
#python Methykit_pipe.py --meta methykit.A123_A119.meta --cpu $PBS_NP
#cat A123_A119.Methykit.methykit.R | R --slave
project=A123_A119
#project=A123_A119.methykit
context=CHG

#compareR1=../../Methykit_A123_A119/A123_A119.Methykit.$context\.DMC.Q0.05.txt
#compareR2=../../Methykit_A123_A119R2/A123_A119R2.Methykit.$context\.DMC.Q0.05.txt
#compare0=../../Methykit_A119_A119/A119_A119.Methykit.$context\.DMC.Q0.05.txt

#compareR1=../../Fisher_A123_A119/MSU7.cytosine_table.fisher_test.$context\.P_adjusted_BH.txt.gz
#compareR2=../../Fisher_A123_A119R2/MSU7.cytosine_table.fisher_test.$context\.P_adjusted_BH.txt.gz
#compare0=../../Fisher_A119_A119/MSU7.cytosine_table.fisher_test.$context\.P_adjusted_BH.Q0.05.txt

zcat ../../Fisher_A123_A119/MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.txt.gz | awk '$21 > 0.1' > $project\.$context\.NDMC1.bed
zcat ../../Fisher_A123_A119R2/MSU7.cytosine_table.fisher_test.CG.P_adjusted_BH.txt.gz | awk '$21 > 0.1' > $project\.$context\.NDMC2.bed

bedtools intersect -f 1 -a $project\.$context\.NDMC1.bed -b $project\.$context\.NDMC1.bed | awk '$22 > 0' | cut -f1-3,8-9,18-19,22 > $project\.$context\.NDMC.bed
#bedtools intersect -f 1 -a $project\.overlap.bed -b $compare0 -v > $project\.DMC_Q0.05.core.bed
#bedtools intersect -f 1 -a $project\.overlap.bed -b $compare0 -v | cut -f1-3,8-9,18-19,22 > $project\.DMC_Q0.05.core.bed 

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

