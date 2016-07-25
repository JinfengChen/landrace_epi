echo "working DMC dataset for A123_A119"
bash DMC_merge.sh
bash NDMC_merge.sh

echo "SNP distribution"
awk '$4<$5' MSU_r7.all.final.full.utr.gff > MSU_r7.all.final.full.utr.corrected.gff
cat MSU_r7.all.final.full.utr.corrected.gff MSU_r7.fa.RepeatMasker.out.gff > MSU_r7.gene_TE.gff 
#CG
python SNP_distr.py --bed A123_A119.CG.DMC_Q0.05.core.bed --gff ../input/MSU_r7.gene_TE.gff --output A123_A119.CG.DMC
python SNP_distr.py --bed A123_A119.CG.NDMC.bed --gff ../input/MSU_r7.gene_TE.gff --output A123_A119.CG.NDMC
#CHG
python SNP_distr.py --bed A123_A119.CHG.DMC_Q0.05.core.bed --gff ../input/MSU_r7.gene_TE.gff --output A123_A119.CHG.DMC
python SNP_distr.py --bed A123_A119.CHG.NDMC.bed --gff ../input/MSU_r7.gene_TE.gff --output A123_A119.CHG.NDMC
#CHH
python SNP_distr.py --bed A123_A119.CHH.DMC_Q0.05.core.bed --gff ../input/MSU_r7.gene_TE.gff --output A123_A119.CHH.DMC
python SNP_distr.py --bed A123_A119.CHH.NDMC.bed --gff ../input/MSU_r7.gene_TE.gff --output A123_A119.CHH.NDMC
