#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l mem=60gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

module load bedtools
python Methykit_pipe.py --meta methykit.A123_A119.meta --cpu $PBS_NP

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

