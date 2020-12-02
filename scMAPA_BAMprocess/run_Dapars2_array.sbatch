#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH -J mb_DaPars2-%A-%a
#SBATCH --output=mb_DaPars2-%A-%a.out
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.
#SBATCH --mem=16g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero
#SBATCH --array=0-21

module load gcc/8.2.0
module load python/anaconda2.7-4.4.0_genomics
module load bedtools/2.27.1
module load r/3.5.1

echo $SLURM_ARRAY_TASK_ID
chrSeqs=($(cat chrIDs.txt))
echo $chrSeqs
chrID=${chrSeqs[${SLURM_ARRAY_TASK_ID}]}
#sampleID=$sampleSeqs
echo $chrID

configure=/ihome/hpark/yub20/NAR_mousebrain/BAMdedup/config.txt
echo $configure
python Dapars2_Multi_Sample_abd.py $configure $chrID
