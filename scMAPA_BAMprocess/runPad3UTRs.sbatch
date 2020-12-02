#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yub20@pitt.edu
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J pad3UTRs-%A-%a
#SBATCH --output=pad3UTRs-%A-%a.out
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.
#SBATCH --mem=8g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-4
##array should start from zero

#echo $SLURM_ARRAY_TASK_ID
chrSeqs=($(cat fileIDs.txt))
#echo $chrSeqs
fileID=${chrSeqs[${SLURM_ARRAY_TASK_ID}]}
#sampleID=$sampleSeqs
echo ${fileID}

module load gcc/8.2.0
module load bedtools/2.27.1
module load pybedtools/0.7.10

python pad3UTRs.py -r mm10_m18_refseq_gene_ensembl_extracted_3UTR.sorted.bed -i ${fileID}
bedtools genomecov -i ${fileID}.trans.bed -g GRCm38.p6.genome.chrom.sizes -bg > ${fileID}.trans.bedgraph
