chrSeqs=($(cat fileIDs.txt))

fileID=${chrSeqs[${SLURM_ARRAY_TASK_ID}]}

echo ${fileID}

module load gcc/8.2.0
module load bedtools/2.27.1
module load pybedtools/0.7.10

python pad3UTRs.py -r mm10_m18_refseq_gene_ensembl_extracted_3UTR.sorted.bed -i ${fileID}
bedtools genomecov -i ${fileID}.trans.bed -g GRCm38.p6.genome.chrom.sizes -bg > ${fileID}.trans.bedgraph
