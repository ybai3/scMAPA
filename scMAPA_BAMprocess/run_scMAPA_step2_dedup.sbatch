#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J PCRdup_rem
#SBATCH --output=PCRdup_rem.out
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.
#SBATCH --mem=64g # Memory pool for all cores (see also --mem-per-cpu)
##array should start from zero

module load gcc/8.2.0
module load samtools
module load bedtools/2.29.0
module load umi-tools/1.0.0

cd /ihome/hpark/yub20/NAR_mousebrain/BAMdedup
bam_files=./*.bam
for file in $bam_files; do
    echo "Startting indexing sorted bam file $file"
        BAM=$file

        bam_file_name=${BAM##*/}
        bam_name=$(echo "$bam_file_name" | cut -f 1 -d '.')

    SORT=$bam_name.sorted.bam
	DEDUP=$bam_name.dedup.bam
    BEDGRAPH=$bam_name.temp.bedgraph
    NEWBEDGRAPH=$bam_name.bedgraph

        echo $BAM
        echo $bam_name
        echo $SORT
		echo $DEDUP
        echo $BEDGRAPH
        echo $NEWBEDGRAPH

        samtools sort $BAM -o $SORT
		samtools index $SORT
		umi_tools dedup -I $SORT -S $DEDUP --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB
        bedtools genomecov -ibam $DEDUP -bg > $BEDGRAPH
        awk '{print "chr" $0}' $BEDGRAPH > $NEWBEDGRAPH
        rm $BEDGRAPH
done

