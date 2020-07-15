#!/bin/bash

module load gcc/8.2.0
module load samtools/1.9
module load bedtools/2.29.0

# ./step2.sh path_to_bam_directory path_to_output_bedgraph_directory

bam_files=$1/*.bam
format=*.bam
for file in $bam_files; do 
    echo "Starting processing clustered bam file $file"
	BAM=$file
	
	bam_file_name=${BAM##*/}
	bam_name=$(echo "$bam_file_name" | cut -f 1 -d '.')
	
    SORT=$bam_name.sorted.bam
    BEDGRAPH=$bam_name.bedgraph
    NEWBEDGRAPH=$2/$bam_name.bedgraph
	
	echo $BAM
	echo $bam_name
	echo $SORT
	echo $BEDGRAPH
	echo $NEWBEDGRAPH
	
	samtools sort $BAM -o $SORT
    bedtools genomecov -ibam $SORT -bg > $BEDGRAPH
    awk '{print "chr" $0}' $BEDGRAPH > $NEWBEDGRAPH
    rm $BEDGRAPH
done



