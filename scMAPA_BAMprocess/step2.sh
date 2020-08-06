#!/bin/bash

# Step2.sh is a shell script that converts splitted BAM to BEDGRAPh.
# Zhenjiang Fan modified it to a shell script.
# Please load samtools and bedtools first.
# ./step2.sh path_to_bam_directory 

bam_files=$1/*.bam
format=*.bam
for file in $bam_files; do 
    echo "Starting processing clustered bam file $file"
	BAM=$file
	
	bam_file_name=${BAM##*/}
	bam_name=$(echo "$bam_file_name" | cut -f 1 -d '.')
	
    SORT=$bam_name.sorted.bam
    BEDGRAPH=$bam_name.temp.bedgraph
    NEWBEDGRAPH=$bam_name.bedgraph
	
	echo $BAM
	echo $bam_name
	echo $SORT
	echo $BEDGRAPH
	echo $NEWBEDGRAPH
	
	samtools sort $BAM -o $SORT
        bedtools genomecov -ibam $SORT -bg > $BEDGRAPH
        awk '{print "chr" $0}' $BEDGRAPH > $NEWBEDGRAPH
        rm $BEDGRAPH
	rm $SORT
done



