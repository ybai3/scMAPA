#!/bin/bash

# step2_2.sh is a shell script that pads 3'tag reads, and converts BAM to BEDGRAPh.
# Please load following tools first.
# module load samtools
# module load bedtools
# module load pybedtools
# Please make sure following script/files are in the same directory.
# mm10_m18_refseq_gene_ensembl_extracted_3UTR.sorted.bed and GRCm38.p6.genome.chrom.sizes if sequencing samples are mouse.
# hg38_v27_refseq_gene_extracted_3UTR_sorted.bed and hg38.chrom.sizes if sequencing samples are human.
# example command
# ./step2.sh path_to_bam_directory 

bam_files=$1/*.bam
format=*.bam
for file in $bam_files; do
    echo "Startting indexing sorted bam file $file"
        BAM=$file

        bam_file_name=${BAM##*/}
        bam_name=$(echo "$bam_file_name" | cut -f 1 -d '.')

	SORT=$bam_name.sorted.bam
	DEDUP=$bam_name.dedup.bam
	DEDUPchr=$bam_name.dedup.chr.bam

    	echo $BAM
        echo $bam_name
        echo $SORT
	echo $DEDUP
	echo $DEDUPchr

        #samtools sort $BAM -o $SORT
	#samtools index $SORT
	#umi_tools dedup -I $SORT -S $DEDUP --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB
	#samtools view -H $DEDUP | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - $DEDUP > $DEDUPchr

	python pad3UTRs.py -r mm10_m18_refseq_gene_ensembl_extracted_3UTR.sorted.bed -i $DEDUPchr
	bedtools genomecov -i $DEDUPchr.trans.bed -g GRCm38.p6.genome.chrom.sizes -bg > $DEDUPchr.trans.bedgraph
done
