module load gcc/8.2.0
module load samtools
module load bedtools/2.29.0
module load umi-tools/1.0.0

bam_files=./*.bam
for file in $bam_files; do
    echo "Startting indexing sorted bam file $file"
        BAM=$file

        bam_file_name=${BAM##*/}
        bam_name=$(echo "$bam_file_name" | cut -f 1 -d '.')

    SORT=$bam_name.sorted.bam
    DEDUP=$bam_name.dedup.bam

    	echo $BAM
        echo $bam_name
        echo $SORT
	echo $DEDUP

        samtools sort $BAM -o $SORT
	samtools index $SORT
	umi_tools dedup -I $SORT -S $DEDUP --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB
done

