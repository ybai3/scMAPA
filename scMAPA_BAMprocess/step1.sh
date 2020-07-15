#!/bin/bash

module load gcc/8.2.0
#module load python/bioconda-3.7-2019.03
#conda acticate fanenv

#./step1.sh path_to_cluster_file path_to_bam_file

python split_bam.py $1 $2

echo "Script executed from the following directory:"
echo ${PWD}
echo "Please copy it and use it in the second step."
