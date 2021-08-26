#!/bin/bash

# Python script called in step1.sh file was coded by Ming Tang, shared on https://hackmd.io/X2K8TpKeRJyb2ibK8PLzew?view
# Zhenjiang Fan modified it to a shell script.
# To run step1.sh, please first load Python with dependent modules.
# ./step1.sh path_to_cluster_file path_to_bam_file

python split_bam.py $1 $2

echo "Script executed from the following directory:"
echo ${PWD}
echo "Please copy it and use it in the second step."
