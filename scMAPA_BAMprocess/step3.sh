#!/bin/bash

#step3.sh will call modified DaPars2 to estimate the proximal PA site and the abundance of long/short isoforms for each cluster.
#The location of input files should be specified in the configuration.txt. The desired chromosome ID should be indicated in chrID.txt.
#Please load the following modules before running step3.sh
#python 3.7x
#./step3.sh configuration.txt chrID.txt

cfg_files=$1
chrID=$2
cat $chrID | while read chr; do
    echo "Starting processing chromosome $chr"
    echo $chr
    echo $cfg_files
    python Dapars2_Multi_Sample_abd.py $cfg_files $chr

done
