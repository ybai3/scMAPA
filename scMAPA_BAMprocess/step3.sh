#!/bin/bash

module load gcc/8.2.0
module load python/anaconda3.7-5.3.1_genomics
module load bedtools/2.27.1
module load r/3.5.1

#./step3.sh path_to_configuration_file_directory

cfg_files=$1/*.txt
for configure in $cfg_files; do 
    echo "Starting processing configuration file $configure"
	echo $configure
    python DaPars_main_3-5-19.py $configure

done
