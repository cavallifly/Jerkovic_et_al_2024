#!/bin/bash

#SBATCH --job-name correlation
#SBATCH -t 1-00:00              # Runtime in D-HH:MM
#SBATCH -o SCC.out # File to which STDOUT will be written
#SBATCH -e SCC.out # File to which STDERR will be written 

for cellType in XXXcellTypeXXX ;
do
    for dir in $(ls -1 | grep XXXdirXXX | head -1);    
    do
	if [[ ! -d ${dir} ]];
	then
	    echo "${dir} is not a directory!"
	    continue
	fi

	bash ../scripts/01c_compute_correlation_between_HiC_and_models_contact_maps.sh ${dir} ${cellType}
    done # Close cycle over $dir
done # Close cycle over $cellType
wait
