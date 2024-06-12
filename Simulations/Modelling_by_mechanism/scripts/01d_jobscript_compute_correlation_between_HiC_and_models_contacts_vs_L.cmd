#!/bin/bash

#SBATCH --job-name corrPs
#SBATCH -t 1-00:00              # Runtime in D-HH:MM
#SBATCH -o 01d_compute_correlation_Ps.out # File to which STDOUT will be written
#SBATCH -e 01d_compute_correlation_Ps.out # File to which STDERR will be written 



for cellType in XXXcellTypeXXX ;
do
    for dir in $(ls -1 | grep XXXdirXXX | head -1);    
    do
	if [[ ! -d ${dir} ]];
	then
	    echo "${dir} is not a directory!"
	    continue
	fi

	check=$(ls -1 ${dir}/C*.txt | wc -l)
	if [[ $check -lt 1 ]];
	then
	    echo "${dir} has not P(s) computed. Please run step 01c first to compute it!"
	    continue
	fi
	
	echo $dir	
	bash ../scripts/01d_compute_correlation_between_HiC_and_models_contacts_vs_L.sh ${dir} ${cellType}
    done
done
wait
