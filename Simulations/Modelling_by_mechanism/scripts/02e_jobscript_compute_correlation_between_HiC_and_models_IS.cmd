#!/bin/bash

#SBATCH --job-name correlation
#SBATCH --mem=32Gb
#SBATCH -n 4                    # Number of cores.
#SBATCH -t 1-00:00              # Runtime in D-HH:MM
#SBATCH -o 02eSCC.out # File to which STDOUT will be written
#SBATCH -e 02e_SCC.out # File to which STDERR will be written 

cellType=XXXcellTypeXXX

for dir in $(ls -1 | grep XXXdirXXX);
do
    if [[ ! -d ${dir} ]];
    then
	echo "${dir} is not a directory!"
	continue
    fi
    
    check=$(ls -1 ${dir}/*.bed | wc -l)
    if [[ $check -lt 1 ]];
    then
	echo "${dir} has not P(s) computed. Please run step 03 first to compute it!"
	continue
    fi
    echo $dir
    
    bash ../scripts/02e_compute_correlation_between_HiC_and_models_IS.sh ${dir} ${cellType}
done
wait
