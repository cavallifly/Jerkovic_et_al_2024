#!/bin/bash

#SBATCH --job-name loop
#SBATCH -n 4                    # Number of cores.
#SBATCH -t 1-00:00              # Runtime in D-HH:MM
#SBATCH -o 03_loop_analysis.out # File to which STDOUT will be written
#SBATCH -e 03_loop_analysis.out # File to which STDERR will be written 

for condition in XXXcellTypeXXX ;
do
    for dir in $(ls -1 | grep XXXdirXXX);    
    do
	if [[ ! -d ${dir} ]];
	then
	    echo "${dir} is not a directory!"
	    continue
	fi
	echo $dir
	
	bash ../scripts/03a_loop_analysis.sh ${dir} ${condition}
	
    done # Close cycle over $dir
done # Close cycle over $condition
wait
