#!/bin/bash

#SBATCH --job-name compFig
#SBATCH -n 1                # Number of cores. 
#SBATCH -t 1-00:00          # Runtime in D-HH:MM
#SBATCH -o 05_create_composite_figure.out # File of STDOUT
#SBATCH -e 05_create_composite_figure.out # File of STDERR

for dir in $1 ;
do
    
    bash ../scripts/05_create_composite_figure.sh ${dir}

done
wait

