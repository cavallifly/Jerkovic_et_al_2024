#!/bin/bash
#SBATCH --job-name shuffling
#SBATCH -n 2                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem=250G             # allocated memory: yes, we need a lot of memory to generate .hic at 200bp resolution!
##SBATCH --mem=128G              # allocated memory: yes, we need a lot of memory to generate .hic at 200bp resolution!
#SBATCH --partition=computepart # specify queue partiton
#SBATCH -t 10-00:00             # Runtime in D-HH:MM
#SBATCH -o 07_shamanNormalization.out   # File to which STDOUT will be written
#SBATCH -e 07_shamanNormalization.out   # File to which STDERR will be written

sample=$1

Rscript ./scripts/07_shamanNormalization.R ${sample} &>> 07_shamanNormalization_${sample}.out
