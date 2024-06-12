#!/bin/bash
#SBATCH --job-name scoring
#SBATCH -n 5                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem=400G              # allocated memory
#SBATCH --partition=computepart # specify queue partiton
#SBATCH -t 10-00:00             # Runtime in D-HH:MM
#SBATCH -o 08_shamanScoring.out   # File to which STDOUT will be written
#SBATCH -e 08_shamanScoring.out   # File to which STDERR will be written

sample=$1

Rscript ./scripts/08_computeScoreTrack_parallel.R ${sample} &>> 08_shamanScoring_${sample}.out


