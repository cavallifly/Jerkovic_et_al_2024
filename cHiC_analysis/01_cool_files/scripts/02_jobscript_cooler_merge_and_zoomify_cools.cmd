#!/bin/bash
#SBATCH --job-name MergeZoomify 
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem=16000             # allocated memory per CPU
#SBATCH --partition=computepart # specify queue partiton
#SBATCH -t 10-00:00              # Runtime in D-HH:MM
#SBATCH -o cooler_merge_and_zoomify_cools.out   # File to which STDOUT will be written
#SBATCH -e cooler_merge_and_zoomify_cools.out   # File to which STDERR will be written


coolFile=$1

#mcoolFile=${coolFile%.cool}_5Mb.mcool
mcoolFile=${coolFile%.cool}.mcool

if [[ ! -e ${mcoolFile} ]];
then
    touch ${outFile}
    cooler zoomify -r 2000,5000,10000,15000,20000,40000,100000 ${coolFile} -o ${coolFile%.cool}.mcool --balance &> 02_cool2mcool_${coolFile%.cool}.out
    #cooler zoomify -r 20000,40000,100000,5000000 ${coolFile} -o ${coolFile%.cool}.mcool --balance &> 02_cool2mcool_${coolFile%.cool}.out
    #cooler zoomify -r 5000000 ${coolFile} -o ${mcoolFile} --balance &> 02_cool2mcool_${coolFile%.cool}.out         
fi

