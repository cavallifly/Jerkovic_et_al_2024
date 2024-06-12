#!/bin/bash

### mustache analysis is DONE on .mcool files
for mcoolFile in $( ls -1 *.mcool);
do

    #For very high resolutions (e.g., 1kb) use:
    #https://github.com/ay-lab/mustache/issues/13
    # mustache -f input.mcool -ch chr3R -r 1000 -o output.tsv -p 16 -pt .05 -cz dm6.chrom.sizes -st 0.7
    # Let me know how it works/worked for 1kb. I would also try -st 08-0.88 (especially when I am getting loops that don't make sense).

    outDir=${mcoolFile%.mcool}_default
    #if [[ -d ${outDir} ]]
    #then
	#	continue
    #fi 
    mkdir -p ${outDir}
    for resolution in 5000 10000 20000 ; 
    do
	if [[ ! -e ${outDir}/${mcoolFile%.mcool}_${resolution}bp.tsv ]];
	then
	    mustache -f ${mcoolFile} -r ${resolution} -p 8 -o ${outDir}/${mcoolFile%.mcool}_${resolution}bp.tsv > mustache_analysis_${mcoolFile%.mcool}.out &
        fi
    done # Close cycle over $resolution
    wait
done # Close cycle over $mcoolFile
