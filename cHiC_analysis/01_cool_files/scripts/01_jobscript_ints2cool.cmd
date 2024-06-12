#!/bin/bash
#SBATCH --job-name ints2cool 
#SBATCH -n 16                   # Number of cores. For now 56 is the number max of core available
#SBATCH --mem=10G              # allocated memory: yes, we need a lot of memory to generate .hic at 200bp resolution!
#SBATCH --partition=computepart # specify queue partiton
#SBATCH -t 1-00:00              # Runtime in D-HH:MM
#SBATCH -o ints2cool.out   # File to which STDOUT will be written
#SBATCH -o ints2cool.out   # File to which STDERR will be written

inFile=$1
outName=$(echo $inFile | sed -e "s/\.noDup\.ints//g" -e "s,files/, ,g" -e "s,/,,g" | awk '{print $2}')
echo $outName

assembly=mm10
author=IJ
outFile=chic_${outName%_ontarget}_${assembly}_${author}.cool
echo $inFile $outFile

if [[ ! -e ${outFile} ]];
then
    touch ${outFile}

    # Capture region
    start=52598845
    end=57599408

    # Capture region only at 100 bp
    cooler cload pairs --assembly ${assembly} --chrom1 1 --pos1 2 --chrom2 4 --pos2 5 ./scripts/chrom_sizes.txt:100 <( awk '{if($1=="chrom1"){next}; for(i=0;i<=$NF;i++){printf("%s\t%d\t%d\t%s\t%d\t%d\t1\n",$1,$2,$3,$4,$5,$6)}}' ${inFile}) ${outFile}
fi

outFile=chic_${outName%_ontarget}_${assembly}_${author}_20kb.cool
echo $inFile $outFile

if [[ ! -e ${outFile} ]];
then
    touch ${outFile}

    # Capture region
    start=52598845
    end=57599408
    
    # Entire chr18 at 20kb
    cooler cload pairs --assembly ${assembly} --chrom1 1 --pos1 2 --chrom2 4 --pos2 5 ./scripts/chrom_sizes.txt:20000 <( awk '{if($1=="chrom1"){next}; for(i=0;i<=$NF;i++){printf("%s\t%d\t%d\t%s\t%d\t%d\t1\n",$1,$2,$3,$4,$5,$6)}}' ${inFile}) ${outFile}
fi
