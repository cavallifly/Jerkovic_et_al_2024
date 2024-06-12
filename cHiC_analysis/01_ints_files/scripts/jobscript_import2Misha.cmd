#!/bin/bash

#SBATCH --job-name import2misha
#SBATCH -n 2                # Number of cores. For now 56 is the number max of core available
#SBATCH --mem=10G           # Memory per cpu ti allocate
#SBATCH -t 0-20:00          # Runtime in D-HH:MM
#SBATCH -o import2misha.log   # File to which STDOUT will be written
#SBATCH -e import2misha.log   # File to which STDERR will be written 

sample=$1
refGenome=mm10

echo "### GENERATE .hic : START ###";

outName=${sample}

if [[ ! -e ${outName}.noDup.ints ]];
then
    # GWNJ-0901:483:GW1907222321st:2:1213:11140:2012	chr10	3100287	-	chr18	29879118	-	286	HIC_chr10_2	HIC_chr18_67228	42	42
    cat ${outName}.allValidPairs | awk 'BEGIN{printf("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tobs\n")}{printf("%s\t%s\t%s\t%s\t%s\t%s\t1\n",$2,$3,$3+1,$5,$6,$6+1)}' > ${outName}.noDup.ints
fi

echo "### IMPORT TO MISHA : START ###";

#dbFolder="/media/data/mishaDB/trackdb/$refGenome/";
#mDBloc="/media/data/mishaDB/trackdb/"
dbFolder="/zdata/data/mishaDB/trackdb/$refGenome/";
mDBloc="/zdata/data/mishaDB/trackdb/"
deDupFile=${sample}.noDup.ints
outName=chic_${sample%_ontarget}_mm10_IJ
#outName=chic_${sample%_ontarget}_merge_mm10_BB
echo $deDupFile $outName

rScript="importToMisha_${sample}.R"

echo "#!/usr/bin/Rscript"  >  ${rScript}
echo >> ${rScript}
echo "mDBloc <- \"${mDBloc}\"" >> ${rScript}
echo >> ${rScript}
echo "require('misha', quietly=TRUE)"  >> ${rScript}
echo "require('shaman', quietly=TRUE)" >> ${rScript}
echo >> ${rScript}
echo "db <- '$dbFolder'" >> ${rScript}
echo "gsetroot(db)"      >> ${rScript}
echo "gdb.reload()"      >> ${rScript}
echo >> ${rScript}
echo "files <- '$deDupFile'"   >> ${rScript}
echo "trackName <- '$outName'" >> ${rScript}
echo "hicName <- 'chic'"     >> ${rScript}
echo "outTrack <- paste0(hicName,'.',trackName)" >> ${rScript}
echo "message(outTrack)"       >> ${rScript}
echo >> ${rScript}
echo "gtrack.2d.import_contacts(outTrack, paste(\"import chic 2d track for\", hicName), files)" >> ${rScript}
#echo "print(paste0('DONE'))" >> ${rScript}

Rscript importToMisha_${sample}.R &>> import2misha_${sample}.log
