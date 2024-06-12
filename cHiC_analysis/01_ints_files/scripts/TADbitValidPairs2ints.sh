chr=chr18
start=52598845
end=57599408

wDir=/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/cHiC_TADbit_mapping

samples1=$(ls -1 ${wDir} | grep Rep | sed "s/_Rep/ /g" | awk '{print $1}' | sort | uniq | grep -v TAD) ; echo $samples1
samples2=$(ls -1 ${wDir} | grep Rep | sed "s/_/ /g" | awk '{print $1"_"$2"_"$3}' | sort | uniq | grep -v TAD) ; echo $samples2 # .ints per replicate
#exit

for sample in ${samples2} ${samples1} ;
do
    echo $sample
    sampleName=${sample}
    outLog=TADbitValidPairs2ints_${sampleName}.out
    nRep=$(echo $sampleName | sed 's/Rep/ /g' | awk '{print $2}')
    if [[ $nRep == "" ]];
    then
	sampleName=${sampleName}_merge
    fi
    echo $sample $sampleName $nRep

    echo "# Capture region only"
    #V300085140L3C003R0690845074:0:0:0:0~2~	chr10	3100130	1	50	3100130	3100670	chr18	37382008	0	150	37381617	37382507
    outFile=${sampleName}_ontarget.noDup.ints
    ls -lrtha ${wDir}/${sample}_*/03_filtered_reads/valid_r1-r2_intersection_*.tsv 2> /dev/null > ${outLog}
    echo "# We consider the starting nucleotide of the restriction fragments at columns 6 and 12" >> ${outLog}
    if [[ ! -e ${outFile} ]];
    then
	touch ${outFile}
	echo "chrom1 start1 end1 chrom2 start2 end2 obs" | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}' > ${outFile}
	awk -v c=${chr} -v s=${start} -v e=${end} '{if($2==c && $8==$2 && (s<=$6 && $6<=e) && (s<=$12 && $12<=e))printf("%s\t%d\t%d\t%s\t%d\t%d\n",$2,$6,$6+1,$8,$12,$12+1)}' ${wDir}/${sample}_*/03_filtered_reads/valid_r1-r2_intersection_*.tsv | awk '{h[$0]++}END{for(i in h){printf("%s\t%d\n",i,h[i])}}' | sort -k 2,2n -k 4,4n >> ${outFile}
    fi
    if [[ -e ${outFile} ]];
    then
	ls -lrtha ${outFile} >> ${outLog}
    fi
    echo ${outFile} $(awk 'BEGIN{c=0}{if($1=="chrom1"){next};c+=$NF}END{printf("cis-contacts %f\n",c)}' ${outFile}) >> ${outLog}

    echo "# Entire chr18"
    outFile=${sampleName}_onchr18.noDup.ints
    if [[ ! -e ${outFile} ]];
    then
	touch ${outFile}
	echo "chrom1 start1 end1 chrom2 start2 end2 obs" | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}' > ${outFile}	
	awk -v c=${chr} '{if($2==c && $8==$2)printf("%s\t%d\t%d\t%s\t%d\t%d\n",$2,$6,$6+1,$8,$12,$12+1)}' ${wDir}/${sample}*/03_filtered_reads/valid_r1-r2_intersection_*.tsv | awk '{h[$0]++}END{for(i in h){printf("%s\t%d\n",i,h[i])}}' >> ${outFile}
    fi
    if [[ -e ${outFile} ]];
    then
	ls -lrtha ${outFile} >> ${outLog}
    fi
    echo ${outFile} $(awk 'BEGIN{c=0}{if($1=="chrom1"){next};c+=$NF}END{printf("cis-contacts %f\n",c)}' ${outFile}) >> ${outLog}    
    continue
    
    echo "# Entire chr18 no Capture"
    outFile=${sampleName}_onchr18NoCapture.noDup.ints
    if [[ ! -e ${outFile} ]];
    then
	touch ${outFile}
	echo "chrom1 start1 end1 chrom2 start2 end2 obs" | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}' > ${outFile}	
	# Get all contacts on chr18 of target
	awk -v c=${chr} -v s=${start} -v e=${end} '{if($2==c && $8==$2){if((s<=$6 && $6<=e) || (s<=$12 && $12<=e)){next}else{printf("%s\t%d\t%d\t%s\t%d\t%d\n",$2,$6,$6+1,$8,$12,$12+1)}}}' ${wDir}/${sample}_*/03_filtered_reads/valid_r1-r2_intersection_*.tsv | awk '{h[$0]++}END{for(i in h){printf("%s\t%d\n",i,h[i])}}' | sort -k 2,2n -k 4,4n >> ${outFile}
    fi
    if [[ -e ${outFile} ]];
    then
	ls -lrtha ${outFile} >> ${outLog}
    fi
    awk 'BEGIN{c=0}{if($1=="chrom1"){next};c+=$NF}END{printf("cis-contacts %f\n",c)}' ${outFile} >> ${outLog}
    
    echo "# Entire chr18 with CaptureDS"
    outFile=${sampleName}_onchr18CaptureDS.noDup.ints
    if [[ ! -e ${outFile} ]];
    then
	touch ${outFile}
	start=50000000
	end=60000000
	grep ${sampleName} ../01_cool_files/onD*
	grep ${sampleName} ../01_cool_files/offD*
	nContactsON=$( grep ${sampleName} ../01_cool_files/onD*  | grep NoC | awk '{print int($4*2)}')
	nContactsOFF=$(grep ${sampleName} ../01_cool_files/offD* | grep NoC | awk '{print int($4*2)}')
	echo ${nContactsON} ${nContactsOFF}
	echo "chrom1 start1 end1 chrom2 start2 end2 obs" | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}' > ${outFile}	
	echo "# Get all contacts on chr18 not in the capture region"
	awk -v c=${chr} -v s=${start} -v e=${end} '{if($2==c && $8==$2){if((s<=$6 && $6<=e) || (s<=$12 && $12<=e)){next}else{printf("%s\t%d\t%d\t%s\t%d\t%d\n",$2,$6,$6+1,$8,$12,$12+1)}}}' ${wDir}/${sample}_*/03_filtered_reads/valid_r1-r2_intersection_*.tsv | awk '{h[$0]++}END{for(i in h){printf("%s\t%d\n",i,h[i])}}' >> ${outFile}
	echo # Down-sample all contacts within the capture region to ${nContactsON}"
	awk -v c=${chr} -v s=${start} -v e=${end} '{if($2==c && $8==$2 && ((s<=$6 && $6<=e) || (s<=$12 && $12<=e))){printf("%s\t%d\t%d\t%s\t%d\t%d\n",$2,$6,$6+1,$8,$12,$12+1)}}' ${wDir}/${sample}_*/03_filtered_reads/valid_r1-r2_intersection_*.tsv | shuf -n ${nContactsON} | awk '{h[$0]++}END{for(i in h){printf("%s\t%d\n",i,h[i])}}' >> ${outFile}
	echo # Down-sample all contacts of the capture region with other parts of chromosome 18 to ${nContactsOFF}"
	awk -v c=${chr} -v s=${start} -v e=${end} '{if($2==c && $8==$2 && ((s<=$6 && $6<=e) || (s<=$12 && $12<=e))){next}else{if($2==c && $8==$2){if((s<=$6 && $6<=e)){printf("%s\t%d\t%d\t%s\t%d\t%d\n",$2,$6,$6+1,$8,$12,$12+1)}; if((s<=$12 && $12<=e)){printf("%s\t%d\t%d\t%s\t%d\t%d\n",$2,$6,$6+1,$8,$12,$12+1)}}}}' ${wDir}/${sample}_*/03_filtered_reads/valid_r1-r2_intersection_*.tsv | shuf -n ${nContactsOFF} | awk '{h[$0]++}END{for(i in h){printf("%s\t%d\n",i,h[i])}}' >> ${outFile}	
	cat ${outFile} | sort -k 2,2n -k 4,4n > _tmp_ints_${outName} ; mv _tmp_ints_${outName} ${outFile}
    fi
    if [[ -e ${outFile} ]];
    then
	ls -lrtha ${outFile} >> ${outLog}
    fi
    awk 'BEGIN{c=0}{if($1=="chrom1"){next};c+=$NF}END{printf("cis-contacts %f\n",c)}' ${outFile} >> ${outLog}

done
