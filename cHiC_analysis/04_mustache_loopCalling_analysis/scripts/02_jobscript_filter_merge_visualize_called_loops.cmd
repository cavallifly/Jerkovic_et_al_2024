#!/bin/bash

# Capture region
s1=52775574
s2=57266687
resolutions="5000 10000 20000"

for dir in $(ls -1 | grep default)
do
    if [[ ! -d $dir ]];
    then
	continue
    fi
    echo $dir
    cd $dir

    outFile=Loop_chr18_${s1}_${s2}_chr18_${s1}_${s2}_filtered_merged.tsv
    outLog=filter_merge_visualize_called_loops.out
    
    sample=$(echo $dir | sed -e "s/chic_//g" -e "s/_merge_mm10_IJ_default//g")
    echo $sample > ${outLog}

    loopFile=_loops.tsv
    ls -lrtha *.tsv >> ${outLog}
    
    if [[ ! -e ${outFile} ]];
    then
	tooCloseBad=40000
	echo "Removing too close (<${tooCloseBad}bp) to bad columns" >> ${outLog}
	rm -fr _loopFile
	for resolution in ${resolutions} ;
	do
	    badColumnsFile=../../../01_cool_files/badColumns/MADmax_badColumns_chic_${sample}_merge_mm10_IJ_at_${resolution}bp_resolution.txt
	    awk -v m=${tooCloseBad} '{if(NF==4){for(i=$2-m;i<=$3+m;i++){bad[i]=1;}}else{for(i in bad){if(($2<=i && i<=$3) || ($5<=i && i<=$6)){next}}; print $0}}' <(sed -e "s,(,,g" -e "s,),,g" -e "s/,//g" ${badColumnsFile}) <(awk '{if(NR>1) printf("%s\t%d\t%d\t%s\t%d\t%d\n",$1,$2,$3,$4,$5,$6)}' *_${resolution}bp.tsv) >> _loopFile
	done
	echo "Remaining loops $(wc -l _loopFile | awk '{print $1}')" >> ${outLog}
	
	fromLimits=40000
	echo "Remove loops too close (<${fromLimits}bp) to the limits of the capture regions" >> ${outLog}
	awk -v m=${fromLimits} -v s1=${s1} -v s2=$s2 '{if(NR==1){print $0; next}; l1=s1+m; l2=s2-m; if((l1<$2 && $2<l2) && (l1<$3 && $3<l2) && (l1<$5 && $5<l2) && (l1<$6 && $6<l2)){print $0}}' _loopFile > _tmp ; mv _tmp _loopFile   
	echo "Remaining loops $(wc -l _loopFile | awk '{print $1}')" >> ${outLog}
	
	tooClose=40000
	echo "Remove loops with too close (<${tooClose}bp) anchors" >> ${outLog}
	awk -v m=${tooClose} '{if(NR==1){print $0; next}; a1=($2+$3)/2; a2=($5+$6)/2; dp=a1-a2; d=sqrt(dp*dp); if(d>m){print $0};}' _loopFile > _tmp ; mv _tmp _loopFile
	echo "Remaining loops $(wc -l _loopFile | awk '{print $1}')" >> ${outLog}
	
	echo "Merge overlapping loops" >> ${outLog}
	sort -k 2,2n -k 5,5n _loopFile | awk 'BEGIN{printf("%s\t%s\t%s\t%s\t%s\t%s\n","chrom1","start1","end1","chrom2","start2","end2")}{if(NR==1){n=0; c1[n]=$1; s1[n]=$2; e1[n]=$3; c2[n]=$4; s2[n]=$5; e2[n]=$6}else{for(i=0;i<=n;i++){a11=0; a22=0; for(j=$2;j<=$3;j++){if(s1[i]<=j && j<=e1[i]){a11=1;break}}; for(k=$5;k<=$6;k++){if(s2[i]<=k && k<=e2[i]){a22=1;break}}; if(a11==1 && a22==1){if($2<s1[i]){s1[i]=$2}; if(e1[i]<$3){e1[i]=$3}; if($5<s2[i]){s2[i]=$5}; if(e2[i]<$6){e2[i]=$6}; break}}; if(a11==0 || a22==0){n++; c1[n]=$1; s1[n]=$2; e1[n]=$3; c2[n]=$4; s2[n]=$5; e2[n]=$6}}}END{for(i in c1){printf("%s\t%d\t%d\t%s\t%d\t%d\n",c1[i],s1[i],e1[i],c2[i],s2[i],e2[i])}}' | sort -k 2,2n -k 5,5n > _tmp ; mv _tmp _loopFile
	echo "Remaining loops $(grep -v chrom _loopFile | wc -l | awk '{print $1}')" >> ${outLog}
	cat _loopFile >> ${outLog}
	
	cp _loopFile ${outFile}
    fi
    
    for resolution in ${resolutions} ;
    do
	outFile=Loop_chr18_${s1}_${s2}_chr18_${s1}_${s2}_merged_with_loops_at_${resolution}bp.png
	if [[ -e ${outFile} ]];
	then
  	    ls -lrtha $outFile
	    continue
	fi

	if [[ ${resolution} -eq 5000 ]];
	then
	    vmax=0.010
	fi	
	if [[ ${resolution} -eq 10000 ]];
	then
	    vmax=0.010
	fi
	if [[ ${resolution} -eq 20000 ]];
	then
	    vmax=0.020
	fi
	
	python ../scripts/02_visualize_called_loops_v1.py ../$(echo ${dir} | sed -e "s/_default//g" -e "s/_pt/ /g" | awk '{print $1}').mcool _loopFile ../captureRegion.bed ${resolution} ${vmax} &> /dev/null		    
	
	mv -v Loop_chr18_${s1}_${s2}_chr18_${s1}_${s2}_at_${resolution}bp_with_loops.png ${outFile} >> ${outLog}
    done # Close cycle over $resolution			
    rm -v _*
    cd ..
done # Close cycle over $dir
