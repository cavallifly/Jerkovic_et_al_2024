for condition in "mESC_DAB" "mESC_DA" "mESC_DB" "mESC_dCTCF" "mESC_Dprom" "mESC_polyA" "mESC_UPPax6" "mESC_WT" "NPC_DAB" "NPC_DA" "NPC_DB" "NPC_dCTCF" "NPC_Dprom" "NPC_polyA" "NPC_WT";
do
    echo $condition
    for resolution in 2000 5000 ;
    do
	echo $resolution
	inFile=insulationProfiles_${condition}_r${resolution}bp.tab	
	for window in $(seq 150000 50000 500001)
	do
	    echo $window
	    w=$(echo $window     | awk '{print int($1/1000)}')
	    r=$(echo $resolution | awk '{print int($1/1000)}')	    
	    outFile=chic_${condition}_merge_mm10_IJ_w${w}kb_r${r}kb_complete.tab
	    tmpFile=_tmp_${outFile}
	    if [[ ! -e ${outFile} ]];
	    then
		touch ${outFile}		
		#head -1 ${inFile}
		head -1 ${inFile} | awk '{for(i=1 ; i<=NF; i++){print $i,i}}' | grep $window | awk 'BEGIN{for(i=1 ; i<=3; i++) print i; print 5}{print $2}' > ${tmpFile}
		#grep -w chr18 ${inFile} | grep -w 52708000
		awk '{if(NF==1){c[$1]=1}else{for(i=1;i<=NF;i++){if(c[i]==1)printf("%s\t",$i)};printf("\n")}}' ${tmpFile} ${inFile} > ${outFile}	    
		#| grep -v NaN | grep -w 52708000 
		rm ${tmpFile}
	    fi
	    outFile=chic_${condition}_merge_mm10_IJ_w${w}kb_r${r}kb.tab
	    tmpFile=_tmp_${outFile}	    
	    if [[ ! -e ${outFile} ]];
	    then
		touch ${outFile}		
		head -1 ${inFile} | awk '{for(i=1 ; i<=NF; i++){print $i,i}}' | grep $window | grep log | awk 'BEGIN{for(i=1 ; i<=3; i++) print i}{print $2}' > ${tmpFile}
		echo "chrom start end ins" | awk '{for(i=1 ; i<NF; i++){printf("%s\t",$i)}; printf("%s\n",$NF)}' > ${outFile}
		awk '{if(NF==1){c[$1]=1}else{for(i=1;i<=NF;i++){if(c[i]==1)printf("%s\t",$i)};printf("\n")}}' ${tmpFile} <(grep -v log ${inFile}) >> ${outFile}	    
		rm ${tmpFile}
	    fi
	done
    done
done
