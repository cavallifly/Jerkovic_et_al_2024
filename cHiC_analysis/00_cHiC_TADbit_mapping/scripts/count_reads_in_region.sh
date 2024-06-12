
outFile=stats_TADbit_mapping_per_sample.tsv
if [[ ! -e ${outFile} ]];
#if [[ -e ${outFile} ]];
then
    echo "#sample total_reads uniquely_mapped valid_pairs" | awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' > ${outFile}    

    for validPairs in $(ls -1 ./*/03_*/valid*tsv) ;
    do
        check=$(grep $validPairs ${outFile} | wc -l | awk '{print $1}')
        if [[ $check -eq 1 ]];
	then
	    continue
	fi

	wDir=$(echo ${validPairs} | sed "s/03_/ /g" | awk '{print $1}')
	echo $validPairs
	echo $wDir
	if [[ ! -e ./${wDir}/TADbit_describe.log ]];
	then
	    tadbit describe -w ${wDir} > ./${wDir}/TADbit_describe.log 2> /dev/null
	fi
	#|  1 |     19 |        self-circle |      5,714 |    True |     4 |
	#|  2 |     20 |       dangling-end |    909,984 |    True |     4 |
	#|  3 |     21 |              error |      4,941 |    True |     4 |
	#|  4 |     22 | extra dangling-end | 11,944,046 |    True |     4 |
	#|  5 |     23 | too close from RES | 13,676,445 |   False |     4 |
	#|  6 |     24 |          too short |  2,439,216 |    True |     4 |
	#|  7 |     25 |          too large |        931 |    True |     4 |
	#|  8 |     26 |   over-represented | 17,263,633 |   False |     4 |
	#|  9 |     27 |         duplicated |  8,856,707 |    True |     4 |
	#| 10 |     28 |      random breaks |     94,332 |    True |     4 |
	#| 11 |     14 |        valid-pairs | 25,959,153 |         |     4 |
	totalReads=$(grep -w full ./${wDir}/TADbit_describe.log | head -1 | awk '{print $6}')
	uniquelyMapped=$(grep -A 4 INTERSECTION_OUTPUTs ./${wDir}/TADbit_describe.log | tail -1 | awk '{print $6}')
	vPairs=$(grep valid-pairs ./${wDir}/TADbit_describe.log | tail -1 | awk '{print $8}')    
	echo $validPairs ${totalReads} $uniquelyMapped ${vPairs} | awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' >> ${outFile}
    done
    
    wc -l stats_TADbit_mapping_*tsv
    sort stats_TADbit_mapping_per_sample.tsv | uniq | grep -v uniq > _tmp
    echo "#sample total_reads uniquely_mapped valid_pairs" | awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' > stats_TADbit_mapping_per_sample.tsv
    cat _tmp >> stats_TADbit_mapping_per_sample.tsv ; rm _tmp
    
    outFile=stats_TADbit_mapping_per_condition_per_replicate.tsv
    echo "#sample total_reads uniquely_mapped valid_pairs" | awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' > ${outFile}    
    sed -e "s/Rep1/Rep1 /g" -e "s/Rep2/Rep2 /g" -e "s/,//g" -e "s,\./,,g" -e "s,/,,g" stats_TADbit_mapping_per_sample.tsv | grep -v uniq | awk '{print $1,$3,$4,$5,$6,$7}' | awk '{tr[$1]+=$2; um[$1]+=$3; vp[$1]+=$4;}END{for(i in tr){printf("%s\t%s\t%s\t%s\n",i,tr[i],um[i],vp[i])}}' | sort -k 4,4n >> stats_TADbit_mapping_per_condition_per_replicate.tsv
    
    outFile=stats_TADbit_mapping_per_condition_merged.tsv
    echo "#sample total_reads uniquely_mapped valid_pairs" | awk '{printf("%s\t%s\t%s\t%s\n",$1,$2,$3,$4)}' > ${outFile}    
    sed -e "s/_Rep1/ /g" -e "s/_Rep2/ /g" stats_TADbit_mapping_per_condition_per_replicate.tsv | grep -v uniq | awk '{tr[$1]+=$2; um[$1]+=$3; vp[$1]+=$4;}END{for(i in tr){printf("%s\t%s\t%s\t%s\n",i,tr[i],um[i],vp[i])}}' | sort -k 4,4n >> stats_TADbit_mapping_per_condition_merged.tsv
fi
wc -l stats_TADbit_mapping_*tsv


for start in 52598845 ;
do	     
    for end in 57599408 ;
    do
	if [[ ${end} -le ${start} ]];
	then
	   continue
	fi
	#outFile=stats_contacts_TADbit_mapping_in_region_chr18_${start}_${end}.tsv	
	outFile=stats_contacts_TADbit_mapping_in_region_chr18_${start}_${end}_new.tsv	
	if [[ ! -e ${outFile} ]];
	then
	    echo "#sample valid_pairs cis_contacts trans_contacts cis_trans_contacts out_contacts" | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6)}' > ${outFile}
	fi
	sort ${outFile} | uniq | grep -v uniq > _tmp ; mv _tmp ${outFile}
	echo $outFile
	
	for validPairs in $(ls -1 ./*/03_*/valid*tsv | grep NPC_DA_) ;
	do
	    echo $validPairs
	    nlines=$(grep -w ${validPairs} ${outFile} | wc -l | awk '{print $1}')
	    if [[ $nlines -gt 0 ]];
	    then
		grep -w ${validPairs} ${outFile}
		continue
	    fi
	    
	    #GWNJ-0901:483:GW1907222321st:2:2124:9181:51746	chr10	3100130	1	150	3100130	3100670	chr1	142332698	0	150	142332695	142333992
	    #echo $ ch18:$start-$end $(awk -v s=${start} -v e=${end} 'BEGIN{c=0}{if($1=="chr18" && $4==$1 && (s<=$2 && $2<=e) && (s<=$5 && $5<=e)) c++}END{print c}' ${validPairs}) >> ${outFile}

	    allValidPairs=$(grep -v CRM ${validPairs} | wc -l | awk '{print $1}')	    

	    cisContacts=$(grep -v CRM ${validPairs} | awk -v s=${start} -v e=${end} 'BEGIN{c=0}{if($2=="chr18" && $8==$2 && (s<=$3 && $3<=e) && (s<=$9 && $9<=e)) c++}END{print c}')
	    transContacts=$(grep -v CRM ${validPairs} | awk -v s=${start} -v e=${end} 'BEGIN{c=0}{if($2=="chr18" && s<=$3 && $3<=e){if($8!="chr18" || s>$9 || $9>e) c++}; if($8=="chr18" && s<=$9 && $9<=e){if($2!="chr18" || s>$3 || $3>e) c++}}END{print c}')
	    cisTransContacts=$(grep -v CRM ${validPairs} | awk -v s=${start} -v e=${end} 'BEGIN{c=0}{if($2=="chr18" && s<=$3 && $3<=e){c++; next}; if($8=="chr18" && s<=$9 && $9<=e){c++; next}}END{print c}')
	    outContacts=$(grep -v CRM ${validPairs} | awk -v s=${start} -v e=${end} 'BEGIN{c=0}{if($2=="chr18" && s<=$3 && $3<=e){next}; if($8=="chr18" && s<=$9 && $9<=e){next}; c++}END{print c}')
	    echo $validPairs ${allValidPairs} ${cisContacts} ${transContacts} ${cisTransContacts} ${outContacts} >> ${outFile}

	done # Close cycle over $condition
	
	outFile=stats_contacts_TADbit_mapping_in_region_chr18_${start}_${end}_per_condition_per_replicate_new.tsv
	echo "#sample valid_pairs cis_contacts trans_contacts cis_trans_contacts out_contacts" | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6)}' > ${outFile}
	sed -e "s/_10/ /g"  -e "s/_2/ /g" -e "s/_3/ /g" -e "s/_4/ /g" -e "s/_5/ /g" -e "s/_6/ /g" -e "s/_7/ /g" -e "s/_8/ /g" -e "s/_9/ /g" -e "s/_1/ /g" -e "s/_Rep1/_Rep1 /g" -e "s/_Rep2/_Rep2 /g" -e "s,\./,,g" -e "s,/,,g" stats_contacts_TADbit_mapping_in_region_chr18_${start}_${end}_new.tsv | grep -v "#" | awk '{print $1,$4,$5,$6,$7,$8}' | awk '{all[$1]+=$2; cis[$1]+=$3; trans[$1]+=$4; cistrans[$1]+=$5; out[$1]+=$6}END{for(i in all){printf("%s\t%s\t%s\t%s\t%s\t%s\n",i,all[i],cis[i],trans[i],cistrans[i],out[i])}}' | sort -k 3,3n >> ${outFile}
	
	outFile=stats_contacts_TADbit_mapping_in_region_chr18_${start}_${end}_per_condition_merged_new.tsv
	echo "#sample valid_pairs cis_contacts trans_contacts cis_trans_contacts out_contacts" | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6)}' > ${outFile}
	sed -e "s/_Rep1/ /g" -e "s/_Rep2/ /g" -e "s,\./,,g" -e "s,/,,g" stats_contacts_TADbit_mapping_in_region_chr18_${start}_${end}_new.tsv | grep -v "#" | awk '{print $1,$3,$4,$5,$6,$7}' | awk '{all[$1]+=$2; cis[$1]+=$3; trans[$1]+=$4; cistrans[$1]+=$5; out[$1]+=$6}END{for(i in all){printf("%s\t%s\t%s\t%s\t%s\t%s\n",i,all[i],cis[i],trans[i],cistrans[i],out[i])}}' | sort -k 3,3n >> ${outFile}

    done # Close cycle over $end
done # Close cycle over $start    
tail stats_contacts*tsv    
wc -l stats_contacts*tsv
