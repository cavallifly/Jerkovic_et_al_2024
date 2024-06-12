ls -1 ${PWD}/../01_cool_files/*Rep*20kb.cool > _coolFiles
nFiles=$(wc -l _coolFiles | awk '{print $1}')

n=0
for resolution in 20000 ; #10000 5000 ;
do    
    for lbr in 0 1000000 2000000 3000000 4000000 ;
    do
	for ubr in 0 1000000 2000000 3000000 4000000 ;    
	do
	    if [[ $lbr -ge $ubr ]];
	    then
		continue
	    fi
	    outFile=01_hicRep_analysis_captureRegion_from_${lbr}bp_to_${ubr}bp_at_${resolution}bp_coolFiles.out
	    if [[ -e ${outFile} ]];
	    then
		sed -e "s,/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/02_hicRep_analysis/../01_cool_files/chic_,,g" -e "s,_mm10_IJ.cool,,g" ${outFile} > _tmp ; mv _tmp ${outFile}
		continue
	    fi
	    for nFile1 in $(seq 1 1 ${nFiles})
	    do
		file1=$(awk -v n=${nFile1} '{if(NR==n) print $1}' _coolFiles)
		for nFile2 in $(seq $((${nFile1}+1)) 1 ${nFiles})
		do
		    file2=$(awk -v n=${nFile2} '{if(NR==n) print $1}' _coolFiles)
		    n=$((${n}+1))
		    echo $file1 $file2 $n
		    
		    sed -e "s,XXXhic1XXX,${file1},g" -e "s,XXXhic2XXX,${file2},g" -e "s,XXXresolutionXXX,${resolution},g" -e "s,XXXlbrXXX,${lbr},g" -e "s,XXXubrXXX,${ubr},g" ./scripts/01_hicRep_analysis_from_coolFiles.R > _tmp${n}.R
		    Rscript _tmp${n}.R &
		    
		    if [[ ${n} -eq 60 ]];
		    then
			n=0
			wait
			rm -fvr _tmp*.R
		    fi		    
		done
	    done
	    wait	    
	    sed -e "s,/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/02_hicRep_analysis/../01_cool_files/chic_,,g" -e "s,_mm10_IJ.cool,,g" ${outFile} > _tmp ; mv _tmp ${outFile}
	done
    done
done
rm -fvr _tmp*.R    
exit
for file in $(ls -1 *out | grep coolFiles);
do
    sed -e "s,/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/02_hicRep_analysis/../01_cool_files/chic_,,g" -e "s,_mm10_IJ.cool,,g" ${file} > _tmp ; mv _tmp ${file}    
done
