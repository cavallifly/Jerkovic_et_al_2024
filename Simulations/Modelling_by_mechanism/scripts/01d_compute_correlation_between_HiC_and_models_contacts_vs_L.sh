scriptsdir=../../scripts/
mapResolution=5000

cd $1
pwd

timestep=NA

outDiff=differences_contacts_vs_L.txt

for condition in $2 ;
do
    for modMatrix in $(ls -1 C_*rc_*rep*.txt 2> /dev/null);
    do
	modMatrix0=${modMatrix}

	n=$(grep ${modMatrix0} ${outDiff} 2> /dev/null | wc -l | awk '{print $1}')
	if [[ $n -gt 0 ]];
	then
	    echo "Difference for ${modMatrix} DONE!"		
	    continue
	fi
	
	echo $modMatrix | sed -e "s/rc_/ /g" -e "s/nm/ /g" | awk '{print $2}'
	
	size=$(tail -1 ${modMatrix} | awk '{print $1+1}')
	
	# Compare with balanced
	cHiCMatrix="../../cHiCMatrices/average_contacts_chic_${condition}_merge_mm10_IJ_modelledRegion_chr18_53700000_56700000_at_${mapResolution}bp_comparison.tab"
	ls -lrtha $cHiCMatrix
	head $cHiCMatrix

	name=$(echo $modMatrix | sed -e "s/\.txt//g" -e "s/NextrudersPerMb/N/g" -e "s/extrusionSpeed/eS/g" -e "s/loopyingBy/lB/g" -e "s/C_vs_L/Ps/g")
	
	awk '{if($1!=0) print $1,0,$2}' ${modMatrix} > _tmp_mod
	modMatrix=_tmp_mod
	    
	awk '{if($1!=0) print $1,0,$2}' ${cHiCMatrix} > _tmp_cHiC	
	expMatrix=_tmp_cHiC
	expFactor=1
	    
	modFactor=1
	echo $size $modFactor $expFactor
	wc -l $expMatrix $modMatrix	
	    
	# Model
	diagOFF=1
	paste ${modMatrix} ${expMatrix} | grep -v Nan | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$3/fM; if($6=="NaN"){next}; if(sqrt(($1-$2)*($1-$2))<diagOFF){next}; print $1,$2,v}}' > _model_${timestep}
	# Experiment
	paste ${modMatrix} ${expMatrix} | grep -v Nan | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$6/fE; if($6=="NaN"){next}; if(sqrt(($1-$2)*($1-$2))<diagOFF){next}; print $1,$2,v}}' > _exp_${timestep}
	paste _model_${timestep} _exp_${timestep} | grep -v NaN | head
	    
	### Diff ###
	n=$(grep ${modMatrix0} ${outDiff} 2> /dev/null | wc -l | awk '{print $1}')
	if [[ $n -eq 0 ]];
	then
	    echo $modMatrix0 >> ${outDiff}
	    paste _model_${timestep} _exp_${timestep} | grep -v NaN | awk '{printf("%d\t%.10f\t%d\t%.10f\n",$1,$3,$4,$6)}'
	    paste _model_${timestep} _exp_${timestep} | grep -v NaN | awk '{printf("%d\t%.10f\t%d\t%.10f\n",$1,$3,$4,$6)}' | awk '{if($1==$3 && $1>4){dp=($2-$4); d+=sqrt(dp*dp); t+=sqrt($4*$4)}}END{print (t-d)/t}' | awk -v c=$condition '{print "All",$0,c,"balanced"}' >> ${outDiff} 
	fi

	rm -fr _tmp* _exp_${timestep} _model_${timestep}
	
    done # Close cycle over $modMatrix
done
cd .. # Exit ${condition}

