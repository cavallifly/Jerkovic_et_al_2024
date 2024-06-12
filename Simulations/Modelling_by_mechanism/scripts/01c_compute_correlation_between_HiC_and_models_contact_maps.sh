scriptsdir=../../scripts/
mapResolution=5000
start=53700000

#pwd
condition=$2 
#pwd
cd $1
pwd
echo "Condition $condition"

rc=100
timestep=0
genDist=1.0

outCorr=correlations_contact_matrix.txt
#cat ${outCorr}
for condition in $2 ;
do
    for modMatrix in $(ls -1 rc_*.tab rc_*_10_rep*.tab rc_*_100_rep*.tab 2> /dev/null | grep -v Obs);
    do
	modMatrix0=${modMatrix}
	n=$(grep ${modMatrix0} ${outCorr} 2> /dev/null | wc -l | awk '{print $1}')
	if [[ $n -gt 0 ]];
	then
	    nlines=$(wc -l ${outCorr} | awk '{print $1}')
	    if [[ $nlines -gt 1 ]];
	    then
		echo "Correlation for ${modMatrix} DONE!"
		continue
	    fi
	fi
	
	rc=100 #$(echo $modMatrix | sed -e "s/rc_/ /g" -e "s/nm/ /g" | awk '{print $2}')
	from=$(echo $modMatrix | sed -e "s/_from_/ /g" -e "s/_to_/ /g" | awk '{print $2}')
	to=$(echo $modMatrix | sed -e "s/_every_/ /g" -e "s/_to_/ /g" | awk '{print $2}')
	nreplicas=$(echo $modMatrix | sed -e "s/_/ /g" | awk '{print $(NF-1)}')
	size=$(tail -1 ${modMatrix} | awk '{print $1}')
	
	# Compare with balanced
	cHiCMatrix="../../cHiCMatrices/chic_${condition}_merge_mm10_IJ_modelledRegion_chr18_53700000_56700000_at_5000bp_comparison.tab"
	ls -lrtha $cHiCMatrix
	
	name=$(echo $modMatrix | sed -e "s/\.tab//g" -e "s/NextrudersPerMb/N/g" -e "s/extrusionSpeed/eS/g" -e "s/loopyingBy/lB/g")

	awk -v size=$size '{if($1!=$2 && ($1!=size && $2!=size)) print $1-40,$2-40,$3}' ${modMatrix} | awk -v np=${nparticles} '{if(0<=$1 && $1<600 && 0<=$2 && $2<600) print $0}' > _tmp_mod
	modMatrix=_tmp_mod    
	
	echo $modMatrix0
	echo $modMatrix0  >> ${outCorr}

	awk -v size=$size '{if($1!=$2 && ($1!=size && $2!=size)) print $0}' ${cHiCMatrix} > _tmp_cHiC	
	expMatrix=_tmp_cHiC
	expFactor=$(awk '{d=sqrt(($1-$2)*($1-$2)); if(d==10){sum+=$3;cnt++}}END{print sum/cnt}' ${expMatrix})
	
	modFactor=$(awk '{d=sqrt(($1-$2)*($1-$2)); if(d==10){sum+=$3;cnt++}}END{print sum/cnt}' ${modMatrix})
	echo $size $modFactor $expFactor
	wc -l $expMatrix $modMatrix	
	    
	# Model
	diagOFF=1
	paste ${modMatrix} ${expMatrix} | grep -v Nan | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$3/fM; if($6=="NaN" || $6==0){next}; if(sqrt(($1-$2)*($1-$2))<diagOFF){next}; if($1< $2) print $1,$2,v}}' > _model_${timestep}
	# Experiment
	paste ${modMatrix} ${expMatrix} | grep -v Nan | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$6/fE; if($6=="NaN" || $6==0){next}; if(sqrt(($1-$2)*($1-$2))<diagOFF){next}; if($1< $2) print $1,$2,v}}' > _exp_${timestep}
	
        python ${scriptsdir}/01c_compute_correlation_between_HiC_and_models_contact_maps.py $((${size}-80)) ${rc} ${timestep} ${mapResolution} ${genDist} | awk -v c=$condition -v f=${from} -v t=${to} -v nr=${nreplicas} '{print $0,c,f,t,nr,"bal"}' >> ${outCorr}

	rm -fr _tmp* _exp_${timestep} _model_${timestep}
	
    done # Close cycle over $modMatrix
done # Close cycle over $condition
cd .. # Exit ${condition}

