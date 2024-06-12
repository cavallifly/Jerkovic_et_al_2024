scriptsdir=../../../scripts/
mapResolution=5000
start=53700000

cellTypes=$2
cd $1

rc=100
timestep=27200000
genDist=0.5

outCorr=correlations_loop_contacts.txt

for modMatrix in $(ls -1 rc*.tab 2> /dev/null) ;
do
    echo $modMatrix
    size=$(tail -1 ${modMatrix} | awk '{print $1}')
    echo "Size $size"

    modMatrix0=${modMatrix}

    name=$(echo $modMatrix0 | sed -e "s/\.tab//g")
    outFile=loops_${name}.tab

    tmpDir=_tmp_${name}_loop
    if [[ -d ${tmpDir} ]];
    then
	continue
    fi
    mkdir $tmpDir
    cd ${tmpDir}
    
    awk -v size=$size '{if($1!=$2 && ($1!=size && $2!=size)) print $1-40,$2-40,$3}' ../${modMatrix} | awk '{if(0<=$1 && $1<600 && 0<=$2 && $2<600) print $0}' > _tmp_mod
    modMatrix=_tmp_mod

    for cellType in ${cellTypes};
    do
	cellType=$(echo $cellType | sed -e "s/_/ /g" | awk '{print $1"_"$2}')
	
	if [[ ! -e ../${outFile} ]];
	then

	    cHiCMatrix="../../../cHiCMatrices/chic_${cellType}_merge_mm10_IJ_modelledRegion_chr18_53700000_56700000_at_5000bp_comparison.tab"
	    awk -v size=$size '{if($1!=$2 && ($1!=size && $2!=size)) print $0}' ${cHiCMatrix} > _tmp_cHiC
	    expMatrix=_tmp_cHiC
	    expFactor=$(awk '{d=sqrt(($1-$2)*($1-$2)); if(d==10){sum+=$3;cnt++}}END{print sum/cnt}' ${expMatrix})
	    
	    modFactor=$(awk '{d=sqrt(($1-$2)*($1-$2)); if(d==10){sum+=$3;cnt++}}END{print sum/cnt}' ${modMatrix})
	    echo $size $modFactor $expFactor
	    maxRank=$(wc -l $expMatrix | awk '{print $1}')
	    wc -l $expMatrix $modMatrix
	    
	    diagOff=5
	    ### Compute ranks
	    # Model vs Hi-C
	    #Observed -> general ranking -> check which ranks have the same expected and assign to all of those the same ranking!
	    paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk '{print $1,$2,$3,NR,$4,$5,$6}' | awk '{if(r[$3]==""){r[$3]=$4}; print $1,$2,r[$3],$5,$6,$7}' | awk -v m=$maxRank '{if($6=="NaN" || $6==0){print $1,$2,-100,$4,$5,$6}else{print $1,$2,($3-1)/(m-1)*200-100,$4,$5,$6}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$3/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fM}; if($1<$2){print $1,$2,v}}' >   contact_map.tab
	    paste ${modMatrix} ${expMatrix} | sort -k 6,6n | awk '{print $1,$2,$3,$4,$5,$6,NR}' | awk '{if(r[$6]==""){r[$6]=$7}; print $1,$2,$3,$4,$5,r[$6]}' | awk -v m=$maxRank '{if($6=="NaN" || $6==0){print $1,$2,-100,$4,$5,$6}else{print $1,$2,$3,$4,$5,($6-1)/(m-1)*200-100}}' | sort -k 1,1n -k 2,2n | awk -v fM=1 -v fE=1 -v diagOFF=${diagOFF} '{if($1!=$4 || $2!=$5){next}; v=$6/fM; if(sqrt(($1-$2)*($1-$2))<diagOFF || $6=="NaN"){v=0/fM}; if($1<$2){print $2,$1,v}}' >>  contact_map.tab	

	    awk 'BEGIN{for(i=1;i<=600;i++){for(j=1;j<=600;j++){m[i,j]=0}}}{ind1=int($1); ind2=int($2); if(NF==2){m[ind1,ind2]=1;m[ind2,ind1]=1;}else{tag="Models";if(ind1>ind2){tag="CHi-C"};if(m[ind1,ind2]==1){printf("%s\t%s\t%s\t%s\n",$1,$2,$3,tag)}}}' ../../scripts/loops_${cellType}_mustache_and_manualAnnot.txt contact_map.tab > ${outFile}
	    cp ${outFile} ../
	    
	    #n=$(grep ${modMatrix0} ../${outCorr} 2> /dev/null | wc -l | awk '{print $1}')
            #if [[ $n -eq 0 ]];
            #then

		
		#grep Models ${outFile} | awk '{print $1,$2,$3}' > _tmp_mod
		#modMatrix=_tmp_mod
		#modFactor=1		
		#echo $modMatrix0

		#grep CHi ${outFile} | awk '{print $2,$1,$3}' > _tmp_cHiC
		#expMatrix=_tmp_cHiC

		#expFactor=1

		#echo $size $modFactor $expFactor
		#wc -l $expMatrix $modMatrix
		
		#diagOFF=1
		#paste ${modMatrix} ${expMatrix} | grep -iv Nan | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$3/fM; if($6=="NaN"){next}; if(sqrt(($1-$2)*($1-$2))<diagOFF){next}; print $1,$2,v}}' > _model_${timestep}
		# Experiment
		#paste ${modMatrix} ${expMatrix} | grep -iv Nan | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$6/fE; if($6=="NaN"){next}; if(sqrt(($1-$2)*($1-$2))<diagOFF){next}; print $1,$2,v}}' > _exp_${timestep}
		#paste _model_${timestep} _exp_${timestep} | grep -iv NaN | head
		
		#echo $modMatrix0 >> ${outCorr}
		#python ${scriptsdir}/04_compute_correlation_between_HiC_and_models_contact_maps_genotoul.py ${size} ${rc} ${timestep} ${mapResolution} ${genDist} | awk -v c=$condition '{if($1=="All") print $0,c," balanced"}' >> ${outCorr}
            #fi

	    
	fi
	if [[ ! -e ../${outFile%tab}pdf ]];
	then
	    awk '{if($4=="Exp"){print $2+41,$1+41,$3,$4}else{print $1+41,$2+41,$3,$4}}' ../${outFile} > _tmp.tab
	    head _tmp.tab
	    Rscript ${scriptsdir}/03a_loop_analysis.R -100
	    mv clipboard.tab ../${outFile%tab}txt
	    dos2unix ${outFile%tab}txt
	    mv _tmp.pdf ../${outFile%tab}pdf
	fi

	#cat ${outCorr} >> ../${outCorr}
	
    done # Close cycle over $modMatrix

    cd ../ # Exit $tmpDir
    rm -fvr $tmpDir
    
done # Close cycle over $cellType
rm -fr contact_map.tab contact_map_diff.tab _tmp_mod _tmp_cHiC _exp_rank _mod_rank distance.txt _tmp.tab _model_${timestep} _exp_${timestep} data.txt
	
cd .. # Exit ${cellTypes}
