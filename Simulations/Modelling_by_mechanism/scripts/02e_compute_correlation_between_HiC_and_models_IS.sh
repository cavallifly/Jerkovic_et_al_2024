scriptsdir=../../../scripts/
mapResolution=5000
start=53700000
end=56700000

rc=100

cd $1

outCorr=differences_IS.txt

tag=$1

for modMatrix in $(ls -1 *rc_*.bed 2> /dev/null | grep ${dir}_ | grep 400kb | grep -v Obs);
do      

    timestepF=$(echo $modMatrix | sed "s/_/ /g" | awk '{print $(NF-5)}')
    timestepI=$(echo $modMatrix | sed "s/_/ /g" | awk '{print $(NF-7)}')    

    condition=$2

    echo $modMatrix
    echo $modMatrix | sed -e "s/rc_/ /g" -e "s/nm/ /g" | awk '{print $2}'
    
    size=$(wc -l ${modMatrix} | awk '{print $1+1}')
    
    # Compare with balanced
    cHiCMatrix=../../insulation_profiles/insulationProfiles_${condition}_w400000bp_r5000bp_for_modelling.tab

    name=$(echo $modMatrix | sed -e "s/\.bed//g" -e "s/NextrudersPerMb/N/g" -e "s/extrusionSpeed/eS/g" -e "s/loopyingBy/lB/g")
    check=$(grep ${name} ${outCorr} | awk -v tsI=${timestepI} -v tsF=${timestepF} '{if(($(NF-2)==tsI) && ($(NF-1)==tsF)){print $0}}' | wc -l | awk '{print $1}')
    if [[ $check -gt 0 ]];
    then
	continue
    fi    
    
    ### Difference ###
    border=0
    awk -v s=${start} -v e=${end} -v b=${border} '{if($2<s+b){next}; if($2>e-b){next}; printf("%s\t%d\t%d\t%f\n",$1,$2,$3,-$6)}' ${cHiCMatrix} > _tmp_${tag}_exp
    head _tmp_${tag}_exp
    minExp=$(sort -k 4,4n _tmp_${tag}_exp | head -1 | awk '{print $4}')
    maxExp=$(sort -k 4,4n _tmp_${tag}_exp | tail -1 | awk '{print $4}')
    echo $minExp $maxExp
    awk -v max=${maxExp} -v min=${minExp} '{print $1,$2,$3,($4-min)/(max-min)}' _tmp_${tag}_exp > _tmp_${tag} ; mv _tmp_${tag} _tmp_${tag}_exp
    minExp=$(sort -k 4,4n _tmp_${tag}_exp | head -1 | awk '{print $4}')
    maxExp=$(sort -k 4,4n _tmp_${tag}_exp | tail -1 | awk '{print $4}')
    echo $minExp $maxExp
    awk -v s=${start} -v e=${end} -v b=${border} '{if($2<s+b){next}; if($2>e-b){next}; print $0}' _tmp_${tag}_exp | awk '{if($2-pp>10000){print ""}; print $0; pp=$2}' > _tmp_${tag} ; mv _tmp_${tag} _tmp_${tag}_exp
    head _tmp_${tag}_exp
    
    awk -v s=${start} -v e=${end} -v b=${border} '{if($2<s+b){next}; if($2>e-b){next}; print $0}' ${modMatrix} > _tmp_${tag}_mod
    minMod=$(sort -k 4,4n _tmp_${tag}_mod | head -1 | awk '{print $4}')
    maxMod=$(sort -k 4,4n _tmp_${tag}_mod | tail -1 | awk '{print $4}')
    echo $minMod $maxMod
    awk -v max=${maxMod} -v min=${minMod} '{print $1,$2,$3,($4-min)/(max-min)}' _tmp_${tag}_mod > _tmp_${tag} ; mv _tmp_${tag} _tmp_${tag}_mod
    minMod=$(sort -k 4,4n _tmp_${tag}_mod | head -1 | awk '{print $4}')
    maxMod=$(sort -k 4,4n _tmp_${tag}_mod | tail -1 | awk '{print $4}')
    echo $minMod $maxMod
    awk -v s=${start} -v e=${end} -v b=${border} '{if($2<s+b){next}; if($2>e-b){next}; print $0}' _tmp_${tag}_mod | awk '{if($2-pp>10000){print ""}; print $0; pp=$2}' > _tmp_${tag} ; mv _tmp_${tag} _tmp_${tag}_mod
    
    # Compare only around the peaks
    echo "9 53 243 258 370 529 573" | awk '{for(i=1;i<=NF;i++) print $i}' > _barriers_${tag}
    if [[ $condition == "mESC_WT" ]];
    then
	echo "10 55 148 323 529 573" | awk '{for(i=1;i<=NF;i++) print $i}' > _barriers_${tag}
    fi	
    cat _barriers_${tag}
    paste _tmp_${tag}_mod _tmp_${tag}_exp > _data
    cat _barriers_${tag} _data | awk '{if(NF==1){for(i=0;i<3;i++){nr[$1-i]=1;nr[$1+i]=1;}}else{n++; if(nr[n]==1){print $0}}}' | awk -v c=$condition -v tsF=${timestepF} -v tsI=${timestepI} -v rc=${rc} -v n=${name} 'BEGIN{d=0; t=0}{dd=$4-$8; d+=sqrt(dd*dd); t++}END{print n,"matrix dcutoff",rc,"Diff",(t-d)/t,"npoints",tsI,tsF,c}' >> ${outCorr}

    rm -fvr _barriers_${tag} _data _tmp_${tag}_mod _tmp_${tag}_exp _exp_${timestep} _model_${timestep}
    
done # Close cycle over $modMatrix
cd .. # Exit ${condition}
