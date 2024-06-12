condition=$(cat _condition)
inFile=$(ls -1 correlation_analysis_Ea_${condition}.log)
nLines=$(cat ${inFile} | wc -l | awk '{print $1}')
echo $nLines
ts=15
condition=$(cat _condition)
outFile=tested_parameters_${condition}.txt
rm -fvr ${outFile}

for nl in $(seq 4 1 ${nLines});
do
    name=$(awk -v nl=${nl} '{if(NR==nl){print $1}}' ${inFile} | sed "s/_from/ /g"| awk '{print $1}')
    Le=$(paste <(awk '{print $1}' compute_lifetime_distribution.out | sed "s/\./_/g") <(awk '{for(i=2;i<=NF;i++){printf(" %s",$i)}; printf("\n")}' compute_lifetime_distribution.out) | grep $name | grep -i aver | grep -i real | awk '{printf("%.0f",$5/60)}')
    if [[ $Le == "" ]];
    then
	Le=NA
    fi
    name=$(awk -v nl=${nl} '{if(NR==nl){print $1}}' ${inFile})    
    corrs=$(grep ${name} ${inFile} | awk '{for(i=2;i<=NF;i++){printf "%s ",$i}}')

    echo $name | sed "s/_/ /g" | sed -e "s/AA//g" -e "s/BB//g" -e "s/AB//g" -e "s/fc//g"
    echo $(echo $name | sed "s/_/ /g" | sed -e "s/AA//g" -e "s/BB//g" -e "s/AB//g" -e "s/fc//g" | awk -v ts=$ts -v Le=${Le} '{print $2"."$3,$4"."$5,$6"."$7,$9"."$10,1,$14,Le,$18"."$19,$20"."$21,$22"."$23,$24"."$25,$26"."$27,$28"."$29,$44,$45,$46,$47,$48,$49,$51,$52,$53,$54,$55,$56,$NF*0.006/ts/60}') $corrs >> ${outFile}
    tail -1 ${outFile}

done
	  
