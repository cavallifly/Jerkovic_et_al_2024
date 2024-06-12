# modeledRegion
chr=chr18
start=53700000
end=56700000
resolution=5000

size=$(awk -v start=${start} -v end=$end -v resolution=$resolution 'BEGIN{print (end-start)/resolution}')
echo "Size $size"

conditions="mESC_DAB mESC_DA mESC_DB mESC_dCTCF mESC_Dprom mESC_polyA mESC_UPPax6 mESC_WT NPC_DAB NPC_DA NPC_DB NPC_dCTCF NPC_Dprom NPC_polyA NPC_WT"

for condition in ${conditions} ;
do
    inFile=./bal*/chic_${condition}_merge_mm10_IJ_modelledRegion_${chr}_${start}_${end}_at_${resolution}bp.tab
    head $inFile
    wc -l $inFile
    outFile=chic_${condition}_merge_mm10_IJ_modelledRegion_${chr}_${start}_${end}_at_${resolution}bp_comparison.tab

    awk -v start=${start} -v resolution=$resolution -v size=$size 'BEGIN{for(i=0;i<size;i++){for(j=0;j<size;j++){m[i,j]}}}{if($7!=-1){v=$7}else{v="NaN"}; m[int(($2-start)/resolution),int(($5-start)/resolution)]=v}END{for(i=0;i<size;i++){for(j=i;j<size;j++){print i,j,m[i,j]}}}' ${inFile} > ${outFile}
done # Close cycle over $condition

conditions="mESC_DABDStomESCWT mESC_DADStomESCWT mESC_DBDStomESCWT mESC_dCTCFDStomESCWT mESC_DpromDStomESCWT mESC_polyADStomESCWT mESC_UPPax6DStomESCWT mESC_WTDStomESCWT NPC_DABDStoNPCWT NPC_DADStoNPCWT NPC_DBDStoNPCWT NPC_dCTCFDStoNPCWT NPC_DpromDStoNPCWT NPC_polyADStoNPCWT NPC_WTDStomESCWT NPC_WTDStoNPCWT"

for condition in ${conditions} ;
do
    inFile=./bal*/chic_${condition}_mm10_IJ_modelledRegion_${chr}_${start}_${end}_at_${resolution}bp.tab
    head $inFile
    wc -l $inFile
    outFile=chic_${condition}_merge_mm10_IJ_modelledRegion_${chr}_${start}_${end}_at_${resolution}bp_comparison.tab

    awk -v start=${start} -v resolution=$resolution -v size=$size 'BEGIN{for(i=0;i<size;i++){for(j=0;j<size;j++){m[i,j]}}}{if($7!=-1){v=$7}else{v="NaN"}; m[int(($2-start)/resolution),int(($5-start)/resolution)]=v}END{for(i=0;i<size;i++){for(j=i;j<size;j++){print i,j,m[i,j]}}}' ${inFile} > ${outFile}
done # Close cycle over $condition
