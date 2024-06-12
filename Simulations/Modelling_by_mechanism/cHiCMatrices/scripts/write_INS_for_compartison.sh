# modeledRegion
chr=chr18
start=53905000
end=56050000
resolution=5000

size=$(awk -v start=${start} -v end=$end -v resolution=$resolution 'BEGIN{print (end-start)/resolution}')

for condition in NPC_WT mESC_WT ;
do
    inFile=Insulation_chic_${condition}_merge_mm10_IJ_INS_Scale_200kb_Res_5000bp.bedGraph
    outFile=Insulation_chic_${condition}_merge_mm10_IJ_INS_Scale_200kb_Res_5000bp_comparison.bedGraph
    awk -v c=$chr -v start=${start} -v end=$end '{if((start<=$2 && $2<=end) && (start<=$3 && $3<=end)) print $2,$NF}' ${inFile} > ${outFile}

done # Close cycle over $condition
