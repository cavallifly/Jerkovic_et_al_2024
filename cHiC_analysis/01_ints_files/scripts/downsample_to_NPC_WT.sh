#TADbitValidPairs2ints_NPC_WT.out:cis-contacts 6596808.000000
nContacts=6596808 # Cis-contacts of NPC_WT

start=52598845
end=57599408

sample=NPC_WT
inFile=${sample}_merge_ontarget.noDup.ints
outFile=${sample}DStoNPCWT_ontarget.noDup.ints
cp ${inFile} ${outFile}
awk -v s=${sample} '{c+=$NF}END{print "Cis-contacts "s"DStoNPCWT_ontarget",c}' ${outFile}
head ${outFile}

for sample in NPC_DA NPC_DB NPC_dCTCF NPC_Dprom NPC_polyA NPC_DAB mESC_Pax6recr ;
do
    (
	inFile=${sample}_merge_ontarget.noDup.ints
	outFile=${sample}DStoNPCWT_ontarget.noDup.ints
	if [[ ! -e ${outFile} ]];
	then
	    head -1 ${inFile} > ${outFile}

            awk -v s=${start} -v e=${end} '{if($1=="chr18" && $4==$1 && (s<=$2 && $2<=e) && (s<=$5 && $5<=e)){for(i=0;i<$7;i++){print $1"_"$2"_"$3"_"$4"_"$5"_"$6,1}}}' ${inFile} | shuf -n ${nContacts} | awk '{h[$1]+=$2}END{for(i in h) print i,h[i]}' | sed "s/_/ /g" | sort -k 2,2n -k 5,5n | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}' >> ${outFile}
	    awk -v s=${sample} '{c+=$NF}END{print "Cis-contacts "s"DStoNPCWT_ontarget",c}' ${outFile}
	fi
    ) &
	
done # Close cycle over $sample

wait
head *DStoNPCWT_ontarget.noDup.ints
