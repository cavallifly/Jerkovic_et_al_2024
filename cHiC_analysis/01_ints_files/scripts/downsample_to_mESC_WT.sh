#TADbitValidPairs2ints_mESC_WT.out:cis-contacts 4043665.000000
nContacts=4043665 # Cis-contacts of mESC_WT

start=52598845
end=57599408


sample=mESC_WT
inFile=${sample}_ontarget.noDup.ints
outFile=${sample}DStomESCWT_ontarget.noDup.ints
cp ${inFile} ${outFile}
awk -v s=${sample} '{c+=$NF}END{print "Cis-contacts "s"DStomESCWT_ontarget",c}' ${outFile}
head ${outFile}

for sample in mESC_DA mESC_DB mESC_dCTCF mESC_Dprom mESC_DAB mESC_polyA NPC_WT ;
do
    (
	inFile=${sample}_ontarget.noDup.ints
	outFile=${sample}DStomESCWT_ontarget.noDup.ints

	if [[ ! -e ${outFile} ]];
	then
	    head -1 ${inFile} > ${outFile}

	    cat ${inFile} | awk -v s=${start} -v e=${end} '{if($1=="chr18" && $4==$1 && (s<=$2 && $2<=e) && (s<=$5 && $5<=e)){for(i=0;i<int($7);i++){print $1"_"$2"_"$3"_"$4"_"$5"_"$6,1}}}' | shuf -n ${nContacts} | awk '{h[$1]+=$2}END{for(i in h) print i,h[i]}' | sed "s/_/ /g" | sort -k 2,2n -k 5,5n | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7)}' >> ${outFile}
	    awk -v s=${sample} '{c+=$NF}END{print "Cis-contacts "s"DStomESCWT_ontarget",c}' ${outFile}
	    head ${outFile}
	fi
    ) &
	
done # Close cycle over $sample

wait
#head *DStomESCWT_ontarget.noDup.ints
