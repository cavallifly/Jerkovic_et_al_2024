python scripts/pull_TADbit_stats.py

inFile=TADbit_mapping_stats_per_sample.txt
grep Input ${inFile} > _tmp
sort -k 1,1d ${inFile} | grep -v Input >> _tmp ; mv _tmp ${inFile}
echo $inFile
head $inFile
echo 

if [[ ! -e ../01_ints_files/TADbitValidPairs2ints_mESC_Pax6recr.out ]];
then
    echo "Please compute the number of cis-contacts in the capture regions in ../01_ints_files before proceeding with this analysis!"
    exit
fi

cisFile=TADbit_mapping_stats_cisContacts.txt
grep cis ../01_ints_files/*out | grep Rep | sed -e "s,../01_ints_files/TADbitValidPairs2ints_,,g" -e "s/\.out/ /g" | uniq | grep -v chr18 | awk '{if($NF!=0) print $1,int($NF)}' | sort -k 1,1d | uniq > ${cisFile}
cat TADbit_mapping_stats_cisContacts.txt

outFile=TADbit_mapping_stats_per_replicate.txt
sed "s/_/ /g" ${inFile} | grep -v Input | awk '{printf("%s_%s_%s\t",$1,$2,$3); for(i=5;i<=NF;i++){printf("%-15d\t",$i)}; printf("\n")}' | awk '{names[$1]=$1;nf=0;for(i=2;i<=NF;i++){m[$1,i-1]+=$i;nf++}}END{for(name in names){printf("%-20s\t",name);for(i=1;i<=nf;i++){printf("%-20s\t",m[name,i])};printf("\n")}}' | sort -k 1,1d > _tmp
echo "Sample	Input_read_pairs	Mapped_full_r1	Mapped_frag_r1	Mapped_full_r2 Mapped_frag_r2	Uniquely_mapped_pairs	Valid-pairs Cis-contacts" | awk '{for(i=1;i<=NF;i++){printf("%-18s\t",$i)};printf("\n")}' > ${outFile}
paste _tmp ${cisFile} | awk '{if($1==$(NF-1)){for(i=1;i<NF;i++){if(i==(NF-1)){continue};printf("%-18s\t",$i)};printf("%-18s\n",$NF)}}' | awk '{if(NF!=0) print $0}' >> ${outFile}
echo $outFile
head $outFile
echo

inFile=TADbit_mapping_stats_per_replicate.txt
outFile=TADbit_mapping_stats_per_condition.txt
head -1 ${inFile} > ${outFile}
sed "s/_/ /g" ${inFile} | grep -v Input | awk '{printf("%s_%s\t",$1,$2); for(i=4;i<=NF;i++){printf("%-15d\t",$i)}; printf("\n")}' | awk '{names[$1]=$1;nf=0;for(i=2;i<=NF;i++){m[$1,i-1]+=$i;nf++}}END{for(name in names){printf("%-20s\t",name);for(i=1;i<=nf;i++){printf("%-20s\t",m[name,i])};printf("\n")}}' | sort -k 1,1d | awk '{if(NF!=0) print $0}' >> ${outFile}
echo $outFile
head $outFile
echo
