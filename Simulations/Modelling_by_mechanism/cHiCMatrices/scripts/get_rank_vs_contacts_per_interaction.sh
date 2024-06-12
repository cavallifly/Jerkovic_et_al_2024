#for file in $(ls -1 *NPC_WT*balanced*masked* | grep -v sorted);
for file in $(ls -1 intScores_Zfp608_NPC_WT_reg2_5kb.txt | grep -v sorted);
do
    echo $file
    head ${file}

    sort -k 3n,3n <(awk '{if(sqrt(($2-$5)*($2-$5))>15000 && $2<$5) print $1"_"$2"_"$3,$4"_"$5"_"$6,$7}' $file) | tac | awk '{print($0,NR)}' > ${file%.txt}_sorted.txt
    head ${file%.txt}_sorted.txt

    for annFile in $(ls -1 ../Modelling_by_mechanisms/filtering_based_on_background_level/Data/annotated_*NPC_WT*);
    do
	echo $annFile
	wc -l $annFile

	mark=XXX

	awk '{if(NF==1){mark[$1]=1}else{if(mark[$1]==1 && mark[$2]==1){print $0}}}' <(awk '{print $1"_"$2"_"$3}' ${annFile}) ${file%.txt}_sorted.txt | head
	awk '{if(NF==1){mark[$1]=1}else{if(mark[$1]==1 && mark[$2]==1){print $0}}}' <(awk '{print $1"_"$2"_"$3}' ${annFile}) ${file%.txt}_sorted.txt | tail
	#grep chr18_56335000_56340000 <(awk '{print $1"_"$2"_"$3}' ${annFile})
	grep chr18_53985000_53990000 <(awk '{print $1"_"$2"_"$3}' ${annFile})

    done # Close cycle over $annFile
done # Close cycle over $file
