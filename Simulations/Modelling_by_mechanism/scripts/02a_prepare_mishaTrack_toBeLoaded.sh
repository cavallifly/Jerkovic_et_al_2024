chrom=chr18
start=53500000
end=56900000
res=5000

mkdir -p virtual_cHiC_matrices

for file in $(ls -1 *XXXdirXXX/rc*.tab | grep -v diff); 
do
    name=$(echo ${file} | sed -e "s,\.tab,,g" -e "s,/, ,g" | awk '{print $NF}')
    dir=$(echo $file | sed "s,rc, ,g" | awk '{print $1}' | sed "s,/,,g")
    part1=$(echo $name | sed -e "s/bp/ /g" | awk '{print $1"bp"}')
    part2=$(echo $name | sed -e "s/bp_/ /g" | awk '{print $2}')
    dirName=$(echo $dir | sed -e "s/\./_/g")
    echo $file $name ${part1} $dirName ${part2}
    
    outFile=${part1}_${dirName}_${part2}.txt
    if [[ -e virtual_cHiC_matrices/${outFile} ]];
    then
	continue
    fi
    
    awk -v chrom=${chrom} -v start=${start} -v res=${res} 'BEGIN{printf("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tscore\n")}{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",chrom,$1*res+start,($1+1)*res+start,chrom,$2*res+start,($2+1)*res+start,$3)}' <(grep -v chrom1 $file) > ${outFile}

    mv -v ${outFile} virtual_cHiC_matrices
done
