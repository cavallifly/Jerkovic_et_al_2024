for file in $(ls -1 obs* | grep -v masked | grep KR);
do
    echo $file
    awk '{h[$2]+=$7; h[$5]+=$7}END{for(i in h){print i,h[i]}}' ${file} | grep -v start | awk '{if($2==0){print $1}}'

    awk '{if(NF==1){masked[$1]=1}else{v=""; if(masked[$2]==1 || masked[$5]==1){v="MASKED"}; print $0,v}}' <( awk '{h[$2]+=$7; h[$5]+=$7}END{for(i in h){print i,h[i]}}' ${file} | grep -v start | awk '{if($2==0){print $1}}' ) ${file} | awk '{if($2<=$5) print $0}' > ${file%.txt}_masked.txt
done
