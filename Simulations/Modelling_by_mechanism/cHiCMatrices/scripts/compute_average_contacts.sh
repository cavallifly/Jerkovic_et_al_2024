#inFile=$1
resolution=5000
start=53700000

#awk -v start=${start} -v res=${resolution} '{print ($2-start)/res,($5-start)/res,$7}' <(cat ${inFile} | grep -v MASKED) | awk '{d=sqrt(($1-$2)*($1-$2)); h[d]+=$3; cnt[d]++}END{f=h[10]/cnt[10]; for(i in h){print i,h[i]/cnt[i]/f}}' | sort -k 1n > average_contacts_${inFile}

for inFile in $(ls -1 chic*comparison*);
do
    echo $inFile
    cat ${inFile} | grep -v NaN | awk '{d=sqrt(($1-$2)*($1-$2)); h[d]+=$3; cnt[d]++}END{f=h[10]/cnt[10]; for(i in h){print i,h[i]/cnt[i]/f}}' | sort -k 1n > average_contacts_${inFile}
done
