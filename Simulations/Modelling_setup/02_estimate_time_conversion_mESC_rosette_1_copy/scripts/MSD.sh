mkdir -p MSD_analysis
b=$(grep "b=" scripts/02*sh | sed "s/b=//g" | awk '{print $1}')
t=$(grep "timestep=" scripts/02*py | sed -e "s/timestep=//g" -e "s/,//g" | awk '{print $1}')
echo "Monomer diameter = $b nm" 
echo "Tau_LJ = $t dt"

awk -v b=${b} -v t=${t} '{v=$2*b*b*0.001*0.001; h[$1]+=v; h2[$1]+=v*v; cnt[$1]++}END{for(i in h){avg=h[i]/cnt[i];avg2=h2[i]/cnt[i];stddev=sqrt(avg2-avg*avg); print i*t,avg,stddev,cnt[i]}}' <( cat replica_*/MSD.txt | grep -v "#" ) | sort -k 1n > ./MSD_analysis/MSD.txt
