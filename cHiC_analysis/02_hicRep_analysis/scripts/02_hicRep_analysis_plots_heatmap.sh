
for inFile in $(ls -1 01_hicRep_analysis_captureRegion_*coolFiles*.out | grep -v "+") ;
do
    echo $inFile
    #head $inFile
    rm -fr _tmp*


    awk '{print $6,$7,$(NF-1)}' ${inFile} | sed -e "s,_mm10_IJ_20kb.cool,,g" | grep Rep > _tmp1
    nSamples=$(awk '{print $1; print $2}' _tmp1 | sort | uniq | wc -l | awk '{print $1}')
    awk '{print $1; print $2}' _tmp1 | sort | uniq | awk '{print $1,NR}' > _tmp2
    min=0 #$(awk 'BEGIN{min=1.0}{if($NF!=0); if($NF<min){min=$NF}}END{printf("%.20f",min)}' _tmp1)
    max=1 #$(awk '{if($NF!=1.0); if($NF>max){max=$NF}}END{print max}' _tmp1)
    echo "Nsamples ${nSamples} Min ${min} Max ${max}"
    #head _tmp1 _tmp2
    awk -v n=$nSamples -v min=${min} -v max=${max} 'BEGIN{for(a=1;a<=n;a++){for(b=1;b<=n;b++){m[a,b]=0.}}}{if(NF==2){i[$1]=int($2)}else{v=(max-$3)/(max-min);m[i[$1],i[$2]]=v;m[i[$2],i[$1]]=v}}END{for(a=1;a<=n;a++){for(c in i){if(i[c]==a){if(a<n){printf("%s\t",c)}else{printf("%s\n",c)}}}};for(a=1;a<=n;a++){for(c in i){if(i[c]==a){printf("%s\t",c)}};for(b=1;b<=n;b++){if(b<n){printf("%.8f\t",m[a,b])}else{printf("%.8f\n",m[a,b])}}}}' _tmp2 _tmp1 > _tmp
    ymin=$(awk '{if(NR>1) for(i=2;i<=NF;i++) print $i}' _tmp | sort -k 1n | head -1)
    ymax=$(awk '{if(NR>1) for(i=2;i<=NF;i++) print $i}' _tmp | sort -k 1n | tail -1)
    echo $ymin $ymax

    sed -e "s/XXXminXXX/${ymin}/g" -e "s/XXXmaxXXX/${ymax}/g" ./scripts/02_hicRep_analysis_plots_heatmap.R > _tmp.R

    Rscript _tmp.R #&> /dev/null
    mv Rplots.png ${inFile%out}png
    rm -fr _tmp* Rplots.pdf
    #done # Close cycle over $resolution
done # Close cycle over $inFile
