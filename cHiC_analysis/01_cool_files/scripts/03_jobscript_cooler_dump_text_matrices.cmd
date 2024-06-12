#!/bin/bash
#SBATCH --job-name DumpText
#SBATCH -n 1                    # Number of cores. For now 56 is the number max of core available
#SBATCH --mem=16000             # allocated memory per CPU
#SBATCH --partition=computepart # specify queue partiton
#SBATCH -t 10-00:00              # Runtime in D-HH:MM
#SBATCH -o cooler_dump_text_matrices.out   # File to which STDOUT will be written
#SBATCH -e cooler_dump_text_matrices.out   # File to which STDERR will be written

#outDir=balanced_text_matrices
outDir=observed_text_matrices
mkdir -p ${outDir} 
rm -i onDiagonal_average* offDiagonal_average*


for resolution in 5000000 ; #20000 10000 5000 ;
do
   for coolFile in $(ls -1 *NoC*.mcool *CaptureDS*.mcool h*BB.mcool );
   do
        for region in chr18 ; #plotRegion captureRegion modelledRegion ;
        do 

	    if [[ $region == "chr18" ]];
            then
                coordinates="chr18:1-90,702,639"
                chrom=chr18
		check=$(echo $coolFile | grep BB | wc -l | awk '{print $1}')
		if [[ $check -eq 1 ]];
		then
                    coordinates="18:1-90,702,639"		    
		    chrom=18
		fi
		start=1
		end=90702639
            fi
	    if [[ $region == "plotRegion" ]];
            then
                coordinates="chr18:53,859,000-56,456,400"
                chrom=chr18
		start=53859000
		end=56456400
            fi	    
	    if [[ $region == "captureRegion" ]];
            then
                coordinates="chr18:52,598,845-57,599,408"
                chrom=chr18
                start=52598845
                end=57599408
            fi
            if [[ $region == "modelledRegion" ]];
            then
                coordinates="chr18:53,700,000-56,700,000"
                chrom=chr18
                start=53700000
                end=56700000
            fi
            outFile=${outDir}/${coolFile%.mcool}_${region}_${chrom}_${start}_${end}_at_${resolution}bp.tab
	    if [[ ! -e ${outFile} ]];
	    then
		#echo $outFile
	        touch ${outFile}
		cooler dump --balanced --join -r ${coordinates} ${coolFile}::/resolutions/${resolution} > ${outFile}
		awk '{if($2==$5) print $0}' ${outFile}		
		awk -v r=$resolution '{if(($2==50000000 || $5==50000000) || ($2==55000000 || $5==55000000)){next}else{if($2==$5 && ($3-$2)==r && NF==8)print $0}}' ${outFile} | awk -v r=$resolution -v f=$coolFile '{s+=$7; s2+=$7*$7; cnt++}END{avg=s/cnt; avg2=s2/cnt; stddev=sqrt(avg2-avg*avg); print "onDiagonal",f,r,avg,stddev}' >> onDiagonal_average_contacts_NoCapture.txt
		#awk -v r=$resolution '{if(($2==50000000 || $5==50000000) || ($2==55000000 || $5==55000000)){next}else{if($2!=$5 && ($3-$2)==r && ($6-$5)==r && NF==8)print $0}}' ${outFile} | awk -v r=$resolution -v f=$coolFile '{s[$2]+=$7; s2[$2]+=$7*$7; cnt[$2]++; s[$5]+=$7; s2[$5]+=$7*$7; cnt[$5]++;}END{for(i in cnt){avg=s[i]/cnt[i]; avg2=s2[i]/cnt[i]; stddev=sqrt(avg2-avg*avg); print "offDiagonal",f,r,avg,stddev}}'
		awk -v r=$resolution '{if(($2==50000000 || $5==50000000) || ($2==55000000 || $5==55000000)){next}else{if($2!=$5 && ($3-$2)==r && ($6-$5)==r && NF==8)print $0}}' ${outFile} | awk '{s[$2]+=$7; cnt[$2]++; s[$5]+=$7; cnt[$5]++;}END{for(i in cnt){avg=s[i]/cnt[i]; print i,avg}}' | awk -v r=$resolution -v f=$coolFile '{s+=$2; s2+=$2*$2; cnt++}END{avg=s/cnt; avg2=s2/cnt; stddev=sqrt(avg2-avg*avg); print "offDiagonal",f,r,avg,stddev}' >> offDiagonal_average_contacts_NoCapture.txt

		#if [[ ! -e ${outFile%.tab}.mat ]];
		#then
		    #awk -v s1=${start} -v r=${resolution} '{ind1=int(($2-s1)/r); ind2=int(($5-s1)/r); if(NF==8){v=$8}else{v=0}; if(v=="NA"){v=0} ; m[ind1,ind2]=v; m[ind2,ind1]=v; if(ind1>max){max=ind1}; if(ind2>max){max=ind2}}END{for(i=0;i<=max;i++){for(j=0;j<=max;j++){v=m[i,j];if(m[i,j]==""){v=0};printf("%s\t",v)};printf("\n")}}' ${outFile} > ${outFile%.tab}.mat
		    #awk -v s1=${start} -v r=${resolution} 'BEGIN{s=int(s1/r)*r}{for(i=1;i<=NF;i++){printf("chr18\t%d\t%d\tchr18\t%d\t%d\t%s\n",s+(i-1)*r ,s+i*r,s+(NR-1)*r,s+NR*r,$i)}}' ${outFile%.tab}.mat > ${outFile}
		#    awk -v s1=${start} -v r=${resolution} '{ind1=int(($2-s1)/r); ind2=int(($5-s1)/r); if(NF==8){v=$7}else{v="NA"} ; m[ind1,ind2]=v; m[ind2,ind1]=v; if(ind1>max){max=ind1}; if(ind2>max){max=ind2}}END{for(i=0;i<=max;i++){for(j=0;j<=max;j++){v=m[i,j];if(m[i,j]==""){v=0};printf("%s\t",v)};printf("\n")}}' ${outFile} > ${outFile%.tab}.mat
		#    awk -v s1=${start} -v r=${resolution} 'BEGIN{s=int(s1/r)*r}{for(i=1;i<=NF;i++){printf("chr18\t%d\t%d\tchr18\t%d\t%d\t%s\n",s+(i-1)*r ,s+i*r,s+(NR-1)*r,s+NR*r,$i)}}' ${outFile%.tab}.mat > ${outFile}		    
	        #fi
	        #echo "Generated ${outFile%.tab}.mat a $(awk '{print NF}' ${outFile%.tab}.mat | uniq)x$(wc -l ${outFile%.tab}.mat | awk '{print $1}') squared matrix"
            fi
	    #if [[ ! -e ${outFile%.tab}.mat ]];
	    #then
		#awk -v s1=${start} -v r=${resolution} '{ind1=int(($2-s1)/r); ind2=int(($5-s1)/r);  m[ind1,ind2]=v; m[ind2,ind1]=v; if(ind1>max){max=ind1}; if(ind2>max){max=ind2}}END{for(i=0;i<=max;i++){for(j=0;j<=max;j++){printf("%d\t",m[i,j])};printf("\n")}}' ${outFile} > ${outFile%.tab}.mat
	    #fi
        done # Close cycle over $region
    done # Close cycle over $coolFile
done # Close cycle over $resolution
