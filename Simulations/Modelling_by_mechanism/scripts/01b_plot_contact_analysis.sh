scriptsdir=../../../scripts/
mapResolution=5000

cd $1
condition=$2

pwd

for modMatrix in $(ls -1 rc*000*_*.tab 2> /dev/null) ;
do

    name=$(echo $modMatrix | sed -e "s/\.tab//g")
    outFile=${name}.png	


    tmpDir=_tmp_${modMatrix}
    if [[ -e ${tmpDir} ]];
    then
	continue
    fi
    mkdir -p $tmpDir
    cd $tmpDir
    
    echo $modMatrix
    size=$(tail -1 ../${modMatrix} | awk '{print $1}')
    echo "Size $size"
    
    modMatrix0=${modMatrix}
    
    awk -v size=$size '{if($1!=$2 && ($1!=size && $2!=size)) print $1-40,$2-40,$3}' ../${modMatrix} | awk '{if(0<=$1 && $1<600 && 0<=$2 && $2<600) print $0}' > _tmp_mod
    modMatrix=_tmp_mod
    
    for exp in ${condition};
    do
	exp=$(echo $condition | sed -e "s/_/ /g" | awk '{print $1"_"$2}')
	
	name=$(echo $modMatrix0 | sed -e "s/\.tab//g")
	outFile=${name}.png	
	if [[ -e ../${outFile} ]];
	then
	    if [[ -e ../C_vs_L_${name}_${exp}.txt ]];
	    then
		if [[ -e ../C_vs_L_${name}_${exp}.png ]];
		then
		    ls -lrtha ../${outFile} ../C_vs_L_${name}_${exp}.txt ../C_vs_L_${name}_${exp}.png
		    continue
		fi
	    fi
	fi
	
	cHiCMatrix="../../../cHiCMatrices/chic_${exp}_merge_mm10_IJ_modelledRegion_chr18_53700000_56700000_at_${mapResolution}bp_comparison.tab"
	awk -v size=$size '{if($1!=$2 && ($1!=size && $2!=size)) print $0}' ${cHiCMatrix} > _tmp_cHiC

	expMatrix=_tmp_cHiC
	expFactor=$(awk '{d=sqrt(($1-$2)*($1-$2)); if(d==10){sum+=$3;cnt++}}END{print sum/cnt}' ${expMatrix})
	
	modFactor=$(awk '{d=sqrt(($1-$2)*($1-$2)); if(d==10){sum+=$3;cnt++}}END{print sum/cnt}' ${modMatrix})
	#head ${modMatrix}
	echo $size $modFactor $expFactor
	wc -l $expMatrix $modMatrix

	diagOff=5
	# Model
	paste ${modMatrix} ${expMatrix} | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$3/fM; if($7=="MASKED"){v=0/fM}; if(sqrt(($1-$2)*($1-$2))<diagOFF){v=0/fM}; if($1<=$2) print $1,$2,v}}' >  contact_map.tab
	# Hi-C
	paste ${modMatrix} ${expMatrix} | awk -v fM=${modFactor} -v fE=${expFactor} -v diagOFF=${diagOFF} '{if($1==$4 && $2==$5){v=$6/fE; if($7=="MASKED"){v=0/fE}; if(sqrt(($1-$2)*($1-$2))<diagOFF){v=0/fE}; if($1< $2) print $2,$1,v}}' >> contact_map.tab 

	# Compute and plot P(s)
	if [[ ! -e ../C_vs_L_${name}_${exp}.txt ]];
	then
	    paste ${modMatrix} ${expMatrix} | sort -k 3,3n | awk '{if($6!="NaN"){print $1,$2,$3}}' | grep -v NaN | awk '{d=sqrt(($1-$2)*($1-$2)); h[d]+=$3; cnt[d]++}END{f=h[10]/cnt[10]; for(i in h){print i,h[i]/cnt[i]/f}}' | sort -k 1n > C_vs_L_${name}_${exp}.txt
	fi	
	cp ../C_vs_L_${name}_${exp}.txt .
	distance=$(sort -k 2,2n C_vs_L_${name}_${exp}.txt | head -1 | awk '{print $1}')	    
	if [[ ! -e ../C_vs_L_${name}_${exp}.png ]];
	then
	    cp ../../../cHiCMatrices/average_contacts_chic_${exp}_merge_mm10_IJ_modelledRegion_chr18_53700000_56700000_at_5000bp_comparison.tab exp.txt
	    cp C_vs_L_${name}_${exp}.txt data.txt
	    echo $distance 0.01 >  distance.txt
	    echo $distance 1 >> distance.txt
	    /usr/bin/gnuplot ${scriptsdir}/01b_plot_contacts_vs_L.gp
	    mv data.ps C_vs_L_${name}_${exp}.ps
	    head data.txt exp.txt
	    rm -fr data.txt exp.txt
	fi
	mv C_vs_L_${name}_${exp}.txt ../
	
	echo ${distance} 0 >  distance.txt # Lower line
	echo 600 $(echo $distance | awk '{print 600-$1}') >> distance.txt # Lower line
	echo >> distance.txt
	echo 0 ${distance} >> distance.txt # Upper line
	echo $(echo $distance | awk '{print 600-$1}') 600 >>  distance.txt # Upper line
	
	awk -v size=${size} 'BEGIN{for(i=0;i<size-80;i++){print i,i,0}}' >> contact_map.tab
	sort -k 1,1n -k 2,2n contact_map.tab > _a.tab ; mv _a.tab contact_map.tab
	
	cbmin=$(awk '{if(NR==1){min=$3}; if($3<min) min=$3}END{print min}' contact_map.tab)
	cbmax=$(awk '{if($3>max) max=$3}END{print max}' contact_map.tab)
	echo $cbmin $cbmax
	for maxFactor in 0.1 ;
	do
	    for minFactor in 100 ;
	    do
		for scale in 0.75 ;
		do
		    for factor in 1 ;
		    do
			if [[ ! -e ${outFile} ]];
			then
			    wc -l contact_map.tab
			    sed -e "s/XXXmaxFactorXXX/${maxFactor}/g" -e "s/XXXminFactorXXX/${minFactor}/g" -e "s/XXXcbmaxXXX/${cbmax}/g" -e "s/XXXcbminXXX/${cbmin}/g" -e "s/XXXsizeXXX/$((${size}-80))/g" -e "s/XXXfactorXXX/${factor}/g" -e "s/XXXscaleXXX/${scale}/g" -e "s/XXXconditionXXX/${exp}/g" ${scriptsdir}/01b_plot_contact_matrix.gp | /usr/bin/gnuplot #2> /dev/null			    
			    mv contact_map.ps ${outFile%.png}.ps

			    bash ${scriptsdir}/ps2pdf.sh #&> /dev/null
			    mv ${outFile} ${outFile%png}pdf ../
			    mv C_vs_L_${name}_${exp}.p?? ../	   
			fi			
			ls -lrtha ../${outFile}
		    done
		done
	    done
	done
    done # Close cycle over $modMatrix
    cd ../ # Exit ${tmpDir}
    rm -fvr ${tmpDir}

done # Close cycle over $exp
rm -fr contact_map.tab _tmp_mod _tmp_cHiC distance.txt
	
cd .. # Exit ${condition}
