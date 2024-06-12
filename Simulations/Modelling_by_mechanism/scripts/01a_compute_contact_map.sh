
condition=$1
cd ${condition}

modelResolution=5000
mapResolution=5000

res=$(awk -v cres=${modelResolution} -v mres=${mapResolution} 'BEGIN{print int(mres/cres)}')
nparticles=600 # All the particles in the system
ncopies=1


b=$(grep "b=" 01*py | grep -v "#" | sed "s/b=//g" | awk '{print $1}')
echo ${res} ${nparticles} ${rc} ${b}

tdelta=$5

dir=$(echo ${PWD} | sed "s,/, ,g" | awk '{print $NF}')
mainDir=${PWD}

for rc in 100 ;
do
    for nreplicas in ${2} ;
    do
	for tmin in ${3} ;
	do
	    for tmax in ${4} ;
	    do
		for tdelta in ${5} ;
		do
		    for copies in intra ; #inter ;
		    do
			dirName=$(echo $dir | sed -e "s/\./_/g")
			outMatrix=rc_${rc}nm_res_${mapResolution}bp_from_${tmin}_to_${tmax}_every_${tdelta}_${nreplicas}_rep.tab			
			i=1
			echo "Copies $i $copies"

			if [[ -e ${outMatrix} ]];
			then
			    continue
			fi
			
			if [[ ${tmin} -gt ${tmax} ]];
			then	    
			    continue
			fi	    

			tmpDir=_tmp_${tmin}_${tdelta}_${tmax}
			if [[ -d ${tmpDir} ]];
			then
			    continue
			fi			      
			mkdir ${tmpDir}
			cd ${tmpDir}
			
			rm -fvr DATA_FILES_INPUT.txt			
			echo ${tmin} ${tdelta} ${tmax}
			for t in $(seq ${tmin} ${tdelta} ${tmax});
			do
			    #echo $t
			    for r in $(seq 1 1 ${nreplicas});
			    do
				ls -1 ${mainDir}/replica_${r}/*_${t}.lammpstrj ${mainDir}/replica_${r}/*_${t}.XYZ 2> /dev/null >> DATA_FILES_INPUT.txt
			    done # Close cycle over ${r}
			done # Close cycle over ${t}		       
			#cat DATA_FILES_INPUT.txt

			nlines=$(cat DATA_FILES_INPUT.txt | wc -l | awk '{print $1}')
			if [[ $nlines -eq 0 ]];
			then
			    rm -fr ${outMatrix} ${outMatrix%.tab}.png output.log distances.txt DATA_FILES_INPUT.txt	    
			    continue
			fi
			
			echo $outMatrix #$outProb
			touch ../${outMatrix}
			#touch ${outFile}
			
			# -r:		 Resolution of the map in monomers
			# -p:		 Number of particles in the system
			# -d:		 Contact distance cutoff between particles (nm)
			# -b:		 Monomer diameter (nm)
			# -k:		 Ploidy of the system
			# -i:                Consider (0) or not (1) inter copies contacts
			echo ${res} ${nparticles} ${rc} ${b} ${ncopies} ${i}
			# Add the handles of 400kb
			np=$(echo ${nparticles} | awk '{print $1+80}')
			../../../scripts/01a_compute_contact_map -r ${res} -p ${np} -d ${rc} -b ${b} -k ${ncopies} -i ${i}
			size=$(awk -v np=${nparticles} -v r=${res} 'BEGIN{print int(np/r)}')
			echo "Size of the matrix ${size}"
			
			# Remove the 200kb of handles
			mv contacts.tab ../${outMatrix}
			rm -fr output.log distances.txt DATA_FILES_INPUT.txt contacts.tab
			cd ../ # Exit $tmpDir
			rm -fvr ${tmpDir}
			
		    done # Close cycle over ${copies}
		done # Close cycle over ${tmin}
	    done # Close cycle over ${tmax}
	done # Close cycle over ${tdelta}
    done # Close cycle over ${nreplicas}
done # Close cycle over ${rc}
