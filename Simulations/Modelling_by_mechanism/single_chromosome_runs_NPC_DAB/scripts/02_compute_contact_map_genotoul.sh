# Arguments
chr=chr18
start=53700000
end=56700000

condition=$1
cd ${condition}

modelResolution=5000
mapResolution=5000

res=$(awk -v cres=${modelResolution} -v mres=${mapResolution} 'BEGIN{print int(mres/cres)}')
nparticles=600 # All the particles in the system
ncopies=1

b=35.1459
echo ${res} ${nparticles} ${rc} ${b}

tdelta=$5

dir=$(echo ${PWD} | sed "s,/, ,g" | awk '{print $NF}')

for rc in 100 ; # From tests on loop-extrusion only 50 100 150 ; # in nm
do
    for nreplicas in $2 ; #From tests on loop-extrusion only $(seq 1 1 10);
    #for nreplicas in $(seq 1 1 10);
    do
	for tdelta in $5 ;
	do
	    #for tmin in 3400000 ; # 1 SMC lifetime
	    for tmin in ${3} ; # 1 SMC lifetime
	    do
		#for tmax in 27200000 ; #$(seq ${tmin} $((${tdelta})) 27200000);
		for tmax in ${4} ;
		do
		    for copies in intra ; #inter ;
		    do
			dirName=$(echo $dir | sed -e "s/\./_/g")
			#outMatrix=rc_${rc}nm_res_${mapResolution}bp_${dirName}_from_${tmin}_to_${tmax}_every_${tdelta}_${nreplicas}_rep.tab
			outMatrix=rc_${rc}nm_res_${mapResolution}bp_from_${tmin}_to_${tmax}_every_${tdelta}_${nreplicas}_rep.tab			
			i=1
			if [[ ${copies} == "inter" ]];
			then
			    outMatrix=rc_${rc}nm_res_${mapResolution}bp_${dirName}_from_${tmin}_to_${tmax}_every_${tdelta}_${nreplicas}_rep_inter.tab			    
			    i=0
			fi
			echo "Copies $i $copies"

			if [[ -e ${outMatrix} ]];
			then
			    continue
			fi
			
			if [[ ${tmin} -ge ${tmax} ]];
			then	    
			    continue
			fi	    
			
			rm -fvr DATA_FILES_INPUT.txt			
			echo ${tmin} ${tdelta} ${tmax}
			for t in $(seq ${tmin} ${tdelta} ${tmax});
			do
			    #echo $t
			    for r in $(seq 1 1 ${nreplicas});
			    do
				ls -1 ${PWD}/replica_${r}/*_${t}.lammpstrj ${PWD}/replica_${r}/*_${t}.XYZ 2> /dev/null >> DATA_FILES_INPUT.txt
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
			touch ${outMatrix}
			#touch ${outFile}
			
			# -r:		 Resolution of the map in monomers
			# -p:		 Number of particles in the system
			# -d:		 Contact distance cutoff between particles (nm)
			# -b:		 Monomer diameter (nm)
			# -k:		 Ploidy of the system
			# -i:                Consider (0) or not (1) inter copies contacts
			echo ${res} ${nparticles} ${rc} ${b} ${ncopies} 1
			#~/TOOLS/compute_contact_map -r ${res} -p ${nparticles} -d ${rc} -b ${b} -k ${ncopies} -i 1
			# Add the handles of 400kb
			np=$(echo ${nparticles} | awk '{print $1+80}')
			~/TOOLS/compute_contact_map -r ${res} -p ${np} -d ${rc} -b ${b} -k ${ncopies} -i ${i}
			size=$(awk -v np=${nparticles} -v r=${res} 'BEGIN{print int(np/r)}')
			echo "Size of the matrix ${size}"
			
			# Remove the 200kb of handles
			mv contacts.tab ${outMatrix}
			#awk '{print $1-40,$2-40,$3}' contacts.tab | awk -v np=${nparticles} '{if(0<=$1 && $1<np && 0<=$2 && $2<np) print $0}' > ${outMatrix}
			rm -fr output.log distances.txt DATA_FILES_INPUT.txt contacts.tab
			
		    done # Close cycle over ${copies}
		done # Close cycle over ${tmin}
	    done # Close cycle over ${tmax}
	done # Close cycle over ${tdelta}
    done # Close cycle over ${nreplicas}
done # Close cycle over ${rc}
