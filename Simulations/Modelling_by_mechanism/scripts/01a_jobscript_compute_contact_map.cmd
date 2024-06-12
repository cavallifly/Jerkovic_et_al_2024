#!/bin/bash

for dir in $(ls -1 | grep XXXdirXXX | head -1 | grep -v "png\|scripts\|out\|log");
do

    if [[ ! -d ${dir} ]];
    then
	echo "${dir} is not a directory!"
	continue
    fi
    echo $dir
    rm -fvr ./${dir}/*/c*_?.*Z ./${dir}/*/c*_??.*Z ./${dir}/*/c*_???.*Z ./${dir}/*/c*_????.*Z

    for minReplicates in 50 ;
    do
	nFiles=$(ls -1 ${dir}/replica_1/*XYZ ${dir}/replica_1/*lammpstrj 2> /dev/null | wc -l | awk '{print int($1)}')
	echo $nFiles
	
	for nF in $(echo ${nFiles} | awk '{print int(($1-1)/2+1)}') ${nFiles} ;
	do
	    echo "Nfiles "$nF
	    tmax=$(ls -lrtha ${dir}/replica_1/*XYZ ${dir}/replica_1/*lammpstrj 2> /dev/null | awk '{print $NF}' | head -${nF} | tail -1 | sed -e "s/_/ /g" -e "s/\.XYZ/ /g" -e "s/\.lammpstrj/ /g" | awk '{print $NF}')
	    
	    for tmin in 0 ;
	    do
		tdelta=$(ls -lrtha ${dir}/replica_1/*XYZ ${dir}/replica_1/*lammpstrj 2> /dev/null | awk '{print $NF}' | head -2 | tail -1 | sed -e "s/_/ /g" -e "s/\.XYZ/ /g" -e "s/\.lammpstrj/ /g" | awk '{print $NF}')
		echo $tmin $tdelta $tmax
		outMatrix=rc_100nm_res_5000bp_from_${tmin}_to_${tmax}_every_${tdelta}_${minReplicates}_rep.tab			
		if [[ -e ./${dir}/${outMatrix} ]];
		then
		    echo "${outmatrix} DONE!"
		    continue
		fi
		
		tmpDir=_tmp_${tmin}_${tdelta}_${tmax}
		if [[ -d ./${dir}/${tmpDir} ]];
		then
		    continue
		fi			      

		check=$(ls -1 ${dir}/rep*/*_${tmax}.lammpstrj ${dir}/rep*/*_${tmax}.XYZ 2> /dev/null | wc -l)
		if [[ $check -lt ${minReplicates} ]];
		then
		    echo "${dir} has ${check} replicate, less than ${minReplicates} done. Not worth analysing it!"
		    continue
		fi		

		bash ../scripts/01a_compute_contact_map.sh ${dir} ${minReplicates} ${tmin} ${tmax} ${tdelta}
	    done
	done
    done
done
