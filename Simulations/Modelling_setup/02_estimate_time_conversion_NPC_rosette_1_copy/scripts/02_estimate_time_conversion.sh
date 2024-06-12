#!/bin/bash

#SBATCH --job-name estTime
#SBATCH -n 1                   # Number of cores.
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 02_estimate_time_conversion.out # File to which STDOUT will be written
#SBATCH -e 02_estimate_time_conversion.out # File to which STDERR will be written 

tag=Zfp608
resolution=5000 
chr=chr18 
s1=53500000 # Region of interest to model
s2=56900000 # Region of interest to model

#NPC
#nuCG (DNA content of a monomer in b) = 5000 bp
#bCG (Diameter of a bead in nm) = 35.1459 nm
#lkCG (Kuhn length of the chain in nm) = 55.7773 nm
#NCG (Number of monomers to represent the chromosome) = 680
#sideCG (Size of the cubic simulation box for the CG model in monomer diameters) = 19.2377 # 9.61885

nparticles=680
b=35.1459
lk=55.7773
side=19.2377 # 9.61885

chrlength=3400000
size=680
ncopies=1
echo "Length of the chain ${chrlength} bp - ${size} ${resolution}bp-beads"

pressure=0.44
prevDir=01_compression_to_NPC_phi_rosette_1_copy

for replica in $(seq 1 1 100) ;
do
    # Decompaction
    replicadir=replica_${replica}_decompaction
    if [[ ! -e ../${prevDir}/replica_${replica}_pressure_${pressure}/compressed_conformation.txt ]];
    then
	exit
    fi
    
    if [[ ! -d ${replicadir} ]];
    then
	
	echo ${replicadir}
	mkdir -p ${replicadir}
	cd ${replicadir}
	
	# 5,000,000 relaxation to collect MSD
	run=5000000
	pyfile=02_estimate_time_conversion.py
	ln -s ../../${prevDir}/replica_${replica}_pressure_${pressure}/compressed_conformation.txt .
	sed -e "s/XXXbXXX/${b}/g" -e "s/XXXlkXXX/${lk}/g" -e "s/XXXrunXXX/${run}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXnparticlesXXX/${nparticles}/g" ../scripts/${pyfile} > ${pyfile%.py}_${replica}.py

	python ${pyfile%.py}_${replica}.py &>> replica_${replica}.out
	rm -fvr mini*XYZ

	# Restart at $timestep
	timestep=100000 # 6 tau_LJ ~ 1s -> 10 min of simulations ~ 7 / 0.006 * 60 * 10 dt = 700,000 dt
	outfile=relaxed_conformation_${timestep}.txt
	
	awk 'BEGIN{flag=0}{if(flag==0) print $0; if($1=="Atoms"){flag=1; print ""}}' relaxed_conformation.txt > ${outfile}
	#cat ${outfile}
	awk -v np=${nparticles} -v t=${timestep} '{if(NR%(np+9)==2){if($1!=t){f=0}}; if(NR%(np+9)==2){if($1==t){f=1}}; if(f==1 && NF==4 && $1!="ITEM:"){print $1,1,1,$2,$3,$4};}' initial_relaxation.XYZ >> ${outfile}
	wc -l ${outfile}
	
	awk 'BEGIN{flag=1}{if($1=="Bonds"){ flag=0; print ""}; if(flag==0){print $0}}' relaxed_conformation.txt >> ${outfile}
	wc -l ${outfile}
	
	cd .. # Exit ${replicadir}
    fi
done # Close cycle over replica
