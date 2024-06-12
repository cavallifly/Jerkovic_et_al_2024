#!/bin/bash                                                                                                                                                                                                


#SBATCH --job-name tarDens
#SBATCH -n 1                   # Number of cores.
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 01_compression_to_the_target_density.out # File to which STDOUT will be written
#SBATCH -e 01_compression_to_the_target_density.out # File to which STDERR will be written 

tag=Zfp608
resolution=5000 
chr=chr18 
s1=53700000 # Region of interest to model
s2=56700000 # Region of interest to model

#mESC
#nuCG (DNA content of a monomer in b) = 5000 bp
#bCG (Diameter of a bead in nm) = 36.2806 nm
#lkCG (Kuhn length of the chain in nm) = 54.0335 nm
#NCG (Number of monomers to represent the chromosome) = 680
#sideCG (Size of the cubic simulation box for the CG model in monomer diameters) = 22.8087 # 11.4044

nparticles=680
b=36.2806
lk=54.0335
side=22.8087 # 11.4044

chrlength=3400000
size=680
ncopies=1
echo "Length of the chain ${chrlength} bp - ${size} ${resolution}bp-beads"

prevDir=00_generate_initial_conformation_rosette_1_copy

for replica in $(seq 1 1 50);
do

    if [[ ! -e ../${prevDir}/Initial_rosette_conformation_with_pbc_replica_${replica}.dat ]];
    then
	echo "The file ../${prevDir}/Initial_rosette_conformation_with_pbc_replica_${replica}.dat is not generated!"
	continue
    fi
    #for pressure in $(seq 0.1 0.1 1.0);
    #for pressure in $(seq 0.2 0.01 0.3)
    for pressure in 0.23
    do
	replicadir=replica_${replica}_pressure_${pressure}
	if [[ -d ${replicadir} ]];
	then
	    continue
	fi
	echo $replicadir
	mkdir -p ${replicadir}
	cd ${replicadir}
	
	seed=$(od -An -N3 -i /dev/urandom)
	pyfile=01_compression_to_the_target_density.py
	initial_conformation=../../${prevDir}/Initial_rosette_conformation_with_pbc_replica_${replica}.dat
	sed -e "s/XXXreplicaXXX/${replica}/g" -e "s,XXXwdirXXX,${PWD},g" -e "s,XXXinitial_conformationXXX,${initial_conformation},g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXpressureXXX/${pressure}/g" -e "s/XXXlkXXX/${lk}/g" -e "s/XXXbXXX/${b}/g" -e "s/XXXnparticlesXXX/${nparticles}/g" -e "s/XXXsideXXX/${side}/g"	../scripts/${pyfile} > ${pyfile%.py}_${replica}.py

	python ${pyfile%.py}_${replica}.py
	    
	cd .. # Exit ${replicadir}
    done # Close cycle over $pressure
done # Close cycle over $replica
