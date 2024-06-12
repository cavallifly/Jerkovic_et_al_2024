#!/bin/bash

#SBATCH --job-name initRos
#SBATCH -n 1                   # Number of cores.
#SBATCH -t 4-00:00              # Runtime in D-HH:MM
#SBATCH -o 00_generate_initial_conformation_rosette.out # File to which STDOUT will be written
#SBATCH -e 00_generate_initial_conformation_rosette.out # File to which STDERR will be written 

tag=Zfp608
resolution=5000 
chr=chr18 
s1=53500000 # Region of interest to model
s2=56900000 # Region of interest to model

#mESC
#nuCG (DNA content of a monomer in b) = 5000 bp
#bCG (Diameter of a bead in nm) = 36.2806 nm
#NCG (Number of monomers to represent the chromosome) = 3000
#rhoCG (Genome density in bp/nm³) = 0.006 nm

#NPC
#nuCG (DNA content of a monomer in b) = 5000 bp
#bCG (Diameter of a bead in nm) = 34.1407 nm
#lkCG (Kuhn length of the chain in nm) = 57.4171 nm
#NCG (Number of monomers to represent the chromosome) = 3000
#sideCG (Size of the cubic simulation box for the CG model in monomer diameters) = 31.5523 # 15.7761
#rhoCG (Genome density in bp/nm³) = 0.012 nm
side=100.

chrlength=3400000
size=680
ncopies=1
echo "Length of the chain ${chrlength} bp - ${size} ${resolution}bp-beads"

# Generate the initial rod-like conformation using TADphys
for replica in $(seq 1 1 50);
do
    if [[ -e Initial_rosette_conformation_with_pbc_replica_${replica}.dat ]];
    then
	continue
    fi
    
    seed=$(od -An -N3 -i /dev/random | awk '{print $1}')

    echo $replica $seed >> used_seeds.txt

    sed -e "s/XXXresolutionXXX/${resolution}/g" -e "s/XXXsizeXXX/${size}/g" -e "s/XXXchrlengthXXX/${chrlength}/g" -e "s/XXXseedXXX/${seed}/g" -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXncopiesXXX/${ncopies}/g" -e "s/XXXsideXXX/${side}/g" ./scripts/Initial_conformation.py | python
    
done
