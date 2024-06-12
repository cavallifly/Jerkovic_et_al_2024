
fastaFile="/zdata/data/DBs/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index//genome.fa"
indexDir="/zdata/data/DBs/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
binaryBowtie2="/usr/bin/bowtie2"
ls -lrtha /zdata/data/DBs/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/
ls -lrtha ${fastaFile} 
ls -lrtha ${binaryBowtie2}

wdir=$PWD
enzyme=DpnII


samples="mESC_DAB_Rep1 mESC_DAB_Rep2 mESC_DA_Rep1 mESC_DA_Rep2 mESC_DB_Rep1 mESC_DB_Rep2 mESC_dCTCF_Rep1 mESC_dCTCF_Rep2 mESC_Dprom_Rep1 mESC_Dprom_Rep2 mESC_Pax6recr_Rep1 mESC_Pax6recr_Rep2 mESC_polyA_Rep1 mESC_polyA_Rep2 mESC_UPPax6_Rep1 mESC_UPPax6_Rep2 mESC_WT_Rep1 mESC_WT_Rep2 NPC_DAB_Rep1 NPC_DAB_Rep2 NPC_DA_Rep1 NPC_DA_Rep2 NPC_DB_Rep1 NPC_DB_Rep2 NPC_dCTCF_Rep1 NPC_dCTCF_Rep2 NPC_Dprom_Rep1 NPC_Dprom_Rep2 NPC_polyA_Rep1 NPC_polyA_Rep2 NPC_WT_Rep1 NPC_WT_Rep2"


for sample in ${samples}
do
    n=0
    for fastqFile1 in $(ls -1 ./raw_data/${sample}/*_1.fq 2> /dev/null) ;    
    do
	n=$((${n}+1))

	experiment=${sample}_${n}
	echo $experiment $enzyme
	
	if [[ -d ./${experiment}/03_filtered_reads ]];
	then
	    echo "Analysis done on ${experiment}. Going to the next sample!"
	    ls -lrtha ./${experiment}/03_filtered_reads/valid*	    
	    continue
	fi
	mkdir -p ${experiment}
	mkdir -p ./${experiment}/03_filtered_reads

	if [[ ! -d ./${experiment}/02_parsed_reads ]];
	then
	    # Mapping	    
	    fastqFile2=$(echo ${fastqFile1} | sed -e "s/_1\.fq\.gz/_2\.fq\.gz/g" -e "s/_1\.fq/_2\.fq/g" -e "s/_R1_001\.fastq\.gz/_R2_001\.fastq\.gz/g")
	    echo $fastqFile1 $fastqFile2

	    if [[ ! -d ./${experiment}/01_mapped_r1 ]];
	    then
		mkdir -p ./${experiment}/01_mapped_r1		
		tadbit map -w ${experiment} --mapper bowtie2 --genome ${fastaFile} --index ${indexDir} --read 1 --fastq ${fastqFile1} -C 8 --renz DpnII --mapper_binary ${binaryBowtie2} 
	    fi
	    if [[ ! -d ./${experiment}/01_mapped_r2 ]];
	    then
		mkdir -p ./${experiment}/01_mapped_r2		
		tadbit map -w ${experiment} --mapper bowtie2 --genome ${fastaFile} --index ${indexDir} --read 2 --fastq ${fastqFile2} -C 8 --renz DpnII --mapper_binary ${binaryBowtie2}
	    fi
	fi
	    
	if [[ ! -d ./${experiment}/02_parsed_reads ]];
	then
	    # Parsing
	    tadbit parse -w ./${experiment} --genome ${fastaFile} --filter_chrom "chr18*" #&>> output_${experiment}.log
	fi
	
	# Filtering
	tadbit filter -w ./${experiment} --cpus 8 #&>> output_${experiment}.log

    done # Close cycle over $experiment
done # Close cycle over $sample
