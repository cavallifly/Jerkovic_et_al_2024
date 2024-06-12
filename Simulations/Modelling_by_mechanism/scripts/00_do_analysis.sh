#!/bin/bash

#SBATCH --job-name analysis
#SBATCH -t 1-00:00              # Runtime in D-HH:MM
#SBATCH -o 00_do_analysis.out # File to which STDOUT will be written
#SBATCH -e 00_do_analysis.out # File to which STDERR will be written 

cellType=$(echo ${PWD} | sed "s/_/ /g" | awk '{print $(NF-1)"_"$NF}')
echo $cellType

mainDir=${PWD}
echo $mainDir

scriptsDir=../scripts/

mkdir -p ../mishaDB/
mkdir -p ../mishaDB/trackdb/
mkdir -p ../mishaDB/trackdb/mm10/
mkdir -p ../mishaDB/trackdb/mm10/tracks/chic/
mkdir -p ../mishaDB/trackdb/mm10/tracks/chic/models_${cellType}/
mkdir -p ../mishaDB/trackdb/mm10/tracks/insulation/models_${cellType}/
rsync -avz ${scriptsDir}/chrom_sizes.txt ../mishaDB/trackdb/mm10/

outFile=correlation_analysis_${cellType}.log

for dir in $(ls -1 | grep -v "log\|out\|_v")
do
    echo $dir
    echo
    if [[ ${dir} == "" ]];
    then
	continue
    fi
    if [[ ! -d ${dir} ]];
    then
	continue
    fi

    echo "01 - Models' contact map analysis"
    echo "01a - Compute models' contact map"
    sed -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/01a_jobscript_compute_contact_map.cmd   > _tmp_${dir}.cmd
    bash _tmp_${dir}.cmd &>> 01a_compute_contact_map.out   ; rm -fr _tmp_${dir}.cmd
    echo "01b - Plot contact analysis"
    sed -e "s/XXXcellTypeXXX/${cellType}/g" -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/01b_jobscript_plot_contact_analysis.cmd > _tmp_${dir}.cmd
    bash _tmp_${dir}.cmd &>> 01b_plot_contact_analysis.out ; rm -fr _tmp_${dir}.cmd    
    echo "01c - Compute correlation between cHi-C and models contact maps"
    sed -e "s/XXXcellTypeXXX/${cellType}/g" -e "s/XXXnXXX/${n}/g" -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/01c_jobscript_compute_correlation_between_HiC_and_models_contact_maps.cmd > _tmp_${dir}.cmd
    bash _tmp_${dir}.cmd &>> 01c_compute_correlation_between_HiC_and_models_contact_maps.out  ; rm -fr _tmp_${dir}.cmd
    echo "01d - Compute difference between cHi-C and models P(s)"
    sed -e "s/XXXcellTypeXXX/${cellType}/g" -e "s/XXXnXXX/${n}/g" -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/01d_jobscript_compute_correlation_between_HiC_and_models_contacts_vs_L.cmd > _tmp_${dir}.cmd
    bash _tmp_${dir}.cmd &>> 01d_compute_correlation_between_HiC_and_models_contacts_vs_L.out ; rm -fr _tmp_${dir}.cmd
    echo
    
    echo "02  - Insulation score analysis analysis"
    echo "02a - Load contacts in misha database"
    sed -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/02a_prepare_mishaTrack_toBeLoaded.sh | bash                  &>>  02a_prepare_mishaTrack.out
    sed -e "s/XXXcellTypeXXX/${cellType}/g" -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/02a_prepare_mishaTrack.R > _tmp_${dir}.R ; Rscript _tmp_${dir}.R             &>> 02a_prepare_mishaTrack.out
    echo "02b - Compute insulation score on models' matrices"            
    sed -e "s/XXXcellTypeXXX/${cellType}/g" -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/02b_computeInsulationTracks_models.R > _tmp_${dir}.R ; Rscript _tmp_${dir}.R &>> 02b_compute_IS.out
    echo "02c - Write the models' insulation score profiles"
    mkdir -p insulation_profiles_models
    sed -e "s/XXXcellTypeXXX/${cellType}/g" -e "s/XXXdirXXX/${dir}/g" -e "s,XXXdirXXX,insulation_profiles_models,g" ${scriptsDir}/02c_writeInsulationTracks_models.R > _tmp_${dir}.R  ; Rscript _tmp_${dir}.R  &> 02c_write_IS.out
    echo "02d - Plot the models' and cHi-C insulation score profiles"
    sed -e "s/XXXcellTypeXXX/${cellType}/g" -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/02d_plotInsulationTracks_models.sh > _tmp_${dir}.sh ; bash _tmp_${dir}.sh &>> 02d_plot_IS.out
    echo "02e - Compute models' and Hi-C insulation score difference"
    sed -e "s/XXXcellTypeXXX/${cellType}/g" -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/02e_jobscript_compute_correlation_between_HiC_and_models_IS.cmd | bash &> 02e_compute_SCC_IS.out ; rm -fr _tmp_${dir}.sh _tmp_${dir}.R
    echo
    
    echo "03 - Loop analysis"
    echo "03a - Plot and compare loop strengths' distribution"
    sed -e "s/XXXdirXXX/${dir}/g" -e "s/XXXcellTypeXXX/${cellType}/g" ${scriptsDir}/03a_jobscript_loop_analysis.cmd > _tmp_${dir}.cmd
    bash _tmp_${dir}.cmd &>> 03a_loop_analysis.out ; rm -fr _tmp_${dir}.cmd    
    echo
    
    echo "04 - Compute the average effective SMC lifetime"
    if [[ ! -e ${dir}/average_SMClifetime.txt ]];
    then
	eT=$(head -22 ${dir}/replica_1/replica_1.out | grep -w s | grep "Extrusion time" | awk '{print $4}')    
	grep Relocating ${dir}/replica_*/replica_*.out | sed -e "s,\[,,g" -e "s,\],,g" > _tmp_lifetime
	Textruders=$(wc -l _tmp_lifetime | awk '{print $1}')
	echo "${dir} Average real   lifetime $(cat _tmp_lifetime | awk -v eT=$eT '{s+=$5 ; t++}END{print s/t*eT}') s (${Textruders})"           >> ${dir}/average_SMClifetime.txt
	echo "${dir} Maximum real   lifetime $(cat _tmp_lifetime | sort -k 5,5n | tail -1 | awk -v eT=$eT '{print $5*eT}') s (${Textruders})"   >> ${dir}/average_SMClifetime.txt
	echo "${dir} Average target lifetime $(cat _tmp_lifetime | awk -v eT=$eT '{s+=$NF; t++}END{print s/t*eT}') s (${Textruders})"           >> ${dir}/average_SMClifetime.txt
	echo "${dir} Maximum target lifetime $(cat _tmp_lifetime | sort -k 7,7n | tail -1 | awk -v eT=$eT '{print $NF*eT}') s (${Textruders})"  >> ${dir}/average_SMClifetime.txt
	rm _tmp_lifetime
    fi
    head ${dir}/average_SMClifetime.txt
    echo
    
    echo "05 - Create final composite figure"
    sed -e "s/XXXdirXXX/${dir}/g" ${scriptsDir}/05_jobscript_create_composite_figure.cmd > _tmp_${dir}.cmd
    bash _tmp_${dir}.cmd ${dir} &>> 05_create_composite_figure.out ; rm -fr _tmp_${dir}.cmd
    echo
    continue
    echo "06 - Zip the directories for each replica's folder to save storage space"
    cd $dir
    for replica in $(seq 1 1 100);
    do
	if [[ -d replica_${replica} ]];
	then
	    if [[ ! -e replica_${replica}.zip ]];
	    then
		zip -rm replica_${replica}.zip replica_${replica}
	    else
		rm -fr replica_${replica}
	    fi
	fi
    done
    cd ..
    
done

echo "07 - Write the table of the correlations"
rm -fr correlation_analysis_${cellType}.log
sed -e "s/XXXcellTypeXXX/${cellType}/g" ${scriptsDir}/07_write_correlations.sh | bash &> 07_write_correlations.out
echo "#$cellType" &>> ${outFile}
echo "#Simulations sorted by the mean of the 3 scores" &>> ${outFile}
echo "#system SCC_matrix Diff_Ps Diff_IS Loop_pv" &>> ${outFile}
awk '{if(NF==5) print $0,$2+$3+$4+$5}' _tmp | sort -k 6,6n | awk '{printf("%s %.2f %.2f %.2f %.2f\n",$1,$2,$3,$4,$5)}' &>> ${outFile}
echo "All analysis done on $(grep -v "#" ${outFile} | wc -l) files" &>> 07_write_correlations.out
echo "Total number of directories $(ls -1 | grep Ea | grep -v log | grep -v png | wc -l)" &>> 07_write_correlations.out
rm -fr _tmp
