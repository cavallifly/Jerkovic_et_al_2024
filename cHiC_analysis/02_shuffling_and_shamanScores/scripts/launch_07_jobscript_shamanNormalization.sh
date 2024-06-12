for condition in $(ls -1 *ints | grep -v Rep | grep DS | sed -e "s/_ontarget\.noDup\.ints//g" | grep mESC_Pax);
do
    check=$(echo $condition | grep DS | wc -l | awk '{print $1}')
    if [[ $check -eq 1 ]];
    then
	condition=chic_${condition}_merge_mm10_IJ
    else
	condition=chic_${condition}_mm10_IJ
    fi
    echo $condition

    outDir=/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/02_shuffling_and_shamanScores/chic_${condition}shuffle_500Small_1000High.track
    if [[ -d ${outDir} ]];
    then
        echo "Condition ${condition} already treated."
	ls -lrtha ${outDir}
        echo "Consider removing this directory if you wish to repeat the analysis!"
	continue
    fi

    #sbatch scripts/07_jobscript_shamanNormalization.cmd ${condition}
    bash scripts/07_jobscript_shamanNormalization.cmd ${condition} 
done
