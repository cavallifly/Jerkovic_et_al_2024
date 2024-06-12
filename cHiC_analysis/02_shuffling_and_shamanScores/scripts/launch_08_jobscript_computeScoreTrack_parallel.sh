for condition in $(ls -1 *ints | grep -v Rep | grep -v DS | sed -e "s/_ontarget\.noDup\.ints//g" | grep mESC_Pax);
do
    check=$(echo $condition | grep DS | wc -l | awk '{print $1}')
    if [[ $check -eq 1 ]];
    then
	condition=${condition}_merge
    else
	condition=${condition}
    fi
    echo $condition

    #sbatch scripts/08_jobscript_computeScoreTrack_parallel.cmd ${condition}
    bash scripts/08_jobscript_computeScoreTrack_parallel.cmd ${condition}
done
