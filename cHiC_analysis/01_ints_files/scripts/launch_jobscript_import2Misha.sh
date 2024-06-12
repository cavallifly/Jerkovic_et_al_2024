for condition in $(ls -1 *ints | sed "s/\.noDup\.ints//g" | grep -v Rep | grep -v chr18 | grep -v DS);
do
    echo $condition

    #sbatch scripts/jobscript_import2Misha.cmd ${condition}
    bash scripts/jobscript_import2Misha.cmd ${condition} 
done
