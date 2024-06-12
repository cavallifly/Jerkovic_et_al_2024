#for coolFile in $(ls -1 c*.cool | grep -v 20kb)
for coolFile in $(ls -1 *.cool | grep -v 20kb | grep -v merge | grep DS); # | grep CaptureDS) # grep NoC)
do
    echo $coolFile

    bash scripts/02_jobscript_cooler_merge_and_zoomify_cools.cmd ${coolFile} &>> cooler_merge_and_zoomify_cools.out
    #sbatch scripts/11_jobscript_cooler_merge_and_zoomify_cools.cmd ${coolFile}

done


