n=0

for job in $(ls -1 ls -1 ./*/rep*/*cmd 2> /dev/null);
do 
    
    if [[ -e $job ]];
    then
	echo $job
    	sbatch $job | awk -v j=${job} '{print j,$0}' &>> ../output_jobsID.log
    	mv ${job} ${job}_launched
    	n=$(($n+1))
    	echo "Launched job number $n"
    	echo
    fi

done # Close cycle over $dir
