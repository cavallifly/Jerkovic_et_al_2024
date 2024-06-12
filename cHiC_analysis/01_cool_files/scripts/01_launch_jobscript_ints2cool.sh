
#inFiles=$(ls -1 ../01_ints_files/*ints | grep -v Rep | grep chr18 | grep -v DS | grep mESC_Pax)
inFiles=$(ls -1 ../01_ints_files/*ints | grep -v chr18 | grep mESC_Pax)
echo $inFiles

for inFile in $inFiles ;
do
    echo $inFile

    bash scripts/01_jobscript_ints2cool.cmd ${inFile} &>> ints2cool.out
    #sbatch scripts/jobscript_ints2cool.cmd ${inFile}     

done


