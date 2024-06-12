cellType=XXXcellTypeXXX
outFile=correlation_analysis.log
rm -fr _tmp

mainDir=${PWD}

nDir=$(ls -1 | grep Ea | grep -v corr | wc -l)
n=0

echo

for dir in $(ls -1 | grep Ea | grep -v corr | grep -v png);
do
    n=$((${n}+1))
    if [[ ! -d ${dir} ]];
    then
	echo "${dir} is not a directory!"
	continue
    fi    
    echo "$dir $n of $nDir"

    cd $dir

    for modMatrix in $(ls -1 rc*000*_*rep.png 2> /dev/null);
    do

	timestepI=$(echo $modMatrix | sed -e "s/_from_/ /g" -e "s/_to_/ /g" | awk '{print $2}')
	timestepF=$(echo $modMatrix | sed -e "s/_to_/ /g" -e "s/_every_/ /g" | awk '{print $2}')
	echo ${timestepI} ${timestepF}

	name=$(echo $modMatrix | sed -e "s/\.png//g" -e "s,\.,_,g")
	dirName=$(echo $dir | sed -e "s/\./_/g")
	echo $name
	
	SCC1=$(grep -A 1 ${name} correlations_contact_matrix.txt 2> /dev/null | grep All | awk '{printf("%.3f\n",$8)}' | uniq)    
	if [[ ${SCC1} == "" ]];
	then
	    echo "Matrix correlation NA for ${dir}"
	    echo
	    cd ..
	    exit
	fi
	echo "Matrix ${SCC1}"

	SCC2=$(grep -A 1 ${name} differences_contacts_vs_L.txt 2> /dev/null  | grep All | grep ${cellType} | awk '{printf("%.3f\n",$2)}' | uniq)    
	if [[ ${SCC2} == "" ]];
	then
	    echo "Ps difference NA for ${dir}"	
	    echo
	    cd ..
	    exit
	fi
	echo "Ps $SCC2"
	
	SCC3=$(grep -A 1 _${dirName}_ differences_IS.txt 2> /dev/null | grep -w ${timestepI} | grep -w ${timestepF} | awk '{printf("%.3f\n",$6)}' | uniq | head -1)
	if [[ ${SCC3} == "" ]];
	then
	    echo "IS correlation NA for ${dir} ${timestepI} ${timestepF}"
	    echo
	    cd ..
	    exit
	fi
	echo "IS $SCC3"
	
	SCC4=$(cat loops_rc_100nm_res_5000bp_from_${timestepI}_to_${timestepF}_*txt 2> /dev/null | grep -v adj | awk '{printf("%.2f\n",$5)}' 2> /dev/null)
	if [[ ${SCC4} == "" ]];
	then
	    echo "Loop SCC is NA for ${dir}"	    
	    echo
	    cd ..
	    exit
	fi
	echo "Loops $SCC4"
	
	echo ${dirName}_from_${timestepI}_to_${timestepF} ${SCC1} ${SCC2} ${SCC3} ${SCC4} >> ../_tmp

    done

    cd ../

done
