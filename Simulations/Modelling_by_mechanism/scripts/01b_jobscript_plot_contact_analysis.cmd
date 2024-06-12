#!/bin/bash

condition=XXXcellTypeXXX

for dir in $(ls -1 | grep XXXdirXXX | grep -v png | grep -v log | head -1);    
do
    if [[ ! -d ${dir} ]];
    then
	echo "${dir} is not a directory!"
	continue
    fi
    
    check=$(ls -1 ${dir}/*.tab 2> /dev/null | wc -l)
    if [[ $check -lt 1 ]];
    then
	echo "${dir} has not matrix computed. Please run step 01a first to compute it!"
	continue
    fi    
    check=$(ls -1 ${dir}/cont*100nm*excluding*.png 2> /dev/null | wc -l)
    if [[ $check -eq 121 ]];
    then
	echo "${dir} has already the plot of the matrix. Please remove if you want to re-do the plot!"
	continue
    fi    
    echo $dir
    pwd
    
    bash ../scripts/01b_plot_contact_analysis.sh ${dir} ${condition}
	
done
