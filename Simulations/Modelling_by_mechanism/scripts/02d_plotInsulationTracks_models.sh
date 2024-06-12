scriptsdir=../../scripts/
mapResolution=5000
start=53700000
end=56700000

cd XXXdirXXX

for modProfile in $(ls -1 *bed | grep _XXXdirXXX_) ;
do
    name=$(echo $modProfile | sed -e "s/\.bed//g")
    echo "Model Profile "${modProfile}    
    
    exp=XXXcellTypeXXX

    expProfile=../../insulation_profiles/insulationProfiles_${exp}_w400000bp_r5000bp_for_modelling.tab
    echo $expProfile

    if [[ ! -e ${name}.png ]];
    then
	borders=0
	awk -v s=${start} -v e=${end} -v b=${borders} '{if($2>=s+b && $3<=e-b){printf("%s\t%d\t%d\t%f\n",$1,$2,$3,-$6)}}' ${expProfile} > _tmp_XXXdirXXX_exp       
	minExp=$(sort -k 4,4n _tmp_XXXdirXXX_exp | head -1 | grep -v NA | awk '{print $4}')    
	maxExp=$(sort -k 4,4n _tmp_XXXdirXXX_exp | tail -1 | grep -v NA | awk '{print $4}')        
	echo $minExp $maxExp
	awk -v max=${maxExp} -v min=${minExp} '{if($4=="NA"){print $0}else{print $1,$2,$3,($4-min)/(max-min)}}' _tmp_XXXdirXXX_exp > _tmp_XXXdirXXX ; mv _tmp_XXXdirXXX _tmp_XXXdirXXX_exp
	minExp=$(sort -k 4,4n _tmp_XXXdirXXX_exp | head -1 | awk '{print $4}')    
	maxExp=$(sort -k 4,4n _tmp_XXXdirXXX_exp | tail -1 | awk '{print $4}')        
	echo $minExp $maxExp

	awk -v s=${start} -v e=${end} -v b=${borders} '{if($2>=s+b && $3<=e-b) print $0}' ${modProfile} > _tmp_XXXdirXXX_mod
	minMod=$(sort -k 4,4n _tmp_XXXdirXXX_mod | head -1 | awk '{print $4}')    
	maxMod=$(sort -k 4,4n _tmp_XXXdirXXX_mod | tail -1 | awk '{print $4}')        
	echo $minMod $maxMod
	awk -v max=${maxMod} -v min=${minMod} '{print $1,$2,$3,($4-min)/(max-min)}' _tmp_XXXdirXXX_mod > _tmp ; mv _tmp _tmp_XXXdirXXX_mod
	minMod=$(sort -k 4,4n _tmp_XXXdirXXX_mod | head -1 | awk '{print $4}')    
	maxMod=$(sort -k 4,4n _tmp_XXXdirXXX_mod | tail -1 | awk '{print $4}')        
	echo $minMod $maxMod
	wc -l _tmp_XXXdirXXX_mod _tmp_XXXdirXXX_exp

	awk '{print $2,$NF}' _tmp_XXXdirXXX_mod > data_XXXdirXXX.txt
	awk '{print $2,$NF}' _tmp_XXXdirXXX_exp > exp_XXXdirXXX.txt
	head exp_XXXdirXXX.txt data_XXXdirXXX.txt	
	sed -e "s/XXXtagXXX/XXXdirXXX/g" ${scriptsdir}/02d_plotInsulationTracks_models.gp | /usr/bin/gnuplot 
	mv data_XXXdirXXX.ps ${name}.ps
	bash ${scriptsdir}/ps2pdf.sh ${name}.ps #&> /dev/null	
	rm _tmp_XXXdirXXX_mod _tmp_XXXdirXXX_exp
	rm -fr data_XXXdirXXX.txt exp_XXXdirXXX.txt	
    fi
	
done # Close cycle over $modProfile

cd .. # Exit $dir

