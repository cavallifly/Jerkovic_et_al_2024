pwd
dir=$1
cd $dir

dirName=$(echo $dir | sed -e "s/\./_/g")
echo $dir
echo $dirName

for modMatrix in $(ls -1 rc*.png 2> /dev/null | grep -v diff) ; 
do    
    
    name=$(echo ${modMatrix} | sed -e "s/\.png//g")
    echo $name
    tmpDir=_tmp_${name}_Fig
    if [[ -d ${tmpDir} ]];
    then
	continue
    fi
    
    part1=$(echo $name | sed -e "s/bp/ /g" | awk '{print $1"bp"}')
    part2=$(echo $name | sed -e "s/bp_/ /g" | awk '{print $2}')

    #from_0_to_59520000_every_62000_50_rep.tab
    timestepI=$(echo $part2 | sed "s/_/ /g" | awk '{print $2}')
    timestepF=$(echo $part2 | sed "s/_/ /g" | awk '{print $4}')    
    echo ${timestepI} ${timestepF}
    outFile=Fig_${part1}_${dirName}_${part2}.png
    echo $outFile   
    if [[ -e ${outFile} ]];
    then
	continue
    fi
    
    mkdir ${tmpDir}
    cd ${tmpDir}

    # Get the estimated extruder lifetime in seconds
    SMClifetime=$(grep Average ../average_SMClifetime.txt | grep real | awk '{print int($5/60)}')
    name1=$(echo $name | sed -e "s/_rep_${cellType}/_rep/g")

    # Get the time mapping to write the simuation time in the title
    td=18.6
    check=$(pwd | grep mESC | wc -l)
    if [[ ${check} -eq 1 ]];
    then
	td=15
    fi   
    t=$(echo ${timestepF} | awk -v td=${td} '{printf("%.1lf",$1*0.006/td/60)}')
    echo $t $td

    # Prepare the frame 1: contact map
    # Get the SCC of the contact matrix
    grep -A 1 ${name} ../correlations_contact_matrix.txt | grep All
    SCC=$(grep -A 1 ${name} ../correlations_contact_matrix.txt | grep All | awk '{printf("%.2f\n",$8)}' | uniq) ; echo $SCC
    # Prepare the .png of the frame 1
    convert -resize 1000x1000 -density 300 -units PixelsPerInch ../${name}.png -gravity NorthWest -fill black -pointsize 7 -annotate +120+60 "SCC = ${SCC}" frame_1.png
    convert -density 300 -units PixelsPerInch frame_1.png -gravity NorthWest -fill black -pointsize 6 -annotate +120+30 "Models at time ${t} min" frame_1a.png
    convert -density 300 -units PixelsPerInch frame_1a.png -gravity NorthWest -fill black -pointsize 6 -annotate +700+670 "CHi-C" frame_1.png

    # Prepare the frame 2 upper part: P(s)
    # Get the difference of the P(s)
    grep -A 1 ${name} ../differences_contacts_vs_L.txt | grep All    
    SCC=$(grep -A 1 ${name} ../differences_contacts_vs_L.txt | grep All | awk '{printf("%.2f\n",$2)}' | uniq) ; echo $SCC    
    # Prepare the .png of the frame 2 upper part
    convert -density 300 -units PixelsPerInch ../C_vs_L_${name}*.png                          -gravity NorthWest -pointsize 50 -annotate +900+2000 "D = ${SCC}" frame_2a.png    
    # Prepare the frame 2 lower part: Insulation profiles
    # Get the difference of the insulation score profiles
    head ../differences_IS.txt
    grep ${dirName}_from ../differences_IS.txt | grep -w "${timestepI} ${timestepF}"
    SCC=$(grep ${dirName}_from_${timestepI}_to_${timestepF} ../differences_IS.txt | awk '{printf("%.2f\n",$6)}' | uniq) ; echo $SCC
    # Prepare the .png of the frame 2 lower part
    convert -density 300 -units PixelsPerInch ../*${dirName}_from_${timestepI}_to_${timestepF}_*.png -gravity NorthWest -pointsize 50 -annotate +1500+700 "D = ${SCC} - L = ${SMClifetime} min" frame_2b.png
    # Merge the .png of the frame 2 upper part with the frame 2 lower part
    convert -resize 1000x500 -density 300 -units PixelsPerInch -append frame_2a.png frame_2b.png frame_2.png

    # Prepare the frame 3: Loops-strength distribututions
    # Get the p-value of the Loops-strength distribututions
    SCCl=$(cat ../loops_${name}.txt 2> /dev/null | grep -v adj | awk '{printf("%.2f\n",$5)}' 2> /dev/null)
    # Prepare the .png of the frame 3
    convert -resize 500x500 -density 300 -units PixelsPerInch ../loops_${name}.pdf -gravity NorthWest -pointsize 5 -annotate +250+300 "p = ${SCCl}" frame_3.png
    #convert -resize 1000x500 -density 300 -units PixelsPerInch -append frame_5.png frame_7.png

    
    convert -density 300 -units PixelsPerInch -gravity center +append frame_1.png frame_2.png frame_3.png ../${outFile}  

    cd ../ # Exit ${tmpDir}
    rm -fvr ${tmpDir}
    
done # Close cycle over modMatrix
rm -fr frame*.png
cd .. # Exit $dir
