pi=3.14159265359
DNAcontent=5451075338 #$(awk '{s+=$2}END{print s*2}' <( grep -v M /media/data/mishaDB/trackdb/mm10/chrom_sizes.txt))

# mESC
cellType=mESC
Vnucleus=1000    ; echo "Vnucleus   (${cellType}) = ${Vnucleus}   um^3"
Vnucleolus=130   ; echo "Vnucleolus (${cellType}) = ${Vnucleolus} um^3"
Dchromatin=0.014 ; echo "Dchromatin (clutches) = ${Dchromatin} um"
Vchromatin=$(awk -v pi=${pi} -v DNA=${DNAcontent} -v Dchrom=${Dchromatin} 'BEGIN{print DNA/1000.*4./3.*pi*Dchrom*Dchrom*Dchrom/8.}') ; echo "Vchromatin (${cellType}) = ${Vchromatin} um^3"
echo $DNAcontent $Vnucleus $Vnucleolus $Vchromatin

outfile=volumetric_density_vs_n_for_${cellType}.txt
rm -fvr ${outfile}
for n in $(seq 0 1 65);
do
    Phi=$(awk -v Vchromatin=${Vchromatin} -v Vnucleus=${Vnucleus} -v Vnucleolus=${Vnucleolus} -v n=${n} 'BEGIN{print (n+1)*Vchromatin/(Vnucleus-Vnucleolus)}')
    rho=$(awk -v DNA=${DNAcontent} -v Vnucleus=${Vnucleus} -v Vnucleolus=${Vnucleolus} -v n=${n} 'BEGIN{print DNA/((Vnucleus-Vnucleolus)*1000*1000*1000)}')
    echo ${n} ${Phi} ${rho} >> ${outfile}
done # Close cycle over ${n}
#cat ${outfile}

echo ""

# NPC
cellType=NPC
Vnucleus=522     ; echo "Vnucleus   (${cellType}) = ${Vnucleus}   um^3"
Vnucleolus=33    ; echo "Vnucleolus (${cellType}) = ${Vnucleolus} um^3"
Dchromatin=0.014 ; echo "Dchromatin (clutches) = ${Dchromatin} um"
Vchromatin=$(awk -v pi=${pi} -v DNA=${DNAcontent} -v Dchrom=${Dchromatin} 'BEGIN{print DNA/1000.*4./3.*pi*Dchrom*Dchrom*Dchrom/8.}') ; echo "Vchromatin (${cellType}) = ${Vchromatin} um^3"
echo $DNAcontent $Vnucleus $Vnucleolus $Vchromatin

outfile=volumetric_density_vs_n_for_${cellType}.txt
rm -fvr ${outfile}
for n in $(seq 0 1 65);
do
    Phi=$(awk -v Vchromatin=${Vchromatin} -v Vnucleus=${Vnucleus} -v Vnucleolus=${Vnucleolus} -v n=${n} 'BEGIN{print (n+1)*Vchromatin/(Vnucleus-Vnucleolus)}')
    rho=$(awk -v DNA=${DNAcontent} -v Vnucleus=${Vnucleus} -v Vnucleolus=${Vnucleolus} -v n=${n} 'BEGIN{print DNA/((Vnucleus-Vnucleolus)*1000*1000*1000)}')
    echo ${n} ${Phi} ${rho} >> ${outfile}
done # Close cycle over ${n}
#cat ${outfile}

