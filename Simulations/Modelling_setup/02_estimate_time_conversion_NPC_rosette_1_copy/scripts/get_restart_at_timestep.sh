timestep=100000 # 6 tau_LJ ~ 1s -> 10 min of simulations ~ 7 / 0.006 * 60 * 10 dt = 700,000 dt


for r in $(seq 1 1 100) ;
do

    replicaDir=replica_${r}_decompaction
    outfile=${replicaDir}/relaxed_conformation_${timestep}.txt
    nparticles=$(grep atoms ${replicaDir}/relaxed_conformation.txt | awk '{print $1}')
    echo $replicaDir ${nparticles}

    awk 'BEGIN{flag=0}{if(flag==0) print $0; if($1=="Atoms"){flag=1; print ""}}' ${replicaDir}/relaxed_conformation.txt > ${outfile}
    #cat ${outfile}
    awk -v np=${nparticles} -v t=${timestep} '{if(NR%(np+9)==2){if($1!=t){f=0}}; if(NR%(np+9)==2){if($1==t){f=1}}; if(f==1 && NF==4 && $1!="ITEM:"){print $1,1,1,$2,$3,$4};}' ${replicaDir}/initial_relaxation.XYZ >> ${outfile}
    wc -l ${outfile}
    
    awk 'BEGIN{flag=1}{if($1=="Bonds"){ flag=0; print ""}; if(flag==0){print $0}}' ${replicaDir}/relaxed_conformation.txt >> ${outfile}
    #awk 'BEGIN{flag=1}{if($1=="Angles"){flag=0; print ""}; if(flag==0){print $0}}' ${replicaDir}/relaxed_conformation.txt >> ${outfile}
    #cat ./scripts/bonds.txt  >> ${outfile}
    #cat ./scripts/angles.txt >> ${outfile}

    wc -l ${outfile}
done
