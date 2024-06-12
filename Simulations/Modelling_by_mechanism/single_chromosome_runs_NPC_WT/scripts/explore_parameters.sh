cellLine=NPC
condition=NPC_WT
scriptsDir=../scripts/

extrusionSpeeds=1000                          
loop_extruders_encounter_rule=crossing        
loop_extruder_barrier_encounter_rule=stalling

# Best set of parameters
#Ea_0.02AA_0.02BB_0.02AB_La_3.00_Sl_4800_Ne_8_eS_1000_P_0.85_0.30_0.80_0.75_0.85_0.30_0.90_Ll_30_30_0_0_30_60_30_Rl_30_300_0_0_30_30_30
attrEnergiesEpi=0.02
attrEnergiesLoop=3.00
SMClifetimes=4800
NextrudersPerMbs=8
# Permabilities of the barriers
Permeabilities1=0.85
Permeabilities2=0.30
Permeabilities3=0.80
Permeabilities4=0.75
Permeabilities5=0.85
Permeabilities6=0.30
Permeabilities7=0.90
# Maximum allowed lifetime of an extruders bumping on a barrier
Llifetimes1=30
Llifetimes2=30
Llifetimes3=0
Llifetimes4=0
Llifetimes5=30
Llifetimes6=60
Llifetimes7=30
Rlifetimes1=30
Rlifetimes2=300
Rlifetimes3=0
Rlifetimes4=0
Rlifetimes5=30
Rlifetimes6=30
Rlifetimes7=30

pythonScript=${scriptsDir}/01_runs_template_${condition}.py

for attrEnergyEpi in ${attrEnergiesEpi} ; do
    attrEnergyBB=${attrEnergyEpi}
    attrEnergyAA=${attrEnergyEpi}
    attrEnergyAB=${attrEnergyEpi}    
    for attrEnergyLoop in ${attrEnergiesLoop} ; do
	for SMClifetime in ${SMClifetimes} ; do
	    for NextrudersPerMb in ${NextrudersPerMbs} ; do
		for extrusionSpeed in ${extrusionSpeeds} ; do
		    for Permeability1 in ${Permeabilities1} ; do 
			for Llifetime1 in ${Llifetimes1} ; do
			    for Rlifetime1 in ${Rlifetimes1} ; do
				for Permeability2 in ${Permeabilities2} ; do 
				    for Llifetime2 in ${Llifetimes2} ; do
					for Rlifetime2 in ${Rlifetimes2} ; do						    
					    for Permeability3 in ${Permeabilities3} ; do 
						for Llifetime3 in ${Llifetimes3} ; do
						    for Rlifetime3 in ${Rlifetimes3} ; do
							for Permeability4 in ${Permeabilities4} ; do 
							    for Llifetime4 in ${Llifetimes4} ; do
								for Rlifetime4 in ${Rlifetimes4} ; do
								    for Permeability5 in ${Permeabilities5} ; do 
									for Llifetime5 in ${Llifetimes5} ; do
									    for Rlifetime5 in ${Rlifetimes5} ; do
										for Permeability6 in ${Permeabilities6} ; do 
										    for Llifetime6 in ${Llifetimes6} ; do
											for Rlifetime6 in ${Rlifetimes6} ; do
											    for Permeability7 in ${Permeabilities7} ; do 
												for Llifetime7 in ${Llifetimes7} ; do
												    for Rlifetime7 in ${Rlifetimes7} ; do
													
													if [[ ${NextrudersPerMb} -eq 0 ]];
													then
													    Permeability1=1.00
													    Permeability2=1.00
													    Permeability3=1.00
													    Permeability4=1.00
													    Permeability5=1.00
													    Permeability6=1.00
													    Permeability7=1.00
													    Rlifetime1=0
													    Rlifetime2=0
													    Rlifetime3=0
													    Rlifetime4=0
													    Rlifetime5=0
													    Rlifetime6=0
													    Rlifetime7=0
													    Llifetime1=0
													    Llifetime2=0
													    Llifetime3=0
													    Llifetime4=0
													    Llifetime5=0
													    Llifetime6=0
													    Llifetime7=0
													    SMClifetime=600
													    extrusionSpeeds=1000
													    lifetime=0
													fi
													
													wDir=Ea_${attrEnergyAA}AA_${attrEnergyBB}BB_${attrEnergyAB}AB_La_${attrEnergyLoop}_Sl_${SMClifetime}_Ne_${NextrudersPerMb}_eS_${extrusionSpeed}_P_${Permeability1}_${Permeability2}_${Permeability3}_${Permeability4}_${Permeability5}_${Permeability6}_${Permeability7}_Ll_${Llifetime1}_${Llifetime2}_${Llifetime3}_${Llifetime4}_${Llifetime5}_${Llifetime6}_${Llifetime7}_Rl_${Rlifetime1}_${Rlifetime2}_${Rlifetime3}_${Rlifetime4}_${Rlifetime5}_${Rlifetime6}_${Rlifetime7}
													echo $wDir
													mkdir -p ${wDir}
													cd ${wDir}
													
													if [[ ! -e 01_runs_${condition}.py ]];
													then
													    awk '{print $0}' ${pythonScript} | sed -e "s/attrEnergyBB/${attrEnergyBB}/g" -e "s/attrEnergyAA/${attrEnergyAA}/g" -e "s/attrEnergyAB/${attrEnergyAB}/g" -e "s/attrEnergyEpi/${attrEnergyEpi}/g" -e "s/attrEnergyLoop/${attrEnergyLoop}/g" -e "s/XXXloop_extruders_encounter_ruleXXX/${loop_extruders_encounter_rule}/g" -e "s/XXXloop_extruder_barrier_encounter_ruleXXX/${loop_extruder_barrier_encounter_rule}/g" -e "s/XXXSMClifetimeXXX/${SMClifetime}/g" -e "s/XXXNextrudersPerMbXXX/${NextrudersPerMb}/g" -e "s/XXXextrusionSpeedXXX/${extrusionSpeed}/g" -e "s/XXXloopingByXXX/${loopingBy}/g" -e "s/Permeability1/${Permeability1}/g" -e "s/Permeability2/${Permeability2}/g" -e "s/Permeability3/${Permeability3}/g" -e "s/Permeability4/${Permeability4}/g" -e "s/Permeability5/${Permeability5}/g" -e "s/Permeability6/${Permeability6}/g" -e "s/Permeability7/${Permeability7}/g" -e "s/Llifetime1/${Llifetime1}/g" -e "s/Llifetime2/${Llifetime2}/g" -e "s/Llifetime3/${Llifetime3}/g" -e "s/Llifetime4/${Llifetime4}/g" -e "s/Llifetime5/${Llifetime5}/g" -e "s/Llifetime6/${Llifetime6}/g" -e "s/Llifetime7/${Llifetime7}/g" -e "s/Rlifetime1/${Rlifetime1}/g" -e "s/Rlifetime2/${Rlifetime2}/g" -e "s/Rlifetime3/${Rlifetime3}/g" -e "s/Rlifetime4/${Rlifetime4}/g" -e "s/Rlifetime5/${Rlifetime5}/g" -e "s/Rlifetime6/${Rlifetime6}/g" -e "s/Rlifetime7/${Rlifetime7}/g" > 01_runs_${condition}.py
													fi
													
													for replica in $(seq 1 1 50);
													do
													    replicaDir=replica_${replica}
													    if [[ -d ${replicaDir} ]];
													    then
														continue
													    fi
													    if [[ -d ${replicaDir}.zip ]];
													    then
														continue
													    fi

													    prevDir=../../../Modelling_setup/02_estimate_time_conversion_${cellLine}_rosette_1_copy/replica_${replica}_decompaction/
													    if [[ ! -e ${prevDir}/relaxed_conformation_100000.txt ]];
													    then
														continue
													    fi
													    echo $replicaDir
													    mkdir -p ${replicaDir}
													    cd ${replicaDir}
													    
													    rsync -avz ../${prevDir}/relaxed_conformation_100000.txt initial_conformation.txt &> /dev/null
													    
													    sed -e "s/XXXreplicaXXX/${replica}/g" ../01_runs_${condition}.py > 01_runs_${condition}_replica_${replica}.py
													    sed -e "s/XXXreplicaXXX/${replica}/g" -e "s/XXXconditionXXX/${condition}/g" -e "s,XXXdirXXX,${PWD},g" ../${scriptsDir}/01_runs_template.cmd > jobscript_${condition}_replica_${replica}.cmd 
													    
													    cd .. # Exit $replicaDir
													done # Close cycle over $replica
													cd .. # Exit $wDir
												    done
												done
											    done
											done
										    done
										done
									    done
									done
								    done
								done
							    done
							done
						    done
						done
					    done
					done
				    done
				done
			    done
			done
		    done
		done
	    done
	done
    done
done
