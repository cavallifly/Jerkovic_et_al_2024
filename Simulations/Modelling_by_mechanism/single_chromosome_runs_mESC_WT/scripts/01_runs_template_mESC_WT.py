from tadphys.modelling.lammps_modelling import run_lammps
import random
from numpy import zeros

# From 52795 to 58989 both included
nparticles  = 680
resolution  = 5000

nchrs = 1
chromosome_particle_numbers=[nparticles]*nchrs
print("Nparticles",nparticles,"Ncopies",nchrs,"Particles per copy",chromosome_particle_numbers)

timestep=0.006


secondInTauLJ = 15 # tau_LJ for 1s in real time
second        = int(secondInTauLJ / 0.012 * int(0.012/timestep)) # tau_LJ for 1s in real time
print("1 Second in simulation time = ",secondInTauLJ,"TauLJ")
print("1 Second in simulation time = ",second,"dt")
print("")

# Exp extrusion speed 50kb/min
extrudedChromatinEachSecond = XXXextrusionSpeedXXX # bp/s
extrudedChromatinEachStep   = 2*resolution # On average
extrusionTimeInSeconds      = extrudedChromatinEachStep / extrudedChromatinEachSecond
extrusionTimeInTauLJ        = extrusionTimeInSeconds * secondInTauLJ
extrusionTime               = int(extrusionTimeInTauLJ / 0.012 * int(0.012/timestep))
print("Extrusion speed = ",extrudedChromatinEachSecond,"bp/s")
print("Extruded chromatin at each extrusion step = ",extrudedChromatinEachStep,"bp/step")
print("Extrusion time = ",extrusionTimeInSeconds,"s")
print("Extrusion time = ",extrusionTimeInTauLJ,"TauLJ")
print("Extrusion time = ",extrusionTime,"dt")
print("")

# Life time of cohesin on chromatin ~20min
SMClifetimeInTauLJ          = XXXSMClifetimeXXX * secondInTauLJ
SMClifetime                 = int(SMClifetimeInTauLJ / 0.012 * int(0.012/timestep))
print("SMC lifetime = ",int(XXXSMClifetimeXXX/60),"min")
print("SMC lifetime = ",XXXSMClifetimeXXX,"s")
print("SMC lifetime = ",SMClifetimeInTauLJ,"TauLJ")
print("SMC lifetime = ",SMClifetime,"dt")
print("SMC lifetime = ",int(SMClifetime/extrusionTime),"extrusionTimes")
print("Processivity (Extruded chromatin in 1 lifetime) = ",SMClifetime/extrusionTime*extrudedChromatinEachStep,"bp")
print("")

if XXXNextrudersPerMbXXX != 0 :
    lifetimes = 4 * (1200/XXXSMClifetimeXXX) * (1000/XXXextrusionSpeedXXX) * (32/XXXNextrudersPerMbXXX)
else:
    lifetimes = 2 * (1200/XXXSMClifetimeXXX) * (1000/XXXextrusionSpeedXXX)
runtime = int(lifetimes * SMClifetime)
print("Run time (",lifetimes," SMC lifetimes) =",lifetimes*XXXSMClifetimeXXX,"s")
print("Run time (",lifetimes," SMC lifetimes) =",runtime,"dt")
dumpTime = int(2*(extrusionTime))
print("Dump time =",dumpTime,"dt")

size=15
loopMatrix = zeros((size,size),dtype=float)
TADsVector = zeros((size),dtype=float)
compsVector = zeros((size),dtype=float)
# Ranked by counts
loopMatrix[4][5] = 0.1670*attrEnergyLoop # loop1
loopMatrix[6][8] = 0.3330*attrEnergyLoop # loop2
loopMatrix[7][9] = 0.5000*attrEnergyLoop # loop3
loopMatrix[10][11] = 0.6670*attrEnergyLoop # loop4
loopMatrix[12][13] = 0.8330*attrEnergyLoop # loop5
loopMatrix[12][14] = 1.0000*attrEnergyLoop # loop6

TADsVector[1] = 1
TADsVector[2] = 2
TADsVector[3] = 3
for i in list(range(4,5+1)):
    TADsVector[i] = 1
for i in list(range(5+1,11+1)):
    TADsVector[i] = 2
for i in list(range(11+1,size)):
    TADsVector[i] = 3
    
compsVector[1] = 1 # "A"
compsVector[2] = -1 # "B"
compsVector[3] = 1
for i in list(range(4,5+1))+list(range(11+1,size)):
    compsVector[i] = 1
for i in list(range(5+1,11+1)):
    compsVector[i] = -1    

attrEnergyTADs = 0
interactions = {}
for i in range(1,size):
    for j in range(1,size):
        if i <= j:
            if compsVector[i] * compsVector[j] == 1 and compsVector[i] == 1:
                compEnergy = attrEnergyAA
            elif compsVector[i] * compsVector[j] == 1 and compsVector[i] == -1:
                compEnergy = attrEnergyBB
            else:
                compEnergy = attrEnergyAB

            if TADsVector[i] == TADsVector[j]:
                TADEnergy = attrEnergyTADs
            else:
                TADEnergy = 0
                
            interactions[(i,j)] = ["attraction",compEnergy+TADEnergy+loopMatrix[i][j]]
print(interactions)

boundary1=96
boundary2=569
compartmentalizationABLoops = {
    "partition" : {1 : range(1,boundary1,1), # TAD1 A
                   2 : range(boundary1,boundary2,1), # TAD2 B
                   3 : range(boundary2,681,1), # TAD3 A                      
                   4 : [50,51], # loop1 TAD1 A
                   5 : [86,87], # loop1 TAD1 A
                   6 : [101,102], # loop2 TAD2 B
                   7 : [102,103], # loop3 TAD2 B
                   8 : [158,159], # loop2 TAD2 B
                   9 : [183,184], # loop3 TAD2 B
                   10 : [258,259], # loop4 TAD2 B
                   11 : [362,363], # loop4 TAD2 B
                   12 : [575,576], # loop5 loop6 TAD3 A
                   13 : [605,606], # loop5 TAD3 A
                   14 : [619,620], # loop6 TAD3 A
    },
    "radii" : 0.5,
    "interactions" : interactions,
}

barriers  = [50,95,188,363,569,613]

Permeability = [Permeability1,Permeability2,Permeability3,Permeability4,Permeability5,Permeability6]
Rlifetime    = [Rlifetime1,Rlifetime2,Rlifetime3,Rlifetime4,Rlifetime5,Rlifetime6]
Llifetime    = [Llifetime1,Llifetime2,Llifetime3,Llifetime4,Llifetime5,Llifetime6]

if XXXNextrudersPerMbXXX != 0 :
    separation = int(1000000/XXXNextrudersPerMbXXX/resolution)
else:
    separation = int(nparticles*resolution+1000000)

# Loop-extrusion genomic distance between extruding cohesins of ~186–372 kb
loop_extrusion_dynamics = { "separation"            : separation,
                            "lifetime"              : int(SMClifetime/extrusionTime),
                            "right_extrusion_rate"  : 1,
                            "left_extrusion_rate"   : 1,
                            "extrusion_time"        : extrusionTime,
                            "barriers"              : barriers,
                            "barriers_right_permeability"  : Permeability,
                            "barriers_left_permeability"   : Permeability,
                            "lifetimeAtBarriersRight"      : Rlifetime,
                            "lifetimeAtBarriersLeft"       : Llifetime,
                            "loop_extruders_encounter_rule" : 'XXXloop_extruders_encounter_ruleXXX', # 'crossing' 'relocating' 'stalling'
                            "loop_extruder_barrier_encounter_rule" : 'XXXloop_extruder_barrier_encounter_ruleXXX', # 'relocating' 'stalling'
                            "chrlength"             : [nparticles]*nchrs,
                            "attraction_strength"   : 10.0,
                            "equilibrium_distance"  : 0.0}

r = XXXreplicaXXX

b=36.2806
lk=54.0335

#At the new coarse-graining, 14nm~1kb: What should be the size of a single bead?
for replica in range(r, r+1,1):
    initial_conformation = "initial_conformation.txt"
    
    run_lammps(initial_conformation=initial_conformation,
               minimize  = False, 
               tethering = False,
               compartmentalization    = compartmentalizationABLoops,               
               loop_extrusion_dynamics = loop_extrusion_dynamics,
               to_dump = dumpTime,
               persistence_length = lk / b / 2,               
               run_time = runtime,
               lammps_folder = "./",
               chromosome_particle_numbers = [nparticles]*nchrs,
               timestep=timestep)
