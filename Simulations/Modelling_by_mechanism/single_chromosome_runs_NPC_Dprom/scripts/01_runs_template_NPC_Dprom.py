from tadphys.modelling.lammps_modelling import run_lammps
from numpy import zeros

# From 52795 to 58989 both included
nparticles = 680
resolution = 5000

nchrs = 1
chromosome_particle_numbers=[nparticles]*nchrs
print("Nparticles",nparticles,"Ncopies",nchrs,"Particles per copy",chromosome_particle_numbers)

timestep=0.006


secondInTauLJ = 18.6 # tau_LJ for 1s in real time
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
print("SMC lifetime = ",int(SMClifetime/extrusionTime),"extrusionTimes")
print("SMC lifetime = ",SMClifetime,"dt")
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

# Compartmentalization
compartmentalizationAll = {
    "partition" : {1 : list(range(1,283,1))+list(range(299,681,1)), 2 : range(283,299,1)},
    "radii" : {1 : 0.5, 2 : 0.5},
    "interactions" : {(1,1) : ["attraction",attrEnergyEpi],(1,2) : ["attraction",attrEnergyEpi],(2,2) : ["attraction",attrEnergyEpi]}
    ###                                                                                                                                                    
}

#barriers  = [49,93,283,298,410,569,613]

size=52
loopMatrix = zeros((size,size),dtype=float)
TADsVector = zeros((size),dtype=float)
compsVector = zeros((size),dtype=float)
# Ranked by counts
LPfactor = XXXLPfactorXXX
loopMatrix[4][8] = 0.8040*attrEnergyLoop # loop1
loopMatrix[5][7] = 0.7060*attrEnergyLoop # loop2
loopMatrix[6][8] = 0.5100*attrEnergyLoop # loop3
loopMatrix[9][50] = 0.0200*attrEnergyLoop # loop4
loopMatrix[10][47] = 0.2160*attrEnergyLoop # loop5
loopMatrix[11][33] = 0.7250*attrEnergyLoop # loop6
loopMatrix[11][39] = 0.8630*attrEnergyLoop # loop7
loopMatrix[11][22] = 0.1760*attrEnergyLoop # loop8
loopMatrix[11][25] = 0.3920*LPfactor*3.00 # loop9
loopMatrix[11][30] = 0.1370*attrEnergyLoop # loop10
loopMatrix[12][14] = 0.7840*attrEnergyLoop # loop11
loopMatrix[13][24] = 0.2550*attrEnergyLoop # loop12
loopMatrix[13][25] = 0.1180*LPfactor*3.00 # loop13
loopMatrix[15][24] = 0.0590*attrEnergyLoop # loop14
loopMatrix[15][25] = 0.0980*LPfactor*3.00 # loop15
loopMatrix[16][22] = 0.2350*attrEnergyLoop # loop16
loopMatrix[16][25] = 0.4120*LPfactor*3.00 # loop17
loopMatrix[17][21] = 0.5490*attrEnergyLoop # loop18
loopMatrix[17][25] = 0.1960*LPfactor*3.00 # loop19
loopMatrix[18][26] = 0.4510*LPfactor*3.00 # loop20
loopMatrix[19][26] = 0.3140*LPfactor*3.00 # loop21
loopMatrix[20][26] = 0.7450*LPfactor*3.00 # loop22
loopMatrix[22][27] = 1.0000*LPfactor*3.00 # loop23
loopMatrix[23][34] = 0.6670*attrEnergyLoop # loop24
loopMatrix[23][31] = 0.5290*attrEnergyLoop # loop25
loopMatrix[23][40] = 0.5690*attrEnergyLoop # loop26
loopMatrix[23][43] = 0.2750*attrEnergyLoop # loop27
loopMatrix[23][46] = 0.2940*attrEnergyLoop # loop28
#loopMatrix[28][35] = 0.9410*LPfactor*attrEnergyLoop # loop29
#loopMatrix[28][36] = 0.7650*LPfactor*attrEnergyLoop # loop30
#loopMatrix[28][37] = 0.8820*LPfactor*attrEnergyLoop # loop31
#loopMatrix[28][40] = 0.9610*LPfactor*attrEnergyLoop # loop32
#loopMatrix[28][42] = 0.8430*LPfactor*attrEnergyLoop # loop33
#loopMatrix[28][46] = 0.8240*LPfactor*attrEnergyLoop # loop34
#loopMatrix[29][31] = 0.9800*LPfactor*attrEnergyLoop # loop35
#loopMatrix[29][38] = 0.4900*LPfactor*attrEnergyLoop # loop36
#loopMatrix[29][47] = 0.3730*LPfactor*attrEnergyLoop # loop37
loopMatrix[28][35] = 0.9410*LPfactor*3.00 # loop29
loopMatrix[28][36] = 0.7650*LPfactor*3.00 # loop30
loopMatrix[28][37] = 0.8820*LPfactor*3.00 # loop31
loopMatrix[28][40] = 0.9610*LPfactor*3.00 # loop32
loopMatrix[28][42] = 0.8430*LPfactor*3.00 # loop33
loopMatrix[28][46] = 0.8240*LPfactor*3.00 # loop34
loopMatrix[29][31] = 0.9800*LPfactor*3.00 # loop35
loopMatrix[29][38] = 0.4900*LPfactor*3.00 # loop36
loopMatrix[29][47] = 0.3730*LPfactor*3.00 # loop37
loopMatrix[32][40] = 0.6470*attrEnergyLoop # loop38
loopMatrix[32][43] = 0.4310*attrEnergyLoop # loop39
loopMatrix[32][45] = 0.3530*attrEnergyLoop # loop40
loopMatrix[32][47] = 0.0390*attrEnergyLoop # loop41
loopMatrix[35][47] = 0.5880*attrEnergyLoop # loop42
loopMatrix[35][40] = 0.3330*attrEnergyLoop # loop43
loopMatrix[35][42] = 0.1570*attrEnergyLoop # loop44
loopMatrix[35][45] = 0.0780*attrEnergyLoop # loop45
loopMatrix[41][45] = 0.6860*attrEnergyLoop # loop46
loopMatrix[41][47] = 0.6080*attrEnergyLoop # loop47
loopMatrix[41][42] = 0.4710*attrEnergyLoop # loop48
loopMatrix[44][47] = 0.6270*attrEnergyLoop # loop49
loopMatrix[48][51] = 0.9020*attrEnergyLoop # loop50
loopMatrix[48][49] = 0.9220*attrEnergyLoop # loop51

TADsVector[1] = 1
TADsVector[2] = 2
TADsVector[3] = 3
for i in list(range(4,8+1)):
    TADsVector[i] = 1
for i in list(range(8+1,48)):
    TADsVector[i] = 2
for i in list(range(48+1,size)):
    TADsVector[i] = 3
    
compsVector[1] = 1 # "A"
compsVector[2] = -1 # "B"
compsVector[3] = 1
for i in list(range(4,8+1))+list(range(48,size)):
    compsVector[i] = 1
for i in list(range(8+1,48)):
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

compartmentalizationABLoops = {
    "partition" : {
        1 : range(1,94,1),    # TAD1 A
        2 : range(94,570,1),  # TAD2 B
        3 : range(570,681,1), # TAD3 A
        4 : [47,48], # loop1 TAD1 A
        5 : [50,51], # loop2 TAD1 A
        6 : [54,55], # loop3 TAD1 A
        7 : [71,72], # loop2 TAD1 A
        8 : [90,91], # loop1 loop3 TAD1 A
        9 : [95,96], # loop4 TAD2 B
        10 : [97,98], # loop5 TAD2 B
        11 : [98,99], # loop10 loop6 loop7 loop8 loop9 TAD2 B
        12 : [100,101], # loop11 TAD2 B
        13 : [154,155], # loop12 loop13 TAD2 B
        14 : [159,160], # loop11 TAD2 B
        15 : [161,162], # loop14 loop15 TAD2 B
        16 : [191,192], # loop16 loop17 TAD2 B
        17 : [221,222], # loop18 loop19 TAD2 B
        18 : [241,242], # loop20 TAD2 B
        19 : [252,253], # loop21 TAD2 B
        20 : [265,266], # loop22 TAD2 B
        21 : [279,280], # loop18 TAD2 B
        22 : [280,281], # loop16 loop23 loop8 TAD2 B
        23 : [281,282], # loop24 loop25 loop26 loop27 loop28 TAD2 B
        24 : [282,283], # loop12 loop14 TAD2 B
        25 : [295,296], # loop13 loop15 loop17 loop19 loop9 TAD2 B
        26 : [296,297], # loop20 loop21 loop22 TAD2 B
        27 : [297,298], # loop23 TAD2 B
        28 : [298,299], # loop29 loop30 loop31 loop32 loop33 loop34 TAD2 B
        29 : [299,300], # loop35 loop36 loop37 TAD2 B
        30 : [325,326], # loop10 TAD2 B
        31 : [326,327], # loop25 loop35 TAD2 B
        32 : [327,328], # loop38 loop39 loop40 loop41 TAD2 B
        33 : [338,339], # loop6 TAD2 B
        34 : [339,340], # loop24 TAD2 B
        35 : [340,341], # loop29 loop42 loop43 loop44 loop45 TAD2 B
        36 : [352,353], # loop30 TAD2 B
        37 : [376,377], # loop31 TAD2 B
        38 : [396,397], # loop36 TAD2 B
        39 : [407,408], # loop7 TAD2 B
        40 : [409,410], # loop26 loop32 loop38 loop43 TAD2 B
        41 : [411,412], # loop46 loop47 loop48 TAD2 B
        42 : [451,452], # loop33 loop44 loop48 TAD2 B
        43 : [452,453], # loop27 loop39 TAD2 B
        44 : [453,454], # loop49 TAD2 B
        45 : [461,462], # loop40 loop45 loop46 TAD2 B
        46 : [463,464], # loop28 loop34 TAD2 B
        47 : [567,568], # loop37 loop41 loop42 loop47 loop49 loop5 TAD2 B
        48 : [575,576], # loop50 loop51 TAD3 A
        49 : [603,604], # loop51 TAD3 A
        50 : [605,606], # loop4 TAD3 A
        51 : [615,616], # loop50 TAD3 A
    },
    "radii" : 0.5,
    "interactions" : interactions,
}

#chr18	49	49
#chr18	93	93  # TAD1
#chr18	283	283 # 3'
#chr18	298	298 # Promoter
#chr18	410	410 # EnhA
#chr18	569	569 # TAD2
#chr18	613	613

barriers  = [49,93,283,298,410,569,613]

Permeability = [Permeability1,Permeability2,Permeability3,Permeability4,Permeability5,Permeability6,Permeability7]
Rlifetime    = [Rlifetime1    ,Rlifetime2    ,Rlifetime3    ,Rlifetime4    ,Rlifetime5    ,Rlifetime6    ,Rlifetime7    ]
Llifetime    = [Llifetime1    ,Llifetime2    ,Llifetime3    ,Llifetime4    ,Llifetime5    ,Llifetime6    ,Llifetime7    ]

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
                            "equilibrium_distance"  : 0.0
}

r = XXXreplicaXXX

b=35.1459
lk=55.7773

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
