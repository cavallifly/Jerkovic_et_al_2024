from os import rename

from tadphys.modelling.lammps_modelling import run_lammps

# From 52795 to 58989 both included
nparticles = XXXnparticlesXXX
nchrs = 1
chromosome_particle_numbers=[nparticles]*nchrs

lk = XXXlkXXX
b  = XXXbXXX
replica = XXXreplicaXXX
side=XXXsideXXX

initial_conformation = "XXXinitial_conformationXXX"

wdir = "XXXwdirXXX/"

run_lammps(initial_conformation=initial_conformation,
           minimize   = True,
           compress_with_pbc = [0.01, XXXpressureXXX, 10000],
           to_dump = 1000,
           persistence_length = lk / b / 2,
           timestep = 0.001,
           run_time = 0,
           lammps_folder = wdir,
           confining_environment = ['cube',side],
           chromosome_particle_numbers = chromosome_particle_numbers)
