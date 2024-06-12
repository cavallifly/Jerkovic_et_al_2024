from os import rename

from tadphys.modelling.lammps_modelling import run_lammps

nparticles = XXXnparticlesXXX

lk=XXXlkXXX
b=XXXbXXX
replica = XXXreplicaXXX


initial_conformation = "compressed_conformation.txt"
initial_relaxation = { "relaxation_time" : XXXrunXXX,
                       "MSD"             : 100,
                      }

run_lammps(initial_conformation=initial_conformation,
           minimize   = False, 
           initial_relaxation = initial_relaxation,
           to_dump = 5000,
           timestep=0.012,
           persistence_length = lk / b / 2,
           run_time = 0,
           lammps_folder = "./",
           chromosome_particle_numbers = [nparticles])
