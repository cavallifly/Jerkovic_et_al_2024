from tadphys.modelling.lammps_modelling import generate_chromosome_rosettes_conformation_with_pbc

nparticles = XXXsizeXXX
nchrs = XXXncopiesXXX
chromosome_particle_numbers=[nparticles]*nchrs
print(nparticles, nchrs, chromosome_particle_numbers)

side = XXXsideXXX

r = XXXreplicaXXX
seed = XXXseedXXX

#At the new coarse-graining, 14nm~1kb: What should be the size of a single bead?
for replica in range(r, r+1,1):
    initial_conformation = "Initial_rosette_conformation_with_pbc_replica_%s.dat" % (replica)
    generate_chromosome_rosettes_conformation_with_pbc ( chromosome_particle_numbers,
                                                         confining_environment=['cube',side],
                                                         rosette_radius=12.0,
                                                         particle_radius=0.5 ,
                                                         seed_of_the_random_number_generator=seed,
                                                         number_of_conformations=1,
                                                         outfile = initial_conformation,
                                                         atom_types = XXXsizeXXX+1)

