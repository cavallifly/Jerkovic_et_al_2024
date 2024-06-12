"""
01 Oct 2020
"""

from string import ascii_uppercase as uc, ascii_lowercase as lc
from os.path import exists
from random import uniform, randint, seed, random, sample, shuffle, choice
from pickle import load, dump
from multiprocessing.dummy import Pool as ThreadPool
from functools import partial
from math import atan2
from itertools import combinations, product, chain
from shutil import copyfile
from operator import itemgetter

import sys
import copy
import os
import shutil
import multiprocessing
import time

from numpy import sin, cos, arccos, sqrt, fabs, pi, zeros, log, exp, array_equal, full_like, ones
from scipy import spatial
import numpy as np
from mpi4py import MPI
from lammps import lammps

from tadphys.modelling import LAMMPS_CONFIG as CONFIG

class InitalConformationError(Exception):
    """
    Exception to handle failed initial conformation.
    """
    pass

### START run_lammps ###
### interactions: chain connectivity (FENE) ; excluded volume (WLC) ; and bending rigidity (KP) ###
def run_lammps(lammps_folder,
               run_time,
               initial_conformation        = None,
               chromosome_particle_numbers = None,
               confining_environment       = None,
               neighbor           = CONFIG.neighbor,
               connectivity       = {1:["FENE", 30.0, 1.5, 1.0, 1.0]},
               excludedVolume     = {(1,1):["LJ", 1.0, 1.0, 1.12246152962189]},
               persistence_length = CONFIG.persistence_length,
               fixed_particles    = None,
               tethering          = False,                              
               minimize           = True,
               initial_relaxation = None,               
               compress_with_pbc    = None,
               compartmentalization = None,
               loop_extrusion_dynamics = None,
               gamma    = CONFIG.gamma,
               timestep = CONFIG.timestep,
               reset_timestep = 0,
               to_dump  = 10000,
               write_restart = False):
    """
    Performs the LAMMPS dynamics.
    :param lammps_folder: Path to the folder to store the output of the simulation.
    :param run_time: Number of timesteps to simulate.
    :param None initial_conformation: path to initial conformation file or None for random walk initial conformation.
    :param None chromosome_particle_numbers: List of the number of particles per simulated chain.
    :param None confining_environment: dictionary with the confining environment of the conformation
      Possible confining environments:
      ['cube',edge_width]
      ['sphere',radius]
    :param CONFIG.neighbor neighbor: see LAMMPS_CONFIG.py.
    :param FENE connectivity: use FENE for a fene bond or harmonic for harmonic
      potential for neighbours (see init_lammps_run)
    :param LJ   excludedVolume: use {(1,1) : ["LJ", 1.0, 1.0, 1.12246152962189]} for purely-repulsive Lennard-Jones interation    
    :param CONFIG.persistence_length persistence_length : see LAMMPS_CONFIG.py.
    :param None fixed_particles: List of particles which don't move during the simulation.
    :param False tethering: whether to apply tethering command or not.
    :param True minimize: whether to apply minimize command or not.    
    :param None  initial_relaxation: whether to apply an initial relation or not.
    :param None compress_with_pbc: List of parameters to apply the compression dynamics in case of a
      system with cubic confinement and pbc. This compression step is used to obtain a system
      with the desired particle density. The input have to be a list 
      of three elements:
      0 - Initial value of the external pressure
      1 - Final   value of the external pressure
      2 - The compression simulation time span (in timesteps).
      e.g. compress_with_pbc=[0.01, 0.01, 100000]
    :param None  compartmentalization: Dictionary with all the parameters to perform the compartmentalization runs:
      XXX
    :param None loop_extrusion_dynamics: Dictionary with all the parameters to perform the loop-extrusion dynamics:
      XXX
    :param CONFIG.gamma    gamma         : see LAMMPS_CONFIG.py.
    :param CONFIG.timestep timestep      : see LAMMPS_CONFIG.py.
    :param 0               reset_timestep: Initial number of integration timesteps.
    :param 10000           to_dump       : Number of integration timesteps between two dumped conformations.
    :param False write_restart: path to file to save a text file with the final LAMMPs configuration.
    :returns: None
    """

    lmp = lammps(cmdargs=['-screen','none','-log',lammps_folder+'log.lammps','-nocite'])
    me = MPI.COMM_WORLD.Get_rank()
    nprocs = MPI.COMM_WORLD.Get_size()

    init_lammps_run(lmp,
                    initial_conformation,
                    neighbor=neighbor,
                    timestep=timestep,
                    reset_timestep=reset_timestep,
                    chromosome_particle_numbers=chromosome_particle_numbers,
                    connectivity=connectivity,
                    excludedVolume=excludedVolume,
                    persistence_length=persistence_length)
        
    lmp.command("dump    1       all    custom    %i   %slangevin_dynamics_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
    lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")

    #######################################################
    # Set up fixes                                        #
    # use NVE ensemble                                    #
    # Langevin integrator Tstart Tstop 1/friction rndseed #
    # => sampling NVT ensamble                            #
    #######################################################
    # Set the group of particles that will be moved during the molecular dynamics
    nparticles = lmp.get_natoms()
    print("Number of atoms in the system =",nparticles)

    # Define fixed and mobile particles
    fixed_group  = ""
    mobile_group = ""
    if not fixed_particles:
        fixed_particles = []
    for particle in range(1,nparticles+1):
        if particle in fixed_particles:
            fixed_group += "%d " % (particle)
            #print("set atom %d vx 0.0 vy 0.0 vz 0.0" % particle)
            lmp.command("set atom %d vx 0.0 vy 0.0 vz 0.0" % particle)
        else:
            mobile_group += "%d " % (particle)
            
    #print("group mobile_particles id %s" % mobile_group)
    lmp.command("group mobile_particles id %s" % mobile_group)

    if fixed_particles != []:
        #print("group fixed_particles id %s" % fixed_group)
        lmp.command("group fixed_particles id %s" % fixed_group)
        #print("fix freeze fixed_particles setforce 0.0 0.0 0.0")
        lmp.command("fix freeze fixed_particles setforce 0.0 0.0 0.0")

    # Define the langevin dynamics integrator
    lmp.command("fix 1 mobile_particles nve")
    lmp.command("fix 2 mobile_particles langevin 1.0  1.0  %f %i" % (gamma,randint(1,100000)))

    # Define confininn environment
    if confining_environment:
        if confining_environment[0] == "sphere":

            radius = confining_environment[1]
            lmp.command("region sphere sphere 0.0 0.0 0.0 %f units box side in" % radius)
            print("region sphere sphere 0.0 0.0 0.0 %f units box side in" % radius)
            
            # Performing the simulation
            lmp.command("fix 5 all  wall/region sphere lj126 1.0 1.0 1.12246152962189")
            print("fix 5 all  wall/region sphere lj126 1.0 1.0 1.12246152962189")
    
    # Define the tethering to the center of the confining environment
    if tethering:
        lmp.command("fix 3 all spring tether 50.0 0.0 0.0 0.0 0.0")

    # Do a minimization step to prevent particles overlaps
    if minimize:

        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %sminimization_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")
        
        print("Performing minimization run...")
        lmp.command("minimize 1.0e-4 1.0e-6 100000 100000")
        
        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %slangevin_dynamics_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")        
        lmp.command("reset_timestep 0") 

    # Do a compression to reach the desired density
    if compress_with_pbc:

        lmp.command("velocity all create 1.0 %s" % randint(1,100000))
        
        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %scompress_with_pbc_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append no")

        # Re-setting the initial timestep to 0
        lmp.command("reset_timestep 0")

        lmp.command("unfix 1")
        lmp.command("unfix 2")
        lmp.command("thermo %d" % (compress_with_pbc[2]/100))
        lmp.command("thermo_style   custom   step temp etotal pxx pyy pzz pxy pxz pyz xlo xhi ylo yhi zlo zhi")
        
        # default as in PLoS Comp Biol Di Stefano et al. 2013 compress_with_pbc = [0.01, 0.01, 100000]
        lmp.command("fix 1 all   nph   iso   %s %s   2.0" % (compress_with_pbc[0], 
                                                             compress_with_pbc[1]))
        lmp.command("fix 2 mobile_particles langevin 1.0  1.0  %f %i" % (gamma,randint(1,100000)))
        print("run %i" % compress_with_pbc[2])
        lmp.command("run %i" % compress_with_pbc[2])

        lmp.command("unfix 1")
        lmp.command("unfix 2")

        lmp.command("fix 1 mobile_particles nve")
        lmp.command("fix 2 mobile_particles langevin 1.0  1.0  %f %i" % (gamma,randint(1,100000)))
                    
        # Here We have to re-define the confining environment
        print("# Previous particle density (nparticles/volume)", lmp.get_natoms()/(confining_environment[1]**3))
        confining_environment[1] = lmp.extract_box()[1][0] - lmp.extract_box()[0][0]
        print("")
        print("# New cubic box dimensions after isotropic compression")
        print(lmp.extract_box()[0][0],lmp.extract_box()[1][0])
        print(lmp.extract_box()[0][1],lmp.extract_box()[1][1])
        print(lmp.extract_box()[0][2],lmp.extract_box()[1][2])
        print("# New confining environment", confining_environment)
        print("# New volumetric density (Vpolymer/Vbox)", lmp.get_natoms()*0.5*0.5*0.5*4./3.*pi/(confining_environment[1]**3))
        print("")

        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %slangevin_dynamics_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            #lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")        
        lmp.command("write_data compressed_conformation.txt nocoeff")
        lmp.command("reset_timestep 0")         
            
    if initial_relaxation:        
        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %sinitial_relaxation.XYZ  id xu yu zu" % (to_dump,lammps_folder))
            lmp.command("dump_modify     1 format line \"%d %.5f %.5f %.5f\" sort id append yes")
        if "MSD" in initial_relaxation:
            lmp.command("compute MSD all msd")
            lmp.command("variable MSD equal c_MSD[4]")
            lmp.command("variable step  equal step")
            lmp.command("fix MSD all print %i \"${step} ${MSD}\" file MSD.txt" % (initial_relaxation["MSD"]))

        if 'distances' in initial_relaxation:
            lmp.command("compute pairs     all property/local patom1 patom2")
            lmp.command("compute distances all pair/local dist")
            lmp.command("dump distances all local %i distances.txt c_pairs[1] c_pairs[2] c_distances" % (initial_relaxation['distances']))
            
        lmp.command("reset_timestep 0")
        lmp.command("run %i" % initial_relaxation["relaxation_time"])
        lmp.command("write_data relaxed_conformation.txt nocoeff")

        if "MSD" in initial_relaxation:
            lmp.command("uncompute MSD")

        if "distances" in initial_relaxation:
            lmp.command("uncompute distances")
            lmp.command("uncompute pairs")
            lmp.command("undump    distances")
            
    # Set interactions for chromosome compartmentalization (block co-polymer model)
    if compartmentalization:

        if 'gyration' in compartmentalization:
            lmp.command("compute RgSquared all gyration")
            lmp.command("variable RgSquared equal c_RgSquared")
            lmp.command("variable step equal step")
            lmp.command("fix RgSquared all print %i \"${step} ${RgSquared}\" file %sRg.txt" % (compartmentalization['gyration'],lammps_folder))

        if 'Ree' in compartmentalization:
            lmp.command("compute xu all property/atom xu")
            lmp.command("compute yu all property/atom yu")
            lmp.command("compute zu all property/atom zu")
            lmp.command("variable Reex equal c_xu[1]-c_xu[%i]" % (chromosome_particle_numbers[0]))
            lmp.command("variable Reey equal c_yu[1]-c_yu[%i]" % (chromosome_particle_numbers[0]))
            lmp.command("variable Reez equal c_zu[1]-c_zu[%i]" % (chromosome_particle_numbers[0]))
            lmp.command("variable step equal step")
            lmp.command("fix Ree all print %i \"${step} ${Reex} ${Reey} ${Reez}\" file %sRee.txt" % (compartmentalization['Ree'],lammps_folder))
            
        if 'distances' in compartmentalization:
            lmp.command("compute pairs     all property/local patom1 patom2")
            lmp.command("compute distances all pair/local dist")
            lmp.command("dump distances all local %i %sdistances_*.txt c_pairs[1] c_pairs[2] c_distances" % (compartmentalization['distances'],lammps_folder))
            lmp.command("dump_modify distances delay 1")
            
        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %scompartmentalization_*.XYZ  id  xu yu zu" % (to_dump,lammps_folder))
            lmp.command("dump_modify 1 format line \"%d %.5f %.5f %.5f\" sort id append no")
            
        # First we have to partition the genome in the defined compartments
        for group in compartmentalization['partition']:
            if isinstance(compartmentalization['partition'][group], (int)):
                compartmentalization['partition'][group] = [compartmentalization['partition'][group]]
            list_of_particles = get_list(compartmentalization['partition'][group])
            
            for atom in list_of_particles:
                #print("set atom %s type %s" % (atom,group+1))
                lmp.command("set atom %s type %s" % (atom,group+1))

        lmp.command("pair_coeff * * lj/cut 1.0 1.0 1.12246152962189")
        
        # Second we have to define the type of interactions
        for pair in compartmentalization['interactions']:

            t1 = pair[0]+1
            t2 = pair[1]+1 
            if t1 > t2:
                t1 = pair[1]+1
                t2 = pair[0]+1 

            epsilon = compartmentalization['interactions'][pair][1]

            try:
                sigma1 = compartmentalization['radii'][pair[0]]
            except:
                sigma1 = 0.5
            try: 
                sigma2 = compartmentalization['radii'][pair[1]]
            except:
                sigma2 = 0.5
            sigma = sigma1 + sigma2
                
            if compartmentalization['interactions'][pair][0] == "attraction":
                rc = sigma * 2.5
            if compartmentalization['interactions'][pair][0] == "repulsion":
                sigma = 3.0
                rc    = sigma * 1.12246152962189

            if epsilon != 0.0:
                lmp.command("pair_coeff %s %s lj/cut %s %s %s" % (t1,t2,epsilon,sigma,rc))

        try:
            lmp.command("minimize 1.0e-4 1.0e-6 100000 100000")
            lmp.command("reset_timestep 0")
            lmp.command("run %s" % (compartmentalization['runtime']))
            lmp.command("write_data %srestart_conformation.txt nocoeff" % lammps_folder)
        except:
            pass
        
        
    # Perform loop-extrusion dynamics
    if loop_extrusion_dynamics:        
        
        np.random.seed(424242)
        
        print("left_extrusion_rate",loop_extrusion_dynamics["left_extrusion_rate"])
        print("right_extrusion_rate",loop_extrusion_dynamics["right_extrusion_rate"])
        
        if 'lifetimeConstant' in loop_extrusion_dynamics:
            lifetimeExponential = 0
            print("Target extruders lifetimes will be always equal to",loop_extrusion_dynamics['lifetime'])
        else:
            lifetimeExponential = 1
            print("Target extruders lifetimes will be drawn from an exponential distribution with average equal to",loop_extrusion_dynamics['lifetime'])

        natoms = lmp.get_natoms()
        print(natoms)
        lmp.command("fix f_unwr all store/state 1 xu yu zu")
        xc_tmp = np.array(lmp.extract_fix("f_unwr",1,2).contents[0:(natoms-2)])
        distances = compute_particles_distance(xc_tmp)

        # Start relaxation step
        try:
            lmp.command("run %i" % loop_extrusion_dynamics['timesteps_relaxation'])
        except:
            pass

        # Start Loop extrusion dynamics
        if to_dump:
            lmp.command("undump 1")
            lmp.command("dump    1       all    custom    %i   %sloop_extrusion_MD_*.lammpstrj  id  xu yu zu" % (to_dump,lammps_folder))
            lmp.command("dump_modify 1 format line \"%d %.5f %.5f %.5f\" sort id append no")
            
        # Get the positions of the fixed extruders and 
        extruders_positions = []
        extruders_lifetimes = []
        extruders_target_lifetimes = []
        right_bumps = []
        left_bumps = []
        
        try:
            print("Defined barriers",loop_extrusion_dynamics['barriers'])
        except:
            pass
        
        try:
            print("Defined topology",loop_extrusion_dynamics['topology'])
        except:
            loop_extrusion_dynamics['topology'] = "Linear"
        
        ### Define active barriers ###
        print("### BEGIN Defined barriers ###")
        print("#Barrier left_permeability right_permeability lifetimeAtBarrierLeft lifetimeAtBarrierRight")
        if 'barriers' in loop_extrusion_dynamics:
            for barrier in loop_extrusion_dynamics['barriers']:
                left_permeability    = "NA"
                right_permeability   = "NA"
                barrierLifetimeLeft  = "NA"
                barrierLifetimeRight = "NA"
                expression_barrier   = "NA"
                if 'barriers_left_permeability' in loop_extrusion_dynamics:
                    left_permeability    = loop_extrusion_dynamics['barriers_left_permeability'][loop_extrusion_dynamics['barriers'].index(barrier)]
                if 'barriers_right_permeability' in loop_extrusion_dynamics:
                    right_permeability   = loop_extrusion_dynamics['barriers_right_permeability'][loop_extrusion_dynamics['barriers'].index(barrier)]
                if 'lifetimeAtBarriersRight' in loop_extrusion_dynamics:
                    barrierLifetimeRight = loop_extrusion_dynamics['lifetimeAtBarriersRight'][loop_extrusion_dynamics['barriers'].index(barrier)]                    
                if 'lifetimeAtBarriersLeft' in loop_extrusion_dynamics:
                    barrierLifetimeLeft  = loop_extrusion_dynamics['lifetimeAtBarriersLeft'][loop_extrusion_dynamics['barriers'].index(barrier)]
                print(barrier, left_permeability, right_permeability, barrierLifetimeLeft, barrierLifetimeRight,"extrusionTimes")
        else:
            print("You didn't define any barriers for loop-extruders")
            loop_extrusion_dynamics['barriers_left_permeability']  = []
            loop_extrusion_dynamics['barriers_right_permeability'] = []
            loop_extrusion_dynamics['barriers']                    = []
        print("### END Defined barriers ###")

        barriersOccupation = []
        extruderOnBarrierLeft   = []
        extruderOnBarrierRight  = []
        if 'barriers' in loop_extrusion_dynamics:
            for barrier in loop_extrusion_dynamics['barriers']:
                extruderOnBarrierLeft.append(-1)
                extruderOnBarrierRight.append(-1)
            

        occupied_positions = barriersOccupation
        print("Barriers' occupied positions",sorted(occupied_positions))
        lifetimeFactor = 1 #int(log(loop_extrusion_dynamics['lifetime'])/log(10))
        print("Lifetime ",loop_extrusion_dynamics['lifetime'])
        print("Lifetime factor ",lifetimeFactor)
                            
        sys.stdout.flush()        

        print("### BEGIN Positioning extruders ###")
        offset = 0 
        print("Positioning extruders on random loading points:")
        for c in range(len(loop_extrusion_dynamics['chrlength'])): 
            nextruders = int(loop_extrusion_dynamics['chrlength'][c]/loop_extrusion_dynamics['separation'])
            print("Number of extruders ",nextruders,"for copy",c+1)       
            for extruder in range(nextruders):
                print("Positioning extruder",extruder+1)

                new_positions = draw_loop_extruder_loading_site(loop_extrusion_dynamics['chrlength'][c],distances)
                new_positions[0] = offset + new_positions[0]
                new_positions[1] = offset + new_positions[1]
                while (new_positions[0] in occupied_positions) or (new_positions[1] in occupied_positions):
                    new_positions = draw_loop_extruder_loading_site(loop_extrusion_dynamics['chrlength'][0],distances)

                extruders_positions.append(new_positions)
                occupied_positions = occupied_positions + new_positions
                print("Occupied positions",sorted(occupied_positions))

                # Initialise the lifetime of each extruder and the flag to mark it they bumped on a barrier on the right or the left
                extruders_lifetimes.append(int(0))
                if lifetimeExponential == 1:
                    lifetime = int(lifetimeFactor*np.random.exponential(loop_extrusion_dynamics['lifetime']/lifetimeFactor, size=1)[0])
                    while lifetime == 0:
                        lifetime = int(lifetimeFactor*np.random.exponential(loop_extrusion_dynamics['lifetime']/lifetimeFactor, size=1)[0])
                else:
                    lifetime = loop_extrusion_dynamics['lifetime']
                extruders_target_lifetimes.append(lifetime)
                right_bumps.append(int(0))
                left_bumps.append(int(0))
            offset += loop_extrusion_dynamics['chrlength'][c]
        print("### END Positioning extruders ###")

        extruders_positions = sorted(extruders_positions, key=itemgetter(0))
        print("Initial extruders' positions (All %d)"      % (len(extruders_positions)),extruders_positions)        
        print("Initial extruders' lifetimes (Variable %d)" % (len(extruders_lifetimes)),extruders_lifetimes)
        print("Initial extruders' target lifetimes (Variable %d)" % (len(extruders_target_lifetimes)),extruders_lifetimes)
            
        lmp.command("compute xu all property/atom xu")
        lmp.command("compute yu all property/atom yu")
        lmp.command("compute zu all property/atom zu")
            
        print("Define the variable extruders")        
        for LES in range(int(run_time/loop_extrusion_dynamics['extrusion_time'])):
            print("### Extrusion round",LES,"###")
            thermo_style="thermo_style   custom   step temp epair emol"
            sys.stdout.flush()
            
            # Update the bond restraint for variable extruders
            variable_extruder_number = 0
            for particle1,particle2 in extruders_positions:                    
                variable_extruder_number += 1
                print("# fix LE%i all restrain bond %i  %i %f %f %f %f" % (variable_extruder_number,
                                                                        particle1,
                                                                        particle2,
                                                                        loop_extrusion_dynamics['attraction_strength'],
                                                                        loop_extrusion_dynamics['attraction_strength'],
                                                                        1.0,
                                                                        loop_extrusion_dynamics['equilibrium_distance']))
                
                lmp.command("fix LE%i all restrain bond %i  %i %f %f %f %f" % (variable_extruder_number,
                                                                                  particle1,
                                                                                  particle2,
                                                                                  loop_extrusion_dynamics['attraction_strength'],
                                                                                  loop_extrusion_dynamics['attraction_strength'],
                                                                                  1.0,
                                                                                  loop_extrusion_dynamics['equilibrium_distance']))
                lmp.command("variable x%i equal c_xu[%i]" % (particle1, particle1))
                lmp.command("variable x%i equal c_xu[%i]" % (particle2, particle2))
                lmp.command("variable y%i equal c_yu[%i]" % (particle1, particle1))
                lmp.command("variable y%i equal c_yu[%i]" % (particle2, particle2))
                lmp.command("variable z%i equal c_zu[%i]" % (particle1, particle1))
                lmp.command("variable z%i equal c_zu[%i]" % (particle2, particle2))
                
                lmp.command("variable xLE%i equal v_x%i-v_x%i" % (variable_extruder_number, particle1, particle2))
                lmp.command("variable yLE%i equal v_y%i-v_y%i" % (variable_extruder_number, particle1, particle2))
                lmp.command("variable zLE%i equal v_z%i-v_z%i" % (variable_extruder_number, particle1, particle2))
                lmp.command("variable dist_%i_%i equal sqrt(v_xLE%i*v_xLE%i+v_yLE%i*v_yLE%i+v_zLE%i*v_zLE%i)" % (particle1,
                                                                                                                 particle2,
                                                                                                                 variable_extruder_number,
                                                                                                                 variable_extruder_number,
                                                                                                                 variable_extruder_number,
                                                                                                                 variable_extruder_number,
                                                                                                                 variable_extruder_number,
                                                                                                                 variable_extruder_number))
                thermo_style += " v_dist_%i_%i" % (particle1, particle2)
            print("Defined",variable_extruder_number,"variable extruders")

            lmp.command("%s" % thermo_style)
            # Doing the LES
            if '1D_run' in loop_extrusion_dynamics:
                lmp.command("run 0")
            else:
                lmp.command("run %i" % loop_extrusion_dynamics['extrusion_time'])

            #exit(1)
            #lmp.command("fix f_unwr all store/state 1 xu yu zu")
            xc_tmp = np.array(lmp.extract_fix("f_unwr",1,2).contents[0:(natoms-2)])
            distances = compute_particles_distance(xc_tmp)
            #print(distances)
                
            # update the lifetime of each extruder
            for extruder in range(len(extruders_positions)):
                extruders_lifetimes[extruder] = extruders_lifetimes[extruder] + 1
            
            # Remove the bond restraints of variable extruders!
            loop_number = 1
            for particle1,particle2 in extruders_positions:
                #print("# unfix LE%i" % (loop_number))
                lmp.command("unfix LE%i" % (loop_number))

                loop_number += 1

            # Update the particles involved in the loop extrusion interaction:
            # decrease the initial_start by one until you get to start
            # increase the initial_stop by one until you get to stop
            extruders_to_relocate = [1]*len(extruders_positions) # 0 if the extruder needs to be relocated and 1 if it doesn't!
            force_extruders_to_relocate = [1]*len(extruders_positions) # 0 if the extruder needs to be relocated and 1 if it doesn't!
            for extruder in range(len(extruders_positions)):
                print("")
                print("#Moving extruder",extruder)
                
                # Keep in memory the current positions
                tmp_extruder = extruders_positions[extruder].copy()
                # Keep in memory the chromosome limits
                total = 0
                start = 1
                nchr  = 0
                for c in loop_extrusion_dynamics['chrlength']:
                    nchr += 1
                    stop  = start + c - 1
                    #print("Chromosome",nchr,"goes from bead",start,"to bead",stop)
                    if start <= extruders_positions[extruder][0] and extruders_positions[extruder][0] <= stop:
                        break                    
                    start = stop + 1                    
                print("Chromosome",nchr,"from bead",start,"to bead",stop,"includes the extruder of position",extruders_positions[extruder])            
                
                # 1. Propose a move of the extruder with probabilities "left_extrusion_rate' or 'right_extrusion_rate'
                # and distinguishing in linear or ring topology
                if loop_extrusion_dynamics['topology'] == "Linear":
                    random_number = uniform(0, 1)
                    if right_bumps[extruder] == 1:
                        random_number = 0
                    if random_number <= float(loop_extrusion_dynamics["left_extrusion_rate"]):       
                        # If the left part reaches the start of the chromosome, put the extruder in another position and re-initialize its lifetime -> Routine to relocate extruder
                        if extruders_positions[extruder][0] > start:
                            extruders_positions[extruder][0] -= 1
                            #print("Propose moving the left arm of the extruder(",random_number,"<=",loop_extrusion_dynamics["left_extrusion_rate"],")",extruder,"from",tmp_extruder[0],"to",extruders_positions[extruder][0])
                        else:
                            force_extruders_to_relocate[extruder] = 0 # 0 if the extruder needs to be relocated and 1 if it is doesn't!
                            print("Relocate the extruder",extruder,"because it reached the start of the chain",tmp_extruder[0])
                            
                    random_number = uniform(0, 1)
                    if left_bumps[extruder] == 1:
                        random_number = 0
                    if random_number <= float(loop_extrusion_dynamics["right_extrusion_rate"]):
                        # If the right part reaches the end of the chromosome, put the extruder in another position and re-initialize its lifetime -> Routine to relocate extruder
                        if extruders_positions[extruder][1] < stop:
                            extruders_positions[extruder][1] += 1
                            #print("Propose moving the right arm of the extruder(",random_number,"<=",loop_extrusion_dynamics["right_extrusion_rate"],")",extruder,"from",tmp_extruder[1],"to",extruders_positions[extruder][1])
                        else:
                            force_extruders_to_relocate[extruder] = 0 # 0 if the extruder needs to be relocated and 1 if it is doesn't!
                            print("Relocate the extruder",extruder,"because it reached the end of the chain",tmp_extruder[1])                            

                # Move the extruder if it doesn't hit the chromosome limits
                if loop_extrusion_dynamics['topology'] == "Ring":
                    # If the left part reaches the start of the chromosome -> Go to the end
                    random_number = uniform(0, 1)
                    if random_number <= loop_extrusion_dynamics['left_extrusion_rate']:
                        if tmp_extruder[0] > start:
                            extruders_positions[extruder][0] -= 1
                        elif tmp_extruder[0] == start:
                            extruders_positions[extruder][0] = stop
                        elif extruders_positions[extruder][0] == tmp_extruder[1]:
                            extruders_positions[extruder][0] = tmp_extruder[0]
                            
                    # If the right part reaches the end of the chromosome -> Go to start
                    random_number = uniform(0, 1)
                    if random_number <= loop_extrusion_dynamics['right_extrusion_rate']:
                        if tmp_extruder[1] <  stop:
                            extruders_positions[extruder][1] += 1
                        if tmp_extruder[1] == stop:
                            extruders_positions[extruder][1] = start
                        if extruders_positions[extruder][1] == tmp_extruder[0]:
                            extruders_positions[extruder][1] = tmp_extruder[1]
                    #print("Propose moving the extruder",extruder,"from",tmp_extruder,"to",extruders_positions[extruder])
                            
                # 2. If the extruder bumps into another extruder (fixed or variable)
                tmp_extruders_positions = [extruders_positions[x] for x in range(len(extruders_positions)) if x != extruder]

                occupied_positions = list(chain(*tmp_extruders_positions))
                #print("Proposed extruder positions",extruders_positions[extruder])
                #print("Occupied_positions (Extruder excluded)",occupied_positions)

                if loop_extrusion_dynamics["loop_extruders_encounter_rule"] == "crossing":
                    #Option extruders can cross each other
                    if extruders_positions[extruder][0] in occupied_positions:
                        print("The left  arm of the extruder",extruder,"bumped into an occupied position, but they can cross each other")
                    if extruders_positions[extruder][1] in occupied_positions:
                        print("The right arm of the extruder",extruder,"bumped into an occupied position, but they can cross each other")

                # 3. if the extruder reaches a barrier it is stop or not depending on the permeability of the barrier
                if   loop_extrusion_dynamics["loop_extruder_barrier_encounter_rule"] == "stalling":

                    if extruders_positions[extruder][0] in loop_extrusion_dynamics['barriers']:
                        perm = loop_extrusion_dynamics['barriers_right_permeability'][loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][0])]
                        if 'lifetimeAtBarriersRight' in loop_extrusion_dynamics:
                            barrierLifetimeRight = loop_extrusion_dynamics['lifetimeAtBarriersRight'][loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][0])]
                        # If the extruder tries to overcome a barrier we stop it with a probability > than the permeability
                        # If the barrier is on the left of the extruders, which is extruding contrary to the chain index, we have to re-put the extruder forwards
                        print("Found a barrier coming from the right at monomer %d with permeability %f" % (extruders_positions[extruder][0],perm))
                        # If the extruder was already stopped right_bumps[extruder] == 1, we have to block it
                        randomProb = uniform(0,1)
                        if right_bumps[extruder] == 1:
                            print("The extruder was already stopped at this barrier, it will stay here until the end of its lifetime ",extruders_target_lifetimes[extruder])
                            if 'stalling' in loop_extrusion_dynamics:
                                extruders_positions[extruder][0] = tmp_extruder[0]
                                extruders_positions[extruder][1] = tmp_extruder[1]
                            randomProb = 1.
                        if randomProb > perm:
                            # Option 1: The extruder stops at the barrier
                            if right_bumps[extruder] != 1:
                                extruders_target_lifetimes[extruder] = barrierLifetimeRight
                            right_bumps[extruder] = 1
                            extruders_to_relocate[extruder]      = 1  # 0 if the extruder needs to be relocated and 1 if it doesn't!
                            # Now, if the extruders can cross each other, they tend to concentrate at barriers, to avoid that we forbid an extruder to stall at an occupied barrier and relocate it
                            if loop_extrusion_dynamics["loop_extruders_encounter_rule"] == "crossing":
                                print("ExtrudeOnbarrier",extruderOnBarrierRight)
                                print("Barriers",loop_extrusion_dynamics['barriers'])
                                print("Index of the barrier",loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][0]))
                                if extruderOnBarrierRight[loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][0])] != extruder and extruderOnBarrierRight[loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][0])] != -1:
                                    print("When the extruders can cross each other, they tend to concentrate at barriers.")
                                    print("To avoid that we forbid an extruder to stall at an occupied barrier and relocate it.")
                                    print("Accordingly I will relocate extruder",extruder)
                                    extruders_to_relocate[extruder]      = 0  # 0 if the extruder needs to be relocated and 1 if it doesn't!
                                else:
                                    extruderOnBarrierRight[loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][0])] = extruder
                            print("The left  part of the extruder",extruder,"bumped in a barrier: the new position",extruders_positions[extruder][0],"is brought to the previous one",tmp_extruder[0])
                            extruders_positions[extruder][0] = tmp_extruder[0]

                    if extruders_positions[extruder][1] in loop_extrusion_dynamics['barriers']:
                        perm = loop_extrusion_dynamics['barriers_left_permeability'][loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][1])]
                        if 'lifetimeAtBarriersLeft' in loop_extrusion_dynamics:
                            barrierLifetimeLeft = loop_extrusion_dynamics['lifetimeAtBarriersLeft'][loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][1])]
                        # If the extruder tries to overcome a barrier we stop it with a probability > than the permeability
                        # If the barrier is on the right of the extruders, which is extruding with the chain index, we have to re-put the extruder backwards
                        print("Found a barrier coming from the left at monomer %d with permeability %f" % (extruders_positions[extruder][1],perm))
                        # If the extruder was already stopped left_bumps[extruder] == 1, we have to block it
                        randomProb = uniform(0,1)
                        if left_bumps[extruder] == 1:
                            print("The extruder was already stopped at this barrier, it will stay here until the end of its lifetime ",extruders_target_lifetimes[extruder])
                            if 'stalling' in loop_extrusion_dynamics:
                                extruders_positions[extruder][0] = tmp_extruder[0]
                                extruders_positions[extruder][1] = tmp_extruder[1]
                            randomProb = 1
                        if randomProb > perm:
                            # Option 1: The extruder stops at the barrier
                            if left_bumps[extruder] != 1:
                                extruders_target_lifetimes[extruder] = barrierLifetimeLeft
                            left_bumps[extruder] = 1
                            # Avoid to relocate the extruder!
                            extruders_to_relocate[extruder] = 1  # 0 if the extruder needs to be relocated and 1 if it doesn't!
                            # Now, if the extruders can cross each other, they tend to concentrate at barriers, to avoid that we forbid an extruder to stall at an occupied barrier and relocate it
                            if loop_extrusion_dynamics["loop_extruders_encounter_rule"] == "crossing":
                                print("ExtrudeOnbarrier",extruderOnBarrierLeft)
                                print("Barriers",loop_extrusion_dynamics['barriers'])
                                print("Index of the barrier",loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][1]))
                                if extruderOnBarrierLeft[loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][1])] != extruder and extruderOnBarrierLeft[loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][1])] != -1:
                                    print("When the extruders can cross each other, they tend to concentrate at barriers.")
                                    print("To avoid that we forbid an extruder to stall at an occupied barrier and relocate it.")
                                    print("Accordingly I will relocate extruder",extruder)
                                    extruders_to_relocate[extruder]      = 0  # 0 if the extruder needs to be relocated and 1 if it doesn't!
                                else:
                                    extruderOnBarrierLeft[loop_extrusion_dynamics['barriers'].index(extruders_positions[extruder][1])] = extruder
                            print("The right part of the extruder",extruder,"bumped in a barrier: the new position",extruders_positions[extruder][1],"is brought to the previous one",tmp_extruder[1])
                            extruders_positions[extruder][1] = tmp_extruder[1]

                # 4. If the extruder reached its lifetime, put it in another position and re-initialize its lifetime -> Routine to relocate extruder
                if extruders_lifetimes[extruder] >= extruders_target_lifetimes[extruder]:
                    extruders_to_relocate[extruder] = 0  # 0 if the extruder needs to be relocated and 1 if it doesn't!      

                # Routine to relocate extruder
                if extruders_to_relocate[extruder] == 0 or force_extruders_to_relocate[extruder] == 0:  # 0 if the extruder needs to be relocated and 1 if it doesn't!
                    print("Relocating extruder ",extruder," lifetime ",extruders_lifetimes[extruder]," Target-lifetime ",extruders_target_lifetimes[extruder])                    
                    tmp_extruders_positions = [extruders_positions[x] for x in range(len(extruders_positions)) if x != extruder]
                    occupied_positions = list(chain(*tmp_extruders_positions))+barriersOccupation
                    print("Occupied_positions (Extruder excluded)",sorted(occupied_positions))

                    extruders_positions[extruder]    = draw_loop_extruder_loading_site(loop_extrusion_dynamics['chrlength'][nchr-1],distances)
                    extruders_positions[extruder][0] = extruders_positions[extruder][0] + start
                    extruders_positions[extruder][1] = extruders_positions[extruder][1] + start
                    print("Proposed random position",extruders_positions[extruder])

                    while (extruders_positions[extruder][0] in occupied_positions) or (extruders_positions[extruder][1] in occupied_positions):
                        print("One of the proposed random positions",extruders_positions[extruder],"is occupied")
                        extruders_positions[extruder]    = draw_loop_extruder_loading_site(loop_extrusion_dynamics['chrlength'][nchr-1],distances)
                        extruders_positions[extruder][0] = extruders_positions[extruder][0] + start
                        extruders_positions[extruder][1] = extruders_positions[extruder][1] + start
                        print("Proposed random position",extruders_positions[extruder])

                    # Re-initialise the lifetime of the extruder
                    extruders_lifetimes[extruder] = 0
                    right_bumps[extruder] = 0
                    left_bumps[extruder] = 0
                    if lifetimeExponential == 1:
                        lifetime = int(lifetimeFactor*np.random.exponential(loop_extrusion_dynamics['lifetime']/lifetimeFactor, size=1)[0])
                        while lifetime == 0:
                            lifetime = int(lifetimeFactor*np.random.exponential(loop_extrusion_dynamics['lifetime']/lifetimeFactor, size=1)[0])
                    else:
                        lifetime = loop_extrusion_dynamics['lifetime']
                    extruders_target_lifetimes[extruder] = lifetime

                    # Re-initialise the occupation of barriers on the left, if the extruder was stalling there
                    for i in range(len(extruderOnBarrierLeft)):
                        if extruderOnBarrierLeft[i] == extruder:                            
                            extruderOnBarrierLeft[i] = -1
                            print("Barrier",loop_extrusion_dynamics['barriers'][i],"is now free of extruders at the left!")
                    # Re-initialise the occupation of barriers on the right, if the extruder was stalling there
                    for i in range(len(extruderOnBarrierRight)):
                        if extruderOnBarrierRight[i] == extruder:                            
                            extruderOnBarrierRight[i] = -1
                            print("Barrier",loop_extrusion_dynamics['barriers'][i],"is now free of extruders at the right!")

            print("Extruders positions at step",LES,extruders_positions)
            print("Extruders lifetimes at step",LES,extruders_lifetimes)
            print("Extruders target lifetimes at step",LES,extruders_target_lifetimes)
            print("Extruder on barriers left", LES,extruderOnBarrierLeft)
            print("Extruder on barriers right",LES,extruderOnBarrierRight)
            print("Left bumps at step",LES,left_bumps)
            print("Right bumps at step",LES,right_bumps)
            sys.stdout.flush()
            
    if not compartmentalization and not loop_extrusion_dynamics and run_time > 0:
        lmp.command("run %i" % run_time)

    if write_restart:
        # Save the final configuration: positions and velocities to restart the trajectory
        lmp.command("write_data final_conformation.txt nocoeff")

### END run_lammps ###

### Initialize lammps_run ###
def init_lammps_run(lmp,
                    initial_conformation,
                    neighbor,
                    timestep,
                    reset_timestep,
                    chromosome_particle_numbers,
                    connectivity,
                    excludedVolume,
                    persistence_length):

    """
    Initialise the parameters for the computation in lammps job

    :param lmp: lammps instance object.
    :param initial_conformation: lammps input data file with the particles initial conformation.
    :param neighbor: see LAMMPS_CONFIG.py.
    :param timestep: see LAMMPS_CONFIG.py.
    :param connectivity  : use a FENE bond or harmonic for harmonic potential for neighbours
    :param excludedVolume: use a LJ for purely-repulsive Lennard-Jones interation
    :param CONFIG.persistence_length persistence_length : see LAMMPS_CONFIG.py.
    """

    #######################################################
    # Box and units  (use LJ units and period boundaries) #
    #######################################################
    lmp.command("units %s" % CONFIG.units)
    lmp.command("atom_style %s" % CONFIG.atom_style) #with stiffness
    lmp.command("boundary %s" % CONFIG.boundary)

    ##########################
    # READ "start" data file #
    ##########################
    if not initial_conformation:    
        initial_conformation = lammps_folder+'initial_conformation.dat'
        generate_chromosome_random_walks_conformation (chromosome_particle_numbers,
                                                       outfile=initial_conformation,
                                                       seed_of_the_random_number_generator=randint(1,100000),
                                                       confining_environment=confining_environment)

    lmp.command("read_data %s" % initial_conformation)        
    lmp.command("mass %s" % CONFIG.mass)

    ##################################################################
    # Pair interactions require lists of neighbours to be calculated #
    ##################################################################
    lmp.command("neighbor %s" % neighbor)
    lmp.command("neigh_modify %s" % CONFIG.neigh_modify)
    
    ##############################################################
    # Sample thermodynamic info  (temperature, energy, pressure) #
    ##############################################################
    lmp.command("thermo %i" % CONFIG.thermo)
    
    ###############################
    # Stiffness term              #
    # E = K * (1+cos(theta)), K>0 #
    ###############################
    lmp.command("angle_style %s" % CONFIG.angle_style)
    lmp.command("angle_coeff * %f" % persistence_length)
    
    ###################################################################
    # Pair interaction between non-bonded atoms                       #
    #                                                                 #
    #  Lennard-Jones 12-6 potential with cutoff:                      #
    #  potential E=4epsilon[ (sigma/r)^12 - (sigma/r)^6]  for r<r_cut #
    #  r_cut =1.12246 = 2^(1/6) is the minimum of the potential       #
    ###################################################################    
    lmp.command("pair_style hybrid/overlay lj/cut %f morse 0.0" % CONFIG.PurelyRepulsiveLJcutoff)

    ################################################################
    #  pair_modify shift yes adds a constant to the potential such #
    #  that E(r_cut)=0. Forces remains unchanged.                  #
    ################################################################
    lmp.command("pair_modify shift yes")
    
    for atomTypePair in excludedVolume:
    
        ######################################
        #  pair_coeff for lj/cut, specify 4: #
        #    * atom type interacting with    #
        #    * atom type                     #
        #    * epsilon (energy units)        #
        #    * sigma (distance units)        #
        ######################################
        lmp.command("pair_coeff %d %d lj/cut %f %f %f" % (atomTypePair[0],
                                                          atomTypePair[1],
                                                          excludedVolume[atomTypePair][1],
                                                          excludedVolume[atomTypePair][2],
                                                          excludedVolume[atomTypePair][3]))
        
    lmp.command("pair_coeff * * morse  %f %f %f" % (0.0, 0.0, 0.0))

    for bondType in connectivity:
        if connectivity[bondType][0] == "FENE":
            #########################################################
            # Pair interaction between bonded atoms                 #
            #                                                       #
            # Fene potential + Lennard Jones 12-6:                  #
            #  E= - 0.5 K R0^2 ln[ 1- (r/R0)^2]                     #
            #     + 4epsilon[ (sigma/r)^12 - (sigma/r)^6] + epsilon #
            #########################################################
            lmp.command("bond_style fene")
            
            ########################################
            # For style fene, specify:             #
            #   * bond type                        #
            #   * K (energy/distance^2)            #
            #   * R0 (distance)                    #
            #   * epsilon (energy)  (LJ component) #
            #   * sigma (distance)  (LJ component) #
            ########################################
            lmp.command("bond_coeff %d %f %f %f %f" % (bondType, connectivity[bondType][1], connectivity[bondType][2], connectivity[bondType][3], connectivity[bondType][4]))
            lmp.command("special_bonds fene") #<=== I M P O R T A N T (new command)
            
        if connectivity[bondType][0] == "harmonic":
            lmp.command("bond_style harmonic")
            lmp.command("bond_coeff * 100.0 1.2")
        if connectivity[bondType][0] == "FENEspecial":
            lmp.command("bond_style fene")        
            lmp.command("bond_coeff %d %f %f %f %f" % (bondType, connectivity[bondType][1], connectivity[bondType][2], connectivity[bondType][3], connectivity[bondType][4]))
            lmp.command("special_bonds fene") #<=== I M P O R T A N T (new command)

    ##############################
    # set timestep of integrator #
    ##############################
    lmp.command("timestep %f" % timestep)

    ######################################
    # reset the counter of the timesteps #
    ######################################
    lmp.command("reset_timestep %d" % (reset_timestep)) 

### Initialize lammps_run ###



########## Part to generate the initial conformation ##########
def generate_chromosome_random_walks_conformation ( chromosome_particle_numbers ,
                                                    confining_environment = None ,
                                                    particle_radius=0.5 ,
                                                    seed_of_the_random_number_generator=1 ,
                                                    number_of_conformations=1,
                                                    outfile="Initial_random_walk_conformation.dat",
                                                    atom_types=None,
                                                    pbc=False,
                                                    center=True):
    """
    Generates lammps initial conformation file by random walks
    
    :param chromosome_particle_numbers: list with the number of particles of each chromosome.    
    :param ['sphere',100.] confining_environment: dictionary with the confining environment of the conformation
            Possible confining environments:
            ['cube',edge_width]
            ['sphere',radius]
    :param 0.5 particle_radius: Radius of each particle.
    :param 1 seed_of_the_random_number_generator: random seed.
    :param 1 number_of_conformations: copies of the conformation.
    :param outfile: file where to store resulting initial conformation file

    """
    seed(seed_of_the_random_number_generator)
    
    # This allows to organize the largest chromosomes first.
    # This is to get a better acceptance of the chromosome positioning.
    chromosome_particle_numbers = [int(x) for x in chromosome_particle_numbers]
    chromosome_particle_numbers.sort(key=int,reverse=True)

    for cnt in range(number_of_conformations):

        final_random_walks = generate_random_walks(chromosome_particle_numbers,
                                                   particle_radius,
                                                   confining_environment,
                                                   pbc=pbc,
                                                   center=center)

        # Writing the final_random_walks conformation
        #print "Succesfully generated conformation number %d\n" % (cnt+1)
        write_initial_conformation_file(final_random_walks,
                                        chromosome_particle_numbers,
                                        confining_environment,
                                        atom_types=atom_types,
                                        out_file=outfile)

##########
        
def generate_chromosome_rosettes_conformation ( chromosome_particle_numbers ,
                                                fractional_radial_positions=None,
                                                confining_environment=['sphere',100.] ,
                                                rosette_radius=12.0 , particle_radius=0.5 ,
                                                seed_of_the_random_number_generator=1 ,
                                                number_of_conformations=1,
                                                outfile = "Initial_rosette_conformation.dat",
                                                atom_types=1):
    """
    Generates lammps initial conformation file by rosettes conformation
    
    :param chromosome_particle_numbers: list with the number of particles of each chromosome.    
    :param None fractional_radial_positions: list with fractional radial positions for all the chromosomes.
    :param ['sphere',100.] confining_environment: dictionary with the confining environment of the conformation
            Possible confining environments:
            ['cube',edge_width]
            ['sphere',radius]
    :param 0.5 particle_radius: Radius of each particle.
    :param 1 seed_of_the_random_number_generator: random seed.
    :param 1 number_of_conformations: copies of the conformation.
    :param outfile: file where to store resulting initial conformation file

    """
    seed(seed_of_the_random_number_generator)
    
    # This allows to organize the largest chromosomes first.
    # This is to get a better acceptance of the chromosome positioning.
    chromosome_particle_numbers = [int(x) for x in chromosome_particle_numbers]
    chromosome_particle_numbers.sort(key=int,reverse=True)    

    initial_rosettes , rosettes_lengths = generate_rosettes(chromosome_particle_numbers,
                                                            rosette_radius,
                                                            particle_radius)
    print(rosettes_lengths)

    
    # Constructing the rosettes conformations
    for cnt in range(number_of_conformations):

        temptative = 0
        particle_inside   = 0 # 0 means a particle is outside
        particles_overlap = 0 # 0 means two particles are overlapping
        while particle_inside == 0 or particles_overlap == 0:
            temptative += 1
            print("Temptative number %d" % temptative)
            particle_inside   = 1
            particles_overlap = 1
            segments_P1 = []
            segments_P0 = []
            side = 0
            init_rosettes = copy.deepcopy(initial_rosettes)
            
            # Guess of the initial segment conformation:
            # 1 - each rod is placed inside the confining evironment
            # in a random position and with random orientation
            # 2 - possible clashes between generated rods are checked
            if fractional_radial_positions:
                if len(fractional_radial_positions) != len(chromosome_particle_numbers):
                    print("Please provide the desired fractional radial positions for all the chromosomes")
                    sys.exit()
                segments_P1 , segments_P0 = generate_rods_biased_conformation(rosettes_lengths, rosette_radius,
                                                                              confining_environment,
                                                                              fractional_radial_positions,
                                                                              max_number_of_temptative=1000)
            else:
                segments_P1 , segments_P0 = generate_rods_random_conformation(rosettes_lengths, rosette_radius,
                                                                              confining_environment,
                                                                              max_number_of_temptative=1000)

            # Roto-translation of the rosettes according to the segment position and orientation 
            final_rosettes = rosettes_rototranslation(init_rosettes, segments_P1, segments_P0)
            
            # Checking that the beads are all inside the confining environment and are not overlapping
            for r in range(len(final_rosettes)):
                molecule0 = list(zip(final_rosettes[r]['x'],final_rosettes[r]['y'],final_rosettes[r]['z']))
                print(len(molecule0),len(molecule0[0]))
                if particle_inside == 0:
                    break
                for i in range(len(molecule0)):
                    # Check if the particle is inside the confining_environment

                    particle_inside = check_point_inside_the_confining_environment(molecule0[i][0],
                                                                                   molecule0[i][1],
                                                                                   molecule0[i][2],
                                                                                   1.0,
                                                                                   confining_environment)
                    if particle_inside == 0:
                        break

            if particle_inside == 1:
                for rosette_pair in list(combinations(final_rosettes,2)):
                    molecule0 = list(zip(rosette_pair[0]['x'],rosette_pair[0]['y'],rosette_pair[0]['z']))
                    molecule1 = list(zip(rosette_pair[1]['x'],rosette_pair[1]['y'],rosette_pair[1]['z']))
                    distances = spatial.distance.cdist(molecule1,molecule0)
                    print("Different chromosomes",len(molecule0),len(molecule0[0]),distances.min())
                    if distances.min() < particle_radius*2.0*0.95:
                        particles_overlap = 0
                        break

                if particles_overlap != 0:
                    for r in range(len(final_rosettes)):
                        molecule0 = list(zip(final_rosettes[r]['x'],final_rosettes[r]['y'],final_rosettes[r]['z']))
                        print(len(molecule0),len(molecule0[0]))

                        distances = spatial.distance.cdist(molecule0,molecule0)
                        print("Same chromosome",distances.min())
                        sys.stdout.flush()
                        for i in range(len(molecule0)):
                            for j in range(i+1,len(molecule0)):
                                if distances[(i,j)] < particle_radius*2.0*0.95:
                                    particles_overlap = 0
                                    print("Particles",i,"and",j,"are contacting",distances[(i,j)])
                                    sys.stdout.flush()
                                if particles_overlap == 0:
                                    break
                            if particles_overlap == 0:
                                break
                        if particles_overlap == 0:
                            break
            
        # Writing the final_rosettes conformation
        print("Succesfully generated conformation number %d\n" % (cnt+1))
        write_initial_conformation_file(final_rosettes,
                                        chromosome_particle_numbers,
                                        confining_environment,
                                        out_file=outfile,
                                        atom_types=atom_types)

##########
        
def generate_chromosome_rosettes_conformation_with_pbc ( chromosome_particle_numbers ,
                                                         fractional_radial_positions=None,
                                                         confining_environment=['cube',100.] ,
                                                         rosette_radius=12.0 , particle_radius=0.5 ,
                                                         seed_of_the_random_number_generator=1 ,
                                                         number_of_conformations=1,
                                                         outfile = "Initial_rosette_conformation_with_pbc.dat",
                                                         atom_types=1, k=6., x=0.38, p=1.0):
    """
    Generates lammps initial conformation file by rosettes conformation
    
    :param chromosome_particle_numbers: list with the number of particles of each chromosome.    
    :param None fractional_radial_positions: list with fractional radial positions for all the chromosomes.
    :param ['cube',100.] confining_environment: dictionary with the confining environment of the conformation
            Possible confining environments:
            ['cube',edge_width]
    :param 0.5 particle_radius: Radius of each particle.
    :param 1 seed_of_the_random_number_generator: random seed.
    :param 1 number_of_conformations: copies of the conformation.
    :param outfile: file where to store resulting initial conformation file

    """
    seed(seed_of_the_random_number_generator)
    
    # This allows to organize the largest chromosomes first.
    # This is to get a better acceptance of the chromosome positioning.
    chromosome_particle_numbers = [int(x) for x in chromosome_particle_numbers]
    chromosome_particle_numbers.sort(key=int,reverse=True)    

    initial_rosettes , rosettes_lengths = generate_rosettes(chromosome_particle_numbers,
                                                            rosette_radius,
                                                            particle_radius,
                                                            k=k, x=x, p=p)
    print(rosettes_lengths)

    
    # Constructing the rosettes conformations
    for cnt in range(number_of_conformations):

        particles_overlap = 0 # 0 means two particles are overlapping
        while particles_overlap == 0:
            particles_overlap = 1
            segments_P1 = []
            segments_P0 = []
            side = 0
            init_rosettes = copy.deepcopy(initial_rosettes)
            
            # Guess of the initial segment conformation:
            # 1 - each rod is placed in a random position and with random orientation
            # 2 - possible clashes between generated rods are checked taking into account pbc
            segments_P1 , segments_P0 = generate_rods_random_conformation_with_pbc (
                rosettes_lengths, 
                rosette_radius,
                confining_environment,
                max_number_of_temptative=100000)

            # Roto-translation of the rosettes according to the segment position and orientation 
            final_rosettes = rosettes_rototranslation(init_rosettes, segments_P1, segments_P0)
            
            # Checking that the beads once folded inside the simulation box (for pbc) are not overlapping
            folded_rosettes = copy.deepcopy(final_rosettes)
            for r in range(len(folded_rosettes)):
                particle = 0
                for x, y, z in zip(folded_rosettes[r]['x'],folded_rosettes[r]['y'],folded_rosettes[r]['z']):
                    #inside_1 = check_point_inside_the_confining_environment(x, y, z,
                    #                                                        particle_radius,
                    #                                                        confining_environment)
                    #if inside_1 == 0:
                    #    print inside_1, r, particle, x, y, z
                    
                    while x >  (confining_environment[1]*0.5):
                        x -=  confining_environment[1]
                    while x < -(confining_environment[1]*0.5):
                        x +=  confining_environment[1]
                            
                    while y >  (confining_environment[1]*0.5):
                        y -=  confining_environment[1]
                    while y < -(confining_environment[1]*0.5):
                        y +=  confining_environment[1]
                                    
                    while z >  (confining_environment[1]*0.5):
                        z -=  confining_environment[1]
                    while z < -(confining_environment[1]*0.5):
                        z +=  confining_environment[1]

                    #inside_2 = check_point_inside_the_confining_environment(x, y, z,
                    #                                                      particle_radius,
                    #                                                      confining_environment)
                    #if inside_2 == 1 and inside_1 == 0:
                    #    print inside_2, r, particle, x, y, z
                    folded_rosettes[r]['x'][particle] = x
                    folded_rosettes[r]['y'][particle] = y
                    folded_rosettes[r]['z'][particle] = z
                    particle += 1

            print(len(folded_rosettes))
            sys.stdout.flush()
            if len(folded_rosettes) > 1:
                for rosette_pair in list(combinations(folded_rosettes,2)):
                    molecule0 = list(zip(rosette_pair[0]['x'],rosette_pair[0]['y'],rosette_pair[0]['z']))
                    molecule1 = list(zip(rosette_pair[1]['x'],rosette_pair[1]['y'],rosette_pair[1]['z']))
                    distances = spatial.distance.cdist(molecule1,molecule0)
                    print(len(molecule0),len(molecule0[0]),"Minimum distance",distances.min())
                    if distances.min() < particle_radius*2.0*0.95:
                        particles_overlap = 0
                        break

            if particles_overlap != 0:
                for r in range(len(folded_rosettes)):
                    molecule0 = list(zip(folded_rosettes[r]['x'],folded_rosettes[r]['y'],folded_rosettes[r]['z']))
                    print(len(molecule0),len(molecule0[0]))

                    distances = spatial.distance.cdist(molecule0,molecule0)
                    print("Minimum distance", distances.min())
                    for i in range(len(molecule0)):
                        for j in range(i+1,len(molecule0)):
                            if distances[(i,j)] < particle_radius*2.0*0.7:
                                particles_overlap = 0
                                print("Particles",i,"and",j,"are contacting",distances[(i,j)])
                                sys.stdout.flush()
                            if particles_overlap == 0:
                                break
                        if particles_overlap == 0:
                            break
                    if particles_overlap == 0:
                        break

            
        # Writing the final_rosettes conformation
        print("Succesfully generated conformation number %d\n" % (cnt+1))
        write_initial_conformation_file(final_rosettes,
                                        chromosome_particle_numbers,
                                        confining_environment,
                                        out_file=outfile,
                                        atom_types=atom_types)

##########

def generate_rosettes(chromosome_particle_numbers, rosette_radius, particle_radius, k=6., x=0.38, p=1.0):
    # Genaration of the rosettes
    # A. Rosa et al. https://doi.org/10.1371/journal.pcbi.1000153
    # List to contain the rosettes and the rosettes lengths
    rosettes = []
    rosettes_lengths  = []

    for number_of_particles in chromosome_particle_numbers:
        
        # Variable to build the chain
        phi = 0.0

        # Dictory of lists to contain the rosette
        rosette      = {}
        rosette['x'] = []
        rosette['y'] = []
        rosette['z'] = []        
        
        # Position of the first particle (x_0, 0.0, 0.0)
        rosette['x'].append(rosette_radius * (x + (1 - x) * cos(k*phi) * cos(k*phi)) * cos(phi))
        rosette['y'].append(rosette_radius * (x + (1 - x) * cos(k*phi) * cos(k*phi)) * sin(phi))
        rosette['z'].append(p * phi / (2.0 * pi))

        # Building the chain: The rosette is growing along the positive z-axes
        for particle in range(1,number_of_particles):

            distance = 0.0
            while distance < (particle_radius*2.0): 
                phi   = phi + 0.001
                x_tmp = rosette_radius * (x + (1 - x) * cos(k*phi) * cos(k*phi)) * cos(phi)
                y_tmp = rosette_radius * (x + (1 - x) * cos(k*phi) * cos(k*phi)) * sin(phi)
                z_tmp = phi / (2.0 * pi)     
                distance  = sqrt((x_tmp - rosette['x'][-1])*(x_tmp - rosette['x'][-1]) +
                                 (y_tmp - rosette['y'][-1])*(y_tmp - rosette['y'][-1]) +
                                 (z_tmp - rosette['z'][-1])*(z_tmp - rosette['z'][-1]))

            rosette['x'].append(x_tmp)
            rosette['y'].append(y_tmp)
            rosette['z'].append(z_tmp)
            if distance > ((particle_radius*2.0)*1.2):
                print("%f %d %d %d" % (distance, particle-1, particle))
            
        rosettes.append(rosette)
        rosettes_lengths.append(rosette['z'][-1]-rosette['z'][0])
        
    return rosettes , rosettes_lengths

##########

def generate_rods_biased_conformation(rosettes_lengths, rosette_radius,
                                      confining_environment,
                                      fractional_radial_positions,
                                      max_number_of_temptative=100000):
    # Construction of the rods initial conformation 
    segments_P0 = []
    segments_P1 = []

    if confining_environment[0] != 'sphere':
        print("ERROR: Biased chromosome positioning is currently implemented")
        print("only for spherical confinement. If you need other shapes, please")
        print("contact the developers")
    
    for length , target_radial_position in zip(rosettes_lengths,fractional_radial_positions):
        tentative            = 0
        clashes              = 0 # 0 means that there is an clash -> PROBLEM
        best_radial_position = 1.0
        best_radial_distance = 1.0
        best_segment_P0      = []
        best_segment_P1      = []
        
        # Positioning the rods
        while tentative < 100000 and best_radial_distance > 0.00005:                

            print("Length = %f" % length)

            print("Trying to position terminus 0")
            segment_P0_tmp = []
            segment_P0_tmp = draw_point_inside_the_confining_environment(confining_environment,
                                                                         rosette_radius)
            print("Successfully positioned terminus 0: %f %f %f" % (segment_P0_tmp[0], segment_P0_tmp[1], segment_P0_tmp[2]))
            
            print("Trying to position terminus 1")
            segment_P1_tmp = []                            
            segment_P1_tmp = draw_second_extreme_of_a_segment_inside_the_confining_environment(segment_P0_tmp[0],
                                                                                               segment_P0_tmp[1],
                                                                                               segment_P0_tmp[2],
                                                                                               length,
                                                                                               rosette_radius,
                                                                                               confining_environment)
            print("Successfully positioned terminus 1: %f %f %f" % (segment_P1_tmp[0], segment_P1_tmp[1], segment_P1_tmp[2]))

            # Check clashes with the previously positioned rods
            clashes = 1
            for segment_P1,segment_P0 in zip(segments_P1,segments_P0):
                clashes = check_segments_clashes(segment_P1,
                                                 segment_P0,
                                                 segment_P1_tmp,
                                                 segment_P0_tmp,
                                                 rosette_radius)
                if clashes == 0:
                    break                

            if clashes == 1:
                # Check whether the midpoint of the segment is close to the target radial position
                segment_midpoint = []
                segment_midpoint.append((segment_P0_tmp[0] + segment_P1_tmp[0])*0.5)
                segment_midpoint.append((segment_P0_tmp[1] + segment_P1_tmp[1])*0.5)
                segment_midpoint.append((segment_P0_tmp[2] + segment_P1_tmp[2])*0.5)

                radial_position = sqrt( ( segment_midpoint[0] * segment_midpoint[0] +
                                          segment_midpoint[1] * segment_midpoint[1] +
                                          segment_midpoint[2] * segment_midpoint[2] ) /
                                        (confining_environment[1]*confining_environment[1]))

                radial_distance = fabs(radial_position-target_radial_position)

                print(radial_position , target_radial_position , radial_distance , best_radial_distance , tentative)
                
                # If the midpoint of the segment is closer to the target radial position than the
                # previous guesses. Store the points as the best guesses!
                if radial_distance < best_radial_distance:                    
                    best_radial_distance = radial_distance
                    best_radial_position = radial_position
                    best_tentative       = tentative+1 # The variable tentative starts from 0

                    best_segment_P0 = []
                    best_segment_P1 = []
                    for component_P0 , component_P1 in zip(segment_P0_tmp,segment_P1_tmp):
                        best_segment_P0.append(component_P0)
                        best_segment_P1.append(component_P1)
                    
                tentative = tentative + 1
                
        if best_segment_P0 == []:
            print("Valid placement not found for chromosome rosette after %d tentatives" % tentative)
            sys.exit()

        print("Successfully positioned chromosome of length %lf at tentative %d of %d tentatives" % (length, best_tentative, tentative))        
        segments_P0.append(best_segment_P0)
        segments_P1.append(best_segment_P1)

    print("Successfully generated rod conformation!")
    return segments_P1 , segments_P0
    
##########

def generate_rods_random_conformation(rosettes_lengths, rosette_radius,
                                      confining_environment,
                                      max_number_of_temptative=1000):
    # Construction of the rods initial conformation 
    segments_P0 = []
    segments_P1 = []
    
    for length in rosettes_lengths:
        tentative = 0
        clashes   = 0
        # Random positioning of the rods
        while tentative < max_number_of_temptative and clashes == 0:                

            tentative += 1
            clashes    = 1
            #print "Length = %f" % length

            print("Trying to position terminus 0")
            #pick uniformly within the confining environment using the rejection method 
            first_point = []
            first_point = draw_point_inside_the_confining_environment(confining_environment,
                                                                      rosette_radius)

            print("Successfully positioned terminus 0: %f %f %f" % (first_point[0], first_point[1], first_point[2]))
            
            print("Trying to position terminus 1")
            #pick from P0 another point one the sphere of radius length inside the confining environment
            last_point = []
            last_point = draw_second_extreme_of_a_segment_inside_the_confining_environment(first_point[0],
                                                                                           first_point[1],
                                                                                           first_point[2],
                                                                                           length,
                                                                                           rosette_radius,
                                                                                           confining_environment)
            
            print("Successfully positioned terminus 1: %f %f %f" % (last_point[0], last_point[1], last_point[2]))
                
            # Check clashes with the previously positioned rods
            clashes = 1 
            for segment_P1,segment_P0 in zip(segments_P1,segments_P0):
                clashes = check_segments_clashes(segment_P1,
                                                 segment_P0,
                                                 last_point,
                                                 first_point,
                                                 rosette_radius)
                if clashes == 0:
                    break                

            #print clashes
        print("Successfully positioned chromosome of length %lf at tentative %d\n" % (length, tentative))        
        segments_P1.append(last_point)
        segments_P0.append(first_point)            

    print("Successfully generated rod conformation!")
    return segments_P1 , segments_P0

##########

def generate_rods_random_conformation_with_pbc(rosettes_lengths, rosette_radius,
                                               confining_environment,
                                               max_number_of_temptative=100000):

    # Construction of the rods initial conformation 
    segments_P0 = []
    segments_P1 = []
    
    for length in rosettes_lengths:
        tentative = 0
        clashes   = 0
        # Random positioning of the rods
        while tentative < 100000 and clashes == 0:                

            tentative += 1
            clashes    = 1
            #print "Length = %f" % length

            print("Trying to position terminus 0")
            #pick uniformly within the confining environment using the rejection method 
            first_point = []
            first_point = draw_point_inside_the_confining_environment(confining_environment,
                                                                      rosette_radius)

            print("Successfully positioned terminus 0: %f %f %f" % (first_point[0], first_point[1], first_point[2]))
            
            print("Trying to position terminus 1")
            #pick from P0 another point one the sphere of radius length inside the confining environment
            last_point = []
            last_point = draw_second_extreme_of_a_segment(first_point[0],
                                                          first_point[1],
                                                          first_point[2],
                                                          length,
                                                          rosette_radius)            
            
            print(last_point)
            # Check clashes with the previously positioned rods
            for segment_P1,segment_P0 in zip(segments_P1,segments_P0):
                clashes = check_segments_clashes_with_pbc(segment_P1,
                                                          segment_P0,
                                                          last_point,
                                                          first_point,
                                                          rosette_radius,
                                                          confining_environment)
                if clashes == 0:
                    break                

            #print clashes
        print("Successfully positioned chromosome of length %lf at tentative %d\n" % (length, tentative))        
        segments_P1.append(last_point)
        segments_P0.append(first_point)            

    print("Successfully generated rod conformation!")
    return segments_P1 , segments_P0

##########

def generate_random_walks(chromosome_particle_numbers,
                          particle_radius,
                          confining_environment,
                          center,
                          pbc=False):
    # Construction of the random walks initial conformation 
    random_walks = []
    
    for number_of_particles in chromosome_particle_numbers:
        #print "Trying to position random walk"

        #print "Positioning first particle"            
        particles_overlap = 0
        while particles_overlap == 0:
            random_walk      = {}
            random_walk['x'] = []
            random_walk['y'] = []
            random_walk['z'] = []        

            particles_overlap = 1

            # Generate the first particle
            first_particle = []
            first_particle = draw_point_inside_the_confining_environment(confining_environment,
                                                                         particle_radius)
            random_walk['x'].append(first_particle[0])        
            random_walk['y'].append(first_particle[1])
            random_walk['z'].append(first_particle[2])

            # Generate all the others N-1 particles
            for particle in range(1,number_of_particles):
                #print("Positioning particle %d" % (particle+1))
                sys.stdout.flush()
                particles_overlap = 0 # 0 means that there is an overlap -> PROBLEM
                overlapCounter = -1
                maxIter = 100000000
                overlapCounter += 1
                if overlapCounter > maxIter:
                    # raise error so log file is created to avoid k_seed
                    overlapCounter = -1
                    #errorName = 'ERROR: Initial conformation non ending loop %s' % confining_environment
                    #raise InitalConformationError(errorName)
                particles_overlap = 1
                new_particle = []
                if pbc:
                    new_particle = draw_second_extreme_of_a_segment(
                        random_walk['x'][-1],
                        random_walk['y'][-1],
                        random_walk['z'][-1],
                        2.0*particle_radius,
                        2.0*particle_radius)
                else:
                    new_particle = draw_second_extreme_of_a_segment_inside_the_confining_environment(
                        random_walk['x'][-1],
                        random_walk['y'][-1],
                        random_walk['z'][-1],
                        2.0*particle_radius,
                        2.0*particle_radius,
                        confining_environment)
                random_walk['x'].append(new_particle[0])        
                random_walk['y'].append(new_particle[1])
                random_walk['z'].append(new_particle[2])

            print("Check if there is a particle overlapping with any other particle in the system")
            print("Within random-walks")
            molecule0 = list(zip(random_walk['x'],random_walk['y'],random_walk['z']))
            print(len(molecule0),len(molecule0[0]))
            sys.stdout.flush()
            
            distances = spatial.distance.cdist(molecule0,molecule0)
            print(distances.min())
            """
            for i in range(len(molecule0)):
                for j in range(i+1,len(molecule0)):
                    if distances[(i,j)] < particle_radius*2.0*0.9:
                        print(i,j,distances[(i,j)])
                        particles_overlap = 0
                        if particles_overlap == 0:
                            break
                    if particles_overlap == 0:
                        break
            """
            print("Between Random-walks")
            if particles_overlap != 0:
                for random_walk_pair in list(combinations(random_walks,2)):
                    molecule0 = list(zip(random_walk_pair[0]['x'],random_walk_pair[0]['y'],random_walk_pair[0]['z']))
                    molecule1 = list(zip(random_walk_pair[1]['x'],random_walk_pair[1]['y'],random_walk_pair[1]['z']))
                    distances = spatial.distance.cdist(molecule1,molecule0)
                    print(len(molecule0),len(molecule0[0]),distances.min())
                    sys.stdout.flush()
                    if distances.min() < particle_radius*2.0*0.9:
                        particles_overlap = 0
                        break
            print(particles_overlap)
                
                    
        print("Successfully positioned random walk of %d particles" % number_of_particles)
        random_walks.append(random_walk)

    print("Successfully generated random walk conformation!")
    if center:
        print("Centering random-walk into the origin")
        for random_walk in random_walks:
            x_com, y_com, z_com = (0.0,0.0,0.0)
            cnt = 0
            for (x,y,z) in zip(random_walk['x'],random_walk['y'],random_walk['z']):
                x_com += x
                y_com += y
                z_com += z
                cnt += 1
            x_com, y_com, z_com = (x_com/cnt,y_com/cnt,z_com/cnt)

            for i in range(len(random_walk['x'])):
                random_walk['x'][i] -= x_com
                random_walk['y'][i] -= y_com
                random_walk['z'][i] -= z_com
            
            x_com, y_com, z_com = (0.0,0.0,0.0)
            cnt = 0
            for (x,y,z) in zip(random_walk['x'],random_walk['y'],random_walk['z']):
                x_com += x
                y_com += y
                z_com += z
                cnt += 1
            x_com, y_com, z_com = (x_com/cnt,y_com/cnt,z_com/cnt)
            
    return random_walks

##########

def check_particle_vs_all_overlap(x,y,z,chromosome,overlap_radius):    
    particle_overlap = 1

    for x0, y0, z0 in zip(chromosome['x'],chromosome['y'],chromosome['z']):
        particle_overlap = check_particles_overlap(x0,y0,z0,x,y,z,overlap_radius)
        if particle_overlap == 0:
            return particle_overlap
        
    return particle_overlap
            
##########

def draw_second_extreme_of_a_segment_inside_the_confining_environment(x0, y0, z0, 
                                                                      segment_length, 
                                                                      object_radius, 
                                                                      confining_environment,
                                                                      max_number_of_temptative=1000):
    inside = 0

    temptative = 0
    while inside == 0 and temptative < max_number_of_temptative:
        temptative += 1

        particle = []
        temp_theta  = arccos(2.0*random()-1.0)
        temp_phi    = 2*pi*random()
        particle.append(x0 + segment_length * cos(temp_phi) * sin(temp_theta))
        particle.append(y0 + segment_length * sin(temp_phi) * sin(temp_theta))
        particle.append(z0 + segment_length * cos(temp_theta))
        # Check if the particle is inside the confining_environment
        inside = check_point_inside_the_confining_environment(particle[0],
                                                              particle[1],
                                                              particle[2],
                                                              object_radius,
                                                              confining_environment)

    return particle

##########

def draw_second_extreme_of_a_segment(x0, y0, z0, 
                                     segment_length, 
                                     object_radius):
    particle = []
    temp_theta  = arccos(2.0*random()-1.0)
    temp_phi    = 2*pi*random()
    particle.append(x0 + segment_length * cos(temp_phi) * sin(temp_theta))
    particle.append(y0 + segment_length * sin(temp_phi) * sin(temp_theta))
    particle.append(z0 + segment_length * cos(temp_theta))
    
    return particle

##########

def draw_point_inside_the_confining_environment(confining_environment, object_radius):
    #pick a point uniformly within the confining environment using the rejection method 

    print(confining_environment)
    if confining_environment[0] == 'cube':
        dimension_x = confining_environment[1] * 0.5
        dimension_y = confining_environment[1] * 0.5
        dimension_z = confining_environment[1] * 0.5        
        if len(confining_environment) > 2:
            print("# WARNING: Defined a cubical confining environment with reduntant paramenters.")
            print("# Only 2 are needed the identifier and the side")
        inside = 0
        while inside == 0:
            particle = []
            particle.append((2.0*random()-1.0)*(dimension_x - object_radius))
            particle.append((2.0*random()-1.0)*(dimension_y - object_radius))
            particle.append((2.0*random()-1.0)*(dimension_z - object_radius))
            # Check if the particle is inside the confining_environment
            inside = check_point_inside_the_confining_environment(particle[0],
                                                                  particle[1],
                                                                  particle[2],
                                                                  object_radius,
                                                                  confining_environment)                                
        
    if confining_environment[0] == 'sphere':
        dimension_x = confining_environment[1]-object_radius
        dimension_y = confining_environment[1]-object_radius
        dimension_z = confining_environment[1]-object_radius
        if len(confining_environment) > 2:
            print("# WARNING: Defined a spherical confining environment with reduntant parameters.")
            print("# Only 2 are needed the identifier and the radius")
        
        inside = 0
        while inside == 0:
            particle = []
            
            particle.append((2.0*random()-1.0)*(dimension_x))
            particle.append((2.0*random()-1.0)*(dimension_y))
            particle.append((2.0*random()-1.0)*(dimension_z))
            # Check if the particle is inside the confining_environment
            inside = check_point_inside_the_confining_environment(particle[0],
                                                                  particle[1],
                                                                  particle[2],
                                                                  object_radius,
                                                                  confining_environment)
            print(particle)
            print(object_radius)
            
    return particle
    
##########        

def check_point_inside_the_confining_environment(Px, Py, Pz,
                                                 object_radius,
                                                 confining_environment):
    # The shapes are all centered in the origin
    # - sphere    : radius r
    # - cube      : side

    if confining_environment[0] == 'sphere':
        radius = confining_environment[1] - object_radius
        if ((Px*Px)/(radius*radius) + (Py*Py)/(radius*radius) + (Pz*Pz)/(radius*radius)) < 1.0 : return 1

    if confining_environment[0] == 'cube':
        hside = confining_environment[1] * 0.5 - object_radius
        if (((Px*Px)/(hside*hside)) < 1.0) and (((Py*Py)/(hside*hside)) < 1.0) and (((Pz*Pz)/(hside*hside)) < 1.0) : return 1

    return 0

##########

def check_segments_clashes(s1_P1, s1_P0, s2_P1, s2_P0, rosette_radius):

    # Check steric clashes without periodic boundary conditions
    if distance_between_segments(s1_P1, s1_P0, s2_P1, s2_P0) < 2.0*rosette_radius:
        # print "Clash between segments",s1_P1,s1_P0,"and",s2_P1_tmp,s2_P0_tmp,"at distance", distance
        return 0

    return 1

##########

def check_segments_clashes_with_pbc(s1_P1, s1_P0, s2_P1, s2_P0, 
                                    rosette_radius,
                                    confining_environment):

    # Check steric clashes with periodic boundary conditions
    if distance_between_segments(s1_P1, s1_P0, s2_P1, s2_P0) < 2.0*rosette_radius:
        # print "Clash between segments",s1_P1,s1_P0,"and",s2_P1_tmp,s2_P0_tmp,"at distance", distance
        return 0

    return 1

##########

def distance_between_segments(s1_P1, s1_P0, s2_P1, s2_P0):

    # Inspiration: http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm 
    # Copyright 2001, softSurfer (www.softsurfer.com)
    # This code may be freely used and modified for any purpose
    # providing that this copyright notice is included with it.
    # SoftSurfer makes no warranty for this code, and cannot be held
    # liable for any real or imagined damage resulting from its use.
    # Users of this code must verify correctness for their application.

    u  = []
    v  = []
    w  = []
    dP = []

    for c_s1_P1,c_s1_P0,c_s2_P1,c_s2_P0 in zip(s1_P1, s1_P0, s2_P1, s2_P0):        
        u.append(c_s1_P1 - c_s1_P0)
        v.append(c_s2_P1 - c_s2_P0)
        w.append(c_s1_P0 - c_s2_P0)
    
    a  = scalar_product(u, u)
    b  = scalar_product(u, v)
    c  = scalar_product(v, v)
    d  = scalar_product(u, w)
    e  = scalar_product(v, w)

    D  = a*c - b*b
    sD = tD = D
        
    if D < (1.0e-7):
        # Segments almost parallel 
        sN = 0.0
        sD = 1.0
        tN = e
        tD = c
    else:
        # Get the closest points on the infinite lines
        sN = (b*e - c*d)
        tN = (a*e - b*d)
        if (sN < 0.0):            
            # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD: # sc > 1 => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c

    if tN < 0.0: # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # Recompute sc for this edge
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a        
    
    elif tN > tD: # tc > 1 => the t=1 edge is visible
        tN = tD
        # Recompute sc for this edge
        if (-d + b) < 0.0:
            sN = 0
        elif (-d + b) > a:
            sN = sD;
        else:
            sN = (-d + b)
            sD = a

    # Finally do the division to get sc and tc
    if abs(sN) < (1.0e-7):
        sc = 0.0
    else:
        sc = sN / sD
        
    if abs(tN) < (1.0e-7):
        tc = 0.0
    else:
        tc = tN / tD
     
    # Get the difference of the two closest points
    for i in range(len(w)):    
        dP.append(w[i] + ( sc * u[i] ) - ( tc * v[i] )) # = S1(sc) - S2(tc)
    
    return norm(dP)   # return the closest distance

##########

def rosettes_rototranslation(rosettes, segments_P1, segments_P0):

    for i in range(len(segments_P1)):
        vector = []
        theta  = []

        for component_P1 , component_P0 in zip(segments_P1[i], segments_P0[i]):
            vector.append(component_P1-component_P0)
            
        # Rotation Angles
        theta.append(atan2(vector[1],vector[2]))
      
        x_temp_2 =  vector[0]
        y_temp_2 =  cos(theta[0]) * vector[1] - sin(theta[0]) * vector[2]                
        z_temp_2 =  sin(theta[0]) * vector[1] + cos(theta[0]) * vector[2]        
        theta.append(atan2(x_temp_2,z_temp_2))
        
        x_temp_1 =  cos(theta[1]) * x_temp_2 - sin(theta[1]) * z_temp_2
        y_temp_1 =  y_temp_2
        z_temp_1 =  sin(theta[1]) * x_temp_2 + cos(theta[1]) * z_temp_2
        
        if(z_temp_1 < 0.0):            
            z_temp_1 = -z_temp_1
            theta.append(pi)
        else:
            theta.append(0.0)
        #print x_temp_1 , y_temp_1 , z_temp_1 
        
        # Chromosome roto-translations
        for particle in range(len(rosettes[i]['x'])):

            x_temp_2 =   rosettes[i]['x'][particle]
            y_temp_2 =   cos(theta[2]) * rosettes[i]['y'][particle] + sin(theta[2]) * rosettes[i]['z'][particle]
            z_temp_2 = - sin(theta[2]) * rosettes[i]['y'][particle] + cos(theta[2]) * rosettes[i]['z'][particle]
            
            x_temp_1 =   cos(theta[1]) * x_temp_2 + sin(theta[1]) * z_temp_2
            y_temp_1 =   y_temp_2
            z_temp_1 = - sin(theta[1]) * x_temp_2 + cos(theta[1]) * z_temp_2

            x =   x_temp_1;
            y =   cos(theta[0]) * y_temp_1 + sin(theta[0]) * z_temp_1;
            z = - sin(theta[0]) * y_temp_1 + cos(theta[0]) * z_temp_1;
            
            # Chromosome translations
            rosettes[i]['x'][particle] = segments_P0[i][0] + x;
            rosettes[i]['y'][particle] = segments_P0[i][1] + y;
            rosettes[i]['z'][particle] = segments_P0[i][2] + z;
    return rosettes

##########

def scalar_product(a, b):

    scalar = 0.0
    for c_a,c_b in zip(a,b):
        scalar = scalar + c_a*c_b 

    return scalar

##########

def norm(a):

    return sqrt(scalar_product(a, a))

##########

def write_initial_conformation_file(chromosomes,
                                    chromosome_particle_numbers,
                                    confining_environment,
                                    out_file="Initial_conformation.dat",
                                    atom_types=1,
                                    angle_types=1,
                                    bond_types=1):
    # Choosing the appropriate xlo, xhi...etc...depending on the confining environment
    xlim = []
    ylim = []
    zlim = []
    if confining_environment[0] == 'sphere':
        radius = confining_environment[1] + 1.0
        xlim.append(-radius)
        xlim.append(radius)
        ylim.append(-radius)
        ylim.append(radius)
        zlim.append(-radius)
        zlim.append(radius)
        
    if confining_environment[0] == 'cube':
        hside = confining_environment[1] * 0.5
        xlim.append(-hside)
        xlim.append(hside)
        ylim.append(-hside)
        ylim.append(hside)
        zlim.append(-hside)
        zlim.append(hside)
            
    fileout = open(out_file,'w')
    n_chr=len(chromosomes)
    n_atoms=0
    for n in chromosome_particle_numbers:
        n_atoms+=n    
        
    fileout.write("LAMMPS input data file \n\n")
    fileout.write("%9d atoms\n" % (n_atoms))
    fileout.write("%9d bonds\n" % (n_atoms-n_chr))
    fileout.write("%9d angles\n\n" % (n_atoms-2*n_chr))
    fileout.write("%9s atom types\n" % atom_types)
    fileout.write("%9s bond types\n" % bond_types)
    fileout.write("%9s angle types\n\n" % angle_types)
    fileout.write("%6.3lf    %6.3lf     xlo xhi\n" % (xlim[0], xlim[1]))
    fileout.write("%6.3lf    %6.3lf     ylo yhi\n" % (ylim[0], ylim[1]))
    fileout.write("%6.3lf    %6.3lf     zlo zhi\n" % (zlim[0], zlim[1]))
  
    fileout.write("\n Atoms \n\n")
    particle_number = 1
    for chromosome in chromosomes:
        for x,y,z in zip(chromosome['x'],chromosome['y'],chromosome['z']):          
            fileout.write("%-8d %s %s %7.4lf %7.4lf %7.4lf\n" % (particle_number, "1", "1", x, y, z))
            particle_number += 1
  
    fileout.write("\n Bonds \n\n")
    bond_number          = 1
    first_particle_index = 1
    for chromosome in chromosomes:
        for i in range(len(chromosome['x'])-1):
            fileout.write("%-4d %s %4d %4d\n" % (bond_number, "1", first_particle_index, first_particle_index+1))
            bond_number          += 1
            first_particle_index += 1
        first_particle_index += 1 # I have to go to the end of the chromosome!
          
    fileout.write("\n Angles \n\n")
    angle_number         = 1
    first_particle_index = 1
    for chromosome in chromosomes:        
        for i in range(len(chromosome['x'])-2):
            fileout.write("%-4d %s %5d %5d %5d\n" % (angle_number, "1", first_particle_index, first_particle_index+1, first_particle_index+2))
            angle_number         += 1    
            first_particle_index += 1
        first_particle_index += 2 # I have to go to the end of the chromosome!

    fileout.close()

##########

def distance(x0,y0,z0,x1,y1,z1):
    return sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))

##########

def check_particles_overlap(x0,y0,z0,x1,y1,z1,overlap_radius):
    if distance(x0,y0,z0,x1,y1,z1) < overlap_radius:
        #print "Particle %f %f %f and particle %f %f %f are overlapping\n" % (x0,y0,z0,x1,y1,z1)
        return 0
    return 1

##### Loop extrusion dynamics functions #####
def read_target_loops_input(input_filename, chromosome_length, percentage):
    # Open input file
    fp_input = open(input_filename, "r")

    loops = []
    target_loops = []
    # Get each loop per line and fill the output list of loops
    for line in fp_input.readlines():

        if line.startswith('#'):
            continue

        splitted = line.strip().split()
        loop = []
        loop.append(int(splitted[1]))
        loop.append(int(splitted[2]))

        loops.append(loop)

    #ntarget_loops = int(len(loops)*percentage/100.)    
    ntarget_loops = int(len(loops))    
    shuffle(loops)
    target_loops = loops[0:ntarget_loops]

    return target_loops

##########

def draw_loop_extruder_loading_site(chromosome_length,distances):

    # draw a starting point for extrusion along the chromosome
    random_particle =  randint(1, chromosome_length-2)

    return [random_particle,random_particle+1]

##########

def compute_particles_distance(xc):
    
    particles = []
    distances = {}

    # Getting the coordinates of the particles
    for i in range(0,len(xc),3):
        x = xc[i]  
        y = xc[i+1]
        z = xc[i+2]
        particles.append((x, y, z))

    # Checking whether the restraints are satisfied
    for pair in combinations(range(len(particles)), 2):
        dist = distance(particles[pair[0]][0],
                        particles[pair[0]][1],
                        particles[pair[0]][2],
                        particles[pair[1]][0],
                        particles[pair[1]][1],
                        particles[pair[1]][2])
        distances[pair] = dist

    return distances

##########

def get_list(input_list):

    output_list = []
    #print(input_list)
    
    for element in input_list:
        #print(type(element))
        if isinstance(element, (int)):
            output_list.append(element)
        if isinstance(element, (list)):
            for subelement in element:
                output_list.append(subelement)
        if isinstance(element, (tuple)):
            for subelement in range(element[0],element[1]+1,element[2]):
                output_list.append(subelement)
    return output_list
