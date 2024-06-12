# Insulation & boundaries
# https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html 

#insulation can be computed with multiple methods. One of the most common methods involves using a diamond-window score to generate an *insulation profile*.

# import standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess

# Import python package for working with cooler files and tools for analysis
import cooler
import cooltools
import cooltools.lib.plotting
from cooltools import insulation
from cooltools.lib.peaks import find_peak_prominence
import bioframe
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import find_peaks

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.1'):
    raise AssertionError("tutorials rely on cooltools version 0.5.1 or higher,"+
                         "please check your cooltools version and update to the latest")

from skimage.filters import threshold_li, threshold_otsu



# Functions to help with plotting
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im

from matplotlib.ticker import EngFormatter
bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)


##IMPORT data set
for resolution in [2000, 5000]:

    #for condition in ["NPC_WT"]:
    for condition in ["mESC_DAB","mESC_DA","mESC_DB","mESC_dCTCF","mESC_Dprom","mESC_polyA","mESC_UPPax6","mESC_WT","NPC_DAB","NPC_DA","NPC_DB","NPC_dCTCF","NPC_Dprom","NPC_polyA","NPC_WT"]:

        print(condition)
        flag = 0
        #We first load the Hi-C data at 10 kbp resolution.
        mcoolFile = "../01_cool_files/chic_%s_merge_mm10_IJ.mcool" % (condition)
        clr = cooler.Cooler(mcoolFile + '::resolutions/'+str(resolution))

        # We then compute insulation genome-wide and write it to a file
        
        windows = range(150000, 500001, 50000)
        for window in windows:
            writePath = "insulationProfiles_%s_w%sbp_r%dbp.tab" % (condition, window, resolution)                
            if os.path.exists(writePath):
                continue
            fp = open(writePath, 'w')
            fp.close()
            print(writePath)
            flag = 1 
        if flag == 0:
            continue
        
        insulation_table = insulation(clr, windows, verbose=True)
        #for window in windows:
        writePath = "insulationProfiles_%s_w%sbp.tab" % (condition, resolution)    
        print(writePath)
        with open(writePath, 'w') as f:
            dfAsString = insulation_table.to_string(index=False)
            f.write(dfAsString)
