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
        for window in windows:
            writePath = "insulationProfiles_%s_w%sbp_r%dbp.tab" % (condition, window, resolution)    
            print(writePath)
            with open(writePath, 'w') as f:
                dfAsString = insulation_table.to_string(index=False)
                f.write(dfAsString)
        continue    
        histkwargs = dict(
                bins=10**np.linspace(-4,1,200),
                histtype='step',
                lw=2,
        )

        # Plot of boundary strength distribution
        f, ax = plt.subplots(figsize=(18, 3))
        thresholds_li = {}
        thresholds_otsu = {}
        w=window
        ax.hist(
            insulation_table[f'boundary_strength_{w}'],
            **histkwargs
        )
        thresholds_li[w] = threshold_li(insulation_table[f'boundary_strength_{w}'].dropna().values)
        thresholds_otsu[w] = threshold_otsu(insulation_table[f'boundary_strength_{w}'].dropna().values)
        n_boundaries_li = (insulation_table[f'boundary_strength_{w}'].dropna()>=thresholds_li[w]).sum()
        n_boundaries_otsu = (insulation_table[f'boundary_strength_{w}'].dropna()>=thresholds_otsu[w]).sum()
        ax.axvline(thresholds_li[w], c='green')
        ax.axvline(thresholds_otsu[w], c='magenta')
        ax.text(0.01, 0.9,
                f'Window {w//1000}kb',
                ha='left',
                va='top',
                transform=ax.transAxes)
        ax.text(0.01, 0.7,
                f'{n_boundaries_otsu} boundaries (Otsu)',
                c='magenta',
                ha='left',
                va='top',
                transform=ax.transAxes)
        ax.text(0.01, 0.5,
                f'{n_boundaries_li} boundaries (Li)',
                c='green',
                ha='left',
                va='top',
                transform=ax.transAxes)
            
        ax.set(
            xscale='log',
            ylabel='# boundaries'
        )

        ax.set(xlabel='Boundary strength')
        plt.savefig('boundary_strength_%s.pdf' % (condition))
            
        ### columns of this insulation dataframe report the insulation score,
        first_window_summary =insulation_table.columns[[ str(window) in i for i in insulation_table.columns]]
    
        insulation_table[['chrom','start','end','region','is_bad_bin']+list(first_window_summary)].iloc[1000:1005]
            
        #start plotting
        plt.rcParams['font.size'] = 12
            
        start = 53700000
        end   = 56700000
        region = ('chr18', start, end)
        norm = LogNorm(vmax=0.1, vmin=0.00001)
        data = clr.matrix(balance=True).fetch(region)
        f, ax = plt.subplots(figsize=(18, 3))
        #im = pcolormesh_45deg(ax, data, start=region[1], resolution=resolution, norm=norm, cmap='fall')
        #ax.set_aspect(0.5)
        #ax.set_ylim(0, 1)
        ax.set_xlim(53700000, 56000000)    
        #format_ticks(ax, rotate=False)
        ax.xaxis.set_visible(False)
            
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="1%", pad=0.1, aspect=2)
        #plt.colorbar(im, cax=cax)
            
        insul_region = bioframe.select(insulation_table, region)
            
        #print(insul_region)
            
        #ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)
        ins_ax = ax
        ins_ax.set_prop_cycle(plt.cycler("color", plt.cm.plasma(np.linspace(0,1,5))))
        ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
                    -insul_region['log2_insulation_score_'+str(window)],
                    label=f'Window {window} bp')
            
        #ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);
        
        #format_ticks(ins_ax, y=False, rotate=False)
        #ax.set_xlim(region[1], region[2])
        
        #for res in windows[1:]:
        #    ins_ax.plot(insul_region[['start', 'end']].mean(axis=1), insul_region[f'log2_insulation_score_{res}'], label=f'Window {res} bp')
        #ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);
        
        #The insulation table also has annotations for valleys of the
        #insulation score, which correspond to highly insulating regions,
        #such as TAD boundaries. All potential boundaries have an assigned
        #boundary_strength_ column. Additionally, this strength is
        #thresholded to find regions that insulate particularly strongly,
        #and this is recorded in the is_boundary_ columns.
        
        print(insul_region[~np.isnan(insul_region['log2_insulation_score_'+str(window)])]['log2_insulation_score_'+str(window)])
        print(find_peaks(-1 * insul_region[~np.isnan(insul_region['log2_insulation_score_'+str(window)])]['log2_insulation_score_'+str(window)]))
        print(find_peak_prominence(-1 * insul_region[~np.isnan(insul_region['log2_insulation_score_'+str(window)])]['log2_insulation_score_'+str(window)]))    
        exit(1)
    
        #Let’s repeat the previous plot and show where we found the boundaries:
        boundaries = insul_region[~np.isnan(insul_region[f'boundary_strength_{window}'])]
        weak_boundaries = boundaries[~boundaries[f'is_boundary_{window}']]
        writePath = "weakBoundaries_%s_w%sbp_r%dbp.tab" % (condition, window, resolution)    
        with open(writePath, 'w') as f:
            dfAsString = weak_boundaries.to_string(index=False)
            f.write(dfAsString)    
            strong_boundaries = boundaries[boundaries[f'is_boundary_{window}']]
            writePath = "strongBoundaries_%s_w%sbp_r%dbp.tab" % (condition, window, resolution)    
            with open(writePath, 'w') as f:
                dfAsString = strong_boundaries.to_string(index=False)
                f.write(dfAsString)    
        #print(weak_boundaries)
        #print(strong_boundaries)    
        ins_ax.scatter(weak_boundaries[['start', 'end']].mean(axis=1),
                       -weak_boundaries[f'log2_insulation_score_{window}'], label='Weak boundaries', c="yellow")
        ins_ax.scatter(strong_boundaries[['start', 'end']].mean(axis=1),
                       -strong_boundaries[f'log2_insulation_score_{window}'], label='Strong boundaries', c="cyan")
        
        ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4);
        
        format_ticks(ins_ax, y=False, rotate=False)
        ax.set_xlim(region[1], region[2])
            
    print("end")
            
    #plt.show()
    plt.savefig('plot_insulation_%s_w%skb_r%skb.pdf' % (condition,int(window/1000),int(resolution/1000)))


