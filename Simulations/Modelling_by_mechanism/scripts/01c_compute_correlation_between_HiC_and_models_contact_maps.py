"""
19 Jul 2013
"""
import matplotlib
matplotlib.use('Agg')

from subprocess                       import Popen, PIPE
from math                             import acos, degrees, pi, sqrt
from warnings                         import warn

from numpy                            import exp as np_exp
from numpy                            import median as np_median
from numpy                            import mean as np_mean
from numpy                            import std as np_std, log2
from numpy                            import array, cross, dot, ma, isnan
from numpy                            import histogram, linspace
from numpy                            import zeros
from numpy.linalg                     import norm

from scipy.optimize                   import curve_fit
from scipy.stats                      import spearmanr, pearsonr, chisquare
from scipy.stats                      import linregress
from scipy.stats                      import normaltest, norm as sc_norm
from scipy.cluster.hierarchy          import linkage, fcluster

import sys

try:
    from matplotlib import pyplot as plt
    from matplotlib.cm import jet, bwr
    matplotlib.use('Agg')
except ImportError:
    warn('matplotlib not found\n')


def correlate_with_real_data(model_contact_matrix, exp_interaction_matrix, model_contact_matrix_OoE, exp_interaction_matrix_OoE, nloci, cutoff, corr, 
                             off_diag, timestep, resolution, LSthreshold = 100000, log_corr=True, axe=None, savefig=None):

    """
    Plots the result of a correlation between a given group of models and
    original Hi-C data.
    
    :param None models: if None (default) the correlation will be computed
    using all the models. A list of numbers corresponding to a given set
    of models can be passed
    :param None cluster: compute the correlation only for the models in the
    cluster number 'cluster'
    :param None cutoff: distance cutoff (nm) to define whether two particles
    are in contact or not, default is 2 times resolution, times scale.
    :param None savefig: path to a file where to save the image generated;
    if None, the image will be shown using matplotlib GUI (the extension
    of the file name will determine the desired format).
    :param False plot: to display the plot
    :param True log_corr: log plot for correlation
    :param None axe: a matplotlib.axes.Axes object to define the plot
    appearance
    :param None contact_matrix: input a contact matrix instead of computing
    it from the models
    
    :returns: correlation coefficient rho, between the two
    matrices. A rho value greater than 0.7 indicates a very good
    correlation
    
    """

    expdata       = []
    moddata       = []
    expdata_OoE   = []
    moddata_OoE   = []
    expdata_long  = []
    moddata_long  = []
    expdata_short = []
    moddata_short = []
    for i in range(len(exp_interaction_matrix)):
        for j in range(i + off_diag, len(exp_interaction_matrix)):
            if exp_interaction_matrix[i][j] == -10.:
                exp_interaction_matrix[i][j] = 0.
                continue            
            expdata.append(exp_interaction_matrix[i][j])
            moddata.append(model_contact_matrix[i][j])
            expdata_OoE.append(exp_interaction_matrix_OoE[i][j])
            moddata_OoE.append(model_contact_matrix_OoE[i][j])
            d = abs((i-j)*resolution) 
            if d > LSthreshold:
                expdata_long.append(exp_interaction_matrix[i][j])
                moddata_long.append(model_contact_matrix[i][j])
            else:
                expdata_short.append(exp_interaction_matrix[i][j])
                moddata_short.append(model_contact_matrix[i][j])

    if corr == 'spearman':
        corr    = spearmanr(moddata, expdata)
        corrOoE = spearmanr(moddata_OoE, expdata_OoE)
        corrL   = spearmanr(moddata_long , expdata_long)
        corrS   = spearmanr(moddata_short, expdata_short)
    #elif corr == 'pearson':
    #    corr = pearsonr(moddata, expdata)
    #elif corr == 'logpearson':
    #    corr = pearsonr(nozero_log_list(moddata), nozero_log_list(expdata))
    #elif corr == 'chi2':
    #    corr = chisquare(array(moddata), array(expdata))
    #    corr = 1. / corr[0], corr[1]
    #else:
    #    raise NotImplementedError('ERROR: %s not implemented, must be one ' +
    #'of spearman, pearson or frobenius\n')
    if not savefig:
        print("All matrix",                               "timestep", timestep, "dcutoff", dcutoff, "SCC", corr[0],  "npoints", len(expdata))
        print("Short (<%dMb))" % (int(LSthreshold/1000000)), "timestep", timestep, "dcutoff", dcutoff, "SCC", corrS[0], "npoints", len(expdata_short))
        print("Long (>%dMb)"   % (int(LSthreshold/1000000)), "timestep", timestep, "dcutoff", dcutoff, "SCC", corrL[0], "npoints", len(expdata_long))
        #print("ObsoverExp",                               "timestep", timestep, "dcutoff", dcutoff, "SCC", corrOoE[0], "npoints", len(expdata_OoE))
        return corr, corrS, corrL, corrOoE
    if not axe:
        fig = plt.figure(figsize=(20, 4.5))
    else:
        fig = axe.get_figure()
    fig.suptitle('Correlation between normalized-real and modeled ' +
                 'contact maps (SCC=%.4f, SCC(<%dMb)=%.4f, SCC(>%dMb)=%.4f, SCC(OoE)=%.4f)' % (corr[0],int(LSthreshold/1000000),corrS[0],int(LSthreshold/1000000),corrL[0],corrOoE[0]),
                 size='x-large')
    ax = fig.add_subplot(131)

    # Visualization put to zero the entries below LSthreshold!!!
    # imshow of the modeled data
    #contact_map(model_contact_matrix, cutoff, nloci, axe=ax)
    cmap = plt.get_cmap('jet')
    cmap.set_bad('darkgrey', 1)   
    exp_max=0.
    exp_min=1000000.
    mod_max=0.
    mod_min=1000000.
    for i in range(len(exp_interaction_matrix)):
        for j in range(len(exp_interaction_matrix)):
            d = abs((i-j)*resolution)
            if d <= LSthreshold:
               exp_interaction_matrix[i][j] = 0.0
               model_contact_matrix[i][j]   = 0.0
               continue
            if model_contact_matrix[i][j] > mod_max:
                mod_max=model_contact_matrix[i][j]
            if model_contact_matrix[i][j] < mod_min:
                mod_min=model_contact_matrix[i][j]

            if exp_interaction_matrix[i][j] > exp_max:
                exp_max=exp_interaction_matrix[i][j]
            if exp_interaction_matrix[i][j] < exp_min:
                exp_min=exp_interaction_matrix[i][j]

    for i in range(len(exp_interaction_matrix)):
        for j in range(len(exp_interaction_matrix)):
            d = abs((i-j)*resolution)
            if d <= LSthreshold:
               exp_interaction_matrix[i][j] = 0.0
               model_contact_matrix[i][j]   = 0.0
            model_contact_matrix[i][j]   = (model_contact_matrix[i][j]  -mod_min) / (mod_max-mod_min)
            exp_interaction_matrix[i][j] = (exp_interaction_matrix[i][j]-exp_min) / (exp_max-exp_min)
            if model_contact_matrix[i][j] < 0.0625:
                model_contact_matrix[i][j] = 0.0625
            if exp_interaction_matrix[i][j] < 0.0625:
                exp_interaction_matrix[i][j] = 0.0625

    ims = ax.imshow(log2(model_contact_matrix), origin='lower',
                    interpolation="nearest", cmap=cmap,
                    extent=(0.5, nloci + 0.5, 0.5, nloci + 0.5))
    ax.set_ylabel('Particles')
    ax.set_xlabel('Particles')
    # From the estimate of the time mapping we have than 17 tau_LJ is about 1 s of trajectory in real time
    ax.set_title('Models contacts count at timestep %s (~%.1lf min)' % (timestep,float(timestep/17/60.)))
    #ax.set_title('Model contacts count')
    cbar = ax.figure.colorbar(ims)
    cbar.ax.set_ylabel('Log2 (Model contacts count)')

    # correlation
    ax = fig.add_subplot(132)
    try:
        if log_corr:
            minmoddata = float(min([m for m in moddata if m]))
            minexpdata = float(min([m for m in expdata if m]))
            moddata, expdata = (log2([(m if m else minmoddata / 2) * 100 for m in moddata]),
                                log2([m if m else minexpdata / 2 for m in expdata]))
    except:
        warn('WARNING: unable to log for correlation with real data...')
    slope, intercept, r_value, p_value, _ = linregress(moddata, expdata)
    # slope, intercept, r_value, p_value, std_err = linregress(moddata, expdata)
    midplot = 'hexbin'
    if midplot == 'classic':
        lnr = ax.plot(moddata, intercept + slope * array (moddata), color='k',
                      ls='--', alpha=.7)
        ax.legend(lnr, ['p-value: %.3f, R: %.3f' % (p_value, r_value)])
        ax.plot(moddata, expdata, 'ro', alpha=0.5)
        ax.set_xlabel('Modelled data')
        ax.set_ylabel('Real data')
    elif midplot == 'hexbin':
        hb = ax.hexbin(moddata, expdata, mincnt=1,
                       gridsize=50, cmap=plt.cm.Spectral_r)
        lnr = ax.plot(moddata, intercept + slope * array (moddata), color='k',
                      ls='--', alpha=.7)
        ax.set_xlabel('Models contacts counts at dcutoff %d nm' % (cutoff))
        ax.set_ylabel('Normalized Hi-C count')
        cbaxes = fig.add_axes([0.41, 0.42, 0.005, 0.45])
        cbar = plt.colorbar(hb, cax=cbaxes)  # orientation='horizontal')
        cbar.set_label('Number of particle pairs')
    elif midplot == 'triple':
        maxval = max(expdata)
        minval = min(expdata)
        ax.set_visible(False)
        axleft = fig.add_axes([0.42, 0.18, 0.1, 0.65])
        axleft.spines['right'].set_color('none')
        axleft.spines['bottom'].set_color('none')
        axleft.spines['left'].set_smart_bounds(True)
        axleft.spines['top'].set_smart_bounds(True)
        axleft.xaxis.set_ticks_position('top')
        axleft.yaxis.set_ticks_position('left')
        axleft.set_ylabel('Normalized Hi-C count')
        axleft.patch.set_visible(False)
        axbott = fig.add_axes([0.44, 0.13, 0.17, 0.5])
        axbott.spines['left'].set_color('none')
        axbott.spines['top'].set_color('none')
        axbott.spines['left'].set_smart_bounds(True)
        axbott.spines['bottom'].set_smart_bounds(True)
        axbott.xaxis.set_ticks_position('bottom')
        axbott.yaxis.set_ticks_position('right')
        axbott.patch.set_visible(False)
        axbott.set_xlabel('Models contacts count')
        axmidl = fig.add_axes([0.44, 0.18, 0.17, 0.65])
        axbott.hist(moddata, bins=20, alpha=.2)
        x, _  = histogram([i if str(i) != '-inf' else 0. for i in expdata],
                          bins=20)
        axleft.barh(linspace(minval, maxval, 20), x,
                    height=(maxval - minval) / 20, alpha=.2)
        axleft.set_ylim((minval -
                         (maxval - minval) / 20, maxval +
                         (maxval - minval) / 20))
        axmidl.plot(moddata, expdata, 'k.', alpha=.3)
        axmidl.plot(moddata, intercept + slope * array (moddata), color='k',
                    ls='--', alpha=.7)
        axmidl.set_ylim(axleft.get_ylim())
        axmidl.set_xlim(axbott.get_xlim())
        axmidl.axis('off')
        # axmidl.patch.set_visible(False)
    ax.set_title('Real versus modelled data')
    ax = fig.add_subplot(133)
    cmap = plt.get_cmap('jet')
    cmap.set_bad('darkgrey', 1)
    ims = ax.imshow(log2(exp_interaction_matrix), origin='lower',
                    interpolation="nearest", cmap=cmap,
                    extent=(0.5, nloci + 0.5, 0.5, nloci + 0.5))
    ax.set_ylabel('Particles')
    ax.set_xlabel('Particles')
    ax.set_title('Normalized Hi-C count')
    cbar = ax.figure.colorbar(ims)
    cbar.ax.set_ylabel('Log2 (normalized Hi-C data)')
    plt.close('all')

    print("All matrix",                               "timestep", timestep, "dcutoff", dcutoff, "SCC", corr[0],  "npoints", len(expdata))
    print("Short (<int(%d/1000000))" % (LSthreshold), "timestep", timestep, "dcutoff", dcutoff, "SCC", corrS[0], "npoints", len(expdata_short))
    print("Long (>int(%d/1000000))" % (LSthreshold),  "timestep", timestep, "dcutoff", dcutoff, "SCC", corrL[0], "npoints", len(expdata_long))
    #print("ObsoverExp",                               "timestep", timestep, "dcutoff", dcutoff, "SCC", corrOoE[0], "npoints", len(expdata_OoE))

    return corr[0], corrS[0], corrL[0], corrOoE[0]

# Getting input data
#print sys.argv
size       = int(sys.argv[1])
dcutoff    = float(sys.argv[2])
dt         = sys.argv[3]
timestep   = int(float(sys.argv[3])*0.006)
resolution = int(sys.argv[4])
genDist    = float(sys.argv[5])

model_matrix = zeros((size,size))
exp_matrix   = zeros((size,size))

model_matrix_OoE = zeros((size,size))
exp_matrix_OoE   = zeros((size,size))

for i in range(size):
    for j in range(size):
        model_matrix[i][j] = -10.
        exp_matrix[i][j]   = -10.
        model_matrix_OoE[i][j] = -10.
        exp_matrix_OoE[i][j]   = -10.

fpmodel = open("_model_%s" % dt, "r")
for line in fpmodel.readlines():
    line  = line.strip()
    split = line.split()
    model_matrix[int(split[0])][int(split[1])] = float(split[2])
    model_matrix[int(split[1])][int(split[0])] = float(split[2])

fpexp = open("_exp_%s" % dt, "r")
for line in fpexp.readlines():
    line  = line.strip()
    split = line.split()
    exp_matrix[int(split[0])][int(split[1])] = float(split[2])
    exp_matrix[int(split[1])][int(split[0])] = float(split[2])

correlate_with_real_data(model_matrix, exp_matrix, model_matrix, exp_matrix, LSthreshold=int(genDist*1000000), nloci=size, cutoff=dcutoff, corr='spearman', off_diag=1, timestep=timestep, resolution=resolution)
