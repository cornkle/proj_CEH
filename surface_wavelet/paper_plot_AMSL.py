# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import multiprocessing
import ipdb
import pandas as pd
from wavelet import util as wutil
from utils import u_arrays, constants as cnst, u_met
from scipy.stats import ttest_ind as ttest
from scipy.interpolate import griddata
import pickle as pkl
from matplotlib.gridspec import GridSpec
import collections

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1


def run_all():
    for h in [15,16,17,18,19,20,21, 22,23,0,1,2,3,4,5,6,7]:  #,22,23,0,1,2,3,4,5,6,7
        plot_map_AMSRE(h)


def plot_map_AMSRE(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    key = '2hOverlap'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_mini_day+1_DRY_" + (key) + ".p", "rb")) #UTC_15000_ALL_-60_5slotSmall

    lsta = (dic['lsta'])[0] / dic['lsta'][1]
    amsr = (dic['amsr'])[0] / dic['amsr'][1]

    cores = dic['cores']

    lcnt = dic['lsta'][1]
    acnt = dic['amsr'][1]

    dist=100
    llevels = np.array(list(np.arange(-0.8, 0, 0.1)) + list(np.arange(0.1, 0.81, 0.1)))#*12000
    alevels = np.array(list(np.arange(-2.5, 0, 0.25)) + list(np.arange(0.25, 2.51, 0.25)))#*12000

    #
    # llevels = np.array(list(np.arange(-1.6, 0, 0.2)) + list(np.arange(0.2, 1.61, 0.2)))#*12000  WET
    # alevels = np.array(list(np.arange(-3, 0, 0.25)) + list(np.arange(0.25, 3.25, 0.25)))#*12000

    # llevels = np.array(list(np.arange(-1, 0, 0.2)) + list(np.arange(0.2, 1.01, 0.2)))#*12000 DRY
    # alevels = np.array(list(np.arange(-4, 0, 0.5)) + list(np.arange(0.5, 4.5, 0.5)))#*12000

    f = plt.figure(figsize=(12, 9), dpi=200)
    ax = f.add_subplot(221)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , lsta , cmap='RdBu_r', extend='both', levels=llevels)
    plt.colorbar(label='K')
    cs = plt.contour((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , amsr , colors='k', linewidths=1, linestyles=['dotted'])
    plt.clabel(cs, inline=1, fontsize=8, fmt="%1.1f")
    ax.plot(0,0, 'bo')
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.axvline(x=0, linestyle='dashed', color='k',linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)

    plt.title('LSTA | Nb cores: ' + str(cores) + '| ' + str(hour).zfill(2) + '00UTC',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , amsr , cmap='RdBu', extend='both', levels=alevels)
    plt.colorbar(label='%')
    # plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, lsta, colors='k',
    #             linewidths=0.8, linestyles=['dashed'])


    ax.plot(0,0, 'bo')
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.axvline(x=0, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
    plt.title('AMSRE| Nb cores: ' + str(cores) + '| ' + str(hour).zfill(2) + '00UTC', fontsize=10)

    ax = f.add_subplot(223)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , lcnt, cmap='RdBu_r', extend='both')
    plt.plot(0,0,'bo')
    plt.colorbar(label='')
    #plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[0,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
    plt.axvline(x=0, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('LSTA: Number valid pixels')


    ax = f.add_subplot(224)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , acnt, cmap='RdBu_r', extend='both')
    plt.plot(0,0,'bo')
    plt.colorbar(label='')
    #plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[0,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
    plt.axvline(x=0, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.title('AMSRE: Number valid pixels')


    plt.tight_layout()
    #plt.savefig(path + '/amsreVSlsta/wcoeff_maps_all_AMSL_SMFINITE_'+str(hour).zfill(2)+'.png')
    #plt.show()
    #plt.close('all')


def plot_amsr_lsta_only(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    key = '2hOverlap'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_mini_day-1_DRY_" + (key) + ".p", "rb")) #UTC_15000_ALL_-60_5slotSmall

    lsta = (dic['lsta'])[0] / dic['lsta'][1]
    amsr = (dic['amsr'])[0] / dic['amsr'][1]

    cores = dic['cores']

    lcnt = dic['lsta'][1]
    acnt = dic['amsr'][1]

    dist=200
    # llevels = np.array(list(np.arange(-0.8, 0, 0.1)) + list(np.arange(0.1, 0.81, 0.1)))#*12000
    # alevels = np.array(list(np.arange(-2.5, 0, 0.25)) + list(np.arange(0.25, 2.51, 0.25)))#*12000

    # llevels = np.array(list(np.arange(-1.6, 0, 0.2)) + list(np.arange(0.2, 1.61, 0.2)))#*12000  WET
    # alevels = np.array(list(np.arange(-3, 0, 0.25)) + list(np.arange(0.25, 3.25, 0.25)))#*12000

    llevels = np.array(list(np.arange(-1.5, 0, 0.2)) + list(np.arange(0.2, 1.51, 0.2)))#*12000 DRY
    alevels = np.array(list(np.arange(-5, 0, 0.5)) + list(np.arange(0.5, 5.5, 0.5)))#*12000

    f = plt.figure(figsize=(10, 4), dpi=200)
    ax = f.add_subplot(121)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , lsta , cmap='RdBu_r', extend='both', levels=llevels)
    plt.colorbar(label='K')
    cs = plt.contour((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , amsr , colors='k', linewidths=1, linestyles=['dotted'])
    plt.clabel(cs, inline=1, fontsize=8, fmt="%1.1f")
    ax.plot(0,0, 'bo')
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.axvline(x=0, linestyle='dashed', color='k',linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)

    plt.title('LSTA | Nb cores: ' + str(cores) + '| ' + str(hour).zfill(2) + '00UTC',
              fontsize=10)

    ax = f.add_subplot(122)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , amsr , cmap='RdBu', extend='both', levels=alevels)
    plt.colorbar(label='%')
    # plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, lsta, colors='k',
    #             linewidths=0.8, linestyles=['dashed'])


    ax.plot(0,0, 'bo')
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.axvline(x=0, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
    plt.title('AMSRE| Nb cores: ' + str(cores) + '| ' + str(hour).zfill(2) + '00UTC', fontsize=10)
    plt.tight_layout()
    #plt.show()
    plt.savefig(path + '2hOverlap/amsreVSlsta/wcoeff_maps_all_AMSL_dry-1_' + str(hour).zfill(2) + '.png')