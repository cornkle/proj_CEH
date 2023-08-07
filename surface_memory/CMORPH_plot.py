# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
from scipy import ndimage
from utils import u_arrays, constants as cnst, u_met
from scipy.stats import ttest_ind as ttest
from scipy.interpolate import griddata
import pickle as pkl
from utils import u_arrays as ua
import os
import collections
import warnings
import ipdb
import salem

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1


def plot_amsr(hour):

    path = cnst.network_data + 'figs/CORE_MEMORY/'

    month = 9
    key = 'trackedDRY'

    # UTC_15000_ALL_-60_5slotSmall
    dic1 = pkl.load(open(path + "tables/dominantCores_" + str(
        hour) + "UTC_15000_2dAMSL_ext_day0_ALLS_minusMean_15-24LT_"+str(month).zfill(2)+"_"+key+".p",
                         "rb"))  # UTC_15000_ALL_-60_5slotSmall
    dic2 = pkl.load(open(path + "tables/dominantCores_" + str(
        hour) + "UTC_15000_2dAMSL_ext_day+1_ALLS_minusMean_15-24LT_"+str(month).zfill(2)+"_"+key+".p", "rb"))
    dic3 = pkl.load(open(path + "tables/dominantCores_" + str(
        hour) + "UTC_15000_2dAMSL_ext_day+2_ALLS_minusMean_15-24LT_"+str(month).zfill(2)+"_"+key+".p",
                         "rb"))  # UTC_15000_ALL_-60_5slotSmall
    dic4 = pkl.load(open(path + "tables/dominantCores_" + str(
        hour) + "UTC_15000_2dAMSL_ext_day+3_ALLS_minusMean_15-24LT_"+str(month).zfill(2)+"_"+key+".p", "rb"))
    dic5 = pkl.load(open(path + "tables/dominantCores_" + str(
        hour) + "UTC_15000_2dAMSL_ext_day+4_ALLS_minusMean_15-24LT_"+str(month).zfill(2)+"_"+key+".p",
                         "rb"))  # UTC_15000_ALL_-60_5slotSmall
    dic6 = pkl.load(open(path + "tables/dominantCores_" + str(
        hour) + "UTC_15000_2dAMSL_ext_day+5_ALLS_minusMean_15-24LT_"+str(month).zfill(2)+"_"+key+".p", "rb"))

    import matplotlib

    cmap = matplotlib.cm.get_cmap('viridis')
    rgba = cmap(0.5)

    names = ['DRY - day0', 'DRY - day+1', 'DRY - day+2', 'DRY - day+3', 'DRY - day+4', 'DRY - day+5']

    pick = [dic4, '', dic2]

    f = plt.figure(figsize=(8.5, 7), dpi=300)
    for ids, dic in enumerate([dic1, dic2, dic3, dic4, dic5, dic6]):

        #lsta = np.sum((dic['lsta'])[0],axis=0) / np.sum(dic['lsta'][1], axis=0)
        amsr = (dic['amsr'])[0] / dic['amsr'][1]
        cmorph = (dic['cmorph'])[0] / dic['cmorph'][1]
        msg = (dic['msg'])[0] / dic['msg'][1]
        cores = dic['cores']

        cmorph = ndimage.gaussian_filter(cmorph, 3, mode='nearest')
        amsr = ndimage.gaussian_filter(amsr, 3, mode='nearest')
        #lsta = ndimage.gaussian_filter(lsta, 3, mode='nearest')
        msg = ndimage.gaussian_filter(msg, 3, mode='nearest')

        # if ids >=4:
        #     ipdb.set_trace()



        print('NUMBER OF CORES', cores)

        dist = 200
        llevels = np.array(list(np.arange(-0.8, 0, 0.1)) + list(np.arange(0.1, 0.81, 0.1)))  # *12000
        alevels = np.array(list(np.arange(-2.5, 0, 0.25)) + list(np.arange(0.25, 2.51, 0.25)))  # *12000

        alevels = [-4, -3, -2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2, 3, 4]

        ax = f.add_subplot(3, 2, ids + 1)

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, amsr,
                     cmap='RdBu', extend='both', levels=alevels)
        plt.colorbar(label='%')
        #
        # plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, lsta,
        #              cmap='RdBu', extend='both', levels=llevels)
        # plt.colorbar(label='%')

        if ids == 3:
            lev = [2, 4, 8, 10, 15, 25, 50, 75]  # np.arange(5, 71, 10)
            # lev = np.arange(10, 71, 5)
            colors = [cmap(0.05), cmap(0.5)]
        else:
            lev = [2, 4, 8, 10, 15, 25, 50, 75]  # np.arange(5, 71, 10)
            # lev = np.arange(10, 71, 5)
            colors = [cmap(0.05), cmap(0.6), cmap(0.99)]
        #
        # cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3,
        #                  cmorph * 100 ,
        #                  linewidths=1.2, linestyles=['solid'], levels=lev, colors='k')
        #
        # plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")

        cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3,
                         msg*100 ,
                         linewidths=1.2, linestyles=['solid'], levels=lev, colors='k')

        plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")
        if ids in [0, 2]:
            ax.set_ylabel('km')
        ax.set_xlabel('km')
        plt.axvline(x=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.axhline(y=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.plot(0, 0, marker='o', color='dimgrey')

        plt.title(names[ids], fontweight='bold', fontname='Ubuntu', fontsize=10)  # + ' | ' + str(cores) + ' cores'

    plt.tight_layout()
    # text = ['a', 'b', 'c', 'd']
    # plt.annotate(text[0], xy=(0.04, 0.96), xytext=(0, 4),xycoords=('figure fraction', 'figure fraction'),
    #              textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    # plt.annotate(text[1], xy=(0.54, 0.96), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
    #              textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    # plt.annotate(text[2], xy=(0.04, 0.49), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
    #              textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    # plt.annotate(text[3], xy=(0.54, 0.49), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
    #              textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)

    plt.savefig(
        path + 'plots/cmorph_day+5_dominantCore_ext_15-24LT_'+key+'_'+str(month).zfill(2)+'_' + str(hour).zfill(
            2) + '.png')
    plt.close('all')

