# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
from utils import u_arrays, constants as cnst, u_met
import pickle as pkl

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1
from scipy import ndimage
import matplotlib.patches as patches


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

    f = plt.figure(figsize=(17, 3), dpi=200)
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






def plot_amsr_lsta_trio(hour):
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    key = '2hOverlap'

    f = plt.figure(figsize=(10.5, 3), dpi=300)

    labels = ['Day-1', 'Day0', 'Day+1']

    left = 0.01
    bottom = 0.1
    width = 0.3
    height=0.8

    spot = [[]]

    for ids, ll in enumerate(['day-1', 'day0', 'day+1']):


        dic = pkl.load(
            open(path + "/coeffs_nans_stdkernel_USE_" + str(hour) + "UTC_15000_2dAMSL_"+ll+"_ALLS_minusMean_INIT_" + (key) + ".p",
                 "rb"))

        lsta = (dic['lsta'])[0] / dic['lsta'][1]
        amsr = (dic['amsr'])[0] / dic['amsr'][1]

        amsr = ndimage.gaussian_filter(amsr, 6, mode='nearest')

        cores = dic['cores']

        lcnt = dic['lsta'][1]
        acnt = dic['amsr'][1]

        dist = 200
        llevels = np.array(list(np.arange(-0.8, 0, 0.1)) + list(np.arange(0.1, 0.81, 0.1)))#*12000
        alevels = np.array(list(np.arange(-2.5, 0, 0.5)) + list(np.arange(0.5, 2.51, 0.5)))#*12000
        alevels = [-2.5,-2,-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5,2,2.5]

        # # llevels = np.array(list(np.arange(-1.6, 0, 0.2)) + list(np.arange(0.2, 1.61, 0.2)))#*12000  WET
        # # alevels = np.array(list(np.arange(-3, 0, 0.25)) + list(np.arange(0.25, 3.25, 0.25)))#*12000
        #
        # llevels = np.array(list(np.arange(-1.5, 0, 0.2)) + list(np.arange(0.2, 1.51, 0.2)))  # *12000 DRY
        # alevels = np.array(list(np.arange(-5, 0, 0.5)) + list(np.arange(0.5, 5.5, 0.5)))  # *12000

        ax = f.add_subplot(1,3,ids+1)

        mp1 = plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, lsta,
                     cmap='RdBu_r', extend='both', levels=llevels)
        #plt.colorbar(label='K')
        cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, amsr,
                         colors='k', linewidths=1, linestyles=['solid'], levels=alevels) #cmap='RdBu'
        plt.clabel(cs, inline=1, fontsize=8, fmt="%1.1f")

        ax.set_xlabel('km')
        if ids == 0:
            ax.set_ylabel('km')

        # if ids > 0:
        #     ax.set_yticklabels('')

        plt.axvline(x=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.axhline(y=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.plot(0, 0, marker='o', color='dimgrey')

        plt.title(labels[ids],fontsize=10)


    plt.tight_layout()
    text = ['a', 'b', 'c']
    plt.annotate(text[0], xy=(0.06, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[1], xy=(0.36, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[2], xy=(0.65, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)

    f.subplots_adjust(right=0.91)
    cax = f.add_axes([0.92, 0.18, 0.015, 0.73])
    cbar = f.colorbar(mp1, cax)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('K', fontsize=10)

    plt.savefig(path + '2hOverlap/amsreVSlsta/MAPS_AMSL_TRIO_ALLS_minusMean_noCore_INIT' + str(hour).zfill(2) + '.pdf')
    plt.close('all')


def plot_amsr_lsta_six(hour):
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    key = 'NEWTRACKING'

    f = plt.figure(figsize=(9, 10.5), dpi=300)

    labels = ['Night-1', 'Day-1', 'Night0', 'Day0', 'Night+1', 'Day+1']

    left = 0.01
    bottom = 0.1
    width = 0.3
    height=0.8

    spot = [[]]

    for ids, ll in enumerate(['night-1', 'day-1', 'night0', 'day0', 'night+1', 'day+1']):


        # dic = pkl.load(
        #     open(path + "/coeffs_nans_stdkernel_USE_" + str(hour) + "UTC_15000_2dAMSL_"+ll+"_ALLS_minusMean_INIT_" + (key) + ".p",
        #          "rb"))


        dic = pkl.load(
            open(path + "/coeffs_nans_stdkernel_USE_" + str(hour) + "UTC_15000_2dAMSL_"+ll+"_ALLS_minusMean_INIT_" + (key) + ".p",
                 "rb"))

        lsta = (dic['lsta'])[0] / dic['lsta'][1]
        amsr = (dic['amsr'])[0] / dic['amsr'][1]

        amsr = ndimage.gaussian_filter(amsr, 6, mode='nearest')

        cores = dic['cores']

        lcnt = dic['lsta'][1]
        acnt = dic['amsr'][1]

        dist = 200
        llevels = np.array(list(np.arange(-0.8, 0, 0.1)) + list(np.arange(0.1, 0.81, 0.1)))#*12000
        alevels = np.array(list(np.arange(-2.5, 0, 0.5)) + list(np.arange(0.5, 2.51, 0.5)))#*12000
        alevels = [-2.5,-2,-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5,2,2.5]

        # # llevels = np.array(list(np.arange(-1.6, 0, 0.2)) + list(np.arange(0.2, 1.61, 0.2)))#*12000  WET
        # # alevels = np.array(list(np.arange(-3, 0, 0.25)) + list(np.arange(0.25, 3.25, 0.25)))#*12000
        #
        # llevels = np.array(list(np.arange(-1.5, 0, 0.2)) + list(np.arange(0.2, 1.51, 0.2)))  # *12000 DRY
        # alevels = np.array(list(np.arange(-5, 0, 0.5)) + list(np.arange(0.5, 5.5, 0.5)))  # *12000

        ax = f.add_subplot(3,2,ids+1)

        mp1 = plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, lsta,
                     cmap='RdBu_r', extend='both', levels=llevels)
        #plt.colorbar(label='K')
        cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, amsr,
                         colors='k', linewidths=1, linestyles=['solid'], levels=alevels) #cmap='RdBu'
        plt.clabel(cs, inline=1, fontsize=8, fmt="%1.1f")

        ax.set_xlabel('km')
        if ids == 0:
            ax.set_ylabel('km')

        # if ids > 0:
        #     ax.set_yticklabels('')

        plt.axvline(x=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.axhline(y=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.plot(0, 0, marker='o', color='dimgrey')

        plt.title(labels[ids],fontsize=10)


    plt.tight_layout()
    text = ['a', 'b', 'c']
    plt.annotate(text[0], xy=(0.06, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[1], xy=(0.36, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[2], xy=(0.65, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)

    f.subplots_adjust(right=0.91)
    cax = f.add_axes([0.92, 0.18, 0.015, 0.73])
    cbar = f.colorbar(mp1, cax)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('K', fontsize=10)

    plt.savefig(path + '2hOverlap/amsreVSlsta/MAPS_AMSL_TRIO_ALLS_minusMean_noCore_INIT_SIX_' + str(hour).zfill(2) + '.png')
    plt.close('all')


def plot_amsr_dry_wet(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    key = '2hOverlap'
    daykey = 'day+1'
    # dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_WET_old_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_WET_old_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic3 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_old_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_old_2hOverlap.p", "rb"))


    # dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_WET_INIT_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_WET_INIT_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic3 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_INIT_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_INIT_2hOverlap.p", "rb"))


    dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_WET_INIT_Q20_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_WET_INIT_Q20_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic3 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_INIT_Q20_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_INIT_Q20_NEWTRACKING.p", "rb"))

    import matplotlib

    cmap = matplotlib.cm.get_cmap('viridis')
    rgba = cmap(0.5)

    names = ['DRY - day0', 'DRY - day+1', 'WET - day0', 'WET - day+1']

    pick = [dic4, '',dic2]

    f = plt.figure(figsize=(8.5, 6), dpi=300)
    for ids, dic in enumerate([dic3,dic4,dic1,dic2]):

        #lsta = (dic['lsta'])[0] / dic['lsta'][1]
        amsr = (dic['amsr'])[0] / dic['amsr'][1]
        cmorph = (dic['cmorph'])[0] / dic['cmorph'][1]

        cmorph = ndimage.gaussian_filter(cmorph, 6, mode='nearest')
        amsr = ndimage.gaussian_filter(amsr, 4, mode='nearest')

        msg = (dic['msg'])[0] / dic['msg'][1]
        cores = dic['cores']

        print('NUMBER OF CORES', cores)


        dist=200
        llevels = np.array(list(np.arange(-0.8, 0, 0.1)) + list(np.arange(0.1, 0.81, 0.1)))#*12000
        alevels = np.array(list(np.arange(-2.5, 0, 0.25)) + list(np.arange(0.25, 2.51, 0.25)))#*12000

        alevels = [-4,-3,-2,-1,-0.5, -0.25, 0.25,0.5,1,2,3,4]

        ax = f.add_subplot(2,2,ids+1)

        plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , amsr , cmap='RdBu', extend='both', levels=alevels)
        plt.colorbar(label='%')

        if ids == 3:
            lev = np.arange(10, 71, 20)
            #lev = np.arange(10, 71, 5)
            colors = [cmap(0.05), cmap(0.5)]
        else:
            lev = np.arange(10, 71, 20)
            #lev = np.arange(10, 71, 5)
            colors = [cmap(0.05), cmap(0.6), cmap(0.99)]

        cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100*1.2,
                    linewidths=1.2, linestyles=['solid'], levels=lev, colors='k')
        plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")

        lev = [-99, 50]
        cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100*1.2,
                    linewidths=1.5, linestyles=['solid'], levels=lev, colors='k')
        plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")


        # if ids in [0,2]:
        #     cmorph2 = ((pick[ids])['cmorph'])[0] / ((pick[ids])['cmorph'])[1]
        #     cmorph2 = ndimage.gaussian_filter(cmorph2, 6, mode='nearest')
        #     lev = [-99,50]
        #
        #     cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3,
        #                      cmorph2 * 100 * 1.2,
        #                      linewidths=0.9, linestyles=['dotted'], levels=lev, colors='b')
        #
        #     plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")


        # cs = plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100, cmap='viridis',
        #            levels=np.arange(0,101,10), extend='both')

        # cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, msg*100, cmap='jet',
        #             linewidths=1, linestyles=['solid'], levels=np.arange(10,91,10))


        #plt.colorbar()

        if ids in [0,2]:
            ax.set_ylabel('km')
        ax.set_xlabel('km')
        plt.axvline(x=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.axhline(y=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.plot(0, 0, marker='o', color='dimgrey')

        plt.title(names[ids] , fontweight='bold', fontname='Ubuntu', fontsize=10) #+ ' | ' + str(cores) + ' cores'



    plt.tight_layout()
    text = ['a', 'b', 'c', 'd']
    plt.annotate(text[0], xy=(0.04, 0.96), xytext=(0, 4),xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[1], xy=(0.54, 0.96), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[2], xy=(0.04, 0.49), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[3], xy=(0.54, 0.49), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)

    plt.savefig(path + '2hOverlap/amsreVSlsta/wcoeff_maps_all_AMSL_DRYWET_CMORPH_NEWTRACKING_Q20_' + str(hour).zfill(2) + '.pdf')
    plt.close('all')


def plot_amsr_day3_six(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    key = '2hOverlap'
    daykey = 'day+1'
    # dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_WET_old_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_WET_old_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic3 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_old_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_old_2hOverlap.p", "rb"))


    # dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_WET_INIT_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_WET_INIT_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic3 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_INIT_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_INIT_2hOverlap.p", "rb"))


 #UTC_15000_ALL_-60_5slotSmall
    dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_INIT_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_INIT_NEWTRACKING.p", "rb"))
    dic3 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+2_ALLS_minusMean_CMORPH_DRY_INIT_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+3_ALLS_minusMean_CMORPH_DRY_INIT_NEWTRACKING.p", "rb"))
    dic5 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+4_ALLS_minusMean_CMORPH_DRY_INIT_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic6 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+5_ALLS_minusMean_CMORPH_DRY_INIT_NEWTRACKING.p", "rb"))

    import matplotlib

    cmap = matplotlib.cm.get_cmap('viridis')
    rgba = cmap(0.5)

    names = ['DRY - day0', 'DRY - day+1', 'DRY - day+2', 'DRY - day+3', 'DRY - day+4', 'DRY - day+5']

    pick = [dic4, '',dic2]

    f = plt.figure(figsize=(8.5, 7), dpi=300)
    for ids, dic in enumerate([dic1,dic2,dic3,dic4, dic5, dic6]):

        #lsta = (dic['lsta'])[0] / dic['lsta'][1]
        amsr = (dic['amsr'])[0] / dic['amsr'][1]
        cmorph = (dic['cmorph'])[0] / dic['cmorph'][1]

        cmorph = ndimage.gaussian_filter(cmorph, 6, mode='nearest')
        amsr = ndimage.gaussian_filter(amsr, 3, mode='nearest')

        # if ids >=4:
        #     ipdb.set_trace()

        msg = (dic['msg'])[0] / dic['msg'][1]
        cores = dic['cores']

        print('NUMBER OF CORES', cores)


        dist=200
        llevels = np.array(list(np.arange(-0.8, 0, 0.1)) + list(np.arange(0.1, 0.81, 0.1)))#*12000
        alevels = np.array(list(np.arange(-2.5, 0, 0.25)) + list(np.arange(0.25, 2.51, 0.25)))#*12000

        alevels = [-4,-3,-2,-1,-0.5, -0.25, 0.25,0.5,1,2,3,4]

        ax = f.add_subplot(3,2,ids+1)

        plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , amsr , cmap='RdBu', extend='both', levels=alevels)
        plt.colorbar(label='%')

        if ids == 3:
            lev = [5,10,15,25,50,75]#np.arange(5, 71, 10)
            #lev = np.arange(10, 71, 5)
            colors = [cmap(0.05), cmap(0.5)]
        else:
            lev = [5,10,15,25,50,75]#np.arange(5, 71, 10)
            #lev = np.arange(10, 71, 5)
            colors = [cmap(0.05), cmap(0.6), cmap(0.99)]

        cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100*1.2,
                    linewidths=1.2, linestyles=['solid'], levels=lev, colors='k')
        plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")


        if ids in [0,2]:
            cmorph2 = ((pick[ids])['cmorph'])[0] / ((pick[ids])['cmorph'])[1]
            cmorph2 = ndimage.gaussian_filter(cmorph2, 6, mode='nearest')
            lev = [-99,50]

            # cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3,
            #                  cmorph2 * 100 * 1.2,
            #                  linewidths=0.9, linestyles=['dotted'], levels=lev, colors='b')
            #
            # plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")


        # cs = plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100, cmap='viridis',
        #            levels=np.arange(0,101,10), extend='both')

        # cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, msg*100, cmap='jet',
        #             linewidths=1, linestyles=['solid'], levels=np.arange(10,91,10))


        #plt.colorbar()

        if ids in [0,2]:
            ax.set_ylabel('km')
        ax.set_xlabel('km')
        plt.axvline(x=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.axhline(y=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.plot(0, 0, marker='o', color='dimgrey')

        plt.title(names[ids] , fontweight='bold', fontname='Ubuntu', fontsize=10) #+ ' | ' + str(cores) + ' cores'



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

    plt.savefig(path + '2hOverlap/amsreVSlsta/wcoeff_maps_all_AMSL_DRYWET_CMORPH_NEWTRACKING_DAY3_SIX_' + str(hour).zfill(2) + '.png')
    plt.close('all')



def plot_amsr_dry_wet_trio(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'

    # dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_WET_old_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_WET_old_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic3 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_old_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_old_2hOverlap.p", "rb"))


    # dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_WET_INIT_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_WET_INIT_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic3 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_INIT_2hOverlap.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    # dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_INIT_2hOverlap.p", "rb"))

    dic0 = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_" + str(
        hour) + "UTC_15000_2dAMSL_night0_ALLS_minusMean_CMORPH_WET_INIT_noQ20_NEWTRACKING.p", "rb"))
    dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_WET_INIT_noQ20_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_WET_INIT_noQ20_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic3 = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_" + str(
        hour) + "UTC_15000_2dAMSL_night0_ALLS_minusMean_CMORPH_DRY_INIT_noQ20_NEWTRACKING.p", "rb"))
    dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_INIT_noQ20_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic5 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_INIT_noQ20_NEWTRACKING.p", "rb"))

    names = ['DRY - day0, 0130UTC','DRY - day0, 1330UTC', 'DRY - day+1, 1330UTC', 'WET - day0, 0130UTC', 'WET - day0, 1330UTC', 'WET - day+1, 1330UTC']

    pick = ['',dic5, '','',dic2]

    f = plt.figure(figsize=(12, 6), dpi=300)
    for ids, dic in enumerate([dic3,dic4,dic5,dic0,dic1,dic2]):

        #lsta = (dic['lsta'])[0] / dic['lsta'][1]
        amsr = (dic['amsr'])[0] / dic['amsr'][1]
        cmorph = (dic['cmorph'])[0] / dic['cmorph'][1]

        cmorph = ndimage.gaussian_filter(cmorph, 8, mode='nearest')
        amsr = ndimage.gaussian_filter(amsr, 8, mode='nearest')

        msg = (dic['msg'])[0] / dic['msg'][1]
        cores = dic['cores']

        print('NUMBER OF CORES', cores)


        dist=200

        # alevels = [-5,-4,-3,-2,-1.5,-1,-0.5,-0.25,0.25,0.5,1,1.5,2,3,4]

        cmap = matplotlib.cm.get_cmap('RdBu')

        vals_low = list(np.linspace(-0.01, 0.5,9))#+list([0.5])
        vals_high = (((np.linspace(0, 0.5,9)-1)[::-1])*-1)[1:-1]#(np.arange(0,0.5,0.05)+0.5+0.05)[0:-2]

        # allcol = list(vals_low)+list(vals_high)

        #ipdb.set_trace()
        #
        #
        # def col_colors(l):
        #     col = []
        #     for ll in l :
        #         col.append(cmap(ll))
        #     return col
        #
        # colours = col_colors(allcol)

        #ipdb.set_trace()
        alevels = [-4, -3, -2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2, 3, 4]
        #alevels = [-3.9,-3.3, -2.7, -2.1, -1.5, -0.9, -0.3, 0.3, 0.9, 1.5, 2.1, 2.7, 3.3,3.9]
        ax = f.add_subplot(2,3,ids+1)

        amsr[np.isnan(amsr)]= 0
        #, norm=matplotlib.colors.Normalize(vmin=-6, vmax=4)
        plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , amsr , cmap='RdBu', extend='both', levels=alevels)
        plt.colorbar(label='%')

        if ids == 3:
            lev = np.arange(10, 71, 20)
            #lev = np.arange(10, 71, 5)

        else:
            lev = np.arange(10, 71, 20)
            #lev = np.arange(10, 71, 5)

        if ids in [1,2,4,5]:

            cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100*1.2,
                        linewidths=1.2, linestyles=['solid'], levels=lev, colors='k')
            plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")

            lev = [-99, 50]
            cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100*1.2,
                        linewidths=1.5, linestyles=['solid'], levels=lev, colors='k')
            plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")


        # if ids in [1,4]:
        #     cmorph2 = ((pick[ids])['cmorph'])[0] / ((pick[ids])['cmorph'])[1]
        #     cmorph2 = ndimage.gaussian_filter(cmorph2, 6, mode='nearest')
        #     lev = [-99,50]
        #
        #     cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3,
        #                      cmorph2 * 100 * 1.2,
        #                      linewidths=0.9, linestyles=['dotted'], levels=lev, colors='b')
        #
        #     plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")


        # cs = plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100, cmap='viridis',
        #            levels=np.arange(0,101,10), extend='both')

        # cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, msg*100, cmap='jet',
        #             linewidths=1, linestyles=['solid'], levels=np.arange(10,91,10))


        #plt.colorbar()

        if ids in [0,3]:
            ax.set_ylabel('km')
            rect = patches.Rectangle((-300, -300), 600, 600, linewidth=1.2, edgecolor='dimgrey', facecolor='none')
            # Add the patch to the Axes
            ax.add_patch(rect)


        ax.set_xlabel('km')
        plt.axvline(x=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.axhline(y=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.plot(0, 0, marker='o', color='dimgrey')

        plt.title(names[ids] , fontweight='bold', fontname='Ubuntu', fontsize=10) #+ ' | ' + str(cores) + ' cores'



    plt.tight_layout()
    text = ['a', 'b', 'c', 'd','e','f']
    plt.annotate(text[0], xy=(0.04, 0.96), xytext=(0, 4),xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[1], xy=(0.34, 0.96), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[2], xy=(0.68, 0.96), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[3], xy=(0.04, 0.49), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[4], xy=(0.34, 0.49), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[5], xy=(0.68, 0.49), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)

    plt.savefig(path + '2hOverlap/amsreVSlsta/wcoeff_maps_all_AMSL_DRYWET_CMORPH_NEWTRACKING_noQ20_TRIO_' + str(hour).zfill(2) + '.png')
    plt.close('all')



def plot_amsr_dry_wet_paper(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'

    dic0 = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_" + str(
        hour) + "UTC_15000_2dAMSL_night0_ALLS_minusMean_CMORPH_WET_INIT_noQ20_NEWTRACKING.p", "rb"))
    dic1 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_WET_INIT_noQ20_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic2 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_WET_INIT_noQ20_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic3 = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_" + str(
        hour) + "UTC_15000_2dAMSL_night0_ALLS_minusMean_CMORPH_DRY_INIT_noQ20_NEWTRACKING.p", "rb"))
    dic4 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_INIT_noQ20_NEWTRACKING.p", "rb")) #UTC_15000_ALL_-60_5slotSmall
    dic5 = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_INIT_noQ20_NEWTRACKING.p", "rb"))

    names = ['DRY - day0', 'DRY - day+1', 'WET - day0', 'WET - day+1']

    f = plt.figure(figsize=(8.5, 6), dpi=300)
    for ids, dic in enumerate([dic4,dic5,dic1,dic2]):

        #lsta = (dic['lsta'])[0] / dic['lsta'][1]
        amsr = (dic['amsr'])[0] / dic['amsr'][1]
        cmorph = (dic['cmorph'])[0] / dic['cmorph'][1]

        cmorph = ndimage.gaussian_filter(cmorph, 8, mode='nearest')
        amsr = ndimage.gaussian_filter(amsr, 8, mode='nearest')

        msg = (dic['msg'])[0] / dic['msg'][1]
        cores = dic['cores']

        print('NUMBER OF CORES', cores)


        dist=200

        alevels = [-4, -3, -2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2, 3, 4]
        #alevels = [-3.9,-3.3, -2.7, -2.1, -1.5, -0.9, -0.3, 0.3, 0.9, 1.5, 2.1, 2.7, 3.3,3.9]
        ax = f.add_subplot(2,2,ids+1)

        amsr[np.isnan(amsr)]= 0
        #, norm=matplotlib.colors.Normalize(vmin=-6, vmax=4)
        plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , amsr , cmap='RdBu', extend='both', levels=alevels)
        plt.colorbar(label='%')

        if ids == 3:
            lev = np.arange(10, 71, 20)
            #lev = np.arange(10, 71, 5)

        else:
            lev = np.arange(10, 71, 20)
            #lev = np.arange(10, 71, 5)

        cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100*1.2,
                    linewidths=1.2, linestyles=['solid'], levels=lev, colors='k')
        plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")

        lev = [-99, 50]
        cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, cmorph*100*1.2,
                    linewidths=1.5, linestyles=['solid'], levels=lev, colors='k')
        plt.clabel(cs, inline=1, fontsize=9, fmt="%1.0f")


        ax.set_xlabel('km')
        plt.axvline(x=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.axhline(y=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.plot(0, 0, marker='o', color='dimgrey')

        plt.title(names[ids] , fontweight='bold', fontname='Ubuntu', fontsize=10) #+ ' | ' + str(cores) + ' cores'


    plt.tight_layout()
    text = ['a', 'b', 'c', 'd']
    plt.annotate(text[0], xy=(0.04, 0.96), xytext=(0, 4),xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[1], xy=(0.54, 0.96), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[2], xy=(0.04, 0.49), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[3], xy=(0.54, 0.49), xytext=(0, 4),  xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)

    plt.savefig(path + '2hOverlap/amsreVSlsta/wcoeff_maps_all_AMSL_DRYWET_CMORPH_NEWTRACKING_noQ20_PAPER_' + str(hour).zfill(2) + '.png')
    plt.close('all')