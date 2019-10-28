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
import collections

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1

def plot_all_cross():

    hours = [14,15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]
    for h in hours:
        plot(h)



def plot(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000.p", "rb"))

    scales = dic['scales']
    nbcores = dic['nbcores']
    nbrcores = dic['nbrcores']
    del dic['scales']
    del dic['nbcores']
    del dic['nbrcores']
    del dic['kernel']


    keys = list(dic.keys())
    cnt = (dic['SN-pos'][0]).shape[0]

    #ipdb.set_trace()

    # for l in keys:
    #     if keys == 'scales':
    #         continue
    #
    #     (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
    #     (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)

    l=0
    dist=100

    snblob = (dic[keys[l]])[0] #/ (dic[keys[l]])[3]).T#-(dic[keys[l+2]])[0]
    snrandom = (dic[keys[l]])[1] #/ (dic[keys[l]])[4]).T#-(dic[keys[l+2]])[1]
    snmask = (dic[keys[l]])[2]#-(dic[keys[l+2]])[2]
    # snblob_std = (dic[keys[l]])[3]
    # snrandom_std = (dic[keys[l]])[4]

    weblob = (dic[keys[l+1]])[0] #/ (dic[keys[l+1]])[3]).T#-(dic[keys[l+3]])[0]
    werandom = (dic[keys[l+1]])[1] #/ (dic[keys[l+1]])[4]).T#-(dic[keys[l+3]])[1]
    wemask = (dic[keys[l+1]])[2]#-(dic[keys[l+3]])[2]
    # weblob_std = (dic[keys[l+1]])[3]
    # werandom_std = (dic[keys[l+1]])[4]

    l=2
    dist=100

    f = plt.figure(figsize=(9, 9))
    ax = f.add_subplot(221)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales , (snblob) - (snrandom) , cmap='RdBu_r', levels = [-0.1, -0.08,-0.06, -0.04, -0.02,0.02,0.04,0.06,0.08,0.1])
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3,scales,   (weblob) - (werandom) , cmap='RdBu_r', levels = [-0.1, -0.08,-0.06, -0.04, -0.02,0.02,0.04,0.06,0.08,0.1], extend='both')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, wemask, colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)


    ax = f.add_subplot(223)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales,  snblob,  cmap='RdBu_r', levels = [-0.1, -0.08,-0.06, -0.04, -0.02,0.02,0.04,0.06,0.08,0.1]) #, vmin = -0.1, vmax=0.1)
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels=[0.5, 1],
                 linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbrcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(224)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, weblob    , cmap='RdBu_r', levels = [-0.1, -0.08,-0.06, -0.04, -0.02,0.02,0.04,0.06,0.08,0.1])
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales,wemask, colors='none', hatches='.', levels=[0.5, 1],
                 linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)

    plt.tight_layout()
    plt.savefig(path + '/USE_plots/lsta_hours_'+keys[l]+'_'+str(hour)+'_15000km_-65.png')
    plt.show()
    #plt.close()

def plot_map(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_75000_coldCore.p", "rb"))

    scales = dic['scales']
    nbcores = dic['nbcores']
    nbrcores = dic['nbrcores']
    del dic['scales']
    del dic['nbcores']
    del dic['nbrcores']


    keys = list(dic.keys())
    cnt = (dic['SN-pos'][0]).shape[0]

    l=0
    dist=100
    pos = np.array([1,4,6,10])
    print(scales[pos])
    #ipdb.set_trace()
    kernel = (dic['kernel'])[0] / dic['cnt'][0]
    random = (dic['kernel'])[1] / dic['cnt'][1]
    lsta = (dic['lsta'])[0] / dic['cnt'][0][0,:,:]
    mask = (dic['kernel'])[2]
    extent = ((dic['lsta'][0]).shape[1] - 1) / 2
    dist=100

    f = plt.figure(figsize=(13, 6))
    ax = f.add_subplot(121)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , lsta , cmap='RdBu_r', extend='both',levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5])
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[3,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
    ax.plot(0,0, 'bo')
    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(122)
    print('averaged: ',scales[0], scales[1])
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , kernel[2:3,:,:].sum(axis=0)-random[2:3,:,:].sum(axis=0), cmap='RdBu_r', extend='both', levels=[-0.06,-0.04,-0.02,0.02,0.04,0.06])
    plt.plot(0,0,'bo')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[0,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)


    plt.tight_layout()
    plt.savefig(path + '/USE_plot/lsta_hours_'+keys[l]+'_'+str(hour)+'.png')
    plt.show()

def plot_all():
    hours =  [14,15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]
    for h in hours:
        plot_map_full(h)

def plot_map_full(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_25000.p", "rb")) #coeffs_nans_stdkernel_USE_"+str(hour)+"UTC.p", "rb"))

    scales = dic['scales']
    nbcores = dic['nbcores']

    del dic['scales']
    del dic['nbcores']
    del dic['nbrcores']

    #ipdb.set_trace()
    keys = list(dic.keys())
    cnt = (dic['SN-pos'][0]).shape[0]

    l=0
    dist=100
    print('scales', scales)
    pos = np.array([1,4,6,8])
    print(scales[pos])

    #ipdb.set_trace()
    kernel = (dic['kernel'])[0] / dic['cnt'][0]
    random = (dic['kernel'])[1] / dic['cnt'][1]
    lsta = (dic['lsta'])[0] / dic['lsta'][1]#dic['lsta'][1]#[0,:,:]
    mask = (dic['kernel'])[2]
    extent = ((dic['lsta'][0]).shape[1] - 1) / 2
    dist=100

    f = plt.figure(figsize=(13, 6))
    ax = f.add_subplot(231)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 ,lsta , cmap='RdBu_r', extend='both',levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5])
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[3,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
    ax.plot(0,0, 'bo', markersize=3)
    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(232)
    print('averaged: ',scales[0], scales[1])
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , kernel.mean(axis=0)-random.mean(axis=0), cmap='RdBu_r', extend='both', levels=[-0.07,-0.06,-0.05,-0.04,-0.03,0.03,0.04,0.05,0.06,0.07])
    plt.plot(0,0,'bo', markersize=3)
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[0,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('Mean of coefficients', fontsize=10)


    for ids, scale in enumerate(scales[pos]):

        ax = f.add_subplot(2,3,ids+3)

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3,   #- random[ids, :, :]
                     kernel[ids, :, :], cmap='RdBu_r', extend='both') # - random[ids, :, :]
        plt.plot(0, 0, 'bo', markersize=3)
        plt.colorbar(label='Power difference (Blob-random)')
        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, mask[0, :, :],
                     colors='none', hatches='.', levels=[0.5, 1], linewidth=0.25)

        ax.set_xlabel('km')
        ax.set_ylabel('km')

        plt.title('Scale: '+str(np.round(scale,1)), fontsize=10)




    plt.tight_layout()
    #plt.savefig(path + '/USE_plots/coefficients_both_'+str(hour)+'_75000km_-60.png')
    plt.show()


def plot_gewex2():
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'

    f = plt.figure(figsize=(12,7))


    for id, h in enumerate([14,15,16,17]):

        dic = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_"+str(h)+"UTC_25000.p", "rb"))

        scales = dic['scales']
        nbcores = dic['nbcores']
        nbrcores = dic['nbrcores']
        del dic['scales']
        del dic['nbcores']
        del dic['nbrcores']
        del dic['kernel']

        keys = list(dic.keys())
        cnt = (dic['SN-pos'][0]).shape[0]

        # ipdb.set_trace()

        # for l in keys:
        #     if keys == 'scales':
        #         continue
        #
        #     (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
        #     (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)

        l = 0
        dist = 100

        snblob = (dic[keys[l]])[0]  # / (dic[keys[l]])[3]).T#-(dic[keys[l+2]])[0]
        snrandom = (dic[keys[l]])[1]  # / (dic[keys[l]])[4]).T#-(dic[keys[l+2]])[1]
        snmask = (dic[keys[l]])[2]  # -(dic[keys[l+2]])[2]
        # snblob_std = (dic[keys[l]])[3]
        # snrandom_std = (dic[keys[l]])[4]

        weblob = (dic[keys[l + 1]])[0]  # / (dic[keys[l+1]])[3]).T#-(dic[keys[l+3]])[0]
        werandom = (dic[keys[l + 1]])[1]  # / (dic[keys[l+1]])[4]).T#-(dic[keys[l+3]])[1]
        wemask = (dic[keys[l + 1]])[2]  # -(dic[keys[l+3]])[2]
        # weblob_std = (dic[keys[l+1]])[3]
        # werandom_std = (dic[keys[l+1]])[4]

        l = 2
        dist = 100

        ax = f.add_subplot(2,4, id+1)

        rand=False
        if rand:
            a = (snblob-snrandom)

            b= (weblob - werandom)
        else:
            a = (snblob )

            b = (weblob )
            # b[b < 0] -= 0.03
            # we_mask[b<-0.05]=1

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, a , cmap='RdBu_r', levels=[-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04, -0.03,-0.01, 0.01, 0.03,0.04,0.05,0.06,0.07,0.08,0.09, 0.1], extend='both')
       # plt.colorbar(label='Power difference (Dry-wet)')
        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels=[0.5, 1],
                     linewidth=0.25)
        ax.set_xticks(np.array([-3, -2, -1, 0, 1, 2, 3]) * 100)


        #ax.set_xlabel('Cross-section (km)')
        if id ==0:
            ax.set_ylabel('Scales (km)')

        plt.title('South-North | ' + str(h).zfill(2) + '00UTC',
                  fontsize=10)

        ax = f.add_subplot(2,4,id+4+1)

        mp1 = plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales,
                    b , cmap='RdBu_r',
                     levels=[-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04, -0.03,-0.01, 0.01, 0.03,0.04,0.05,0.06,0.07,0.08,0.09, 0.1], extend='both')
       # plt.colorbar(label='Power difference (Dry-wet)')

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, wemask, colors='none', hatches='.',
                     levels=[0.5, 1], linewidth=0.1)

        ax.set_xlabel('Cross-section (km)')
        ax.set_xticks(np.array([-3, -2, -1,  0, 1,2, 3])*100)
        #ax.set_xticklabels([-3, -2, -1, -0.5, 0,  0.5,1,2, 3])
        if id == 0:
            ax.set_ylabel('Scales (km)')

        plt.title('West-East | ' + str(h).zfill(2) + '00UTC', fontsize=10)

    plt.tight_layout()

    f.subplots_adjust(right=0.87)
    cax = f.add_axes([0.89, 0.57, 0.02, 0.38])
    cbar = f.colorbar(mp1, cax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('Coefficient mean', fontsize=12)

    cax = f.add_axes([0.89, 0.08, 0.02, 0.38])
    cbar = f.colorbar(mp1, cax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('Coefficient mean', fontsize=12)


    plt.show()
