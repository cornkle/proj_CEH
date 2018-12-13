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

def run_hours():

    l = [6, 21, 0, 3, 18]
    for ll in l:
        composite(ll)



def plot(hour):

    path = cnst.network_data + 'figs/LSTA-bullshit/AGU'
    dic = pkl.load(open(path+"/coeffs_test_nans_stdall"+str(hour)+"UTC.p", "rb"))

    scales = dic['scales']
    nbcores = dic['nbcores']
    nbrcores = dic['nbrcores']
    del dic['scales']
    del dic['nbcores']
    del dic['nbrcores']
    del dic['kernel']


    keys = list(dic.keys())
    cnt = (dic['SN-pos'][0]).shape[0]

    # for l in keys:
    #     if keys == 'scales':
    #         continue
    #
    #     (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
    #     (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)

    l=0
    dist=100

    snblob = ((dic[keys[l]])[0].T / (dic[keys[l]])[3]).T#-(dic[keys[l+2]])[0]
    snrandom = ((dic[keys[l]])[1].T / (dic[keys[l]])[4]).T#-(dic[keys[l+2]])[1]
    snmask = (dic[keys[l]])[2]#-(dic[keys[l+2]])[2]
    snblob_std = (dic[keys[l]])[3]
    snrandom_std = (dic[keys[l]])[4]

    weblob = ((dic[keys[l+1]])[0].T / (dic[keys[l+1]])[3]).T#-(dic[keys[l+3]])[0]
    werandom = ((dic[keys[l+1]])[1].T / (dic[keys[l+1]])[4]).T#-(dic[keys[l+3]])[1]
    wemask = (dic[keys[l+1]])[2]#-(dic[keys[l+3]])[2]
    weblob_std = (dic[keys[l+1]])[3]
    werandom_std = (dic[keys[l+1]])[4]

    l=2
    dist=100

    f = plt.figure(figsize=(9, 9))
    ax = f.add_subplot(221)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales , (snblob) - (snrandom) , cmap='RdBu_r', vmin = -0.3, vmax=0.3)
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3,scales,   (weblob) - (werandom) , cmap='RdBu_r', levels = [-0.3,-0.2,-0.1, -0.05,-0.025, 0.025, 0.05,0.1, 0.2,0.3], extend='both')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, wemask, colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)


    ax = f.add_subplot(223)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales,  snrandom,  cmap='RdBu_r', vmin = -0.4, vmax=0.4) #, vmin = -0.1, vmax=0.1)
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels=[0.5, 1],
                 linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbrcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(224)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, werandom    , cmap='RdBu_r', vmin = -0.4, vmax=0.4) # vmin = -0.1, vmax=0.1)
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales,wemask, colors='none', hatches='.', levels=[0.5, 1],
                 linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)

    plt.tight_layout()
    #plt.savefig(path + '/lsta_hours_'+keys[l]+'_'+str(hour)+'.png')
    plt.show()

def plot_map(hour):

    path = cnst.network_data + 'figs/LSTA-bullshit/AGU'
    dic = pkl.load(open(path+"/coeffs_test_nans_stdall"+str(hour)+"UTC.p", "rb"))

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

    dist=100

    plt.figure()
    plt.imshow(dic['cnt'][0][1,:,:]-dic['cnt'][0][2,:,:], origin='lower')
    plt.show()

    f = plt.figure(figsize=(15, 8))
    ax = f.add_subplot(121)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , lsta , cmap='RdBu_r', extend='both')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[3,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(122)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , kernel[0::,:,:].mean(axis=0)-random[0::,:,:].mean(axis=0), cmap='RdBu_r', extend='both', vmin=-15, vmax=15)
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[0,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)


    plt.tight_layout()
    #plt.savefig(path + '/lsta_hours_'+keys[l]+'_'+str(hour)+'.png')
    plt.show()


def plot_gewex2():
    path = '/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/wavelet'

    f = plt.figure(figsize=(12,7))


    for id, h in enumerate([18,21,0,3]):

        dic = pkl.load(open(path + "/c_wet_dry_withzero" + str(h) + "UTC.p", "rb"))

        scales = dic['scales']
        sn_mask = dic['SN-dw_mask']
        we_mask = dic['WE-dw_mask']
        del dic['scales']
        del dic['SN-dw_mask']
        del dic['WE-dw_mask']

        keys = list(dic.keys())
        cnt = (dic['SN-pos'][0]).shape[0]


        l = 0
        dist = 100

        snblob = (dic[keys[l]])[0]  # -(dic[keys[l+2]])[0]
        snrandom = (dic[keys[l]])[1]  # -(dic[keys[l+2]])[1]
        snmask = (dic[keys[l]])[2]  # -(dic[keys[l+2]])[2]

        weblob = (dic[keys[l + 1]])[0]  # -(dic[keys[l+3]])[0]
        werandom = (dic[keys[l + 1]])[1]  # -(dic[keys[l+3]])[1]
        wemask = (dic[keys[l + 1]])[2]  # -(dic[keys[l+3]])[2]

        snmask_r = ~snmask
        wemask_r = ~wemask
        we_mask_r = ~we_mask
        sn_mask_r = ~sn_mask
        l = 2
        dist = 100

        wet_snblob = (dic[keys[l]])[0]  # -(dic[keys[l+2]])[0]
        wet_snrandom = (dic[keys[l]])[1]  # -(dic[keys[l+2]])[1]
        wet_snmask = (dic[keys[l]])[2]  # -(dic[keys[l+2]])[2]

        wet_weblob = (dic[keys[l + 1]])[0]  # -(dic[keys[l+3]])[0]
        wet_werandom = (dic[keys[l + 1]])[1]  # -(dic[keys[l+3]])[1]
        wet_wemask = (dic[keys[l + 1]])[2]  # -(dic[keys[l+3]])[2]

        ax = f.add_subplot(2,4, id+1)

        rand=False
        if rand:
            a = (snblob-snrandom) -(wet_snblob-wet_snrandom)

            b= (weblob - werandom) - (wet_weblob-wet_werandom)
        else:
            a = (snblob  - wet_snblob )

            b = (weblob - wet_weblob )
            # b[b < 0] -= 0.03
            # we_mask[b<-0.05]=1

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, a , cmap='RdBu_r', levels=[-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04, -0.03,-0.01, 0.01, 0.03,0.04,0.05,0.06,0.07,0.08,0.09, 0.1], extend='both')
       # plt.colorbar(label='Power difference (Dry-wet)')
        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, sn_mask, colors='none', hatches='.', levels=[0.5, 1],
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

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, we_mask, colors='none', hatches='.',
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
    cbar.set_label('Wavelet power (Dry-Wet)', fontsize=12)

    cax = f.add_axes([0.89, 0.08, 0.02, 0.38])
    cbar = f.colorbar(mp1, cax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('Wavelet power (Dry-Wet)', fontsize=12)


    plt.show()
