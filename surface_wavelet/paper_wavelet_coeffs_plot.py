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

def plot_all_cross():

    hours = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7] #15,16,
    for h in hours:
        plot(h)

def merge():

        hours = [15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]  # 15,16,
        for h in hours:
            path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
            dic = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_" + str(h) + "UTC_15000_-60.p", "rb"))

            for l in dic.keys():
                if (l == 'scales') | ('core' in l):
                    continue
                print(l)
                if 'pos' in l:
                    (dic[l])[0] = np.nanmedian((dic[l])[0], axis=0)
                    (dic[l])[1] = np.nanmedian((dic[l])[1], axis=0)
                else:
                    (dic[l])[0] = np.nanmedian((dic[l])[0], axis=0)
                    try:
                        (dic[l])[1] = np.nanmedian((dic[l])[1], axis=0)
                    except IndexError:
                        continue

            pkl.dump(dic, open(path + "/coeffs_nans_stdkernel_USE_" + str(h) + "UTC_15000_-60_merge_median.p", "wb"))

def plot(hour):
    name = 'UTC_15000_-60_AMSRE'#'UTC_15000_-60_AMSRE'
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+name+".p", "rb"))

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
    #     if l == 'scales':
    #         continue
    #     if 'pos' in l:
    #         (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
    #         (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)
    #     else:
    #         (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
    #         try:
    #             (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)
    #         except IndexError:
    #             continue

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
    levels= np.linspace(-0.3,0.3,20)#[-0.14, -0.12,-0.1, -0.08,-0.06, -0.04, -0.02,0.02,0.04,0.06,0.08,0.1, 0.12, 0.14]
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales ,  (snrandom) , cmap='RdBu_r', levels = levels, extend='both')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels = [0.5,1])

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3,scales,   (werandom) , cmap='RdBu_r', levels =levels, extend='both')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, wemask, colors='none', hatches='.', levels = [0.5,1])

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)


    ax = f.add_subplot(223)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales,  snblob,  cmap='RdBu_r', levels =levels, extend='both') #, vmin = -0.1, vmax=0.1)
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels=[0.5, 1],
                 linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbrcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(224)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, weblob    , cmap='RdBu_r', levels =levels, extend='both')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales,wemask, colors='none', hatches='.', levels=[0.5, 1],
                 linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)

    plt.tight_layout()
    plt.savefig(path + '/initVSprop/'+name+'_cross_'+str(hour)+'.png')
    plt.show()
    #plt.close()


def plot_diurnal():

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    name = 'UTC_15000_ALL_-60_5slotSmall'
    hours = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]#,7,8,9]#,10,11,12,13]
    dist = 100

    f = plt.figure(figsize=(12,5))
    dnames = ['WE-pos', 'SN-pos']
    pcoll = []
    for ids, s in enumerate([109, 109]):
    #for ids, s in enumerate([58, 58]):

        ax = f.add_subplot(1,2,ids+1)

        arr = np.zeros((len(hours), 2*100+1))

        for hids, h in enumerate(hours):

            print('Doing ', h)

            dic = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_" + str(h) + name + ".p", "rb"))

            scales = dic['scales']
            nbcores = dic['nbcores']
            nbrcores = dic['nbrcores']
            del dic['scales']
            del dic['nbcores']
            del dic['nbrcores']
            del dic['kernel']

            pos = np.where(np.isclose(scales, s, atol=4))
            spos = int(pos[0])
            print(scales)
            print('Scale', scales[spos])

            keys = list(dic.keys())


            # for l in keys:
            #     if l == 'scales':
            #         continue
            #     if 'pos' in l:
            #         (dic[l])[0] = np.nanmedian((dic[l])[0], axis=0)
            #         (dic[l])[1] = np.nanmedian((dic[l])[1], axis=0)
            #     else:
            #         (dic[l])[0] = np.nanmedian((dic[l])[0], axis=0)
            #         try:
            #             (dic[l])[1] = np.nanmedian((dic[l])[1], axis=0)
            #         except IndexError:
            #             continue
            #ipdb.set_trace()
            arr[hids,:] = (dic[dnames[ids]])[0][spos,:]

        levels = list(np.arange(-0.3, 0, 0.05)) +  list(np.arange(0.05, 0.31, 0.05))
        plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, np.arange(len(hours)), arr , cmap='RdBu_r', levels=levels,extend='both') # [-0.14,-0.12,-0.1, -0.08,-0.06, -0.04, -0.03,0.03,0.04,0.06,0.08,0.1,0.12,0.14]
        plt.colorbar(label='Wavelet coefficients')
        plt.contour((np.arange(0, 2*dist+1) - dist) * 3, np.arange(len(hours)), arr, levels = levels, colors='k', linewidths=0.5)

        ax.set_xlabel('km')
        ax.set_ylabel('Time of day')
        ax.set_yticklabels([15,17,19,21,23,1,3,5,7,9])
        plt.axvline(x=0, linestyle='dashed', color='blue')

        plt.title(dnames[ids]+' '+str(int(np.round(scales[spos]))), fontsize=10)

        pcoll.append(arr)
    dout = {
        'plot' : pcoll,
        'hours' : hours,
        'dist' : dist
    }
    plt.tight_layout()
    plt.savefig(path + 'initVSprop/' + name + '_HOV_DIURNAL.png')
    plt.show()
    plt.close()

    #return dout


def plot_map(hour):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_-60_AMSRE.p", "rb")) #UTC_15000_ALL_-60_5slotSmall

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
    pos = np.array([2,4,6,8])#np.array([1,4,6,10])
    print(scales[pos])
    #ipdb.set_trace()
    kernel = (dic['kernel'])[0] / dic['cnt'][0]
    random = (dic['kernel'])[1] / dic['cnt'][1]
    lsta = (dic['lsta'])[0] / dic['cnt'][0][0,:,:]
    mask = (dic['kernel'])[2]
    extent = ((dic['lsta'][0]).shape[1] - 1) / 2
    dist=100
    levels = list(np.arange(-0.3, 0, 0.05)) + list(np.arange(0.05, 0.31, 0.05))

    f = plt.figure(figsize=(10, 8), dpi=200)
    ax = f.add_subplot(221)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , lsta , cmap='RdBu_r', extend='both',levels=list(np.arange(-1, 0, 0.2)) + list(np.arange(0.2, 1.2, 0.2)))
    plt.colorbar(label='K')
    #plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[3,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
    ax.plot(0,0, 'bo')
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.axvline(x=0, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)

    plt.title('LSTA | Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC',
              fontsize=10)

    ax = f.add_subplot(222)
    print('averaged: ',scales[0], scales[1])
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , kernel[1,:,:], cmap='RdBu_r', extend='both', levels=levels)
    plt.plot(0,0,'bo')
    plt.colorbar(label='Wavelet coefficients')
    #plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[0,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
    plt.axvline(x=0, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('30km', fontsize=10)

    ax = f.add_subplot(223)
    print('averaged: ',scales[0], scales[1])
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , kernel[2,:,:], cmap='RdBu_r', extend='both', levels=levels)
    plt.plot(0,0,'bo')
    plt.colorbar(label='Wavelet coefficients')
    #plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[0,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
    plt.axvline(x=0, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('59km', fontsize=10)

    ax = f.add_subplot(224)
    print('averaged: ',scales[0], scales[1])
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , kernel[3,:,:], cmap='RdBu_r', extend='both', levels=levels)
    plt.plot(0,0,'bo')
    plt.colorbar(label='Wavelet coefficients')
    #plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[0,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
    plt.axvline(x=0, linestyle='dashed', color='k', linewidth=1)
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('109km', fontsize=10)


    plt.tight_layout()
    plt.savefig(path + '/initVSprop/wcoeff_maps_all_AMSRE_'+str(hour)+'.png')
    plt.show()

def plot_all():
    hours =  [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]
    for h in hours:
        plot_map_full(h)

def plot_map_full(hour, amsre=False):

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    name = 'UTC_15000_ALL_-60_5slotSmall'
    if amsre:
        tag = '_AMSRE'
    else:
        tag = ''
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+name+tag+".p", "rb")) #coeffs_nans_stdkernel_USE_"+str(hour)+"UTC.p", "rb"))

    scales = dic['scales']
    nbcores = dic['nbcores']

    del dic['scales']
    del dic['nbcores']
    del dic['nbrcores']

    keys = list(dic.keys())

    #ipdb.set_trace()
    keys = list(dic.keys())
    cnt = (dic['SN-pos'][0]).shape[0]

    l=0
    dist=100
    print('scales', scales)
    pos = np.array([2,4,6,8])
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

    if amsre:
        cmap = 'RdBu'
    else:
        cmap = 'RdBu_r'

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 ,lsta , cmap=cmap, extend='both',levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5])
    plt.colorbar(label='Wavelet coefficient')
    #plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[3,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
    ax.plot(0,0, 'bo', markersize=3)
    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('Spatial composite, Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(232)
    print('averaged: ',scales[0], scales[1])
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , kernel.mean(axis=0)-random.mean(axis=0), cmap=cmap, extend='both', levels= [-0.1, -0.08,-0.06, -0.04, -0.02,0.02,0.04,0.06,0.08,0.1])
    plt.plot(0,0,'bo', markersize=3)
    plt.colorbar(label='Wavelet coefficient')
    #plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3, mask[0,:,:], colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('Mean of coefficients', fontsize=10)


    for ids, scale in enumerate(scales[pos]):

        ax = f.add_subplot(2,3,ids+3)

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3,   #- random[ids, :, :]
                     kernel[ids, :, :], cmap=cmap, extend='both',levels= [-0.1, -0.08,-0.06, -0.04, -0.02,0.02,0.04,0.06,0.08,0.1]) # - random[ids, :, :]
        plt.plot(0, 0, 'bo', markersize=3)
        plt.colorbar(label='Wavelet coefficient')
        #plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, mask[0, :, :],
        #             colors='none', hatches='.', levels=[0.5, 1], linewidth=0.25)

        ax.set_xlabel('km')
        ax.set_ylabel('km')

        plt.title('Scale: '+str(np.round(scale,1)), fontsize=10)




    plt.tight_layout()
    plt.savefig(path + 'initVSprop/'+name+'_map_'+str(hour)+'.png')
    plt.show()


def plot_gewex2():

    f = plt.figure(figsize=(12,7))

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    name = 'UTC_15000_ALL_-60_5slotSmall'

    for id, h in enumerate([17,20,23,2]):

        dic = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_" + str(h) + name + ".p", "rb"))

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
        #     if l == 'scales':
        #         continue
        #     if 'pos' in l:
        #         (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
        #         (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)
        #     else:
        #         (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
        #         try:
        #             (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)
        #         except IndexError:
        #             continue

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
        levels = list(np.arange(-0.3, 0, 0.05)) + list(np.arange(0.05, 0.31, 0.05))

        ax = f.add_subplot(2,4, id+1)

        rand=True
        if rand:
            a = (snblob-snrandom)

            b= (weblob - werandom)
        else:
            a = (snblob )

            b = (weblob )
            # b[b < 0] -= 0.03
            # we_mask[b<-0.05]=1

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, a , cmap='RdBu_r', levels =levels, extend='both')
       # plt.colorbar(label='Power difference (Dry-wet)')
        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels=[0.5, 1],
                     linewidth=0.25)

        #plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, scales, a, colors='k', levels=levels, extend='both', linewidths=0.5)
        plt.axvline(x=0, linestyle='dashed', color='blue')
        ax.set_xticks(np.array([-3, -2, -1, 0, 1, 2, 3]) * 100)


        #ax.set_xlabel('Cross-section (km)')
        if id ==0:
            ax.set_ylabel('Scales (km)')

        plt.title('South-North | ' + str(h).zfill(2) + '00UTC',
                  fontsize=10)

        ax = f.add_subplot(2,4,id+4+1)

        mp1 = plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales,
                    b , cmap='RdBu_r',
                     levels = levels, extend='both')
       # plt.colorbar(label='Power difference (Dry-wet)')

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, wemask, colors='none', hatches='.',
                     levels=[0.5, 1], linewidth=0.1)

        #plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, scales, b, colors='k', levels=levels, extend='both',
        #            linewidths=0.5)

        plt.axvline(x=0, linestyle='dashed', color='blue')

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
    plt.savefig(path + 'initVSprop/' + name + '_CROSS_DIURNAL_nocontour.png')


def paper_plot():
    hour= 17
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_-60_merge.p", "rb"))

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
    pos = np.array([2,4,6,8])#np.array([1,4,6,10])
    print(scales[pos])
    #ipdb.set_trace()
    kernel = (dic['kernel'])[0] / dic['cnt'][0]
    random = (dic['kernel'])[1] / dic['cnt'][1]
    lsta = (dic['lsta'])[0] / dic['cnt'][0][0,:,:]
    mask = (dic['kernel'])[2]
    extent = ((dic['lsta'][0]).shape[1] - 1) / 2
    dist=100

    f = plt.figure(figsize=(13,13))
    gs = GridSpec(4,4, figure=f)

    ax1 = f.add_subplot(gs[0,0])

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, (np.arange(0, 2*dist+1) - dist) * 3 , lsta , cmap='RdBu_r', extend='both',levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5])
    plt.colorbar(label='K')
    ax1.plot(0,0, 'bo')
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')
    plt.title('LSTA, Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC',
              fontsize=10)
    coords = [(0, 0), (0, 1), (1, 0), (1,1)]
    for ids, scale in enumerate(scales[pos]):

        if ids == 0:
            continue

        ax = f.add_subplot(gs[coords[ids][0],coords[ids][1]])

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3,   #- random[ids, :, :]
                     kernel[ids, :, :], cmap='RdBu_r', extend='both',levels= [-0.1, -0.08,-0.06, -0.04, -0.02,0.02,0.04,0.06,0.08,0.1]) # - random[ids, :, :]
        plt.plot(0, 0, 'bo', markersize=3)
        plt.colorbar(label='')
        plt.title('Length scale: '+str(int(np.round(scale,1))))
        ax.set_xlabel('km')
        ax.set_ylabel('km')


    dout = plot_diurnal()
    titles = ['East-West transect', 'South-North transect']
    coords = [(0,2), (2,4)]
    for ids, p in enumerate(dout['plot']):

        ax = f.add_subplot(gs[ids, 2:4])

        plt.contourf((np.arange(0, 2 * dout['dist'] + 1) - dout['dist']) * 3, np.arange(len(dout['hours'])), p, cmap='RdBu_r',
                     levels=[-0.15, -0.12, -0.09, -0.06, -0.03, 0.03, 0.06, 0.09, 0.12,
                             0.15])  # [-0.14,-0.12,-0.1, -0.08,-0.06, -0.04, -0.03,0.03,0.04,0.06,0.08,0.1,0.12,0.14]
        plt.colorbar(label='Wavelet coefficients')
        plt.contour((np.arange(0, 2 * dout['dist'] + 1) - dout['dist']) * 3, np.arange(len(dout['hours'])), p,
                    levels=[-0.15, -0.12, -0.09, -0.06, -0.03, 0.03, 0.06, 0.09, 0.12, 0.15], colors='k',
                    linewidths=0.5)

        ax.set_xlabel('km')
        ax.set_ylabel('Time of day')
        ax.set_yticklabels([15, 17, 19, 21, 23, 1, 3, 5, 7, 9])
        plt.vlines(0, ymin=0, ymax=15, linestyle='dashed')

        plt.title(titles[ids] + ' ' + '109km', fontsize=10)


    for id, h in enumerate([16,17,19,21]):

        path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
        dic = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_" + str(h) + "UTC_25000.p", "rb"))

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
        #     if l == 'scales':
        #         continue
        #     if 'pos' in l:
        #         (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
        #         (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)
        #     else:
        #         (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
        #         try:
        #             (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)
        #         except IndexError:
        #             continue

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

        ax = f.add_subplot(gs[2,id])

        rand=True
        if rand:
            a = (snblob-snrandom)

            b= (weblob - werandom)
        else:
            a = (snblob )

            b = (weblob )
            # b[b < 0] -= 0.03
            # we_mask[b<-0.05]=1

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, a , cmap='RdBu_r', levels = [-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1], extend='both')
       # plt.colorbar(label='Power difference (Dry-wet)')
        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels=[0.5, 1],
                     linewidth=0.25)
        ax.set_xticks(np.array([-3, -2, -1, 0, 1, 2, 3]) * 100)


        #ax.set_xlabel('Cross-section (km)')
        if id ==0:
            ax.set_ylabel('Scales (km)')

        plt.title('South-North | ' + str(h).zfill(2) + '00UTC',
                  fontsize=10)

        ax = f.add_subplot(gs[3,id])

        mp1 = plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales,
                    b , cmap='RdBu_r',
                     levels = [-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1], extend='both')
        plt.colorbar(label='Wavelet coefficients')

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, wemask, colors='none', hatches='.',
                     levels=[0.5, 1], linewidth=0.1)

        ax.set_xlabel('Cross-section (km)')
        ax.set_xticks(np.array([-3, -2, -1,  0, 1,2, 3])*100)
        #ax.set_xticklabels([-3, -2, -1, -0.5, 0,  0.5,1,2, 3])
        if id == 0:
            ax.set_ylabel('Scales (km)')

        plt.title('West-East | ' + str(h).zfill(2) + '00UTC', fontsize=10)

    plt.tight_layout()

    # f.subplots_adjust(right=0.87)
    # cax = f.add_axes([0.89, 0.57, 0.02, 0.38])
    # cbar = f.colorbar(mp1, cax)
    # cbar.ax.tick_params(labelsize=12)
    # cbar.set_label('Coefficient mean', fontsize=12)
    #
    # cax = f.add_axes([0.89, 0.08, 0.02, 0.38])
    # cbar = f.colorbar(mp1, cax)
    # cbar.ax.tick_params(labelsize=12)
    # cbar.set_label('Coefficient mean', fontsize=12)
    #plt.savefig(path + '/paper/wcoeff_hours_'+keys[l]+'_'+str(hour)+'.png')
    plt.show()