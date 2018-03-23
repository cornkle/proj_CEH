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
import pdb
import pandas as pd
from scipy import ndimage
from cold_cloud_trend import era_geop_t3d as era_geop
from utils import u_gis

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)

def cut_kernel(zpos, ypos, xpos, da):
    mdist = 5

    try:
        mini_mean = da.isel(time=zpos, lon=slice(xpos - mdist, xpos + mdist + 1),
                               lat=slice(ypos - mdist, ypos + mdist + 1)).values
    except ValueError:
        return

    mini_mmean = np.nanmean(mini_mean)

    if np.isnan(mini_mmean):
        # pdb.set_trace()
        return

    dist = 100

    try:
        kernel = da.isel(time = zpos, lon=slice(xpos - dist, xpos + dist + 1),
                            lat=slice(ypos - dist, ypos + dist + 1)).values
    except ValueError:
        kernel = np.array([0])

    if kernel.shape != (1, 201, 201):
        #return

        kernel = np.zeros([1, 201, 201]) *np.nan

        if xpos - dist >= 0:
            xmin = 0
            xmindist = dist
        else:
            xmin = (xpos - dist) * -1
            xmindist = dist + (xpos - dist)

        if ypos - dist >= 0:
            ymin = 0
            ymindist = dist
        else:
            ymin = (ypos - dist) * -1
            ymindist = dist + (ypos - dist)

        if xpos + dist < da.shape[2]:
            xmax = kernel.shape[2]
            xmaxdist = dist + 1
        else:
            xmax = dist - (xpos - da.shape[2])
            xmaxdist = dist - (xpos + dist - da.shape[2])

        if ypos + dist < da.shape[1]:
            ymax = kernel.shape[1]
            ymaxdist = dist + 1
        else:
            ymax = dist - (ypos - da.shape[1])
            ymaxdist = dist - (ypos + dist - da.shape[1])

        cutk = da[zpos, ypos - ymindist: ypos + ymaxdist, xpos - xmindist: xpos + xmaxdist]

        kernel[:, ymin: ymax, xmin:xmax] = cutk

    ano = kernel - mini_mmean
    kernel3 = kernel - np.nanmean(kernel)

    cnt = kernel.copy() * 0
    cnt[np.isfinite(kernel)] = 1

    return ano, kernel, kernel3, cnt


def composite():
    hour = 17
    mds = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/lsta_daily_*.nc')
    mds = mds.sel(lat=slice(10.5,18), lon=slice(-9.5,9.5))
    
    pos = np.where(mds['cell'].values==hour)

    clocal = np.zeros([1, 201, 201])
    cregional = np.zeros([1, 201, 201])
    coriginal = np.zeros([1, 201, 201])
    ccnt = np.zeros([1, 201, 201])

    for z, y, x in zip(pos[0], pos[1], pos[2]):

        try:
            local, original, regional, cnt = cut_kernel(z, y, x, mds['LSTA'])
        except TypeError:
            continue

        clocal = np.nansum(np.vstack((clocal, local)), axis=0)[np.newaxis,...]
        cregional  = np.nansum(np.vstack((cregional, regional)), axis=0)[np.newaxis,...]
        coriginal = np.nansum(np.vstack((coriginal, original)), axis=0)[np.newaxis,...]
        ccnt = np.nansum(np.vstack((ccnt, cnt)), axis=0)[np.newaxis,...]

    # xfi = mds['cell'].values.shape[2]
    # yfi = mds['cell'].values.shape[1]
    # zfi = mds['cell'].values.shape[0]
    #
    # rlocal = np.zeros([1, 201, 201])
    # rregional = np.zeros([1, 201, 201])
    # roriginal = np.zeros([1, 201, 201])
    # rcnt = np.zeros([1, 201, 201])
    # for z in np.arange(zfi):
    #     randx = np.random.randint(0, xfi, 250)
    #     randy = np.random.randint(0, yfi, 250)
    #     posr = (randy, randx)
    #     for y, x in zip(posr[0], posr[1]):
    #
    #         try:
    #             local, original, regional, cnt = cut_kernel(z, y, x, mds['LSTA'])
    #         except TypeError:
    #             continue
    #
    #         rlocal = np.nansum(np.vstack((rlocal, local)), axis=0)[np.newaxis, ...]
    #         rregional = np.nansum(np.vstack((rregional, regional)), axis=0)[np.newaxis, ...]
    #         roriginal = np.nansum(np.vstack((roriginal, original)), axis=0)[np.newaxis, ...]
    #         rcnt = np.nansum(np.vstack((rcnt, cnt)), axis=0)[np.newaxis, ...]



    plocal = np.squeeze(clocal/ccnt) #- np.squeeze(rlocal/rcnt) #
    pregional =  np.squeeze(cregional/ccnt) #- np.squeeze(rregional/rcnt) #
    poriginal =  np.squeeze(coriginal/ccnt) #- np.squeeze(roriginal/rcnt) #
    pcnt = np.squeeze(ccnt)

    f = plt.figure(figsize=(9, 7))
    ax = f.add_subplot(221)

    plt.contourf( plocal, cmap='RdBu_r',   levels=[-0.5,-0.4,-0.3,-0.2,-0.1, -0.05, 0, 0.05,0.1, 0.2,0.3,0.4,0.5],  extend='both')
    plt.plot(100, 100, 'bo')

    ax.set_xticklabels(np.array((np.linspace(0, 200, 5) - 100) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, 200, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Local Anomaly, Nb cores: ' + str(np.max(ccnt)) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf(pregional , cmap='RdBu_r',  levels=[-0.5,-0.4,-0.3,-0.2,-0.1, -0.05, 0, 0.05,0.1, 0.2,0.3,0.4,0.5], extend='both' )
    plt.plot(100, 100, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, 200, 5) - 100) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, 200, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional Anomaly', fontsize=10)

    ax = f.add_subplot(223)

    plt.contourf(pcnt, cmap='viridis',  extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(100, 100, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, 200, 5) - 100) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, 200, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal Anomaly',
              fontsize=10)

    ax = f.add_subplot(224)

    plt.contourf(poriginal, cmap='RdBu_r',  levels=[-0.5,-0.4,-0.3,-0.2,-0.1, -0.05, 0, 0.05,0.1, 0.2,0.3,0.4,0.5], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(100,100, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, 200, 5) - 100) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, 200, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal Anomaly - Random',
              fontsize=10)

    plt.tight_layout()

   # plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/scales/new/composites_lsta/'+str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
   # plt.close()




        
        