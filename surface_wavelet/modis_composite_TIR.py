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


def loop():

    for l in np.arange(0,24):
        print('Doing '+str(l))
        composite(l)

def composite(h):
    pool = multiprocessing.Pool(processes=8)


    file = '/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'


    # nightp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_0-3UTC_centrePoint.nc'
    # dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'

    hour = h

    # if hour > 6:
    #     file = dayp
    # else:
    #     file = nightp

    msg = xr.open_dataarray(file)

    msg = msg[(msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) ]


    msg = msg.sel(lat=slice(10,18), lon=slice(-10,10))

    #
    res = pool.map(file_loop, (msg))
    pool.close()

    res_list = []
    cnt_list = []
    # pdb.set_trace()
    # for m in msg[0:2]:
    #     res, cnt = file_loop(m)
    #     res_list.append(res)
    #     cnt_list.append(cnt)

    res = [x for x in res if x is not None]

    res_list = []
    res2_list = []
    res3_list = []
    cnt_list = []


    for r in res:
        res_list.append(np.squeeze(r[0]))
        res2_list.append(np.squeeze(r[1]))
        res3_list.append(np.squeeze(r[2]))
        cnt_list.append(np.squeeze(r[3]))


    if res_list  == []:
        return

    kernel_sum = np.nansum(np.squeeze(np.stack(res_list, axis=0)), axis=0)
    kernel2_sum = np.nansum(np.squeeze(np.stack(res2_list, axis=0)), axis=0)
    kernel3_sum = np.nansum(np.squeeze(np.stack(res3_list, axis=0)), axis=0)
    cnt_sum = np.nansum(np.squeeze(np.stack(cnt_list, axis=0)), axis=0)

    if kernel_sum.ndim != 2:
        return
    f = plt.figure(figsize=(9, 7))
    ax = f.add_subplot(221)
    xdist = 60
    plt.contourf(cnt_sum, cmap='viridis',   extend='both')
    plt.plot(xdist, xdist, 'bo')

    # ax.set_xticklabels(np.array((np.linspace(0, 121, 5) - 100) * 3, dtype=int))
    # ax.set_yticklabels(np.array((np.linspace(0, 121, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Local Anomaly, Nb cores: ' + str(np.max(cnt_sum)) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf(kernel3_sum / cnt_sum  ,  levels=np.arange(-75,-64), cmap='viridis',  extend='both' )
    plt.plot(xdist, xdist, 'bo')
    # ax.set_xticklabels(np.array((np.linspace(0, 121, 5) - 100) * 3, dtype=int))
    # ax.set_yticklabels(np.array((np.linspace(0, 121, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional Anomaly', fontsize=10)

    ax = f.add_subplot(223)

    plt.contourf(kernel2_sum / cnt_sum, levels=np.arange(-75,-64),cmap='viridis',   extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(100, 100, 'bo')
    # ax.set_xticklabels(np.array((np.linspace(0, 121, 5) - 100) * 3, dtype=int))
    # ax.set_yticklabels(np.array((np.linspace(0, 121, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal Anomaly',
              fontsize=10)

    ax = f.add_subplot(224)

    plt.contourf(kernel2_sum / cnt_sum, cmap='viridis', extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(100,100, 'bo')
    # ax.set_xticklabels(np.array((np.linspace(0, 121, 5) - 100) * 3, dtype=int))
    # ax.set_yticklabels(np.array((np.linspace(0, 121, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal Anomaly - Random',
              fontsize=10)

    plt.tight_layout()

    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/scales/new/composites_lsta/'+'TIRR_'+str(hour).zfill(2)+'00UTC')
    #plt.close()

    # f = plt.figure(figsize=(9, 4))
    # ax = f.add_subplot(121)
    #
    # plt.hist(((kernel_sum / cnt_sum) - (rkernel_sum / rcnt_sum))[50,50])
    # plt.title('Local Anomaly, Nb cores: ' + str(np.max(cnt_sum)) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
    #           fontsize=10)
    #
    # ax = f.add_subplot(122)
    #
    # plt.hist(((kernel3_sum / cnt_sum) - (rkernel3_sum / rcnt_sum)[50,50]))
    # plt.title('Regional Anomaly', fontsize=10)
    #
    #
    # plt.tight_layout()


def cut_kernel(xpos, ypos, lsta_day2):
    mdist = 2

    mini_mean = lsta_day2.isel(lon=slice(xpos - mdist, xpos + mdist + 1),
                               lat=slice(ypos - mdist, ypos + mdist + 1)).values

    mini_mmean = np.nanmean(mini_mean)
    if np.isnan(mini_mmean):
        #pdb.set_trace()
        return

    dist = 60

    # try:
    kernel = lsta_day2.isel(lon=slice(xpos - dist, xpos + dist + 1),
                            lat=slice(ypos - dist, ypos + dist + 1)).values

    # kernel2 = lsta_day2.isel(lon=slice(xpos - dist, xpos + dist + 1),
    #                         lat=slice(ypos - dist, ypos + dist + 1)).values

    if kernel.shape != (1, 121, 121):
        return

        kernel = np.zeros([1, 121, 121]) *np.nan

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

        if xpos + dist < lsta_day2.shape[2]:
            xmax = kernel.shape[2]
            xmaxdist = dist + 1
        else:
            xmax = dist - (xpos - lsta_day2.shape[2])
            xmaxdist = dist - (xpos + dist - lsta_day2.shape[2])

        if ypos + dist < lsta_day2.shape[1]:
            ymax = kernel.shape[1]
            ymaxdist = dist + 1
        else:
            ymax = dist - (ypos - lsta_day2.shape[1])
            ymaxdist = dist - (ypos + dist - lsta_day2.shape[1])

        cutk = lsta_day2[:, ypos - ymindist: ypos + ymaxdist, xpos - xmindist: xpos + xmaxdist].values


        kernel[:, ymin: ymax, xmin:xmax] = cutk

    ano = kernel - mini_mmean
    kernel3 = kernel - np.nanmean(kernel)

    cnt = np.zeros_like(kernel) 
    cnt[np.isfinite(kernel)] = 1


    return ano, kernel, kernel3, cnt



def file_loop(fi):


    file2 = '/users/global/cornkle/MCSfiles/blob_map_MCSs_-50_JJAS_points_dominant_gt15k.nc'
    mcs = xr.open_dataarray(file2)
    mcs = mcs.sel(lat=slice(10, 18), lon=slice(-10, 10))
    try:
        mcs = mcs[mcs['time']==fi['time']]
    except ValueError:
        return
    mcs.values[mcs.values>-50]=-50

    print('Doing day: ', fi['time'])
    pos = np.where( (fi.values >= 5) & (fi.values < 75))#(fi.values >= 1) & (fi.values <= 20)) #<-50)#

    if np.sum(pos) == 0:
        print('No blobs found')
        return

    kernel_list = []
    kernel2_list = []
    kernel3_list = []
    cnt_list = []


    for y, x in zip(pos[0], pos[1]):

        try:
            ano, kernel2, kernel3, cnt = cut_kernel(x, y, mcs)
        except TypeError:
            continue

        kernel_list.append(ano)
        kernel2_list.append(kernel2)
        kernel3_list.append(kernel3)
        cnt_list.append(cnt)

    if kernel_list == []:
        return None
    print(len(kernel_list))
    if len(kernel_list) == 1:
      return None
    else:
        kernel_sum = np.nansum(np.squeeze(np.stack(kernel_list, axis=0)), axis=0)
        kernel2_sum = np.nansum(np.squeeze(np.stack(kernel2_list, axis=0)), axis=0)
        kernel3_sum = np.nansum(np.squeeze(np.stack(kernel3_list, axis=0)), axis=0)
        cnt_sum = np.nansum(np.squeeze(np.stack(cnt_list, axis=0)), axis=0)

    print('Returning')

    return (kernel_sum, kernel2_sum, kernel3_sum, cnt_sum)

if __name__ == "__main__":
    composite()
