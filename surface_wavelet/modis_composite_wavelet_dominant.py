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
from wavelet import util
import itertools
from scipy.stats import ttest_ind as ttest
from scipy.interpolate import griddata

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1


def composite():
    pool = multiprocessing.Pool(processes=8)

    nightp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_0-3UTC_centrePoint.nc'
    dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'

    # nightp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_0-3UTC.nc'
    # dayp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_15-18UTC.nc'

    hour = 18

    if hour > 6:
        file = dayp
    else:
        file = nightp

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 7) ]

    msg = msg.sel(lat=slice(10,20), lon=slice(-10, 10))

    res = pool.map(file_loop, msg)
    pool.close()
    # cnt_list = []
    # res_list = []
    # for m in msg[0:2]:
    #     res, cnt = file_loop(m)
    #     res_list.append(res)
    #     cnt_list.append(cnt)

    res = [x for x in res if x is not None]

    res_list = []
    res2_list = []


    for r in res:
        res_list.append(np.squeeze(r[0]))
        res2_list.append(np.squeeze(r[1]))


    #pdb.set_trace()
    k = np.squeeze(np.vstack(res_list))
    rk = np.squeeze(np.vstack(res2_list))


    k_mean = np.nanmean(k, axis=0)
    rk_mean = np.nanmean(rk, axis=0)

    f = plt.figure(figsize=(9, 7))
    ax = f.add_subplot(121)

    plt.contourf(k_mean, cmap='RdBu_r')
    plt.plot(50, 50, 'bo')
    ax.set_xticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_yticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Mean scale' + str(k.shape[0]) + '| ' + str(hour).zfill(2) + '00UTC, Aug-Sep',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf(k_mean-rk_mean, cmap='RdBu_r')
    plt.plot(50, 50, 'bo')
    ax.set_xticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_yticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Mean scale-random', fontsize=10)



    plt.tight_layout()


def cut_kernel(xpos, ypos, wl):

    dist = 50
    #pdb.set_trace()
    try:
        kernel = wl[ypos - dist: ypos + dist + 1,  xpos - dist: xpos + dist + 1]
       # pdb.set_trace()
    except ValueError:
        # print('Point at boundary, pass')
        return

    if kernel.shape != (101, 101):

        return

    return kernel


def file_loop(fi):

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))
    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 6:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date #+ dayd

    lsta = xr.open_dataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/' \
                           'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    print('Doing '+ 'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    lsta_day = lsta['LSTA']
    lsta_day2 = lsta_day.copy()
    lsta_day2 = lsta_day2.where(lsta_day2 > -900)

    pos = np.where((fi.values >= 1) & (fi.values <= 20))
    if np.sum(pos) == 0:
        print('No blobs found')
        return

    wav_input = lsta_day2-lsta_day2.mean()
    wav_input = np.squeeze(wav_input.values)
    points = np.where(np.isfinite(wav_input))
    inter = np.where(np.isnan(wav_input))

    # interpolate
    wav_input[inter] = 0#griddata(points, np.ravel(wav_input[points]), inter, method='nearest')

    wav = util.waveletLSTA(np.squeeze(wav_input), 3,method='dry')
    wl = wav['dominant']
    wl[inter]=np.nan

    xfi = fi.shape[1]
    yfi = fi.shape[0]
    randx = np.random.randint(0, xfi, 250)
    randy = np.random.randint(0, yfi, 250)
    posr = (randy, randx)

##############################Random loop
    rk_list = []
    for y, x in zip(posr[0], posr[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        point = lsta_day.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_day['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_day['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            kernel = cut_kernel(xpos, ypos, wl)
        except TypeError:
            continue

        rk_list.append(kernel)  # north-south wavelet


###############################Blob loop
    k_list = []

    for y, x in zip(pos[0], pos[1]):


        lat = fi['lat'][y]
        lon = fi['lon'][x]

        point = lsta_day.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_day['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_day['lat'].values == plat)
        ypos = int(ypos[0])

        try:
           kernel = cut_kernel(xpos, ypos, wl)
        except TypeError:
            continue
        #print(plat,plon)
        k_list.append(kernel)

    k_list = [x for x in k_list if x is not None]
    rk_list = [x for x in rk_list if x is not None]

    if k_list == []:
        return None
    print(len(k_list))
    if (len(k_list) <= 1) | (k_list == None):
      return None
    else:
        k_sum = np.squeeze(np.stack(k_list, axis=0))
        rk_sum = np.squeeze(np.stack(rk_list, axis=0))

    scales = wav['scales']

    print('Returning')

    return (k_sum, rk_sum,  scales)


if __name__ == "__main__":
    composite()
