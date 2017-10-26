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

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


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
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) ]

    msg = msg.sel(lat=slice(14,20), lon=slice(-10,10))

    res = pool.map(file_loop, msg)
    pool.close()
    # res_list = []
    # cnt_list = []
    # for m in msg[0:2]:
    #     res, cnt = file_loop(m)
    #     res_list.append(res)
    #     cnt_list.append(cnt)

    res = [x for x in res if x is not None]

    res_list = []
    res2_list = []
    res3_list = []
    cnt_list = []
    rres_list = []
    rres2_list = []
    rres3_list = []
    rcnt_list = []
    for r in res:
        res_list.append(np.squeeze(r[0]))
        res2_list.append(np.squeeze(r[1]))
        res3_list.append(np.squeeze(r[2]))
        cnt_list.append(np.squeeze(r[3]))
        rres_list.append(np.squeeze(r[4]))
        rres2_list.append(np.squeeze(r[5]))
        rres3_list.append(np.squeeze(r[6]))
        rcnt_list.append(np.squeeze(r[7]))

    kernel_sum = np.nansum(np.squeeze(np.stack(res_list, axis=0)), axis=0)
    kernel2_sum = np.nansum(np.squeeze(np.stack(res2_list, axis=0)), axis=0)
    kernel3_sum = np.nansum(np.squeeze(np.stack(res3_list, axis=0)), axis=0)
    cnt_sum = np.nansum(np.squeeze(np.stack(cnt_list, axis=0)), axis=0)
    rkernel_sum = np.nansum(np.squeeze(np.stack(rres_list, axis=0)), axis=0)
    rkernel2_sum = np.nansum(np.squeeze(np.stack(rres2_list, axis=0)), axis=0)
    rkernel3_sum = np.nansum(np.squeeze(np.stack(rres3_list, axis=0)), axis=0)
    rcnt_sum = np.nansum(np.squeeze(np.stack(rcnt_list, axis=0)), axis=0)

    f = plt.figure(figsize=(9, 7))
    ax = f.add_subplot(221)

    plt.contourf((kernel_sum / cnt_sum)-(rkernel_sum / rcnt_sum), cmap='RdBu_r', vmin=-0.3, vmax=0.3)
    plt.plot(50, 50, 'bo')
    ax.set_xticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_yticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Local Anomaly, Nb cores: ' + str(np.max(cnt_sum)) + '| ' + str(hour).zfill(2) + '00UTC, Aug-Sep',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf((kernel3_sum / cnt_sum) - (rkernel3_sum / rcnt_sum) , cmap='RdBu_r', vmin=-0.3, vmax=0.3)
    plt.plot(50, 50, 'bo')
    ax.set_xticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_yticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional Anomaly', fontsize=10)

    ax = f.add_subplot(223)

    plt.contourf((kernel2_sum / cnt_sum), cmap='RdBu_r', vmin=-2, vmax=2) #-(rkernel2_sum / rcnt_sum)
    plt.plot(50, 50, 'bo')
    ax.set_xticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_yticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal Anomaly',
              fontsize=10)

    ax = f.add_subplot(224)

    plt.contourf((kernel2_sum / cnt_sum)-(rkernel2_sum / rcnt_sum), cmap='RdBu_r', vmin=-0.5, vmax=0.5) #-(rkernel2_sum / rcnt_sum)
    plt.plot(50, 50, 'bo')
    ax.set_xticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_yticklabels((np.linspace(0, 100, 6) - 50) * 3)
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal Anomaly - Random',
              fontsize=10)

    plt.tight_layout()

    f = plt.figure(figsize=(9, 4))
    ax = f.add_subplot(121)

    plt.hist(((kernel_sum / cnt_sum) - (rkernel_sum / rcnt_sum))[50,50])
    plt.title('Local Anomaly, Nb cores: ' + str(np.max(cnt_sum)) + '| ' + str(hour).zfill(2) + '00UTC, Aug-Sep',
              fontsize=10)

    ax = f.add_subplot(122)

    plt.hist(((kernel3_sum / cnt_sum) - (rkernel3_sum / rcnt_sum)[50,50]))
    plt.title('Regional Anomaly', fontsize=10)


    plt.tight_layout()


def cut_kernel(xpos, ypos, lsta_day, lsta_day2):
    mdist = 5

    # if (x<dist) or (x>lsta_day.shape[1]-dist-2):
    try:

        mini_mean = lsta_day.isel(lon=slice(xpos - mdist, xpos + mdist + 1),
                                  lat=slice(ypos - mdist, ypos + mdist + 1)).values
        mini_mean[mini_mean < -900] = np.nan
        mini_mean = np.nanmean(mini_mean)

    except ValueError:
        # print('Point at boundary, pass')
        return

    if np.isnan(mini_mean):
        return

    dist = 50
    try:
        kernel = lsta_day.isel(lon=slice(xpos - dist, xpos + dist + 1),
                               lat=slice(ypos - dist, ypos + dist + 1)).values
        kernel[kernel < -900] = np.nan
        kernel2 = lsta_day2.isel(lon=slice(xpos - dist, xpos + dist + 1),
                                 lat=slice(ypos - dist, ypos + dist + 1)).values
        kernel2[kernel2 < -900] = np.nan

    except ValueError:
        # print('Point at boundary, pass')
        return

    if kernel.shape != (1, 101, 101):
        return

    ano = kernel - mini_mean
    kernel3 = kernel - np.nanmean(kernel)


    cnt = kernel.copy()
    cnt[np.isfinite(cnt)] = 1

    return ano, kernel2, kernel3, cnt



def file_loop(fi):
    print('Doing day: ', fi)

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

    # plt.figure()
    # plt.imshow(lsta_day2[0,:,:])
    # #
    # return

    pos = np.where((fi.values >= 1) & (fi.values <= 20))
    if np.sum(pos) == 0:
        print('No blobs found')
        return

    kernel_list = []
    kernel2_list = []
    kernel3_list = []
    cnt_list = []

    xfi = fi.shape[1]
    yfi = fi.shape[0]
    randx = np.random.randint(0,xfi,300)
    randy = np.random.randint(0,yfi, 300)
    posr = (randy, randx)

    rkernel_list = []
    rkernel2_list = []
    rkernel3_list = []
    rcnt_list = []
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
            ano, kernel2, kernel3, cnt = cut_kernel(xpos, ypos, lsta_day, lsta_day2)
        except TypeError:
            continue

        rkernel_list.append(ano)
        rkernel2_list.append(kernel2)
        rkernel3_list.append(kernel3)
        rcnt_list.append(cnt)

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
            ano, kernel2, kernel3, cnt = cut_kernel(xpos, ypos, lsta_day, lsta_day2)
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

        rkernel_sum = np.nansum(np.squeeze(np.stack(rkernel_list, axis=0)), axis=0)
        rkernel2_sum = np.nansum(np.squeeze(np.stack(rkernel2_list, axis=0)), axis=0)
        rkernel3_sum = np.nansum(np.squeeze(np.stack(rkernel3_list, axis=0)), axis=0)
        rcnt_sum = np.nansum(np.squeeze(np.stack(rcnt_list, axis=0)), axis=0)


    # plt.figure()
    # plt.contourf(kernel_sum/cnt_sum, cmap= 'RdBu', vmin=-3, vmax=3)
    # plt.plot(21,21,'bo')
    # plt.colorbar()

    print('Returning')

    return (kernel_sum, kernel2_sum, kernel3_sum, cnt_sum, rkernel_sum, rkernel2_sum, rkernel3_sum, rcnt_sum )


if __name__ == "__main__":
    composite()
