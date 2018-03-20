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


    file = '/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant_smallstorm.nc'

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

    msg = msg.sel(lat=slice(10.5,18), lon=slice(-9.5,9.5))
    #
    # res = pool.map(file_loop, msg)
    # pool.close()

    res_list = []
    cnt_list = []
    pdb.set_trace()
    for m in msg[0:2]:
        res, cnt = file_loop(m)
        res_list.append(res)
        cnt_list.append(cnt)

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

    if res_list  == []:
        return

    kernel_sum = np.nansum(np.squeeze(np.stack(res_list, axis=0)), axis=0)
    kernel2_sum = np.nansum(np.squeeze(np.stack(res2_list, axis=0)), axis=0)
    kernel3_sum = np.nansum(np.squeeze(np.stack(res3_list, axis=0)), axis=0)
    cnt_sum = np.nansum(np.squeeze(np.stack(cnt_list, axis=0)), axis=0)
    rkernel_sum = np.nansum(np.squeeze(np.stack(rres_list, axis=0)), axis=0)
    rkernel2_sum = np.nansum(np.squeeze(np.stack(rres2_list, axis=0)), axis=0)
    rkernel3_sum = np.nansum(np.squeeze(np.stack(rres3_list, axis=0)), axis=0)
    rcnt_sum = np.nansum(np.squeeze(np.stack(rcnt_list, axis=0)), axis=0)

    if kernel_sum.ndim != 2:
        return

    f = plt.figure(figsize=(9, 7))
    ax = f.add_subplot(221)

    plt.contourf((kernel_sum / cnt_sum)-(rkernel_sum / rcnt_sum), cmap='RdBu_r',  levels=[-0.4,-0.3,-0.2,-0.1, -0.05, 0, 0.05, 0.1, 0.2,0.3,0.4], extend='both')
    plt.plot(100, 100, 'bo')

    ax.set_xticklabels(np.array((np.linspace(0, 200, 5) - 100) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, 200, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Local Anomaly, Nb cores: ' + str(np.max(cnt_sum)) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf((kernel3_sum / cnt_sum) - (rkernel3_sum / rcnt_sum) , cmap='RdBu_r',  levels=[-0.4,-0.3,-0.2,-0.1, -0.05, 0, 0.05,0.1, 0.2,0.3,0.4], extend='both' )
    plt.plot(100, 100, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, 200, 5) - 100) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, 200, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional Anomaly', fontsize=10)

    ax = f.add_subplot(223)

    plt.contourf((kernel2_sum / cnt_sum), cmap='RdBu_r',  levels=[-0.4,-0.3,-0.2,-0.1, -0.05, 0, 0.05, 0.1, 0.2,0.3,0.4], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(100, 100, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, 200, 5) - 100) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, 200, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal Anomaly',
              fontsize=10)

    ax = f.add_subplot(224)

    plt.contourf((kernel2_sum / cnt_sum)-(rkernel2_sum / rcnt_sum), cmap='RdBu_r',  levels=[-0.4,-0.3,-0.2,-0.1, -0.05, 0, 0.05,0.1, 0.2,0.3,0.4], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(100,100, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, 200, 5) - 100) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, 200, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal Anomaly - Random',
              fontsize=10)

    plt.tight_layout()

    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/scales/new/composites_lsta/'+str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()

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


def cut_kernel(xpos, ypos, lsta_day2, month, lon, lat, t, parallax=None):
    mdist = 5

    if parallax:
        height = era_geop.era_Tlapse(month, t, lon, lat)  # height in meters
        km, coords = u_gis.parallax_corr_msg(0, 0, lon, lat, height / 1000)
        lx, ly = km
        print(t, height, lx, ly)
        lx = int(np.round(lx / 3.))
        ly = int(np.round(ly / 3.))  # km into pixels
        xpos = xpos - lx
        ypos = ypos - ly

    mini_mean = lsta_day2.isel(lon=slice(xpos - mdist, xpos + mdist + 1),
                               lat=slice(ypos - mdist, ypos + mdist + 1)).values

    mini_mmean = np.nanmean(mini_mean)

    if np.isnan(mini_mmean):
        # pdb.set_trace()
        return

    dist = 100

    # try:
    kernel = lsta_day2.isel(lon=slice(xpos - dist, xpos + dist + 1),
                            lat=slice(ypos - dist, ypos + dist + 1)).values

    # kernel2 = lsta_day2.isel(lon=slice(xpos - dist, xpos + dist + 1),
    #                         lat=slice(ypos - dist, ypos + dist + 1)).values

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

        cutk = lsta_day2[:, ypos - ymindist: ypos + ymaxdist, xpos - xmindist: xpos + xmaxdist]


        kernel[:, ymin: ymax, xmin:xmax] = cutk

    ano = kernel - mini_mmean
    kernel3 = kernel - np.nanmean(kernel)

    cnt = kernel.copy() * 0
    cnt[np.isfinite(kernel)] = 1

    return ano, kernel, kernel3, cnt



def file_loop(fi):
    print('Doing day: ', fi)

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))
    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 15:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date #+ dayd
    try:
        lsta = xr.open_dataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/' \
                           'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2) + '.nc')
    except OSError:
        return
    print('Doing '+ 'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    lsta_day = lsta['LSTA']
    lsta_day2 = lsta_day.copy()
    lsta_day2 = lsta_day2.where(lsta_day2 > -900)

    # plt.figure()
    # plt.imshow(lsta_day2[0,:,:])
    # #
    # return

    pos = np.where( (fi.values >= 5) & (fi.values < 60))#(fi.values >= 1) & (fi.values <= 20)) #<-50)#
    pdb.set_trace()
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
            ano, kernel2, kernel3, cnt = cut_kernel(xpos, ypos, lsta_day2, daybefore.month, plon, plat, -40, parallax=False)
        except TypeError:
            continue

        rkernel_list.append(ano)
        rkernel2_list.append(kernel2)
        rkernel3_list.append(kernel3)
        rcnt_list.append(cnt)

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        t = fi.sel(lat=lat, lon=lon)

        point = lsta_day2.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_day['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_day['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            ano, kernel2, kernel3, cnt = cut_kernel(xpos, ypos, lsta_day2, daybefore.month, plon, plat, -40, parallax=False)
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

    print('Returning')

    return (kernel_sum, kernel2_sum, kernel3_sum, cnt_sum, rkernel_sum, rkernel2_sum, rkernel3_sum, rcnt_sum )

if __name__ == "__main__":
    composite()
