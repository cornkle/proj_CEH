# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import pdb
import pandas as pd
from utils import u_met, u_parallelise, u_gis, u_arrays, constants as cnst, u_darrays
import ipdb
import pickle as pkl
import itertools


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def diurnal_loop():

    for h in [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5]:  #range(0,24)

        composite(h)


def composite(h):
    #pool = multiprocessing.Pool(processes=8)

    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    file = cnst.MCS_POINTS_DOM

    hour = h

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour ) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) & (msg['time.month'] <= 9)  ]

    msg = msg.sel( lat=slice(10.2,19.3), lon=slice(-9.8, 9.8))

    msg.attrs['refhour'] = h

    dic = u_parallelise.run_flat(3,file_loop,msg,['c30', 's100', 'e100', 'r30', 'rs100', 're100'])

    # res = []
    # for m in msg[15:20]:
    #     out = file_loop(m)
    #     res.append(out)
    # return

    print('Writing pickle')

    pkl.dump(dic, open(path + "/LSTA_histograms_"+str(hour).zfill(2)+"_corrected_SouthBox.p", "wb"))



def cut_kernel(xpos, ypos, arr, dist):

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    kmean = kernel - np.nanmean(kernel)

    if kernel.shape != (2*dist+1, 2*dist+1):
        return
    ycirc30, xcirc30 = u_arrays.draw_circle(dist, dist,5) # 15km radius
    k30 = np.nanmean(kmean[ycirc30, xcirc30])

    if not np.isfinite(np.nansum(k30)):
        return
    #kernel[ycirc30, xcirc30] = 700

    ycirc100, xcirc100 = u_arrays.draw_circle(dist+1, dist-67, 17)  # at - 200km, draw 50km radius circle
    #s100 = np.nanmean(kmean[ycirc100,xcirc100])
    s100 = np.nanmean(kmean[dist-67-17:dist-67+17, dist-50:dist]) #at -200km in box 100km high

    #kernel[ycirc100,xcirc100] = 1000

    ycirc100e, xcirc100e = u_arrays.draw_circle(dist+51, dist+1, 17)  # at - 150km, draw 50km radius circle
    e100 = np.nanmean(kmean[ycirc100e,xcirc100e])
    # kernel[ycirc100e, xcirc100e] = 500
    #
    # f = plt.figure()
    # plt.imshow(kernel)
    #
    # return

    return k30, s100, e100



def file_loop(fi):

    print('Doing day: ', fi.time)

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))

    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 13:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)

    try:
        lsta = xr.open_dataset(cnst.LSTA_NEW + 'lsta_daily_' + fdate + '.nc')
    except OSError:
        return None
    print('Doing ' + 'lsta_daily_' + fdate + '.nc')

    topo = xr.open_dataset(cnst.LSTA_TOPO)
    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    lsta_da = lsta['LSTA'].squeeze()

    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return None


    lsta_da.values[ttopo.values>=450] = np.nan
    lsta_da.values[gradsum>30] = np.nan
    pos = np.where((fi.values >= 5) & (fi.values <= 65) )  #(fi.values >= 5) & (fi.values < 65), fi.values<-40

    if (np.sum(pos) == 0):
        print('No blobs found')
        return None

    c30 = []
    r30 = []
    s100  = []
    rs100 = []
    e100 = []
    re100 = []

    dist = 100

    mcs_hour = xr.open_dataarray(cnst.MCS_HOUR_DAILY)
    mcsimage = xr.open_dataarray(cnst.MCS_15K)
    mcsimage = mcsimage.sel(time=fi.time, lat=slice(10.2,19.3), lon=slice(-9.8, 9.8))

    labels, goodinds = u_arrays.blob_define(mcsimage.values, -50, minmax_area=[600,100000], max_area=None)

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        if (labels[y,x] not in goodinds) | (labels[y,x] == 0):
            print('MCS too small!!')
            continue

        if (mcsimage.values[y,x] > -60):
            print('Core too warm!!')
            continue

        mhour = mcs_hour.sel(time=pd.datetime(int(fdate[0:4]), int(fdate[4:6]), int(fdate[6:8]), 14), lon=lon,
                             lat=lat).values
        if mhour <= 13:
            mhour += 24

        if mhour == 0:
            mhour += 24

        chour = fi['time.hour'].values

        if (chour >= 0) & (chour <= 13):
            chour += 24
        if (mhour < chour) | (np.isnan(mhour)):
            print('Core overlaps: earliest:', mhour, ' core: ', chour)
            continue

        point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            kc30, kcs100, kce100 = cut_kernel(xpos, ypos, lsta_da.values, dist)
        except TypeError:
            continue

        # if np.sum(np.isfinite(kc30)) ==0:
        #     continue

        c30.append(kc30)
        s100.append(kcs100)
        e100.append(kce100)

        ##### random

        rdist = 50
        randy50 = [y - rdist, y - rdist, y - rdist, y, y, y + rdist, y + rdist, y + rdist]
        randx50 = [x - rdist, x, x + rdist, x - rdist, x + rdist, x - rdist, x, x + rdist]
        randy50_100 = [y - rdist, y - rdist, y, y, y + rdist, y + rdist]

        rdist = 100
        randx100 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]

        rdist = 150
        randx150 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]

        randy = np.array(randy50 + randy50_100 + randy50_100)
        randx = np.array(randx50 + randx100 + randx150)
        # ipdb.set_trace()
        # yrand = np.array([y-50,y,y+50])
        # xrand = np.arange(40,fi.shape[1],50)

        #for ry, rx in itertools.product(yrand,xrand):
        for ry, rx in zip(randy,randx):

            if ry < 0:
                continue
            if ry > fi.shape[0] - 1:
                continue

            if rx < 0:
                continue
            if rx > fi.shape[1] - 1:
                continue

            try:
                lat = fi['lat'][ry]
            except IndexError:
                ipdb.set_trace()
            try:
                lon = fi['lon'][rx]
            except IndexError:
                ipdb.set_trace()

            point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
            plat = point['lat'].values
            plon = point['lon'].values

            xpos = np.where(lsta_da['lon'].values == plon)
            xpos = int(xpos[0])
            ypos = np.where(lsta_da['lat'].values == plat)
            ypos = int(ypos[0])

            try:
                rc30, rcs100, rce100 = cut_kernel(xpos, ypos, lsta_da.values, dist)
            except TypeError:
                continue

            r30.append(rc30)
            rs100.append(rcs100)
            re100.append(rce100)


    if c30 == []:
        return None

    # if len(c30) == 1:
    #   return None

    print('Returning with kernel')

    return (c30, s100, e100, r30, rs100, re100)



def plot(h):
    hour=h
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path + "/composite_new_"+str(hour).zfill(2)+".p", "rb"))

    extent = dic['ano'].shape[1]/2-1

    f = plt.figure(figsize=(14, 7))
    ax = f.add_subplot(231)

    plt.contourf(dic['regional'] / dic['cnt'], cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both')
    plt.plot(extent, extent, 'bo')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly, Nb cores: ' + str(np.max(dic['cnt'])) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)

    ax = f.add_subplot(232)

    plt.contourf((dic['regional'] / dic['cnt'])  , cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly', fontsize=10)

    ax = f.add_subplot(233)

    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal anomaly',
              fontsize=10)

    ax = f.add_subplot(234)

    plt.contourf((dic['ano'] / dic['cnt']) , cmap='RdBu_r',  levels=[-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both')
    plt.plot(extent,extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal anomaly - random',
              fontsize=10)

    ax = f.add_subplot(235)

    plt.contourf(dic['cnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Valid count',
              fontsize=10)

    ax = f.add_subplot(236)

    plt.contourf(dic['cnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent,extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Random valid count',
              fontsize=10)

    plt.tight_layout()
    # plt.savefig(cnst.network_data + "/figs/LSTA-bullshit/AGU/" + +str(hour).zfill(2)+'_allplots.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()

def plot_gewex(h):
    hour=h
    path = cnst.network_data + "figs/LSTA/corrected_LSTA/new/LSTA_only_plots/LSTA_old_vs_new_MCSfilter"
    dic = pkl.load(open(path + "/composite_old_"+str(hour).zfill(2)+".p", "rb"))

    extent = (dic['ano'].shape[1]-1)/2

    f = plt.figure(figsize=(7, 5))
    ax = f.add_subplot(111)
    print(dic['ano'].shape)
    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu_r', extend='both',levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]) #-(rkernel2_sum / rcnt_sum)  levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')

    plt.title(str(hour).zfill(2)+'00 UTC | '+str(np.max(dic['cnt']))+' cores', fontsize=17)


    plt.tight_layout()

    plt.savefig(path +'/lsta_old_' + str(hour).zfill(2)+'_single.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_all():

    hours = [14,16,18,20,22,0,2,4,6]#[15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]

    for h in hours:
        plot_gewex(h)
