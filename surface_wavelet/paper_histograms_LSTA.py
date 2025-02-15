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

    pkl.dump(dic, open(path + "/LSTA_histograms_"+str(hour).zfill(2)+"_SlotFilter.p", "wb"))



def cut_kernel(xpos, ypos, arr, dist, nbslot=False):

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    if nbslot is not False:
        nbsl = u_arrays.cut_kernel(nbslot, xpos, ypos, dist)
        if np.sum(nbsl[dist-10:dist+10,dist-10:dist+50]>5) / np.sum(np.isfinite(nbsl[dist-10:dist+10,dist-10:dist+50])) <=0.5:
            print('TOO FEW SLOTS!')
            return

    kmean = kernel - np.nanmean(kernel)

    if kernel.shape != (2*dist+1, 2*dist+1):
        return
    ycirc30, xcirc30 = u_arrays.draw_circle(dist+1, dist+1,6) # 15km radius
    k30 = np.nanmean(kmean[ycirc30, xcirc30])

    if not np.isfinite(np.nansum(k30)):
        return
    #kernel[ycirc30, xcirc30] = 700

    ycirc100, xcirc100 = u_arrays.draw_circle(dist+1, dist-67, 17)  # at - 200km, draw 50km radius circle
    #s100 = np.nanmean(kmean[ycirc100,xcirc100])
    s100 = np.nanmean(kmean[dist-67-17:dist-67+17, dist-50:dist]) #at -200km in box 100km high

    #kernel[ycirc100,xcirc100] = 1000

    ycirc100e, xcirc100e = u_arrays.draw_circle(dist+31, dist+1, 17)  # at - 100/150km, draw 50km radius circle
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
    slot_da = lsta['NbSlot'].squeeze().values

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
        try:
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.04)
        except KeyError:
            continue
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            kc30, kcs100, kce100 = cut_kernel(xpos, ypos, lsta_da.values, dist, nbslot=slot_da)
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
            try:
                point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.04)
            except KeyError:
                continue
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