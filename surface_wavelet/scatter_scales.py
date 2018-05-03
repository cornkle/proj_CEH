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
from utils import u_plot
from utils import constants
import pickle as pkl
import statsmodels.stats.proportion as prop

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def plot():
    dic = pkl.load(
        open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/dominant_scales_save/scatter_scales.p", "rb"))

    cmap = u_plot.discrete_cmap(24, base_cmap='gnuplot2')

    f = plt.figure()
    ax = plt.subplot(111)

    chourly = []
    shourly = []
    sstd = []
    cstd = []
    hh = []
    for h in range(0,24,1):
        c = np.mean(dic['cell'][(dic['hour']==h) & (dic['cell']>=15)  & (dic['surface'] >= 9) & (dic['surface'] < 130)]) #& (dic['surface'] < 130)
        s = np.mean(np.abs(dic['surface'][(dic['hour']==h) & (dic['cell']>=15)& (dic['surface'] >= 9) & (dic['surface'] < 130)]))
        sc = np.std(dic['cell'][(dic['hour']==h) & (dic['cell']>=15)  & (dic['surface'] >= 9) & (dic['surface'] < 130)])
        ss = np.std(np.abs(dic['surface'][(dic['hour']==h) & (dic['cell']>=15)& (dic['surface'] >= 9) & (dic['surface'] < 130)]))

        pdb.set_trace()
        cstd.append(sc)
        sstd.append(ss)

        chourly.append(c)
        shourly.append(s)
        hh.append(h)

    plt.scatter(chourly, shourly, c=hh, cmap=cmap)
   # plt.errorbar(chourly, shourly, cstd, sstd)
    plt.colorbar()

    chourly = []
    hh = []
    for h in range(0,24,1):
        c = np.mean(dic['cell'][(dic['hour']==h) & (dic['cell']>=15)])
        chourly.append(c)
        hh.append(h)
    f = plt.figure()
    ax = plt.subplot(111)
    plt.scatter(hh, chourly, c=hh, cmap=cmap)

    chourly = []
    shourly = []

    for h in np.unique(dic['surface'])[7:-1]:
        c = np.mean((dic['cell'])[(dic['surface'] == h) ])
        s = h

        chourly.append(c)
        shourly.append(s)
        hh.append(h)

    f = plt.figure()
    ax = plt.subplot(111)
    plt.scatter(chourly, shourly)
    #plt.colorbar()


def composite():
    pool = multiprocessing.Pool(processes=4)

    file = constants.MCS_POINTS_DOM

    msg = xr.open_dataarray(file)
    msg = msg[  (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6)] #(msg['time.hour'] >= 17) &

    msg = msg.sel(lat=slice(10.2, 17), lon=slice(-9.5, 9.5))

    res = pool.map(file_loop, msg)
    pool.close()

    # for m in msg[0:10]:
    #     file_loop(m)
    #
    # return

    res = [x for x in res if x is not None]

    cell = []
    surface = []
    hour = []

    for r in res:
        cell.append(r[0])
        surface.append(r[1])
        hour.append(r[2])
    pdb.set_trace()
    cell = [item for sublist in cell for item in sublist]  # flatten list of lists
    surface = [item for sublist in surface for item in sublist]  # flatten list of lists
    hour = [item for sublist in hour for item in sublist]  # flatten list of lists

    cell = np.array(cell, dtype=float)
    cell = cell[np.isfinite(surface)]

    surface = np.array(surface, dtype=float)
    surface = surface[np.isfinite(surface)]

    hour = np.array(hour, dtype=float)
    hour = hour[np.isfinite(surface)]

    dic = {'cell': cell,
           'surface': surface,
           'hour' : hour }

    pkl.dump(dic,
             open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/dominant_scales_save/scatter_scales.p", "wb"))
    print('Successfully written scatter_scales save file')


def cut_kernel(xpos, ypos, lsta_day2):


    kernel = lsta_day2.isel(lon=xpos+4,
                             lat=ypos).values


    return kernel



def file_loop(fi):
    print('Starting day: ', fi)

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))
    dayd = pd.Timedelta('1 days')

    if (fi['time.hour'].values) <= 16:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date

    try:
        lsta = xr.open_dataset(constants.LSTA_DOMSCALE + \
                               'lsta_daily_powerMaxscale_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
            daybefore.day).zfill(2) + '.nc')
    except OSError:
        print('File not found')
        return

    print('Doing ' + 'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    # lsta = lsta.where(np.abs(temp) > 0.2)
    # temp_da = temp_da.where(np.abs(temp) > 0.2)
    topo = xr.open_dataset(constants.LSTA_TOPO)
    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    lsta_da = lsta['LSTA'].squeeze()
    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.40:
        print('Not enough valid')
        return None


    lsta_da.values[ttopo.values >= 450] = np.nan
    lsta_da.values[gradsum > 30] = np.nan

    scale_list = lsta_da.values.flatten()

    # plt.figure()
    # plt.imshow(lsta_day2[0,:,:])
    # #
    # return

    pos = np.where( (fi.values >= 5))#(fi.values >= 1) & (fi.values <= 20)) #<-50)#

    if np.sum(pos) == 0:
        print('No blobs found')
        return

    cell_list = []
    surface_list = []
    hours_list = []

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        t = fi.sel(lat=lat, lon=lon)

        point = lsta.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            ano = cut_kernel(xpos, ypos, lsta_da)
        except TypeError:
            continue

        if ~np.isfinite(ano):
            continue


        cell_list.append(fi.values[y,x])
        surface_list.append(ano)
        hours_list.append(fi['time.hour'].values)

    if surface_list == None:
        return None

    if surface_list == []:
        return None
    print(len(surface_list))
    try:
        if len(surface_list) == 0:
            return None
    except TypeError:
        return None

    surface_list = np.array(surface_list).flatten()
    cell_list = np.array(cell_list).flatten()
    hours_list = np.array(hours_list).flatten()

    print('Returning')

    return (cell_list, surface_list, hours_list)

if __name__ == "__main__":
    composite()
