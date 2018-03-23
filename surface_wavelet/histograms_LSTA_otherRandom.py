# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import pickle as pkl
matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
import pdb
import multiprocessing


def cut_kernel(zpos, ypos, xpos, da):

    pixshift = 40
    xpos = xpos-pixshift

    mdist = 0

    try:
        lpoint = da.isel(time=zpos, lon=slice(xpos - mdist, xpos + mdist + 1),
                               lat=slice(ypos - mdist, ypos + mdist + 1)).values
    except (IndexError, ValueError):
        return

    if lpoint == np.nan:
        return


    try:
        rpoint = da.isel(time=zpos, lon=slice(0, -1),
                               lat=slice(ypos - mdist, ypos + mdist + 1)).values
    except (IndexError, ValueError):
        return

    return lpoint, rpoint.flatten()


def composite(hour):
    #hour = 18
    mds = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/lsta_daily_*.nc', autoclose=True)
    mds = mds.sel(lat=slice(10,19), lon=slice(-10,10))
    #mds = mds.isel(time=slice(0,10))
    pos = np.where((mds['cell'].values==hour) )
    shape = mds['LSTA'].shape

    print('Doing with hour ', hour)

    all = []
    point = []

    for z, y, x in zip(pos[0], pos[1], pos[2]):

        print('Doing pos', z, y, x)

        try:
            lpoint, rpoint = cut_kernel(z, y, x, mds['LSTA'])
        except TypeError:
            continue

        point.extend(lpoint)
        all.extend(rpoint)

    pdb.set_trace()
    all = np.array(all, dtype=float)
    point = np.concatenate(point)

    all = all[np.isfinite(all)]
    point = point[np.isfinite(point)]


    dic = { 'all' : all,
            'point' : point,

    }

    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/histo_other"+str(hour).zfill(2)+".p", "wb"))

    #
    # plt.figure()
    # nball, bins, v = plt.hist(all, bins=np.arange(-10,10,1), normed=True, edgecolor='k', color=None, alpha=0.3)
    # nbpoint, bins, v = plt.hist(point, bins=np.arange(-10,10,1), normed=True,edgecolor='k', color=None, alpha=0.3)


    #return nball, nbpoint, bins

def plot():

    dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/histo_other17.p", "rb"))

    all = dic['all']

    point = dic['point']


    f = plt.figure()
    ax = f.add_subplot(121)

    nball, bins,v = plt.hist(all, bins=np.arange(-10,10,1), normed=True, edgecolor='k', color=None, alpha=0.3)
    nbpoint, bins, v = plt.hist(point, bins=np.arange(-10,10,1), normed=True, edgecolor='k', color=None, alpha=0.3)
    plt.xlabel('Max. LSTA[9x9km box] | 120km East of cell / random point in same LAT ')
    plt.ylabel('Frequency')
    stri = (np.sum(point >= np.percentile(all, 90)) / point.size * 100).round(2)
    plt.title('18-19UTC:'+str(stri)+'% of Cells occur in warmest quartile')

    # dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/histo0.p", "rb"))
    #
    # all = dic['all']
    # allgrad = dic['allgrad']
    # point = dic['point']
    # pointgrad = dic['pointgrad']
    #
    # perc = np.percentile(all, np.arange(0, 101, 10))
    # percg = np.percentile(allgrad, np.arange(0, 101, 10))
    # ax = f.add_subplot(122)
    #
    # nball, bins, v = plt.hist(all, bins=np.arange(-10, 10, 1), normed=True, edgecolor='k', color=None, alpha=0.3)
    # nbpoint, bins, v = plt.hist(point, bins=np.arange(-10, 10, 1), normed=True, edgecolor='k', color=None, alpha=0.3)
    # plt.xlabel('Max. LSTA[9x9km box] | 120km East of cell / random point in same LAT ')
    # plt.ylabel('Frequency')
    # stri = (np.sum(point >= np.percentile(all, 90)) / point.size * 100).round(2)
    # plt.title('0-4UTC:' + str(stri) + '% of Cells occur in warmest quartile')
    #
    #
    # plt.figure()
    # ngall, gbins, v = plt.hist(allgrad, bins=np.arange(-10,10,1), normed=True, edgecolor='k', color=None, alpha=0.3)
    # ngpoint, gbins, v = plt.hist(pointgrad, bins=np.arange(-10,10,1), normed=True, edgecolor='k', color=None, alpha=0.3)
    p = 90
    prob = np.sum(point > np.percentile(all, p)) / point.size
    print(prob)
    print( (prob - (100-p)*0.01) / ((100-p)*0.01)) # percentage of cells in warmest 25% of LSTA

    p = 10
    prob = np.sum(point < np.percentile(all, p)) / point.size
    print(prob)
    print( (prob - p*0.01) / (p*0.01)) # percentage of cells in warmest 25% of LSTA

    # ((np.sum(point >= np.percentile(all, 75)) / point.size) -0.25) / 0.25



    return nball, nbpoint, bins
        