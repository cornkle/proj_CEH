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
from statsmodels.stats.proportion import proportion_confint
from utils import constants as cnst

def diurnal_loop():

    pool = multiprocessing.Pool(processes=8)
    h = np.arange(0,24)
    res = pool.map(composite, h)
    pool.close()

def cut_kernel(zpos, ypos, xpos, da):

    pixshift = 40
    xpos = xpos+pixshift

    mdist = 2

    try:
        mini_mean = da.isel(time=zpos, lon=slice(xpos - mdist, xpos + mdist + 1),
                               lat=slice(ypos - mdist, ypos + mdist + 1)).values
    except (IndexError, ValueError):
        return
    try:
        lpoint = np.nanmean(mini_mean)
    except ValueError:
        return

    if lpoint == 0:
        return

    return lpoint


def composite(hour):
    #hour = 18
    mds = xr.open_mfdataset(cnst.network_data + '/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_*.nc', autoclose=True)
    mds = mds.sel(lat=slice(10,17), lon=slice(-10,10))
    #mds = mds.isel(time=slice(150,300))
    pos = np.where((mds['cell'].values==hour) )
    shape = mds['LSTA'].shape

    print('Doing with hour ', hour)

    all = []
    point = []

    for z, y, x in zip(pos[0], pos[1], pos[2]):

        print('Doing pos', z, y, x)

        try:
            lpoint = cut_kernel(z, y, x, mds['LSTA'])
        except TypeError:
            continue

        print(lpoint)
        point.append(lpoint)

        randx = np.random.randint(0, shape[2], 20)
        randy = np.random.randint(-10, 11, 20)

        for ry, rx in zip(randy, randx):

            if y+ry < 0: ry=0
            if y+ry > shape[1]-1: ry = shape[1]-1

            try:
                apoint = cut_kernel(z, y+ry, rx, mds['LSTA'])
            except TypeError:
                continue

            all.append(apoint)

    all = np.array(all, dtype=float)
    point = np.array(point, dtype=float)

    all = all[np.isfinite(all)]
    point = point[np.isfinite(point)]


    dic = { 'all' : all,
            'point' : point}

    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/scales/P_LSTA_givenCore/histo_test"+str(hour).zfill(2)+".p", "wb"))


def plot():

    dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/histo_shift003.p", "rb"))

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


def plot_diurn():

    percmmax = []
    percmmin = []
    nbmax = []
    nbmin = []
    err90_up = []
    err90_low = []
    err10_up = []
    err10_low = []


    rrange = [16,17,18,19,20,21,22,23,0,1,2,3,4,5, 6]

    for h in rrange:

        dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/histo_noshift"+str(h).zfill(2)+".p", "rb"))
        print('Open')
        all = dic['all']

        point = dic['point']

        p = 90
        pprob = np.sum(point > np.percentile(all, p))
        prob = pprob / point.size

        percmax = (prob - (100-p)*0.01) / ((100-p)*0.01) *100 # percentage of cells in warmest 25% of LSTA
        percmmax.append(percmax)
        nbmax.append(point > np.percentile(all, p))

        low90, upp90 = proportion_confint(pprob, point.size)

        err90_up.append( ((upp90 - (100-p)*0.01) / ((100-p)*0.01) *100) - percmax)
        err90_low.append(percmax -((low90 - (100 - p) * 0.01) / ((100 - p) * 0.01) * 100) )

        p = 10
        pprob = np.sum(point < np.percentile(all, p))
        prob = pprob / point.size
        percmin =  (prob - p*0.01) / (p*0.01) *100 # percentage of cells in warmest 25% of LSTA

        percmmin.append(percmin)
        nbmin.append(len(point))
        low10, upp10 = proportion_confint(pprob, point.size)

        err10_up.append((upp10 - p*0.01) / (p*0.01) *100 - percmin)
        err10_low.append( percmin - (low10 - p*0.01) / (p*0.01) *100 )

    f = plt.figure(figsize=(9,5))
    ax = f.add_subplot(111)
    ax.bar(np.arange(0,15), percmmin,  label='10th centile',yerr=np.vstack((err10_up, err10_low)), edgecolor='k') #
    ax.bar(np.arange(0, 15), percmmax, label='90th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
    ax.set_xticks(np.arange(0, 15))
    ax.set_xticklabels(rrange)

    ax.set_xlabel('Hour')
    #plt.axvline(16, color='k')

    #ax.set_xlim(-2,25)
    #plt.ylabel('NbCells|LSTA gt/lt perc. / NbCells')
    plt.ylabel('Difference in probability (%)')
    plt.legend()

    ax1 = ax.twiny()
    ax1.bar(np.arange(0, 15), percmmin, label='10th centile', yerr=np.vstack((err10_up, err10_low)), edgecolor='k')
    ax1.bar(np.arange(0, 15), percmmax, label='90th centile', yerr=np.vstack((err90_up, err90_low)), edgecolor='k')
    ax1.set_xticks(np.arange(0,15))
    ax1.set_xticklabels(nbmin, rotation=45)
    ax1.set_xlabel('Number of convective cores')

    plt.tight_layout()
    plt.annotate('a)', xy=(0.04, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')  # transform=ax.transAxes,
