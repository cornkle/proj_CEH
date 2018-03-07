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
import pickle as pkl
import statsmodels.stats.proportion as prop

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def loop():
    bins = np.arange(-150,151,15)

    arr_blob = np.zeros([24,len(bins)-1])
    arr_scale = np.zeros([24, len(bins) - 1])
    arr_blobc = np.zeros([24,len(bins)-1])
    arr_scalec = np.zeros([24, len(bins) - 1])
    nblob = []

    for l in np.arange(0,24):
        print('Doing '+str(l))
        blobs, scales, bins, nblobs, blobsc, scalesc = composite(l)

        arr_blob[l,:] = blobs
        arr_scale[l, :] = scales
        arr_blobc[l,:] = blobsc
        arr_scalec[l, :] = scalesc
        nblob.append(nblobs)

    dic = { 'scale': arr_scale,
            'blob' : arr_blob,
            'blobc' : arr_blobc,
            'scalec' : arr_scalec,
            'bin': bins,
            'nblobs' : np.array(nblob)}

    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/scales.p", "wb"))





def plot():


    dic = pkl.load( open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/scales.p", "rb"))

    bin = np.array(dic['bin'])
    center = bin[0:-1] + (bin[1::]-bin[0:-1])

    blobc = dic['blobc']
    scalec = dic['scalec']

    pdb.set_trace()


    data = (dic['blob']-dic['scale'])*100

    hours = np.arange(0, 24)

    pos = np.where((hours > 9) & (hours < 13))
    data[pos[0], :] = 0

    f = plt.figure()
    ax = plt.subplot(111)
    pmap = ax.pcolormesh(data, vmin=-1, vmax=1, cmap='RdBu_r')
    ax.set_xticks(np.arange(dic['blob'].shape[1])+1, minor=False)
    ax.set_xticklabels(center)
    cbar = plt.colorbar(pmap)
    cbar.set_label('Difference in scale frequency | Blobs')

    ax.set_yticks(np.arange(dic['blob'].shape[0]) + 1, minor=False)
    ax.set_yticklabels(np.arange(0,24))
    ax.set_xlabel('Surface Scales of pos/neg deviation to surroundings')
    ax.set_ylabel('Hours')


    ax1 = ax.twinx()
    ax1.set_yticks(np.arange(dic['blob'].shape[0])+1, minor=False)
    ax1.set_yticklabels(dic['nblobs'])
    pdb.set_trace()

    print(np.sum(dic['blobc']>0)/np.sum(dic['nblobs']))

    print(np.sum(np.isfinite(dic['blobc'])))

def plot2():


    dic = pkl.load( open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/scales.p", "rb"))

    bin = np.array(dic['bin'])
    center = bin[0:-1] + (bin[1::]-bin[0:-1])

    data = dic['blob']- (np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc']))#dic['scale']#(np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc']))## (np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc'])) #dic['scale']  #(np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc']))
    db = dic['blobc']
    filler = np.zeros_like(db)
    for i in range(db.shape[0]):
        for j in range(db.shape[1]):
            low , up = prop.proportion_confint(db[i,j], np.nansum(db[i,:]))

            unrange = (db[i,j] / np.nansum(db[i,:])) - low
            filler[i,j] = unrange

    mask = np.zeros_like(db)
    mask[filler>np.abs(data)] = 1
    #data[np.where(mask)] = 0

    f = plt.figure()
    ax = plt.subplot(111)
    pmap = ax.pcolormesh(data*100, vmin=-2, vmax=2, cmap='RdBu_r')
    ax.set_xticks(np.arange(dic['blob'].shape[1])+1, minor=False)
    ax.set_xticklabels(center)
    cbar = plt.colorbar(pmap)
    cbar.set_label('Difference in scale frequency | Blobs')

    ax.set_yticks(np.arange(dic['blob'].shape[0]) + 1, minor=False)
    ax.set_yticklabels(np.arange(0,24))
    ax.set_xlabel('Surface Scales of pos/neg deviation to surroundings')
    ax.set_ylabel('Hours')


    ax1 = ax.twinx()
    ax1.set_yticks(np.arange(dic['blob'].shape[0])+1, minor=False)
    ax1.set_yticklabels(dic['nblobs'])

    print(np.sum(dic['blobc']>0)/np.sum(dic['nblobs']))

    print(np.sum(np.isfinite(dic['blobc'])))

    print(np.sum(data, axis=0))


def composite(h):
    pool = multiprocessing.Pool(processes=8)


    file = '/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'
    hour = h


    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == h) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) ]

    msg = msg.sel(lat=slice(10.5,17.5), lon=slice(-9.5,9.5))


    res = pool.map(file_loop, msg)
    pool.close()


    # for m in msg[0:50]:
    #     file_loop(m)
    #
    # return

    res = [x for x in res if x is not None]

    blobs = []
    scales = []
    sign = []
    signt = []


    for r in res:
        blobs.append(r[0])
        scales.append(r[1])
        sign.append(r[2])

    blobs = [item for sublist in blobs for item in sublist]  # flatten list of lists
    scales = [item for sublist in scales for item in sublist]  # flatten list of lists
    sign = [item for sublist in sign for item in sublist]  # flatten list of lists


    blobs = np.array(blobs, dtype=float)

    blobs = blobs[np.isfinite(blobs)]

    scales = np.array(scales, dtype=float)
    scales = scales[np.isfinite(scales)]

    print(np.unique(blobs), len(np.unique(blobs)))

    weight_blobs = np.ones_like(blobs)/float(len(blobs))
    weight_scales = np.ones_like(scales) / float(len(scales))

    histb, hb = np.histogram(blobs, bins=np.arange(-150,151,15), weights = weight_blobs)
    hists, hs = np.histogram(scales, bins=np.arange(-150,151,15), weights=weight_scales)

    histbc, hb = np.histogram(blobs, bins=np.arange(-150,151,15))
    histsc, hs = np.histogram(scales, bins=np.arange(-150,151,15))

    print('Number of blobs:', blobs.size)

    # f = plt.figure()
    # plt.bar(hb[0:-1], histb, align='edge', width=hb[1::]-hb[0:-1],edgecolor='k')
    #
    # f = plt.figure()
    # plt.bar(hs[0:-1], hists, align='edge', width=hb[1::]-hb[0:-1],edgecolor='k')

    # f = plt.figure()
    # plt.bar(hs[0:-1], histb-hists, align='edge', width=hb[1::]-hb[0:-1],edgecolor='k')

    return histb,hists,hb, blobs.size, histbc, histsc

    # bins = np.percentile(scales,np.arange(0,101,5))
    #
    # #bins = [-140,-112,-92,-63,-43,-24,-13,-9, 9, 16,29,43,63,92,112,140]
    #
    # histb, hb = np.histogram(blobs, bins=bins)
    #
    # f = plt.figure()
    # plt.bar(hb[0:-1], histb , align='edge', width=hb[1::] - hb[0:-1], edgecolor='k')
    #
    # histb, hb = np.histogram(scales, bins=bins)
    # f = plt.figure()
    # plt.bar(hb[0:-1], histb, align='edge', width=hb[1::] - hb[0:-1], edgecolor='k')
    #
    # print('Percentage pixels of similar LSTA and wavelet sign:', np.sum(sign)/len(sign))


    #plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/scales/new/composites_lsta/'+str(hour).zfill(2)+'00UTC_lsta_scaleDist_dominant<20.png')
    #plt.close()


def cut_kernel(xpos, ypos, lsta_day2):


    # kernel = lsta_day2.isel(lon=slice(xpos+3 , xpos + 4 + 1),
    #                         lat=slice(ypos-1 , ypos + 1 + 1)).values
    # #
    # #
    # kernel =  np.mean(kernel)

    kernel = lsta_day2.isel(lon=xpos+5,
                             lat=ypos).values

    # if (kernel == np.nan):
    #     return

    return kernel



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
        lsta = xr.open_dataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/scale_maps/' \
                           'lsta_daily_scale_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2) + '.nc')
    except OSError:
        return
    print('Doing '+ 'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    lsta_check = xr.open_dataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/' \
                           'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')
    try:
        lsta_check = lsta_check.sel(lat=slice(lsta['lat'].values.min(),lsta['lat'].values.max()), lon=slice(lsta['lon'].values.min(),lsta['lon'].values.max()))
    except OSError:
        return

    lsta_day = lsta['LSTA']
    lsta = lsta_check['LSTA']

    lsta_day2 = lsta_day.copy()

    lstav = lsta_check['LSTA'].values

    lsta_day2 = lsta_day2.where(np.abs(lstav) > 0.2)
    lsta = lsta.where(np.abs(lstav) > 0.2)

    if np.sum(lsta_day2) == np.nan:
        return

    #lsta_day2.values[np.abs(lsta_day2.values) < 9] = np.nan

    scale_list = (lsta_day2.values).flatten()

    # plt.figure()
    # plt.imshow(lsta_day2[0,:,:])
    # #
    # return

    pos = np.where( (fi.values >= 5) & (fi.values < 76))#(fi.values >= 1) & (fi.values <= 20)) #<-50)#

    if np.sum(pos) == 0:
        print('No blobs found')
        return

    kernel_list = []
    sign_list = []

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
            ano = cut_kernel(xpos, ypos, lsta_day2)
        except TypeError:
            continue

        if ano == np.nan:
            continue

        try:
            anot = cut_kernel(xpos, ypos, lsta)
        except TypeError:
            continue

        if anot > 0: signt = 1
        else: signt = 0
        if ano > 0: sign = 1
        else: sign = 0

        if signt == sign: bla = 1
        else: bla = 0

        if ano == np.nan: bla = np.nan

        kernel_list.append(ano)
        sign_list.append(bla)

    if kernel_list == None:
        return None

    if kernel_list == []:
        return None
    print(len(kernel_list))
    try:
        if len(kernel_list) == 0:
            return None
    except TypeError:
        return None

    print('Returning')

    return (kernel_list, scale_list, sign_list)

if __name__ == "__main__":
    composite()
