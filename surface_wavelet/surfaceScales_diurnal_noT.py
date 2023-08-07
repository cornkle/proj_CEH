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
from utils import constants as cnst, u_plot as uplot
import numpy.ma as ma
import pickle as pkl
import statsmodels.stats.proportion as prop

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def loop():
    bins = np.arange(-200,201,20)

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

    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/dominant_scales_save/scales_dist.p", "wb"))
    print('Successfully written save file')



def plot():


    dic = pkl.load( open(cnst.network_data + "figs/LSTA-bullshit/scales/new/dominant_scales_save/scales.p", "rb"))

    bin = np.array(dic['bin'])
    bin = bin[3:-3]
    center = bin[0:-1] + (bin[1::]-bin[0:-1])/2

    data = (dic['blob']- dic['scale'])/(dic['scale'])#dic['scale']#(np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc']))## (np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc'])) #dic['scale']  #(np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc']))
    data = data[:,3:-3]*100
    # data[0:3,-5::]= np.abs(data[0:3,-5::])
    # data[3:7, -9:-4] = np.abs(data[3:7, -9:-4])
    f = plt.figure()
    plt.imshow(data, origin='lower')

    db = dic['blobc'][:,3:-3]
    filler = np.zeros_like(db)
    for i in range(db.shape[0]):
        for j in range(db.shape[1]):
            low , up = prop.proportion_confint(db[i,j], np.nansum(db[i,:]))

            unrange = (db[i,j] / np.nansum(db[i,:])) - low
            filler[i,j] = unrange

    mask = np.zeros_like(db)
    #mask[filler>(np.abs(dic['blob'][:,3:-3]- dic['scale'][:,3:-3])+0.004)] = 1
    #
    # mask = mask[:,-3:3]
    # data[np.where(mask)] = 0

    #data[np.where(data > 0)] += 1
    data[np.where(np.abs(data)<8)]=0
    levels=[-40,-30,-20,-10,-5,5,10,20,30,40]
    f = plt.figure(figsize=(7,5))
    ax = plt.subplot(111)

    pmap = ax.contourf(center, np.arange(0,24), data, levels=[-25, -20, -15, -10,-5,5,10,15, 20, 25], cmap='PuOr_r', extend='both')
    plt.axvline(0, c='k', linestyle='--')
    ax.set_xticks(center, minor=False)
    ax.set_xticklabels(np.array(center, dtype=int)) # center
    # cbar = plt.colorbar(pmap)
    # cbar.set_label('Difference in scale frequency | Blobs')

    ax.set_yticks(np.arange(0,24), minor=False)
    ax.set_yticklabels(np.arange(0,24))
    ax.set_xlabel('Scales of pos/neg. amplitude (km)')
    ax.set_ylabel('Time of day')


    ax1 = ax.twinx()
    ax1.set_yticks(np.arange(dic['blob'].shape[0]), minor=False)
    ax1.set_yticklabels(dic['nblobs'])

    plt.tight_layout()
    f.subplots_adjust(right=0.81)
    cax = f.add_axes([0.89, 0.15, 0.02, 0.75])
    cbar = f.colorbar(pmap, cax)
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label('Difference in surface-scale frequency (%)', fontsize=11)

    plt.show()

    print(np.sum(dic['blobc']>0)/np.sum(dic['nblobs']))

    print(np.sum(np.isfinite(dic['blobc'])))

    print(np.sum(data, axis=0))

def plot_colormesh():


    dic = pkl.load( open(cnst.network_data + "figs/LSTA-bullshit/scales/new/dominant_scales_save/scales.p", "rb"))

    bin = np.array(dic['bin'])
    bin = bin[3:-3]
    center = bin[0:-1] + (bin[1::]-bin[0:-1])/2

    data = (dic['blob']- dic['scale'])/(dic['scale'])#dic['scale']#(np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc']))## (np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc'])) #dic['scale']  #(np.sum(dic['blobc'], axis=0)/np.sum(dic['blobc']))
    data = data[:,3:-3]*100
    # data[0:3,-5::]= np.abs(data[0:3,-5::])
    # data[3:7, -9:-4] = np.abs(data[3:7, -9:-4])
    f = plt.figure()
    plt.imshow(data, origin='lower')

    db = dic['blobc']
    filler = np.zeros_like(db)
    for i in range(db.shape[0]):
        for j in range(db.shape[1]):
            low , up = prop.proportion_confint(db[i,j], np.nansum(db[i,:]))

            unrange = (db[i,j] / np.nansum(db[i,:])) - low
            filler[i,j] = unrange

    #mask = np.zeros_like(db)
    #mask[filler[:,3:-3]>(np.abs(dic['blob'][:,3:-3]- dic['scale'][:,3:-3])+0.004)] = 1

    #mask = mask[:,-3:3]
    #data[np.where(mask)] = 0

    #data[np.where(data > 0)] += 1
    data[np.where(np.abs(data*100)<8)]=0
    levels=[-40,-30,-20,-10,-5,5,10,20,30,40]
    f = plt.figure(figsize=(7,5))
    ax = plt.subplot(111)

    pmap = ax.contourf(center, np.arange(0,24), data, levels=[-30, -20, -15, -10,-5,5,10,15, 20, 30], cmap='PuOr_r', extend='both')
    plt.axvline(0, c='k', linestyle='--')
    ax.set_xticks(center, minor=False)
    ax.set_xticklabels(np.array(center, dtype=int)) # center

    X, Y = np.meshgrid(center,np.arange(0,24))
    cmapp = uplot.discrete_cmap(9, base_cmap='RdBu_r')

    Zm = ma.masked_where(np.isnan(data),data)

    mappable = ax.pcolormesh(X, Y, Zm, cmap=cmapp, vmin=-35, vmax=35) # viridis_rvmin=20, vmax=40
    ax.set_ylabel('Time of day')
    ax.set_xlabel('Scale of pos/neg. amplitude (km)')

    cbar = f.colorbar(mappable, ticks=np.linspace(-35,35,10)) # ticks=np.linspace(30,45,11) , ticks=np.linspace(20,40,11)
    cbar.set_label('Difference in surface-scale frequency (%)')


    # cbar = plt.colorbar(pmap)
    # cbar.set_label('Difference in scale frequency | Blobs')

    # ax.set_yticks(np.arange(0,24), minor=False)
    # ax.set_yticklabels(np.arange(0,24))
    # ax.set_xlabel('Scales of pos/neg. amplitude (km)')
    # ax.set_ylabel('Time of day')


    ax1 = ax.twinx()
    ax1.set_yticks(np.arange(dic['blob'].shape[0]), minor=False)
    ax1.set_yticklabels(dic['nblobs'])

    plt.tight_layout()
    f.subplots_adjust(right=0.81)
    cax = f.add_axes([0.89, 0.15, 0.02, 0.75])
    # cbar = f.colorbar(pmap, cax)
    # cbar.ax.tick_params(labelsize=11)
    # cbar.set_label('Difference in surface-scale frequency (%)', fontsize=11)

    plt.show()



def composite(h):
    pool = multiprocessing.Pool(processes=6)


    file = constants.MCS_POINTS_DOM
    hour = h


    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == h) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) ]

    msg = msg.sel(lat=slice(10.5,18), lon=slice(-9.5,9.5))


    res = pool.map(file_loop, msg)
    pool.close()

    # #
    # bla = []
    # for m in msg[0:50]:
    #     res = file_loop(m)
    #     bla.append(res)
    # pdb.set_trace()
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

    histb, hb = np.histogram(blobs, bins=np.arange(-200,201,20), weights = weight_blobs)
    hists, hs = np.histogram(scales, bins=np.arange(-200,201,20), weights=weight_scales)

    histbc, hb = np.histogram(blobs, bins=np.arange(-200,201,20))
    histsc, hs = np.histogram(scales, bins=np.arange(-200,201,20))

    print('Number of blobs:', blobs.size)

    # f = plt.figure()
    # plt.bar(hb[0:-1], histb, align='edge', width=hb[1::]-hb[0:-1],edgecolor='k')
    #
    # f = plt.figure()
    # plt.bar(hs[0:-1], hists, align='edge', width=hb[1::]-hb[0:-1],edgecolor='k')

    # f = plt.figure()
    # plt.bar(hs[0:-1], histb-hists, align='edge', width=hb[1::]-hb[0:-1],edgecolor='k')

    return histb,hists,hb, blobs.size, histbc, histsc


def cut_kernel(xpos, ypos, lsta_day2):


    # kernel = lsta_day2.isel(lon=slice(xpos+3 , xpos + 4 + 1),
    #                         lat=slice(ypos-1 , ypos + 1 + 1)).values
    # #
    # #
    # kernel =  np.mean(kernel)

    kernel = lsta_day2.isel(lon=xpos+10,
                             lat=ypos).values

    # if (kernel == np.nan):
    #     return

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

    # lsta_da.values[np.isnan(lsta_da.values)] = 0

    lsta_da.values[ttopo.values >= 400] = np.nan
    lsta_da.values[gradsum > 30] = np.nan

    #scale_list = lsta_da.values.flatten()

    # plt.figure()
    # plt.imshow(lsta_day2[0,:,:])
    # #
    # return

    pos = np.where( (fi.values >= 5) & (fi.values < 60))#(fi.values >= 1) & (fi.values <= 20)) #<-50)#

    if np.sum(pos) == 0:
        print('No blobs found')
        return

    kernel_list = []
    sign_list = []
    scale_list = []

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

        xfi = fi.shape[1]
        randx = np.random.randint(0, xfi, 100)
        if np.min(pos[0]) == 0:
            ymin = np.min(pos[0])
        else:
            ymin = np.min(pos[0]) - 1
        if np.max(pos[0]) == fi.shape[0] - 1:
            ymax = np.max(pos[0])
        else:
            ymax = np.max(pos[0]) + 1
        randy = np.random.randint(ymin, ymax, 100)
        posr = (randy, randx)

        for y, x in zip(posr[0], posr[1]):

            lat = fi['lat'][y]
            lon = fi['lon'][x]

            point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
            plat = point['lat'].values
            plon = point['lon'].values

            xpos = np.where(lsta_da['lon'].values == plon)
            xpos = int(xpos[0])
            ypos = np.where(lsta_da['lat'].values == plat)
            ypos = int(ypos[0])

            try:
                scale = cut_kernel(xpos, ypos, lsta_da)
            except TypeError:
                print('Kernel random error')
                continue

        if ano > 0:
            bla = 1
        else:
            bla = 0

        kernel_list.append(ano)
        sign_list.append(bla)
        scale_list.append(scale)

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
