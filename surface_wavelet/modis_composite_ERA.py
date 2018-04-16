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


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def loop():

    for l in np.arange(0,24):
        print('Doing '+str(l))
        composite(l)

def composite(h):
    pool = multiprocessing.Pool(processes=8)

    file = '/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'

    hour = h

    msg = xr.open_dataarray(file)

    msg = msg[(msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) ]


    msg = msg.sel(lat=slice(10,18), lon=slice(-10,10))

    #
    res = pool.map(file_loop, (msg))
    pool.close()

    res_list = []
    cnt_list = []
    # pdb.set_trace()
    # for m in msg:
    #     res, cnt = file_loop(m)
    #     res_list.append(res)
    #     cnt_list.append(cnt)

    res = [x for x in res if x is not None]

    shear_ano = []
    shear_kernel = []


    for r in res:
        shear_ano.append(np.squeeze(r[0]))
        shear_kernel.append(np.squeeze(r[1]))

    shear_ano_mean = np.nanmean(np.squeeze(np.stack(shear_ano, axis=0)), axis=0)
    shear_kernel_mean = np.nanmean(np.squeeze(np.stack(shear_kernel, axis=0)), axis=0)



    f = plt.figure(figsize=(9, 7))
    ax = f.add_subplot(121)
    xdist = 60
    plt.pcolormesh(shear_ano_mean, cmap='RdBu',)
    #plt.plot(xdist, xdist, 'bo')

    # ax.set_xticklabels(np.array((np.linspace(0, 121, 5) - 100) * 3, dtype=int))
    # ax.set_yticklabels(np.array((np.linspace(0, 121, 9) - 100) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('ANO',
              fontsize=10)

    ax = f.add_subplot(122)

    plt.pcolormesh(shear_kernel_mean , cmap='viridis' )
    #plt.plot(xdist, xdist, 'bo')

    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('MEAN', fontsize=10)



    plt.tight_layout()

    plt.show()



def cut_kernel(xpos, ypos, lsta_day2, rotate=False):

    dist = 5

    kernel = np.zeros([1, 11, 11]) *np.nan

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

    ano = kernel - np.nanmean(kernel)

    return ano, kernel



def file_loop(fi):
    print('Doing day: ', fi)

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))
    dayd = pd.Timedelta('1 days')

    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 16:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date

    # era = xr.open_dataset('/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCpl.nc')
    # eday = era.sel(time=str(daybefore.year)+'-'+str(daybefore.month)+'-'+str(daybefore.day))
    #
    # u600 = eday['u'].sel(level=600).values
    # u925 = eday['u'].sel(level=925).values
    #
    # shear = u600-u925


    eras = xr.open_dataset('/localscratch/wllf030/cornkle/ERA-I/daily_2006-2010_12UTCsrfc.nc')
    esday = eras.sel(time=str(daybefore.year) + '-' + str(daybefore.month) + '-' + str(daybefore.day))

    u = esday['u10'].values
    v = esday['v10'].values

    div = esday['p84.162'].values

    pos = np.where((fi.values >= 5) & (fi.values < 65))

    if (np.sum(pos) == 0) | (len(pos[0]) < 3):
        print('No blobs found')
        return None

    shear_ano = []
    shear_kernel = []

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        point = esday.sel(latitude=lat, longitude=lon, method='nearest')
        plat = point['latitude'].values
        plon = point['longitude'].values

        xpos = np.where(esday['longitude'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(esday['latitude'].values == plat)
        ypos = int(ypos[0])
        try:
           ano, kernel = cut_kernel(xpos, ypos, u)
        except TypeError:
            continue

        shear_ano.append(ano)
        shear_kernel.append(kernel)



    shear_ano_mean = np.nanmean(np.squeeze(np.stack(shear_ano, axis=0)), axis=0)
    shear_kernel_mean = np.nanmean(np.squeeze(np.stack(shear_kernel, axis=0)), axis=0)


    print('Returning')

    return (shear_ano_mean, shear_kernel_mean)

if __name__ == "__main__":
    composite()
