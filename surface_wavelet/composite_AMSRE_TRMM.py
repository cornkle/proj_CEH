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
from utils import u_met, u_parallelise, u_gis, u_arrays, constants
import pickle as pkl
import salem

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def loop():

    for l in np.arange(0,24):
        print('Doing '+str(l))
        composite(l)

def composite(h):

    file = constants.TRMM5KM_FILE

    hour = h
    # 16UTC - 18UTC
    # 0-3 UTC
    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] >=0) & (msg['time.hour'] <=3) &(
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6)]

    msg = msg.sel(lat=slice(10, 20), lon=slice(-10, 10))

    dic = u_parallelise.run_arrays(7, file_loop, msg, ['ano', 'regional', 'cnt', 'rano', 'rregional', 'rcnt'])

    for k in dic.keys():
       dic[k] = np.nansum(dic[k], axis=0)

    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/amsre_0-"+str(hour).zfill(2)+".p", "wb"))

    extent = dic['ano'].shape[1]/2

    f = plt.figure(figsize=(14, 7))
    ax = f.add_subplot(231)

    plt.contourf(dic['regional'] / dic['cnt'], cmap='RdBu',  vmin=-0.8, vmax=0.8)
    plt.plot(extent, extent, 'bo')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly, Nb cores: ' + str(np.max(dic['cnt'])) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)

    ax = f.add_subplot(232)

    plt.contourf((dic['regional'] / dic['cnt']) - (dic['rregional'] / dic['rcnt']) , cmap='RdBu',  vmin=-0.8, vmax=0.8)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly', fontsize=10)

    ax = f.add_subplot(233)

    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu',  vmin=-2, vmax=2) #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal anomaly',
              fontsize=10)

    ax = f.add_subplot(234)

    plt.contourf((dic['ano'] / dic['cnt']) - (dic['rano'] / dic['rcnt']), cmap='RdBu',  vmin=-2, vmax=2)
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

    plt.contourf(dic['rcnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent,extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Random valid count',
              fontsize=10)

    plt.tight_layout()

    # plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/scales/new/composites_lsta/test/TRMM_AMSR'+str(hour).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()


def cut_kernel(xpos, ypos, arr, date, lon, lat, t, parallax=False, rotate=False):

    if parallax:
        km, coords = u_gis.call_parallax_era(date.month, t, lon, lat, 0, 0)
        lx, ly = km

        lx = int(np.round(lx / 27.5))
        ly = int(np.round(ly / 27.5))  # km into pixels
        xpos = xpos - lx
        ypos = ypos - ly
    #AMSRE 0.25 degrees ~ 27.5 km
    dist = 10

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    if rotate:
        kernel = u_met.era_wind_rotate(kernel,date.month,lat,lon,level=700, ref_angle=90)

    # if (np.sum(np.isfinite(kernel)) < 0.10 * kernel.size):
    #     return

    kernel3 = kernel - np.nanmean(kernel)

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1

    if kernel.shape != (dist*2+1, dist*2+1):
        return None

    return kernel, kernel3, cnt


def file_loop(fi):

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
        daybefore = date

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)
    try:
        lsta = xr.open_dataset(constants.AMSRE_ANO_DAY + 'sma_' + fdate + '.nc')
    except OSError:
        return None
    lsta = lsta.sel(time=str(daybefore.year)+'-'+str(daybefore.month)+'-'+str(daybefore.day))
    lsta = lsta.sel(lon=slice(-11, 11), lat=slice(9, 21))
    print('Doing '+ 'AMSR_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    lsta_da = lsta['SM'].squeeze()

    topo = xr.open_dataset(constants.LSTA_TOPO)
    ttopo = topo['h']
    ttopo = lsta_da.salem.lookup_transform(ttopo)

    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    # if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.50:
    #     print('Not enough valid')
    #     return None

    # lsta_da.values[np.isnan(lsta_da.values)] = 0

    lsta_da.values[ttopo.values >= 450] = np.nan
    lsta_da.values[gradsum > 30] = np.nan

    pos = np.where(fi.values >= 15)

    if (np.sum(pos) == 0) | (len(pos[0]) < 2):
        print('No blobs found')
        return None

    kernel2_list = []
    kernel3_list = []
    cnt_list = []

    xfi = fi.shape[1]

    randx = np.random.randint(0, xfi, 100)
    try:
        randy = np.random.randint(np.min(pos[0]), np.max(pos[0]), 100)
    except ValueError:
        return None
    posr = (randy, randx)

    rkernel2_list = []
    rkernel3_list = []
    rcnt_list = []

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
            rkernel2, rkernel3, rcnt = cut_kernel(xpos, ypos, lsta_da, daybefore, plon, plat, -40, parallax=False,
                                                  rotate=False)
        except TypeError:
            continue

        rkernel2_list.append(rkernel2)
        rkernel3_list.append(rkernel3)
        rcnt_list.append(rcnt)

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        t = fi.sel(lat=lat, lon=lon)

        point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            kernel2, kernel3, cnt = cut_kernel(xpos, ypos, lsta_da, daybefore, plon, plat, -40, parallax=False,
                                               rotate=False)
        except TypeError:
            continue

        kernel2_list.append(kernel2)
        kernel3_list.append(kernel3)
        cnt_list.append(cnt)

    if kernel2_list == []:
        return None

    if len(kernel2_list) == 1:
        return None
    else:

        kernel2_sum = np.nansum(np.stack(kernel2_list, axis=0), axis=0)
        kernel3_sum = np.nansum(np.stack(kernel3_list, axis=0), axis=0)
        cnt_sum = np.nansum(np.stack(cnt_list, axis=0), axis=0)

        rkernel2_sum = np.nansum((np.stack(rkernel2_list, axis=0)), axis=0)
        rkernel3_sum = np.nansum((np.stack(rkernel3_list, axis=0)), axis=0)
        rcnt_sum = np.nansum((np.stack(rcnt_list, axis=0)), axis=0)

    print('Returning')

    return (kernel2_sum, kernel3_sum, cnt_sum,  rkernel2_sum, rkernel3_sum, rcnt_sum)


def plot_gewex():

    dic1 = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/amsre_0-17.p", "rb"))
    dic2 = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/amsre_0-03.p", "rb"))

    extent = (dic1['ano'].shape[1]-1)/2

    f = plt.figure(figsize=(5, 7))
    ax = f.add_subplot(211)

    plt.contourf(dic1['regional'] / dic1['cnt'], cmap='RdBu',  levels=[-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent/2) * 60, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 30, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='Volumetric soil moisture anomaly (%)')
    plt.title('16-1800UTC | '+str(np.max(dic1['cnt']))+' cores', fontsize=12)

    ax = f.add_subplot(212)

    plt.contourf(dic2['regional'] / dic2['cnt'], cmap='RdBu',  levels=[-1,-0.75,-0.5,-0.25,0.25,0.5,0.75,1], extend='both')
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent * 2, 9) - extent/2) * 60, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent * 2, 9) - extent) * 30, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='Volumetric soil moisture anomaly (%)')
    plt.title('00-0300UTC | ' + str(np.max(dic2['cnt'])) + ' cores', fontsize=12)


    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/GEWEX/AMSRE_TRMM.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()