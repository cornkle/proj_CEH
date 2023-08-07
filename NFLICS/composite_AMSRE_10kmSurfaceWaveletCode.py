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
import ipdb
import glob
import pandas as pd
from scipy import ndimage
from utils import u_met, u_parallelise, u_gis, u_arrays, constants as cnst
import pickle as pkl
import salem


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def loop():

    for l in np.arange(0,24):
        print('Doing '+str(l))
        composite(l)

def composite(h):
    file = cnst.MCS_POINTS_DOM  # '/home/ck/DIR/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'
    #and /home/ck/DIR/cornkle/data/OBS/AMSRE/aqua/sma_nc_day_new/'
    hour = h

    # msg = xr.open_dataarray(file)
    # msg = msg[(msg['time.hour'] == h)  & (msg['time.minute'] == 0) & (
    #     msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6)]
    #
    # msg = msg.sel( lat=slice(10,20), lon=slice(-12,15))

    files = glob.glob(cnst.elements_drive + '/Africa/WestAfrica/cores_bigDomain/*.nc')
    msg = xr.open_mfdataset(files)
    hour = h

    # file = cnst.MCS_POINTS_DOM
    # msg = xr.open_dataarray(file)
    # msg = msg[(msg['time.hour'] == h) & (msg['time.minute'] == 0) & (
    #         msg['time.year'] == y) & (msg['time.month'] >= 6)]

    msg = msg.sel(lat=slice(10, 20), lon=slice(-12, 15))

    msg = (msg['dom'])[
        (msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (msg['time.month'] >= 6) &
        (msg['time.month'] <= 9) & (msg['time.year'] >= 2012) & (msg['time.year'] <= 2016)]

    msg = msg.load()

    #ipdb.set_trace()
    dic = u_parallelise.run_arrays(4,file_loop,msg,['ano', 'regional', 'cnt', 'rano', 'rregional', 'rcnt'])

    # res = []
    # #ipdb.set_trace()
    # for m in msg[200:201]:
    #
    #     out = file_loop(m)
    #     res.append(out)
    # return

    for k in dic.keys():
       dic[k] = np.nansum(dic[k], axis=0)

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

    plt.contourf((dic['ano'] / dic['cnt'])- (dic['rano'] / dic['rcnt']), cmap='RdBu',  vmin=-2, vmax=2)
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

    #plt.savefig(cnst.network_data + 'figs/LSTA-bullshit/AGU/'+'AMSRE_' + str(hour).zfill(2)+'.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.show()


def cut_kernel(xpos, ypos, arr, date, lon, lat, t, parallax=False, rotate=False):

    if parallax:
        km, coords = u_gis.call_parallax_era(date.month, t, lon, lat, 0, 0)
        lx, ly = km

        lx = int(np.round(lx / 27.5))
        ly = int(np.round(ly / 27.5))  # km into pixels
        xpos = xpos - lx
        ypos = ypos - ly

    dist = 30
    #dist=100

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

    if (fi['time.hour'].values) <= 16:
        print('Nighttime')
        daybefore = date #- dayd
    else:
        print('Daytime')
        daybefore = date

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)
    #alt_path = cnst.AMSRE_ANO_DAY
    #alt_path = cnst.AMSRE_ANO_DAY_CORR
    try:
        alt_path = glob.glob(cnst.elements_drive + 'global/AMSR2/daily/10km/day_anom/amsr2_*_'+fdate+'.nc')[0]
    except:
        return None
    lsta = xr.open_dataset(alt_path)
    lsta = lsta.sel(time=str(daybefore.year)+'-'+str(daybefore.month)+'-'+str(daybefore.day))
    lsta = lsta.sel(lon=slice(-14, 14), lat=slice(4, 23))
    #ipdb.set_trace()
    print('Doing '+ 'AMSR_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    lsta_da = lsta['ts'].squeeze()

    # topo = xr.open_dataset(cnst.LSTA_TOPO)
    # ttopo = topo['h']
    #ttopo = lsta_da.salem.lookup_transform(ttopo)

    # try:
    #     lsta_da = topo.salem.transform(lsta_da)
    # except RuntimeError:
    #     print('lsta_da on LSTA interpolation problem')
    #     return None

    # grad = np.gradient(ttopo.values)
    # gradsum = abs(grad[0]) + abs(grad[1])

    # if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.50:
    #     print('Not enough valid')
    #     return None

    # lsta_da.values[np.isnan(lsta_da.values)] = 0

    # lsta_da.values[ttopo.values >= 450] = np.nan
    # lsta_da.values[gradsum > 30] = np.nan
    pos = np.where(fi.values<=-5)#np.where((fi.values >= 5) & (fi.values < 65))

    if (np.sum(pos) == 0) | (len(pos[0]) < 3):
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

    # mcs_hour = xr.open_dataarray(cnst.MCS_HOUR_DAILY)
    # mcsimage = xr.open_dataarray(cnst.MCS_15K)
    # mcsimage = mcsimage.sel(time=fi.time, lat=slice(10.2,19), lon=slice(-9.9,9.9))
    # counter = 0
    #
    # labels, goodinds = u_arrays.blob_define(mcsimage.values, -50, minmax_area=[600, 50000], max_area=None)
    #ipdb.set_trace()

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        # if (labels[y,x] not in goodinds) | (labels[y,x] == 0):
        #     print('MCS too small!!')
        #     continue
        #
        # if (mcsimage.values[y,x] > -60):
        #     print('Core too warm!!')
        #     continue
        #
        # mhour = mcs_hour.sel(time=pd.datetime(int(fdate[0:4]), int(fdate[4:6]), int(fdate[6:8]), 14), lon=lon,
        #                      lat=lat).values
        # if mhour <= 13:
        #     mhour += 24
        #
        # chour = fi['time.hour'].values
        #
        # if (chour >= 0) & (chour <= 13):
        #     chour += 24
        # if (mhour+2 < chour) | (np.isnan(mhour)):
        #     print('Core overlaps: earliest:', mhour, ' core: ', chour)
        #     continue

        t = fi.sel(lat=lat, lon=lon)
        try:
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.1)
        except KeyError:
            print('point index error')
            continue
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

        # plt.figure()
        # plt.pcolormesh(kernel2_sum, vmin=-10,vmax=10, cmap='RdBu')
        #
        # plt.figure()
        # plt.pcolormesh(cnt_sum, cmap='viridis')

    print('Returning')

    return (kernel2_sum, kernel3_sum, cnt_sum,  rkernel2_sum, rkernel3_sum, rcnt_sum)
