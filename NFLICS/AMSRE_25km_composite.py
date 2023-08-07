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
from utils import u_met, u_parallelise, u_gis, u_arrays, constants as cnst
import ipdb
import pickle as pkl
import glob


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def diurnal_loop():

    for h in [14,16,18,20,22,0,2,4,6]:  #range(0,24)

        if h <= 15:
            shift = (12 + h) * (-1)
        else:
            shift = 12 - h

        composite(h, shift)


def composite(h, eh):
    #pool = multiprocessing.Pool(processes=8)

    path = cnst.network_data + 'figs/NFLICS/LSTA_stats_study/'

    for y in range(2012,2020):
        files = glob.glob(cnst.elements_drive + '/Africa/WestAfrica/cores_bigDomain/*_'+str(y)+'_*.nc')
        msg = xr.open_mfdataset(files)
        hour = h

        # file = cnst.MCS_POINTS_DOM
        # msg = xr.open_dataarray(file)
        # msg = msg[(msg['time.hour'] == h) & (msg['time.minute'] == 0) & (
        #         msg['time.year'] == y) & (msg['time.month'] >= 6)]

        msg = msg.sel(lat=slice(10, 20), lon=slice(-12, 15))

        msg = (msg['dom'])[(msg['time.hour'] == hour ) & (msg['time.minute'] == 0) & (msg['time.month'] >= 6) & (msg['time.month'] <= 9)  ]

        msg = msg.load()

        msg.attrs['refhour'] = h
        msg.attrs['eh'] = eh
        #
        dic = u_parallelise.run_arrays(4,file_loop,msg,['ano', 'regional', 'cnt'])
        #
        # res = []
        # for m in msg[0:50]:
        #     out = file_loop(m)
        #     res.append(out)
        # ipdb.set_trace()
        # return

        for k in dic.keys():
           dic[k] = np.nansum(dic[k], axis=0)
        print('File written')
        pkl.dump(dic, open(path + "/composite_new_AMSR2_25km_2012-2020_TS_"+str(y)+'_h'+str(hour).zfill(2)+".p", "wb"))



def cut_kernel(xpos, ypos, arr, date, lon, lat, t, parallax=False, rotate=False):

    if parallax:
        km, coords = u_gis.call_parallax_era(date.month, t, lon, lat, 0, 0)
        lx, ly = km

        lx = int(np.round(lx / 27.5))
        ly = int(np.round(ly / 27.5))  # km into pixels
        xpos = xpos - lx
        ypos = ypos - ly

    dist = 11
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
    # alt_path = cnst.AMSRE_ANO_DAY_CORR
    # lsta = xr.open_dataset(alt_path + 'sma_' + fdate + '.nc')

    alt_path = '/media/ck/Elements/global/AMSR2/daily/25km/day_anom/'
    try:
        lsta = xr.open_dataset(alt_path + 'amsr2_25km_anom_' + fdate + '.nc')
    except:
        print('Surface file not found, return')
        return None

    lsta = lsta.sel(time=str(daybefore.year)+'-'+str(daybefore.month)+'-'+str(daybefore.day))
    lsta = lsta.sel(lon=slice(-14, 14), lat=slice(4, 23))
    #ipdb.set_trace()
    print('Doing '+ 'AMSR_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    #lsta_da = lsta['SM'].squeeze()
    lsta_da = lsta['ts'].squeeze()  # soil_moisture_c1

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
    pos = np.where(fi.values<=-5) #np.where((fi.values >= 5) & (fi.values < 65))

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
        try:
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.2)
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
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.2)
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
    del lsta
    return (kernel2_sum, kernel3_sum, cnt_sum,  rkernel2_sum, rkernel3_sum, rcnt_sum)




def plot(h):
    hour=h
    #pin = cnst.network_data + 'figs/NFLICS/LSTA_stats_study/' + "/composite_new_AMSR10km_"
    pin = cnst.network_data + 'figs/NFLICS/LSTA_stats_study/' + "/composite_new_AMSR2_25km_2012-2020_"
    y1 = 2012
    y2 = 2018

    dic = pkl.load(open(
            pin+str(y1)+'_h'+str(h).zfill(2)+".p", "rb"))

    def coll(dic, h, year):
        print(h)
        core = pkl.load(open(
            pin+str(year)+'_h'+str(h).zfill(2)+".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]

    for y in range(y1+1, y2):
        coll(dic, h, y)

    extent = dic['ano'].shape[1]/2-1

    f = plt.figure(figsize=(14, 4))
    ax = f.add_subplot(131)

    plt.contourf(dic['regional'] / dic['cnt'], cmap='RdBu',  levels=np.linspace(-1.2,1.2,10), extend='both')
    plt.plot(extent, extent, 'bo')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='%')
    plt.title('Regional anomaly, Nb cores: ' + str(np.max(dic['cnt'])) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)


    ax = f.add_subplot(132)

    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu',  levels=np.linspace(-2.2,2.2,10), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='%')
    plt.title('Seasonal anomaly',
              fontsize=10)


    ax = f.add_subplot(133)

    plt.contourf(dic['cnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='n')
    plt.title('Valid count',
              fontsize=10)


    plt.tight_layout()
    # plt.savefig(cnst.network_data + "/figs/LSTA-bullshit/AGU/" + +str(hour).zfill(2)+'_allplots.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    # plt.close()


def plot_ts(h):
    hour=h
    #pin = cnst.network_data + 'figs/NFLICS/LSTA_stats_study/' + "/composite_new_AMSR10km_"
    pin = cnst.network_data + 'figs/NFLICS/LSTA_stats_study/' + "/composite_new_AMSR2_25km_2012-2020_TS_"
    y1 = 2012
    y2 = 2018

    dic = pkl.load(open(
            pin+str(y1)+'_h'+str(h).zfill(2)+".p", "rb"))

    def coll(dic, h, year):
        print(h)
        core = pkl.load(open(
            pin+str(year)+'_h'+str(h).zfill(2)+".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]

    for y in range(y1+1, y2):
        coll(dic, h, y)

    extent = dic['ano'].shape[1]/2-1

    f = plt.figure(figsize=(14, 4))
    ax = f.add_subplot(131)

    plt.contourf(dic['regional'] / dic['cnt'], cmap='RdBu_r',  levels=np.linspace(-1.2,1.2,10), extend='both')
    plt.plot(extent, extent, 'bo')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Regional anomaly, Nb cores: ' + str(np.max(dic['cnt'])) + '| ' + str(hour).zfill(2) + '00UTC, Jun-Sep',
              fontsize=10)


    ax = f.add_subplot(132)

    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu_r',  levels=np.linspace(-2.2,2.2,10), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 5) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='K')
    plt.title('Seasonal anomaly',
              fontsize=10)


    ax = f.add_subplot(133)

    plt.contourf(dic['cnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    dictest = pkl.load(
        open('/home/ck/DIR/cornkle/data/GLOBAL_MCS/save_files/WAf_mcs_tracks_extc_20010101_20011231.p', "rb"))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='n')
    plt.title('Valid count',
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
