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
from wavelet import util
import itertools
from scipy.stats import ttest_ind as ttest

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1


def composite():
    pool = multiprocessing.Pool(processes=8)

    nightp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_0-3UTC_centrePoint.nc'
    dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'

    # nightp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_0-3UTC.nc'
    # dayp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_15-18UTC.nc'

    hour = 18

    if hour > 6:
        file = dayp
    else:
        file = nightp

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 7) ]

    msg = msg.sel(lat=slice(10,20), lon=slice(-10, 10))

    res = pool.map(file_loop, msg)
    pool.close()
    # cnt_list = []
    # res_list = []
    # for m in msg[0:2]:
    #     res, cnt = file_loop(m)
    #     res_list.append(res)
    #     cnt_list.append(cnt)

    res = [x for x in res if x is not None]

    res_list = []
    res2_list = []
    res3_list = []
    cnt_list = []

    rres_list = []
    rres2_list = []
    rres3_list = []
    rcnt_list = []

    for r in res:
        res_list.append(np.squeeze(r[0]))
        res2_list.append(np.squeeze(r[1]))
        #res3_list.append(np.squeeze(r[2]))
        #cnt_list.append(np.squeeze(r[3]))

        rres_list.append(np.squeeze(r[2]))
        rres2_list.append(np.squeeze(r[3]))
        #rres3_list.append(np.squeeze(r[6]))
        #rcnt_list.append(np.squeeze(r[7]))

        scales = r[4]


    ns = np.squeeze(np.vstack(res_list))
    ew = np.squeeze(np.vstack(res2_list))
    # cns_sum = np.nansum(np.squeeze(np.stack(res3_list, axis=0)), axis=0)
    # cew_sum = np.nansum(np.squeeze(np.stack(cnt_list, axis=0)), axis=0)

    rns = np.squeeze(np.vstack(rres_list))
    rew = np.squeeze(np.vstack(rres2_list))
    # rcns_sum = np.nansum(np.squeeze(np.stack(rres3_list, axis=0)), axis=0)
    # rcew_sum = np.nansum(np.squeeze(np.stack(rcnt_list, axis=0)), axis=0)

    nsstat, nspvalue = ttest(ns, rns, axis=0, equal_var=False, nan_policy='omit')
    ewstat, ewpvalue = ttest(ew, rew, axis=0, equal_var=False, nan_policy='omit')

    nsmask = nspvalue < 0.05
    ewmask = ewpvalue < 0.05

    ns_mean = np.nanmean(ns, axis=0)
    ew_mean = np.nanmean(ew, axis=0)

    rns_mean = np.nanmean(rns, axis=0)
    rew_mean = np.nanmean(rew, axis=0)

    f = plt.figure(figsize=(9, 5))
    ax = f.add_subplot(121) #- (rns_sum / rcns_sum)
   # pdb.set_trace()
    plt.contourf((np.arange(0, 101) - 50) * 3, scales ,(ns_mean) - (rns_mean)  , cmap='RdBu_r', vmax=0.05, vmin=-0.05)
    plt.colorbar(label='K')
    #plt.contourf((np.arange(0, 101) - 50) * 3, scales, nsmask, colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    #  plt.plot(50, 50, 'bo')
    #ax.set_xticklabels((np.linspace(0, 100, 6) - 50) * 3)

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(ns.shape[0]) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(122) #-(rew_sum / rcew_sum)

    plt.contourf((np.arange(0, 101) - 50) * 3,scales,   (ew_mean)- (rew_mean)  , cmap='RdBu_r', vmax=0.05, vmin=-0.05)
    plt.colorbar(label='K')
   # plt.contourf((np.arange(0, 101) - 50) * 3, scales, ewmask, colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)
   # ax.set_xticklabels((np.linspace(0, 100, 6) - 50) * 3)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)


    plt.tight_layout()




def cut_kernel(xpos, ypos, wl):

    dist = 50
    #pdb.set_trace()
    try:
        kernel = wl[:, ypos - dist: ypos + dist + 1,  xpos - dist: xpos + dist + 1]
       # pdb.set_trace()
    except ValueError:
        # print('Point at boundary, pass')
        return

    if kernel.shape != (np.shape(wl)[0], 101, 101):
        return

    wav_ns = kernel[:,:,50]
    wav_we = kernel[:,50,:]

    cnt = kernel.copy()
    cnt[np.isfinite(cnt)] = 1
    cnt_ns = cnt[:,:,50]
    cnt_we = cnt[:,50,:]

    return wav_ns, wav_we, cnt_ns, cnt_we



def file_loop(fi):



    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))
    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 6:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date #+ dayd

    lsta = xr.open_dataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/' \
                           'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    print('Doing '+ 'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    lsta_day = lsta['LSTA']
    lsta_day2 = lsta_day.copy()
    lsta_day2 = lsta_day2.where(lsta_day2 > -900)

    # plt.figure()
    # plt.imshow(lsta_day2[0,:,:])
    # #
    # return

    pos = np.where((fi.values >= 1) & (fi.values <= 20))
    if np.sum(pos) == 0:
        print('No blobs found')
        return

    wav_input = lsta_day2-lsta_day2.mean()
    #pdb.set_trace()
    wav = util.waveletLSTA(np.squeeze(wav_input.values), 3,method='dry')
    wl = wav['power']



    xfi = fi.shape[1]
    yfi = fi.shape[0]
    randx = np.random.randint(0, xfi, 250)
    randy = np.random.randint(0, yfi, 250)
    posr = (randy, randx)

##############################Random loop
    rns_list = []
    rwe_list = []
    rcns_list = []
    rcwe_list = []
    for y, x in zip(posr[0], posr[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        point = lsta_day.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_day['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_day['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            wav_ns, wav_we, cnt_ns, cnt_we = cut_kernel(xpos, ypos, wl)
        except TypeError:
            continue

        rns_list.append(wav_ns)  # north-south wavelet
        rwe_list.append(wav_we)  # west-east wavelet
        rcns_list.append(cnt_ns)  # north_south valid count
        rcwe_list.append(cnt_we)  # west_east valid count

###############################Blob loop
    ns_list = []
    we_list = []
    cns_list = []
    cwe_list = []

    for y, x in zip(pos[0], pos[1]):


        lat = fi['lat'][y]
        lon = fi['lon'][x]

        point = lsta_day.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_day['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_day['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            wav_ns, wav_we, cnt_ns, cnt_we = cut_kernel(xpos, ypos, wl)
        except TypeError:
            continue
        #print(plat,plon)
        ns_list.append(wav_ns)  #north-south wavelet
        we_list.append(wav_we) # west-east wavelet
        cns_list.append(cnt_ns)  # north_south valid count
        cwe_list.append(cnt_we) # west_east valid count

    # f = plt.figure()
    # f.add_subplot(131)
    # plt.contourf(kernel_list[4][0,], cmap='RdBu_r', vmin=-8, vmax=8)
    # plt.title('Local')
    # plt.colorbar()
    # f.add_subplot(132)
    # plt.contourf(kernel2_list[4][0,], cmap='RdBu_r', vmin=-8, vmax=8)
    # plt.title('Seasonal')
    # plt.colorbar()
    # f.add_subplot(133)
    # plt.contourf(kernel3_list[4][0,], cmap='RdBu_r', vmin=-8, vmax=8)
    # plt.title('Regional (kernel)')
    # plt.colorbar()
    #
    # return

    if ns_list == []:
        return None
    print(len(ns_list))
    if len(ns_list) == 1:
      return None
    else:
        ns_sum = np.squeeze(np.stack(ns_list, axis=0))
        we_sum = np.squeeze(np.stack(we_list, axis=0))

        # cns_sum = np.nansum(np.squeeze(np.stack(cns_list, axis=0)), axis=0)
        # cwe_sum = np.nansum(np.squeeze(np.stack(cwe_list, axis=0)), axis=0)

        rns_sum = np.squeeze(np.stack(rns_list, axis=0))
        rwe_sum = np.squeeze(np.stack(rwe_list, axis=0))
        # rcns_sum = np.nansum(np.squeeze(np.stack(rcns_list, axis=0)), axis=0)
        # rcwe_sum = np.nansum(np.squeeze(np.stack(rcwe_list, axis=0)), axis=0)

    # plt.figure()
    # plt.contourf(kernel_sum/cnt_sum, cmap= 'RdBu', vmin=-3, vmax=3)
    # plt.plot(21,21,'bo')
    # plt.colorbar()
    scales = wav['scales']

    print('Returning')

    return (ns_sum, we_sum,  rns_sum, rwe_sum,  scales)
# rcns_sum, rcwe_sum, cns_sum, cwe_sum,


if __name__ == "__main__":
    composite()
