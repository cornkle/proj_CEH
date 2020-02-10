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
import pandas as pd
from wavelet import util as wutil
from utils import u_arrays, constants as cnst, u_met
from scipy.stats import ttest_ind as ttest
from scipy.interpolate import griddata
import pickle as pkl
from utils import u_arrays as ua
import os
import collections
import warnings
import salem

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1

def run_hours():

    l = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7] #15,16,
    for ll in l:
        composite(ll)

def composite(hour):

    key = '2hOverlap'

    msgopen = pd.read_csv(cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/cores_gt15000km2_table_AMSRE_LSTA_tracking_bigWin_' + key + '_'+str(hour)+'.csv', na_values=-999)
    # msgopen = pd.read_csv(
    #     cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/core_txt/cores_gt15000km2_table_1640_580_' + str(
    #         hour) + '.csv')

    msg = pd.DataFrame.from_dict(msgopen)

    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))


    msg = msg[ (msg['dtime']<=2) & ((msg['SMmean0']<=-4)) & ((msg['SMmean-1']<=-0.5)) ] #  & (np.isfinite(msg['SMmean0'])) & (np.isfinite(msg['SMmean0'])) #& (msg['SMmean0']<-1)
    #msg = msg[(msg['LSTAslotfrac'] >= 0.5) & (msg['dtime'] <= 2) & (np.isfinite(msg['SMmean0']))]
    msgin = msg[(msg['lat']>8.5) & (msg['lat']<20.5) & (msg['topo']<450)]#[msg['initTime']<=3]#[msg['SMwet']==2]

    print('Number of cores', len(msgin))
    #ipdb.set_trace()
    # calculate the chunk size as an integer
    #'chunk_size = int(msg.shape[0] / pnumber)
    msgin.sort_values(by='date')
    #ipdb.set_trace()`-

    msgy = msgin

    chunk, chunk_ind, chunk_count = np.unique(msgy.date, return_index=True, return_counts=True)

    chunks = [msgy.loc[msgy.index[ci:ci + cc]] for ci, cc in zip(chunk_ind, chunk_count)] # daily chunks

    # res = []
    # for m in chunks[::]:
    #     out = file_loop(m)
    #     res.append(out)

    #return
    pool = multiprocessing.Pool(processes=5)

    res = pool.map(file_loop, chunks)
    pool.close()

    print('Returned from parallel')


    res = [x for x in res if x is not None]

    lsta_list = []
    lsta_cnt = []
    amsr_list = []
    amsr_cnt = []
    cores = 0
    for r in res:


        lsta_list.append(r[0])
        lsta_cnt.append(r[1])
        amsr_list.append(r[2])
        amsr_cnt.append(r[3])
        cores += r[4]

    dic = collections.OrderedDict([

                                   ('lsta' , [lsta_list, lsta_cnt]),
                                   ('amsr' , [amsr_list, amsr_cnt]),
                                   ('cores', cores)])

    keys = list(dic.keys())

    for l in keys:

        if l == 'cores':
            continue

        try:
            (dic[l])[0] = np.squeeze(np.vstack((dic[l])[0]))
        except TypeError:
            ipdb.set_trace()
        (dic[l])[1] = np.squeeze(np.vstack((dic[l])[1]))

    for l in keys:

        if l == 'cores':
            continue

        (dic[l])[0] = np.nansum((dic[l])[0], axis=0)
        (dic[l])[1] = np.nansum((dic[l])[1], axis=0)


    outpath = cnst.network_data + '/figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    pkl.dump(dic, open(outpath+"coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMSL_DRYtest_" + key + ".p", "wb"))
    print('Save file written!')




def cut_kernel_lsta(xpos, ypos, arr):

    dist = 100

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    kernel = kernel
    if kernel.shape != (2*dist+1, 2*dist+1):
        print('Kernel wrong shape')
        return

    # plt.figure()
    # plt.imshow(kernel, origin='lower')

    return kernel



def file_loop(df):

    date = df['date'].iloc[0]
    hour = df['hour'].iloc[0]
    print('Doing day: ', date)

    storm_date = date

    dayd = pd.Timedelta('1 days')

    if (hour) <= 13:
        print('Nighttime')
        daybefore = storm_date - dayd
    else:
        print('Daytime')
        daybefore = storm_date

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)

    amsre = xr.open_dataset(cnst.AMSRE_ANO_DAY + 'sma_' + fdate + '.nc')
    amsre = amsre.sel(time=str(daybefore.year)+'-'+str(daybefore.month)+'-'+str(daybefore.day))
    amsre = amsre.sel(lon=slice(-11, 11), lat=slice(8, 21))
    print('Doing '+ 'AMSR_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    amsr_da = amsre['SM'].squeeze()


    try:
        lstal = xr.open_dataset(cnst.LSTA_1330 + 'lsta_daily_' + fdate + '.nc')
    except OSError:
        return None
    print('Doing ' + 'lsta_daily_' + fdate + '.nc')

    lsta_da = lstal['LSTA'].squeeze()
    slot_da = lstal['NbSlot'].squeeze().values

    # if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.01:
    #     print('Not enough valid')
    #     return None

    topo = xr.open_dataset(cnst.LSTA_TOPO)
    #topo = xr.open_dataset(cnst.WA_TOPO_3KM)
    #topo = topo.sel(lat=slice(7,25), lon=slice(-14,14))

    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    # if (np.sum(np.isfinite(amsr_da)) / amsr_da.size) < 0.01:
    #     print('Not enough valid')
    #     return None

    try:
        amsr_da = topo.salem.transform(amsr_da)
    except RuntimeError:
        print('amsr_da on LSTA interpolation problem')
        return None


    # try:
    #     lsta_da = topo.salem.transform(lsta_da)
    # except RuntimeError:
    #     print('lsta_da on LSTA interpolation problem')
    #     return None


    lsta_da.values[ttopo.values >= 450] = np.nan
    lsta_da.values[gradsum > 30] = np.nan

    amsr_da.values[ttopo.values >= 450] = np.nan
    amsr_da.values[gradsum > 30] = np.nan

    del topo

    lsta = []
    amsr = []


    ###############################Blob loop
    cores = 0

    for dids, dit in df.iterrows():

        lat = dit.lat
        lon = dit.lon

        try:
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.04)
        except KeyError:
            print('Nearest point finding error')
            continue

        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            lsta_kernel = cut_kernel_lsta(xpos, ypos, lsta_da.values)
        except TypeError:
            print('LSTA kernel error')
            continue

        try:
            amsre_kernel = cut_kernel_lsta(xpos, ypos, amsr_da.values)
        except TypeError:
            print('AMSR kernel error')
            continue

        lsta.append(lsta_kernel)
        amsr.append(amsre_kernel)
        cores += 1


    del lsta_da
    del amsr_da

    if (np.array(amsr).ndim == 1) | (np.array(lsta).ndim == 1):
        return None


    #try:
    if (len(lsta) <= 1) | (len(amsr) <= 1) | (np.array(amsr).ndim < 3) | (np.nansum(lsta) == 0)| (np.nansum(amsr)==0):
        print('NorthSouth 1')
        return None
    else:

        lsta_sum = np.nansum(np.stack(lsta, axis=0), axis=0)[np.newaxis,...]
        lsta_cnt = np.sum(np.isfinite(np.stack(lsta, axis=0)), axis=0)[np.newaxis,...]
        try:
            amsr_sum = np.nansum(np.stack(amsr, axis=0), axis=0)[np.newaxis,...]
        except ValueError:
            ipdb.set_trace()
        amsr_cnt = np.sum(np.isfinite(np.stack(amsr, axis=0)), axis=0)[np.newaxis,...]
    # except TypeError:
    #     ipdb.set_trace()


    print('Returning with kernel, success!!')

    return (lsta_sum, lsta_cnt, amsr_sum, amsr_cnt, cores)

if __name__ == "__main__":
    run_hours()
