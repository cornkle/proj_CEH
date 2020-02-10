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

    msgopen = pd.read_csv(cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/cores_gt15000km2_table_AMSRE_LSTA_tracking_' + key + '_'+str(hour)+'.csv', na_values=-999)

    msg = pd.DataFrame.from_dict(msgopen)

    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))

    msg = msg[(msg['LSTAslotfrac']>=0.025) & (msg['dtime']<=2) & (np.isfinite(msg['SMmean0']))] #  & (np.isfinite(msg['SMmean0'])) (msg['SMmean0']<-1)

    msgin = msg[(msg['lat']>8.5) & (msg['lat']<20.5) & (msg['topo']<450)]

    print('Number of cores', len(msgin))
    msgin.sort_values(by='date')


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

    amsrk30 = []
    amsre100 = []

    ramsrk30 = []
    ramsre100 = []

    cores = 0
    for r in res:

        amsrk30.extend(r[0])
        amsre100.extend(r[1])


        ramsrk30.extend(r[2])
        ramsre100.extend(r[3])

        cores += r[4]

    dic = collections.OrderedDict([
                                   ('amsr' , [amsrk30, amsre100]),
                                    ('ramsr', [ramsrk30, ramsre100]),
                                   ('cores', cores)])


    outpath = cnst.network_data + '/figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    pkl.dump(dic, open(outpath + "/AMSR_histograms_" + str(hour).zfill(2) + "_SMFINITE.p", "wb"))
    print('Save file written!')




def cut_kernel_lsta(xpos, ypos, arr):

    dist=200

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    kmean = kernel

    if kernel.shape != (2*dist+1, 2*dist+1):
        return
    ycirc30, xcirc30 = u_arrays.draw_circle(dist+1, dist+1,6) # 15km radius
    k30 = np.nanmean(kmean[ycirc30, xcirc30])

    ycirc100e, xcirc100e = u_arrays.draw_circle(dist+51, dist+1, 17)  # at - 150km, draw 50km radius circle


    e100 = np.nanmean(kmean[ycirc100e,xcirc100e])
    if np.sum(np.isfinite(e100)) / e100.size < 0.05:
        return

    return k30, e100



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
    amsre = amsre.sel(lon=slice(-11.5, 11.5), lat=slice(8, 21))
    print('Doing '+ 'AMSR_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    amsr_da = amsre['SM'].squeeze()


    #topo = xr.open_dataset(cnst.LSTA_TOPO)
    topo = xr.open_dataset(cnst.WA_TOPO_3KM)
    topo = topo.sel(lat=slice(7,25), lon=slice(-14,14))

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

    amsr_da.values[ttopo.values >= 450] = np.nan
    amsr_da.values[gradsum > 30] = np.nan

    del topo

    amsrk30 = []
    amsre100 = []

    ramsrk30 = []
    ramsre100 = []


    ###############################Blob loop
    cores = 0

    for dids, dit in df.iterrows():

        lat = dit.lat
        lon = dit.lon

        try:
            point = amsr_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.04)
        except KeyError:
            print('Nearest point finding error')
            continue

        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(amsr_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(amsr_da['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            ak30, ae100 = cut_kernel_lsta(xpos, ypos, amsr_da.values)
        except TypeError:
            print('AMSR kernel error')
            continue


        amsrk30.append(ak30)
        amsre100.append(ae100)

        cores += 1

        ##### random

        y = ypos
        x = xpos

        rdist = 50
        randy50 = [y - rdist, y - rdist, y - rdist, y, y, y + rdist, y + rdist, y + rdist]
        randx50 = [x - rdist, x, x + rdist, x - rdist, x + rdist, x - rdist, x, x + rdist]
        randy50_100 = [y - rdist, y - rdist, y, y, y + rdist, y + rdist]

        rdist = 100
        randx100 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]

        rdist = 150
        randx150 = [x - rdist, x + rdist, x - rdist, x + rdist, x - rdist, x + rdist]

        randy = np.array(randy50 + randy50_100 + randy50_100)
        randx = np.array(randx50 + randx100 + randx150)

        for ry, rx in zip(randy,randx):

            if ry < 0:
                continue
            if ry > amsr_da.shape[0] - 1:
                continue

            if rx < 0:
                continue
            if rx > amsr_da.shape[1] - 1:
                continue

            try:
                lat = amsr_da['lat'][ry]
            except IndexError:
                ipdb.set_trace()
            try:
                lon = amsr_da['lon'][rx]
            except IndexError:
                ipdb.set_trace()

            point = amsr_da.sel(lat=lat, lon=lon, method='nearest')
            plat = point['lat'].values
            plon = point['lon'].values

            xpos = np.where(amsr_da['lon'].values == plon)
            xpos = int(xpos[0])
            ypos = np.where(amsr_da['lat'].values == plat)
            ypos = int(ypos[0])

            try:
                arc30, arce100 = cut_kernel_lsta(xpos, ypos, amsr_da.values)
            except TypeError:
                continue

            ramsrk30.append(arc30)
            ramsre100.append(arce100)

    del amsr_da

    if (len(amsrk30) == 0) :
        return None


    print('Returning with kernel, success!!')

    return (amsrk30, amsre100,
            ramsrk30, ramsre100, cores)

if __name__ == "__main__":
    run_hours()
