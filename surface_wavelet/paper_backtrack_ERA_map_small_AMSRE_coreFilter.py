# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import ipdb
import pandas as pd
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_arrays as ua, constants as cnst, u_darrays
from scipy.interpolate import griddata
import multiprocessing
import os
import glob

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def diurnal_loop():
    afternoon = list(range(14,24))
    night = list(range(0,8))
    all = afternoon + night

    hlist = []
    for hh in all:
        if hh >= 14:
            hlist.append((hh,12-hh))
        else:
            hlist.append((hh, 12-(hh+24)))

    for l in hlist:
        print('Doing '+str(l))
        composite(l[0], l[1])


def eh_loop():

    all = np.arange(-41,1,3)


    #hlist = []
    # for hh in all:
    #     if hh >= 14:
    #         hlist.append((hh,12-hh))
    #     else:
    #         hlist.append((hh, 12-(hh+24)))

    for l in all:
        print('Doing '+str(l))
        composite(20, l)

def rewrite_list(hour):
    path = '/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/ERA5/core_txt/'
    dic = pkl.load(
        open(path + "cores_gt15000km2_table_AMSRE_" + str(hour) + ".p", "rb"))
    # new = dic.copy()
    # for k in new.keys():
    #     new[k] = []
    #
    # for k in dic.keys():
    #     lists = dic[k]
    #     for l in lists:
    #         new[k].extend(l)
    #
    # pkl.dump(new, open(path + "cores_gt15000km2_table_1640_580_" + str(hour) + "_new.p", "wb"))

    df = pd.DataFrame.from_dict(dic)
    df = df.reindex(columns=['year', 'month', 'day', 'hour', 'lon', 'lat', 'xloc', 'yloc', 'area', 'csize', 't', 'storm_id', 'SMmean', 'SMdry', 'SMwet'])
    df.to_csv(path + "cores_gt15000km2_table_AMSRE_tracking_" + str(hour) + ".csv", na_rep=-999, index_label='id')


def composite(h):

    path = '/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/ERA5/core_txt/'
    msgopen = pd.read_csv(
        '/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/wavelet_coefficients/core_txt/cores_gt15000km2_table_1640_580_'+str(h)+'.csv')
    hour = h
    msg = pd.DataFrame.from_dict(msgopen)# &  &

    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))

    # calculate the chunk size as an integer
    #'chunk_size = int(msg.shape[0] / pnumber)
    msg.sort_values(by='date')
    msg['SMmean'] = np.nan
    msg['SMdry'] = np.nan
    msg['SMwet'] = np.nan

    chunk, chunk_ind, chunk_count = np.unique(msg.date, return_index=True, return_counts=True)


    chunks = [msg.ix[msg.index[ci:ci + cc]] for ci, cc in zip(chunk_ind, chunk_count)] # daily chunks

    # res = []
    # for m in chunks[0:100]:
    #     out = file_loop(m)
    #     res.append(out)
    #
    # ipdb.set_trace()
    # return
    pool = multiprocessing.Pool(processes=4)

    res = pool.map(file_loop, chunks)
    pool.close()

    print('Returned from parallel')
    # ipdb.set_trace()
    res = [x for x in res if x is not None]

    df_concat = pd.concat(res)

    dic = df_concat.to_dict()

    pkl.dump(dic, open(path+"/cores_gt15000km2_table_AMSRE_tracking_"+str(hour)+".p", "wb"))  #"+str(hour)+"
    print('Save file written!')
    print('Dumped file')

    rewrite_list(hour)



def cut_kernel(xpos, ypos, arrlist, dist):

    wetflag = 0
    dryflag = 0
    smean = np.nan

    for ids, arr in enumerate(arrlist):


        kernel = ua.cut_kernel(arr,xpos, ypos,dist)
        kernel = kernel - np.nanmean(kernel)
        if kernel.shape != (dist*2+1, dist*2+1):
            print('Kernels shape wrong!')

        if ids == 0:
            if (np.sum(np.isfinite(kernel)) < 2):
                return np.nan, np.nan, np.nan
            # if (np.sum(np.isfinite(kernel)) > 2):
            #     ipdb.set_trace()


            smean = np.nanmean(kernel[dist-30:dist+30, dist:dist+67])

        #ycirc100e, xcirc100e = ua.draw_circle(dist + 100, dist + 1, 100)  # at - 150km, draw 50km radius circle
        wet = np.nansum(kernel[dist-30:dist+30, dist:dist+67]>=1)/np.sum(np.isfinite(kernel[dist-30:dist+30, dist:dist+67]))
        dry = np.nansum(kernel[dist - 30:dist + 30, dist:dist + 67] <= -1) / np.sum(
            np.isfinite(kernel[dist - 30:dist + 30, dist:dist + 67]))


        if wet >= 0.5:
            wetflag +=1

        if dry >= 0.5:
            dryflag +=1

    return smean, wetflag, dryflag



def file_loop(df):

    date = df['date'].iloc[0]
    hour = df['hour'].iloc[0]
    print('Doing day: ', date)

    storm_date = date

    dayd = pd.Timedelta('1 days')

    if (hour) <= 13:
        print('Nighttime')
        lsta_date = storm_date - dayd
    else:
        print('Daytime')
        lsta_date = storm_date

    fdate = str(lsta_date.year) + str(lsta_date.month).zfill(2) + str(lsta_date.day).zfill(2)

    topo = xr.open_dataset(cnst.LSTA_TOPO)
    ttopo = topo['h']
    #
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0])+abs(grad[1])

    smpath = [cnst.AMSRE_ANO_DAY + 'sma_' + fdate + '.nc',
              cnst.AMSRE_ANO_NIGHT + 'sma_' + fdate + '.nc',
              ]

    smlist = []
    dist = 200

    for sid , sp in enumerate(smpath):

        try:
            lsta = xr.open_dataset(sp)
        except OSError:
                return None
        print('Doing '+ sp)

        lsta = lsta.sel(lon=slice(-11, 11), lat=slice(9, 21))

        lsta_da = lsta['SM'].squeeze()

        if sid == 0:
            if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
                print('Not enough valid')
                return None

        try:
            lsta_da = topo.salem.transform(lsta_da)
        except RuntimeError:
            print('lsta_da on LSTA interpolation problem')
            return None


        lsta_da.values[ttopo.values>=450] = np.nan
        lsta_da.values[gradsum>30] = np.nan

        smlist.append(lsta_da)
        del lsta

    if len(smlist)!=2:
        return None

    del topo

    for dids, dit in df.iterrows():

        point = lsta_da.sel(lat=dit.lat, lon=dit.lon, method='nearest')

        plat = point['lat'].values
        plon = point['lon'].values



        xpos = np.where((smlist[0])['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where((smlist[0])['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            smmean, wetflag, dryflag = cut_kernel(xpos, ypos, smlist, dist)
        except TypeError:
            continue

        # if np.isfinite(smmean):
        #     ipdb.set_trace()

        df.loc[dids, 'SMmean'] = smmean
        df.loc[dids,'SMdry'] = dryflag
        df.loc[dids,'SMwet'] = wetflag

    return df