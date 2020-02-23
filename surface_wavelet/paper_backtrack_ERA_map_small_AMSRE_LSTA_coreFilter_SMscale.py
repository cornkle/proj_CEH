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
from wavelet import util as wutil


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


key = '2hOverlap'


def run_hours():

    l = [15,16,17,18,19,20,21,22,23, 0,1,2,3,4,5,6,7]  #
    for ll in l:
        composite(ll)

def rewrite_list(hour):
    path = '/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/ERA5/core_txt/'
    dic = pkl.load(
        open(path + "cores_gt15000km2_table_AMSRE_LSTA_tracking_SMcircle_" + key + '_' + str(hour) + ".p", "rb"))
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
    df = df.reindex(columns=['year', 'month', 'day', 'hour', 'lon', 'lat', 'xloc', 'yloc', 'area', 'csize', 't', 'storm_id', 'topo','SMmean0', 'SMdry0', 'SMwet0','SMmean-1', 'SMdry-1', 'SMwet-1', 'LSTAmean', 'LSTAslotfrac', 'dtime', 'ERAqmean', 'SMscale'])
    df.to_csv(path + "cores_gt15000km2_table_AMSRE_LSTA_tracking_SMcircle_" + key + '_' + str(hour) + ".csv", na_rep=-999, index_label='id')


def composite(h):

    path = '/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/ERA5/core_txt/'
    msgopen = pd.read_csv(
        '/home/ck/DIR/cornkle/figs/LSTA/corrected_LSTA/new/wavelet_coefficients/core_txt/cores_gt15000km2_table_1640_580_'+str(h) + '_' + key +'.csv')
    hour = h
    msg = pd.DataFrame.from_dict(msgopen)# &  &

    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))

    # calculate the chunk size as an integer
    #'chunk_size = int(msg.shape[0] / pnumber)
    msg.sort_values(by='date')
    msg['SMmean0'] = np.nan
    msg['SMdry0'] = np.nan
    msg['SMwet0'] = np.nan
    msg['SMmean-1'] = np.nan
    msg['SMdry-1'] = np.nan
    msg['SMwet-1'] = np.nan
    msg['LSTAmean'] = np.nan
    msg['LSTAslotfrac'] = np.nan
    msg['ERAqmean'] = np.nan
    msg['SMscale'] = np.nan
    #msg['topo'] = np.nan

    chunk, chunk_ind, chunk_count = np.unique(msg.date, return_index=True, return_counts=True)


    chunks = [msg.loc[msg.index[ci:ci + cc]] for ci, cc in zip(chunk_ind, chunk_count)] # daily chunks

    # res = []
    # for m in chunks[150:155]:
    #     out = file_loop(m)
    #     res.append(out)
    #
    # #ipdb.set_trace()
    # return
    pool = multiprocessing.Pool(processes=4)

    res = pool.map(file_loop, chunks)
    pool.close()

    print('Returned from parallel')
    # ipdb.set_trace()
    res = [x for x in res if x is not None]

    df_concat = pd.concat(res)

    dic = df_concat.to_dict()

    #ipdb.set_trace()

    pkl.dump(dic, open(path+"/cores_gt15000km2_table_AMSRE_LSTA_tracking"
                            "_SMcircle_" + key + "_" +str(hour)+".p", "wb"))  #"+str(hour)+"
    print('Save file written!')
    print('Dumped file')

    rewrite_list(hour)



def cut_kernel(xpos, ypos, arrlist, dist, era=False):

    wetflag = np.array([np.nan, np.nan])
    dryflag = np.array([np.nan, np.nan])
    smean = np.array([np.nan, np.nan])
    outscale = 0
    emean = np.nan

    for ids, arr in enumerate(arrlist):


        kernel = ua.cut_kernel(arr,xpos, ypos,dist)
        kernel = kernel #- np.nanmean(kernel)
        if kernel.shape != (dist*2+1, dist*2+1):
            print('Kernels shape wrong!')

        ycirc100e, xcirc100e = ua.draw_circle(dist + 51, dist + 1, 34)  # at - 150km, draw 50km radius circle 17pixel
        outmean = np.nanmean(kernel[ycirc100e, xcirc100e])


        #ipdb.set_trace()
        if ids == 0:

            if (np.sum(np.isfinite(kernel)) < 2):
                return smean, wetflag, dryflag, outscale

            wav_input = kernel.copy()
            wav_input = wav_input - np.nanmean(wav_input)
            # points = np.where(np.isfinite(wav_input))
            inter1 = np.where(np.isnan(wav_input))

            # interpolate over sea from land points
            wav_input[inter1] = 0  # halfway between minus and plus rather than interpolate

            wavpos = wutil.applyHat_pure(wav_input, dataset='METSRFC24')

            wavarr = wavpos['coeffs'].copy()
            orig = wavpos['scales']

            scales = np.round(orig,0)
            #print(scales)


            for inds, wava in enumerate(wavarr):
                wava[inter1] = np.nan
                wavarr[inds, :, :] = (wava) / (np.nanstd(wava))

            try:
                maxx = np.nanmin(wavarr[:, ycirc100e, xcirc100e])

                # if np.isfinite(outmean):
                #     plt.figure()
                #     plt.contourf(np.arange(dist*2+1), scales, wavarr[:, dist+1, :], cmap='RdBu', levels=np.linspace(-3,3,6))
                #     plt.colorbar()
                #
                #     ipdb.set_trace()


                maxpos = np.where(wavarr==maxx)
                if len(maxpos[0]) > 1:
                    print('More than one dominant wavelet scale!')
                outscale = scales[(maxpos[0])[0]]

                if np.abs(maxx) < 1:
                    outscale = 0

            except:
                outscale = 0


            del wavpos
            del wav_input


        #outmean = np.nanmean(kernel[dist-10:dist+10, dist:dist+67])  # dist-30 / dist + 30, dist+100


        smean[ids] = outmean

        #ycirc100e, xcirc100e = ua.draw_circle(dist + 100, dist + 1, 100)  # at - 150km, draw 50km radius circle
        # wet = np.nansum(kernel[dist-10:dist+10, dist:dist+67]>=0.1)/np.sum(np.isfinite(kernel[dist-10:dist+10, dist:dist+67]))
        # dry = np.nansum(kernel[dist - 10:dist + 10, dist:dist + 67] <= -1) / np.sum(
        #     np.isfinite(kernel[dist - 10:dist + 10, dist:dist + 67]))

        wet = np.nansum(kernel[ycirc100e,xcirc100e]>=0.1)/np.sum(np.isfinite(kernel[ycirc100e,xcirc100e]))
        dry = np.nansum(kernel[ycirc100e,xcirc100e] <= -1) / np.sum(
            np.isfinite(kernel[ycirc100e,xcirc100e]))


        wetflag[ids] = wet
        dryflag[ids] = dry

    if era is not False:
        ekernel = ua.cut_kernel(era, xpos, ypos, dist)
        emean = np.round(np.nanmean(ekernel[dist - 10:dist + 10, dist:dist + 67]) * 1000, 2)


    return smean, wetflag, dryflag, outscale, emean


def cut_kernel_lsta(xpos, ypos, arr, nbslot=False):

    dist = 200
    isgood = 1

    kernel = ua.cut_kernel(arr,xpos, ypos,dist)
    kernel = kernel #- np.nanmean(kernel)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    if kernel.shape != (2*dist+1, 2*dist+1):
        return

    if nbslot is not False:
        nbsl = ua.cut_kernel(nbslot, xpos, ypos, dist)
        if np.sum(nbsl[dist-10:dist+10,dist-10:dist+67]>2) / np.sum(np.isfinite(nbsl[dist-10:dist+10,dist-10:dist+67])) <=0.5:
            print('TOO FEW SLOTS!')
            isgood = 0

    outmean = np.nanmean(kernel[dist - 10:dist + 10, dist:dist + 67])
    slot2frac = np.sum(nbsl[dist-10:dist+10,dist-10:dist+67]>2) / np.sum(np.isfinite(nbsl[dist-10:dist+10,dist-10:dist+67]))



    return outmean, slot2frac




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


    topo = xr.open_dataset(cnst.WA_TOPO_3KM)
    topo = topo.sel(lon=slice(-12, 12), lat=slice(8, 21))
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
                print('AMSRE OS ERROR')
                return None
        print('Doing '+ sp)

        lsta = lsta.sel(lon=slice(-12.5, 12.5), lat=slice(7.7, 22))

        sm_da = lsta['SM'].squeeze()

        if sid == 0:
            if (np.sum(np.isfinite(sm_da)) / sm_da.size) < 0.05:
                print('Not enough valid')
                return None

        try:
            sm_da = topo.salem.transform(sm_da)
        except RuntimeError:
            print('lsta_da on LSTA interpolation problem')
            return None

        smlist.append(sm_da)
        del lsta

    if len(smlist)!=2:
        return None

    try:
        lsta = xr.open_dataset(cnst.LSTA_1330 + 'lsta_daily_' + fdate + '.nc')
    except OSError:
        print('LSTA OS ERROR')
        return None
    print('Doing ' + 'lsta_daily_' + fdate + '.nc')

    lsta_da = lsta['LSTA'].squeeze()


    slot_da = lsta['NbSlot'].squeeze().values


    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return None

    file = cnst.ERA5
    edate = storm_date
    edate = edate.replace(hour=12)
    pl_clim = xr.open_dataset(
        file + 'monthly/synop_selfmade/CLIM_2006-2010/ERA5_2006-2010_CLIM_' + str(edate.month).zfill(2) + '-' + str(
            edate.day).zfill(2) + '-' + str(edate.hour).zfill(2) + '_pl.nc').load()
    pl_clim = u_darrays.flip_lat(pl_clim)
    try:
        pl_clim = pl_clim.rename({'lat': 'latitude', 'lon': 'longitude'})
    except:
        pass


    lowq = pl_clim['q'].sel(level=925).squeeze()

    try:
        era_da = topo.salem.transform(lowq)
    except RuntimeError:
        print('lsta_da on LSTA interpolation problem')
        return None

    del topo
    for dids, dit in df.iterrows():

        nolat = 1
        try:
            point = lsta_da.sel(lat=dit.lat, lon=dit.lon, method='nearest', tolerance=0.04)
        except KeyError:
            print('LSTA tolerance error')
            nolat = 0
            lsta_kernel = np.nan
            slot2frac = np.nan

        if nolat:
            plat = point['lat'].values
            plon = point['lon'].values

            xpos = np.where(lsta_da['lon'].values == plon)
            xpos = int(xpos[0])
            ypos = np.where(lsta_da['lat'].values == plat)
            ypos = int(ypos[0])

            try:
                lsta_kernel, slot2frac = cut_kernel_lsta(xpos, ypos, lsta_da, nbslot=slot_da)
            except TypeError:
                print('LSTA Kernel error')
                lsta_kernel = np.nan
                slot2frac = np.nan


        try:
            point = smlist[0].sel(lat=dit.lat, lon=dit.lon, method='nearest', tolerance=0.04)
        except KeyError:
            ipdb.set_trace()
            print('AMSR tolerance error')
            continue

        plat = point['lat'].values
        plon = point['lon'].values

        xxpos = np.where((smlist[0])['lon'].values == plon)
        xxpos = int(xxpos[0])
        yypos = np.where((smlist[0])['lat'].values == plat)
        yypos = int(yypos[0])


        try:

            smmean, wetflag, dryflag, scale, emean = cut_kernel(xxpos, yypos, smlist, dist, era=era_da)
        except:
            print('SM TypeError')
            continue


        # if np.isfinite(smmean):
        #     ipdb.set_trace()
        try:
            df.loc[dids, 'SMmean0'] = smmean[0]
        except TypeError:
            ipdb.set_trace()
        df.loc[dids,'SMdry0'] = dryflag[0]
        df.loc[dids,'SMwet0'] = wetflag[0]
        df.loc[dids, 'SMmean-1'] = smmean[1]
        df.loc[dids, 'SMdry-1'] = dryflag[1]
        df.loc[dids, 'SMwet-1'] = wetflag[1]
        df.loc[dids, 'LSTAmean'] = lsta_kernel
        df.loc[dids, 'LSTAslotfrac'] = slot2frac  # percentage of >2 slots in centred area
       # df.loc[dids, 'topo'] = topo[ypos,xpos]
        df.loc[dids,'ERAqmean'] = emean
        df.loc[dids,'SMscale'] = scale
    #ipdb.set_trace()

    return df