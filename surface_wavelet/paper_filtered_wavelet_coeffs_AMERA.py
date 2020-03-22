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
from utils import u_arrays as ua, u_darrays
import os
import collections
import warnings
from scipy import ndimage
import salem

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1


daykey = 'day-1'


def run_hours():

    l = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7] #15,16,
    for ll in l:
        composite(ll)

def composite(hour):

    key = 'NEWTRACKING'
    #msgopen = pd.read_csv(cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/cores_gt15000km2_table_AMSRE_LSTA_tracking_new_2hOverlap_'+str(h)+'.csv', na_values=-999)

    msgopen = pd.read_csv(
        cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/init_merged2/cores_gt15000km2_table_AMSRE_tracking2_' + str(
            hour) + '_init.csv', na_values=[-999, -99])

    msg = pd.DataFrame.from_dict(msgopen)

    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))

    msgopen = msg

    #basic filter
    msgopen = msgopen[(msgopen['lat']>9.5) & (msgopen['lat']<20.5) & (msgopen['topo']<=450) & (msgopen['dtime']<=2)]
    #propagation filter
    msgopen = msgopen[(msgopen['xdiff']>=100) | (msgopen['initTime'] <= 2.5)]
    #lsta filter
    #msgopen = msgopen[msgopen['LSTAslotfrac']>=0.05]
    #wetness_filter
    msgopen = msgopen[np.isfinite(msgopen['SMmean0'])]# & np.isfinite(msgopen['SMmean-1'])]
    # #eraq_filter
    msgopen = msgopen[(msgopen['ERAqmean'] >= 14)]#& (msgopen['ERAqmean'] <= 15.7)]
    # #dry_filter
    msgopen = msgopen[(msgopen['SMmean0']<=-5.48)]#&(msgopen['SMmean-1']<=-0.01)] #294 cases, with q 312
    # # #wet_filter
    # msgopen = msgopen[(msgopen['SMmean0']>=0.31) & (msgopen['SMmean-1']>=-0.01)] #295 cases, with q 318

    msgin = msgopen
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
    #
    # #return
    pool = multiprocessing.Pool(processes=1)

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
    pkl.dump(dic, open(outpath+"coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_2dAMERA_"+daykey+"_ALLS_minusMean_INIT_" + key + ".p", "wb"))
    print('Save file written!')




def cut_kernel_lsta(xpos, ypos, arr, mean=False):

    dist = 200

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)
    if mean:
        kernel = kernel - np.nanmean(kernel)

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

    prev_day = storm_date - dayd
    next_day = storm_date + dayd

    ctimes = {'day0' : [daybefore, cnst.AMSRE_ANO_DAY],
              'night0': [daybefore, cnst.AMSRE_ANO_NIGHT],
              'day-1' : [prev_day, cnst.AMSRE_ANO_DAY],
              'night-1': [prev_day, cnst.AMSRE_ANO_NIGHT],
              'night+1': [next_day, cnst.AMSRE_ANO_NIGHT],
              'day+1' : [next_day, cnst.AMSRE_ANO_DAY]}

    tag = daykey

    outime = (ctimes[tag])[0]

    fdate = str(outime.year) + str(outime.month).zfill(2) + str(outime.day).zfill(2)


    #topo = xr.open_dataset(cnst.LSTA_TOPO)
    topo = xr.open_dataset(cnst.WA_TOPO_3KM)
    topo = topo.sel(lat=slice(7.5,23), lon=slice(-12,12))

    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    amsre = xr.open_dataset((ctimes[tag])[1] + 'sma_' + fdate + '.nc')
    amsre = amsre.sel(time=str(outime.year)+'-'+str(outime.month)+'-'+str(outime.day))
    amsre = amsre.sel(lon=slice(-12.5, 12.5), lat=slice(7, 23.5))
    print('Doing '+ 'AMSR_' + str(outime.year) + str(outime.month).zfill(2) + str(
        outime.day).zfill(2) + '.nc')

    amsr_da = amsre['SM'].squeeze()

    file = cnst.ERA5
    t1 = outime.replace(hour=12)
    csmm = xr.open_dataset(
        file + 'hourly/surface/ERA5_' + str(t1.year) + '_' + str(t1.month).zfill(2) + '_srfc_SM.nc')
    csmm = u_darrays.flip_lat(csmm)


    sm_clim = xr.open_dataset(file + 'monthly/synop_selfmade/CLIM_2006-2010/ERA5_2006-2010_CLIM_'+str(t1.month).zfill(2)+'-'+str(t1.day).zfill(2)+'-'+str(t1.hour).zfill(2)+'_srfc_SM.nc').load()
    sm_clim = u_darrays.flip_lat(sm_clim)
    try:
        sm_clim = sm_clim.rename({'lat':'latitude', 'lon':'longitude'})
    except:
        pass

    csmm = csmm.sel(longitude=slice(-12.5, 12.5), latitude=slice(7, 23.5))
    sm_clim = sm_clim.sel(longitude=slice(-12.5, 12.5), latitude=slice(7, 23.5))
    ssm = csmm.sel(time=t1)

    lsta_da = ssm['swvl1'] - sm_clim['swvl1'].squeeze()


    print('Doing ' + 'lsta_daily_' + fdate + '.nc')


    try:
        amsr_da = topo.salem.transform(amsr_da)
    except RuntimeError:
        print('amsr_da on LSTA interpolation problem')
        return None


    try:
        lsta_da = topo.salem.transform(lsta_da)
    except RuntimeError:
        print('lsta_da on LSTA interpolation problem')
        return None


    lsta_da.values[ttopo.values >= 450] = np.nan
    lsta_da.values[gradsum > 30] = np.nan

    amsr_da.values[ttopo.values >= 450] = np.nan
    amsr_da.values[gradsum > 30] = np.nan

    del topo
    del sm_clim
    del csmm

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

        if daykey == 'day0':
            sh = 0#90
        else:
            sh = 0

        try:
            lsta_kernel = cut_kernel_lsta(xpos-sh, ypos, lsta_da.values, mean=True)
        except TypeError:
            print('LSTA kernel error')
            continue

        try:
            amsre_kernel = cut_kernel_lsta(xpos, ypos, amsr_da.values, mean=True)
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



def plot_amsr_ERA_trio(hour):
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    key = 'NEWTRACKING'

    f = plt.figure(figsize=(10.5, 3), dpi=300)


    names = ['',

        "/coeffs_nans_stdkernel_USE_" + str(
            hour) + "UTC_15000_2dAMSL_day0_ALLS_minusMean_CMORPH_DRY_INIT_Q20_NEWTRACKING.p",
                    "/coeffs_nans_stdkernel_USE_" + str(
            hour) + "UTC_15000_2dAMSL_day+1_ALLS_minusMean_CMORPH_DRY_INIT_Q20_NEWTRACKING.p"

    ]

    labels = ['Day-1', 'Day0', 'Day+1']

    left = 0.01
    bottom = 0.1
    width = 0.3
    height=0.8

    spot = [[]]

    for ids, ll in enumerate(['day-1', 'day0', 'day+1']):


        dic = pkl.load(
            open(path + "/coeffs_nans_stdkernel_USE_" + str(hour) + "UTC_15000_2dAMERA_"+ll+"_ALLS_minusMean_INIT_" + (key) + ".p",
                 "rb"))



        lsta = (dic['lsta'])[0] / dic['lsta'][1]
        amsr = (dic['amsr'])[0] / dic['amsr'][1]
        #
        if ids > 0:
            dic2 = pkl.load(
                open(path + names[ids],
                     "rb"))
            amsr = (dic2['amsr'])[0] / dic2['amsr'][1]

        amsr = ndimage.gaussian_filter(amsr, 4, mode='nearest')
        lsta = ndimage.gaussian_filter(lsta, 8, mode='nearest')

        cores = dic['cores']

        lcnt = dic['lsta'][1]
        acnt = dic['amsr'][1]

        dist = 200
        llevels = np.array(list(np.arange(-1.5, 0, 0.25)) + list(np.arange(0.25, 1.6, 0.25)))#*12000
        alevels = np.array(list(np.arange(-2.5, 0, 0.5)) + list(np.arange(0.5, 2.51, 0.5)))#*12000
        alevels = [-2.5,-2,-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5,2,2.5]
        alevels = [-4,-3,-2,-1,-0.5, -0.25, 0.25,0.5,1,2,3,4]
        llevels = np.array(list(np.arange(-1.35, 0, 0.5)) + list(np.arange(0.5, 1.36, 0.5)))  # *12000
        llevels = [-1.35,-1,-0.6,0,0.6,1,1.35]

        # # llevels = np.array(list(np.arange(-1.6, 0, 0.2)) + list(np.arange(0.2, 1.61, 0.2)))#*12000  WET
        # # alevels = np.array(list(np.arange(-3, 0, 0.25)) + list(np.arange(0.25, 3.25, 0.25)))#*12000
        #
        # llevels = np.array(list(np.arange(-1.5, 0, 0.2)) + list(np.arange(0.2, 1.51, 0.2)))  # *12000 DRY
        # alevels = np.array(list(np.arange(-5, 0, 0.5)) + list(np.arange(0.5, 5.5, 0.5)))  # *12000

        ax = f.add_subplot(1,3,ids+1)

        mp1 = plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, amsr,
                          levels=alevels, cmap='RdBu', extend='both') #cmap='RdBu'


        cs = plt.contour((np.arange(0, 2 * dist + 1) - dist) * 3, (np.arange(0, 2 * dist + 1) - dist) * 3, lsta*100,
                      extend='both', levels=llevels, colors='k', linewidths=1, linestyles=['solid'])

        plt.clabel(cs, inline=1, fontsize=8, fmt="%1.1f")
        #plt.colorbar(label='K')


        ax.set_xlabel('km')
        if ids == 0:
            ax.set_ylabel('km')

        # if ids > 0:
        #     ax.set_yticklabels('')

        plt.axvline(x=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.axhline(y=0, linestyle='dashed', color='dimgrey', linewidth=1.2)
        plt.plot(0, 0, marker='o', color='dimgrey')

        plt.title(labels[ids],fontsize=10)


    plt.tight_layout()
    text = ['a', 'b', 'c']
    plt.annotate(text[0], xy=(0.06, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[1], xy=(0.36, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)
    plt.annotate(text[2], xy=(0.65, 0.92), xytext=(0, 4), xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points', fontweight='bold', fontname='Ubuntu', fontsize=13)

    f.subplots_adjust(right=0.91)
    cax = f.add_axes([0.92, 0.18, 0.015, 0.73])
    cbar = f.colorbar(mp1, cax)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('%', fontsize=10)

    plt.savefig(path + '2hOverlap/amsreVSlsta/MAPS_AMERA_TRIO_ALLS_minusMean_noCore_INIT2_' + str(hour).zfill(2) + '.png')
    plt.close('all')