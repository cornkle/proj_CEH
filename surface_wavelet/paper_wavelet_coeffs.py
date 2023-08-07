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
import collections
import warnings

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1

def run_hours():

    l = [17, 15,16,18,19,20,21,22,23,0,1,2,3,4,5,6,7,8,9,10,11,12,13] #15,16,
    for ll in l:
        composite(ll)

def composite(hour):
    pool = multiprocessing.Pool(processes=3)

    file = cnst.MCS_POINTS_DOM #MCS_TMIN #
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients' #corrected_LSTA/wavelet/large_scale

    hour = hour

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) & (msg['time.month'] <= 9) ]

    msg = msg.sel(lat=slice(10.2,19.3), lon=slice(-9.8, 9.8))
    res = pool.map(file_loop, msg)
    pool.close()
    print('return parallel')
    # res = []
    # for m in msg[0:30]:
    #     out = file_loop(m)
    #     res.append(out)

    res = [x for x in res if x is not None]

    snpos_list = []
    wepos_list = []

    rsnpos_list = []
    rwepos_list = []

    vkernel_list = []
    rkernel_list = []

    vkernel_cnt = []
    rkernel_cnt = []

    lsta_list = []
    lsta_cnt = []

    for r in res:
        snpos_list.append(np.squeeze(r[0]))
        wepos_list.append(np.squeeze(r[1]))

        rsnpos_list.append(np.squeeze(r[2]))
        rwepos_list.append(np.squeeze(r[3]))

        vkernel_list.append(r[4])
        rkernel_list.append(r[5])

        scales = r[6]

        vkernel_cnt.append(r[7])
        rkernel_cnt.append(r[8])

        lsta_list.append(r[9])
        lsta_cnt.append(r[10])

    dic = collections.OrderedDict([('SN-pos' , [snpos_list, rsnpos_list]),
                                   ('WE-pos' , [wepos_list, rwepos_list]),
                                   ('kernel' , [vkernel_list, rkernel_list]),
                                   ('lsta' , [lsta_list, lsta_cnt]),
                                   ('cnt' , [vkernel_cnt, rkernel_cnt]),
                                   ('scales' , scales)])

    keys = list(dic.keys())

    for l in keys:
        if l == 'scales':
            continue
        (dic[l])[0] = np.squeeze(np.vstack((dic[l])[0]))
        try:
            (dic[l])[1] = np.squeeze(np.vstack((dic[l])[1]))
        except IndexError:
            continue

    dic['nbcores'] = dic['SN-pos'][0].shape[0]
    dic['nbrcores'] = dic['SN-pos'][1].shape[0]

    for l in keys:
        if (l == 'scales') | (l == 'lsta'):
            continue

        a =  dic[l][0]
        b =  dic[l][1]

        nsstat, nspvalue = ttest(a, b, axis=0, equal_var=False, nan_policy='omit')
        mask = nspvalue < 0.05
        dic[l].append(mask)

        # if 'pos' in l:
        #     dic[l].append(np.nanstd(dic[l][0], axis=(0,2)))
        #     dic[l].append(np.nanstd(dic[l][1], axis=(0,2)))

    for l in keys:
        if l == 'scales':
            continue
        if 'pos' in l:
            (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
            (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)
        else:
            (dic[l])[0] = np.nansum((dic[l])[0], axis=0)
            try:
                (dic[l])[1] = np.nansum((dic[l])[1], axis=0)
            except IndexError:
                continue


    pkl.dump(dic, open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_-60_ALL_3slot.p", "wb"))
    print('Save file written!')



def cut_kernel(xpos, ypos, arr, date, lat, lon, rotate=False):

    dist = 100

    kernel = u_arrays.cut_kernel_3d(arr,xpos, ypos,dist)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    if kernel.shape != (arr.shape[0], 2*dist+1, 2*dist+1):
        return

    if rotate:
        kernel = u_met.era_wind_rotate3d(kernel, date, lat, lon, level=850, ref_angle=270)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        wav_ns = np.nanmean(kernel[:,:, 100:103], axis=2) #kernel[:,:,50] #

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        wav_we = np.nanmean(kernel[:,100:103,:], axis=1) #kernel[:,50,:] #



    # filler = kernel[[1,2,4,6,10], :, :]
    # slices = [(1,3), (3,5), (5,7), (7,9), (9,11)]
    # for ids, sl in enumerate(slices):
    #     filler[ids,]
    #

    return wav_ns, wav_we, kernel[[2,4,6,8], :, :]  #[1,4,6,10]


def cut_kernel_lsta(xpos, ypos, arr, date, lat, lon, nbslot=False):

    dist = 100

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)
    kernel = kernel - np.nanmean(kernel)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    if kernel.shape != (2*dist+1, 2*dist+1):
        return

    if nbslot is not False:
        nbsl = u_arrays.cut_kernel(nbslot, xpos, ypos, dist)
        if np.sum(nbsl[dist-10:dist+10,dist-10:dist+67]>5) / np.sum(np.isfinite(nbsl[dist-10:dist+10,dist-10:dist+67])) <=0.5:
            print('TOO FEW SLOTS!')
            return

    return kernel


def file_loop(fi):


    print('Doing day: ', fi.time)

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))

    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 13:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)

    try:
        lsta = xr.open_dataset(cnst.LSTA_NEW + 'lsta_daily_' + fdate + '.nc')
    except OSError:
        return None
    print('Doing ' + 'lsta_daily_' + fdate + '.nc')

    topo = xr.open_dataset(cnst.LSTA_TOPO)
    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    lsta_da = lsta['LSTA'].squeeze()
    slot_da = lsta['NbSlot'].squeeze().values

    #### remove mean from LSTA

    # f = plt.figure()
    # plt.imshow(lsta_da, origin='lower')

    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return None

    lsta_da.values[np.isnan(lsta_da.values)] = 0

    lsta_da.values[ttopo.values >= 450] = np.nan
    lsta_da.values[gradsum > 30] = np.nan
    pos = np.where((fi.values >= 5) & (fi.values <= 65)) # (fi.values >= 5) & (fi.values < 65) #(fi.values >= 5) & (fi.values < 65)
    del topo
    if (np.sum(pos) == 0):
        print('No blobs found')
        return None

    wav_input = lsta_da.values.copy()
    #points = np.where(np.isfinite(wav_input))
    inter1 = np.where(np.isnan(wav_input))

    # interpolate over sea from land points
    wav_input[inter1] = 0  #halfway between minus and plus rather than interpolate
    #
    # try:
    #      wav_input[inter1] = griddata(points, np.ravel(wav_input[points]), inter1, method='linear')
    # except ValueError:
    #     pass
    #
    # inter = np.where(np.isnan(wav_input))
    # try:
    #      wav_input[inter] = griddata(points, np.ravel(wav_input[points]), inter, method='nearest')
    # except ValueError:
    #     wav_input[inter]=0

    wavpos = wutil.applyHat_pure(wav_input, dataset='METSRFC')

    wavarr = wavpos['coeffs'].copy()
    scales = wavpos['scales'].copy()

    for inds, wava in enumerate(wavarr):
        wava[inter1] = np.nan
        #wavarr[inds,:,:] = wava/np.nanstd(wava)
        #wavarr[inds, :, :] = (wava-np.nanmean(wava)) / np.nanstd(wava)
        #input[np.abs(input)<1] = np.nan # minimum of 1*std, significance test
        wavarr[inds, :, :] = (wava) / (np.nanstd(wava))

        #wavarr[inds, inter1[0], inter1[1]] = np.nan
    del wavpos
    wav_input[inter1] = np.nan

    # xfi = fi.shape[1]
    # randx = np.random.randint(0,xfi,20)
    # if np.min(pos[0]) == 0:
    #     ymin = np.min(pos[0])
    # else:
    #     ymin = np.min(pos[0])-1
    # if np.max(pos[0]) == fi.shape[0]-1:
    #     ymax = np.max(pos[0])
    # else:
    #     ymax = np.max(pos[0])+1
    # randy = np.random.randint(-10,11,20)
    #
    # posr = (randy, randx)

    rnspos = []
    rwepos = []
    rkernel = []
    lsta = []


    ###############################Blob loop

    nspos = []
    wepos = []
    vkernel = []
    mcs_hour = xr.open_dataarray(cnst.MCS_HOUR_DAILY)
    mcsimage = xr.open_dataarray(cnst.MCS_15K)
    mcsimage = mcsimage.sel(time=fi.time, lat=slice(10.2,19.3), lon=slice(-9.8, 9.8))

    labels, goodinds = ua.blob_define(mcsimage.values, -50, minmax_area=[600,100000], max_area=None)

    for y, x in zip(pos[0], pos[1]):

        # print('yx', y ,x)
        # print('labels', labels)

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        if (labels[y,x] not in goodinds) | (labels[y,x] == 0):
            print('MCS too small!!')
            continue

        if (mcsimage.values[y,x] > -60):
            print('Core too warm!!')
            continue

        mhour = mcs_hour.sel(time=pd.datetime(int(fdate[0:4]), int(fdate[4:6]), int(fdate[6:8]), 14), lon=lon,
                             lat=lat).values
        if mhour <= 13:
            mhour += 24

        chour = fi['time.hour'].values

        if (chour >= 0) & (chour <= 13):
            chour += 24
        if (mhour < chour) | (np.isnan(mhour)):
            print('Core overlaps: earliest:', mhour, ' core: ', chour)
            continue
        try:
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.04)
        except KeyError:
            continue
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            lsta_kernel = cut_kernel_lsta(xpos, ypos, wav_input, daybefore, plon, plat, nbslot=slot_da)
        except TypeError:
            print('Kernel error')
            continue

        try:
            wavpos_ns, wavpos_we, vvkernel = cut_kernel(xpos, ypos, wavarr, daybefore, plon, plat, rotate=False)
        except TypeError:
            print('Kernel error')
            continue



        nspos.append(wavpos_ns)  # north-south wavelet
        wepos.append(wavpos_we)  # west-east wavelet
        vkernel.append(vvkernel)
        lsta.append(lsta_kernel)

        ##### random

        rdist = 50
        randy50 = [y-rdist, y-rdist, y-rdist, y, y, y+rdist, y+rdist, y+rdist]
        randx50 = [x-rdist, x, x+rdist, x-rdist, x+rdist, x-rdist, x, x + rdist]
        randy50_100 = [y - rdist,  y - rdist, y, y, y + rdist, y + rdist]

        rdist = 100
        randx100 = [x-rdist,  x+rdist, x-rdist, x+rdist, x-rdist, x + rdist]

        rdist = 150
        randx150 = [x-rdist,  x+rdist, x-rdist, x+rdist, x-rdist, x + rdist]


        randy = np.array(randy50 + randy50_100 + randy50_100)
        randx = np.array(randx50 + randx100 + randx150)
        #ipdb.set_trace()
        for ry, rx in zip(randy, randx):

            if ry < 0:
                continue
            if ry > fi.shape[0]-1:
                continue

            if rx < 0:
                continue
            if rx > fi.shape[1]-1:
                continue

            try:
                lat = fi['lat'][ry]
            except IndexError:
                ipdb.set_trace()
            try:
                lon = fi['lon'][rx]
            except IndexError:
                ipdb.set_trace()
            try:
                point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.04)
            except KeyError:
                continue
            plat = point['lat'].values
            plon = point['lon'].values

            xpos = np.where(lsta_da['lon'].values == plon)
            xpos = int(xpos[0])
            ypos = np.where(lsta_da['lat'].values == plat)
            ypos = int(ypos[0])

            try:
                wavpos_ns, wavpos_we, rrkernel = cut_kernel(xpos, ypos, wavarr, daybefore, plon, plat, rotate=False)
            except TypeError:
                print('Kernel random error')
                continue

            rnspos.append(wavpos_ns)  # north-south wavelet
            rwepos.append(wavpos_we)  # west-east wavelet
            rkernel.append(rrkernel)

    del lsta_da
    del mcs_hour
    del mcsimage

    if nspos == []:
        print('NorthSouth empty')
        return None

    if (len(nspos) == 1) | (len(rnspos) == 1):
        print('NorthSouth 1')
        return None
    else:

        nspos_sum = np.squeeze(np.stack(nspos, axis=0))
        wepos_sum = np.squeeze(np.stack(wepos, axis=0))
        vkernel_sum = np.nansum(np.stack(vkernel, axis=0), axis=0)[np.newaxis,...]
        vkernel_cnt = np.sum(np.isfinite(np.stack(vkernel, axis=0)), axis=0)[np.newaxis,...]
        lsta_sum = np.nansum(np.stack(lsta, axis=0), axis=0)[np.newaxis,...]
        lsta_cnt = np.sum(np.isfinite(np.stack(lsta, axis=0)), axis=0)[np.newaxis,...]

        rnspos_sum = np.squeeze(np.stack(rnspos, axis=0))
        rwepos_sum = np.squeeze(np.stack(rwepos, axis=0))
        rkernel_sum = np.nansum(np.stack(rkernel, axis=0), axis=0)[np.newaxis,...]
        rkernel_cnt = np.sum(np.isfinite(np.stack(rkernel, axis=0)), axis=0)[np.newaxis,...]




    print('Returning with kernel, success!!')

    return (nspos_sum, wepos_sum,
            rnspos_sum, rwepos_sum, vkernel_sum, rkernel_sum, scales, vkernel_cnt, rkernel_cnt, lsta_sum, lsta_cnt)
# rcns_sum, rcwe_sum, cns_sum, cwe_sum,


if __name__ == "__main__":
    run_hours()
