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
import collections

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1

def run_hours():

    l = [6, 21, 0, 3, 18]
    for ll in l:
        composite(ll)


def composite(hour):
    pool = multiprocessing.Pool(processes=5)

    file = cnst.MCS_POINTS_DOM #MCS_TMIN #
    path = cnst.network_data + '/figs/LSTA-bullshit/AGU' #corrected_LSTA/wavelet/large_scale

    hour = hour

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) ]

    msg = msg.sel(lat=slice(10.2,19.3), lon=slice(-9.7, 9.7))
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

    dic = collections.OrderedDict([('SN-pos' , [snpos_list, rsnpos_list]),
                                   ('WE-pos' , [wepos_list, rwepos_list]),
                                   ('kernel' , [vkernel_list, rkernel_list]),
                                   ('lsta' , [lsta_list]),
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
        sa = a.shape
        sb = b.shape
        #if 'pos' in l:
        #    a = (a.swapaxes(0,1).reshape(sa[1], sa[0]*sa[2]).T/np.nanstd(dic[l][0], axis=(0,2))).T.reshape(sa[1], sa[0],sa[2]).swapaxes(0,1)
        #    b = (b.swapaxes(0,1).reshape(sb[1], sb[0]*sb[2]).T/np.nanstd(dic[l][1], axis=(0,2))).T.reshape(sb[1], sb[0], sb[2]).swapaxes(0,1)

        nsstat, nspvalue = ttest(a, b, axis=0, equal_var=False, nan_policy='omit')
        mask = nspvalue < 0.05
        dic[l].append(mask)

        if 'pos' in l:
            dic[l].append(np.nanstd(dic[l][0], axis=(0,2)))
            dic[l].append(np.nanstd(dic[l][1], axis=(0,2)))

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


    pkl.dump(dic, open(path+"/coeffs_test_nans_stdkernel"+str(hour)+"UTC.p", "wb"))
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

    with np.errstate(invalid='ignore'):
        wav_ns = np.nanmean(kernel[:,:, 99:102], axis=2) #kernel[:,:,50] #
    # if np.sum(np.isfinite(wav_ns)) < 0.2 * wav_ns.size:
    #     return
        wav_we = np.nanmean(kernel[:,99:102,:], axis=1) #kernel[:,50,:] #



    return wav_ns, wav_we, kernel[[1,4,6,10], :, :]


def cut_kernel_lsta(xpos, ypos, arr, date, lat, lon):

    dist = 100

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    if kernel.shape != (2*dist+1, 2*dist+1):
        return

    return kernel


def file_loop(fi):


    print('Doing day: ', fi.time)

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))

    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 16:
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

    #### remove mean from LSTA

    # f = plt.figure()
    # plt.imshow(lsta_da, origin='lower')

    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return None

    lsta_da.values[np.isnan(lsta_da.values)] = 0

    lsta_da.values[ttopo.values >= 450] = np.nan
    lsta_da.values[gradsum > 30] = np.nan
    pos = np.where((fi.values >= 5) & (fi.values < 65)) # (fi.values >= 5) & (fi.values < 65) #(fi.values >= 5) & (fi.values < 65)

    if (np.sum(pos) == 0):
        print('No blobs found')
        return None

    wav_input = lsta_da.values.copy()
    points = np.where(np.isfinite(wav_input))
    inter1 = np.where(np.isnan(wav_input))

    # interpolate over sea from land points
    wav_input[inter1] = 0  #halfway between minus and plus rather than interpolate

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

    wavpos = wutil.applyHat(wav_input, dataset='METSRFC')

    wavarr = wavpos['coeffs'].copy()

    for inds, wava in enumerate(wavarr):
        wava[inter1] = np.nan
        wavarr[inds,:,:] = wava/np.nanstd(wava)
        #wavarr[inds, inter1[0], inter1[1]] = np.nan

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

    for y, x in zip(pos[0], pos[1]):

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
            wavpos_ns, wavpos_we, vvkernel = cut_kernel(xpos, ypos, wavarr, daybefore, plon, plat, rotate=False)
        except TypeError:
            print('Kernel error')
            continue

        try:
            lsta_kernel = cut_kernel_lsta(xpos, ypos, wav_input, daybefore, plon, plat)
        except TypeError:
            print('Kernel error')
            continue

        nspos.append(wavpos_ns)  # north-south wavelet
        wepos.append(wavpos_we)  # west-east wavelet
        vkernel.append(vvkernel)
        lsta.append(lsta_kernel)

        ##### random

        randx = np.array(np.linspace(50,fi.shape[1]-50, 5), dtype='int16') #np.random.randint(0, fi.shape[1], 10)
        randy = np.random.randint(-10, 11, 9)
        #ipdb.set_trace()
        for ry, rx in zip(randy, randx):

            if y+ry < 0: ry= ry*-1
            if y+ry > fi.shape[0]-1: ry = ry*-1
            try:
                lat = fi['lat'][ry+y]
            except IndexError:
                ipdb.set_trace()
            lon = fi['lon'][rx]

            point = lsta_da.sel(lat=lat, lon=lon, method='nearest')
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

        # nspos_sum = ((nspos_sum - np.nanmin(nspos_sum, axis=2)[..., None]) /
        #         (np.nanmax(nspos_sum, axis=2) - np.nanmin(nspos_sum, axis=2))[..., None])
        #
        # wepos_sum = ((wepos_sum - np.nanmin(wepos_sum, axis=2)[..., None]) /
        #              (np.nanmax(wepos_sum, axis=2) - np.nanmin(wepos_sum, axis=2))[..., None])

        rnspos_sum = np.squeeze(np.stack(rnspos, axis=0))
        rwepos_sum = np.squeeze(np.stack(rwepos, axis=0))
        rkernel_sum = np.nansum(np.stack(rkernel, axis=0), axis=0)[np.newaxis,...]
        rkernel_cnt = np.sum(np.isfinite(np.stack(rkernel, axis=0)), axis=0)[np.newaxis,...]


    scales = wavpos['scales']

    print('Returning')

    return (nspos_sum, wepos_sum,
            rnspos_sum, rwepos_sum, vkernel_sum, rkernel_sum, scales, vkernel_cnt, rkernel_cnt, lsta_sum)
# rcns_sum, rcwe_sum, cns_sum, cwe_sum,


if __name__ == "__main__":
    run_hours()
