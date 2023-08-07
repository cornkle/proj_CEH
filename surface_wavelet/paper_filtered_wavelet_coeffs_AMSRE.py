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
from scipy import ndimage
import collections
import warnings
import salem

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1

def run_hours():

    l = [15,16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]#,8,9,10,11,12,13] #15,16,
    for ll in l:
        composite(ll)

def composite(h):

    #key = '2hOverlap'

    #msgopen = pd.read_csv(cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/cores_gt15000km2_table_AMSRE_LSTA_tracking_new_' + key + '_'+str(hour)+'.csv', na_values=-999)
    key = 'NEWTRACKING'
    #msgopen = pd.read_csv(cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/cores_gt15000km2_table_AMSRE_LSTA_tracking_new_2hOverlap_'+str(h)+'.csv', na_values=-999)

    msgopen = pd.read_csv(
        cnst.network_data + 'figs/LSTA/corrected_LSTA/new/ERA5/core_txt/init_merged2/cores_gt15000km2_table_AMSRE_tracking2_' + str(
            h) + '_init.csv', na_values=[-999, -99])

    hour = h
    msg = pd.DataFrame.from_dict(msgopen)# &  &

    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))

    msgopen = msg

    # #basic filter
    # msgopen = msgopen[(msgopen['lat']>9.5) & (msgopen['lat']<20.5) & (msgopen['topo']<=450) & (msgopen['dtime']<=2)]
    # #propagation filter
    # msgopen = msgopen[(msgopen['xdiff']>=100) | (msgopen['initTime'] <= 2.5)]
    # #wetness_filter
    # msgopen = msgopen[np.isfinite(msgopen['SMmean0'])]# & (msgopen['LSTAslotfrac']>=0.05)]# & np.isfinite(msgopen['SMmean-1'])]
    # # #eraq_filter
    # msgopen = msgopen[(msgopen['ERAqmean'] >= 14)]#& (msgopen['ERAqmean'] <= 15.7)]
    # # #dry_filter
    # msgopen = msgopen[(msgopen['SMmean0']<=-5.23)]#&(msgopen['SMmean-1']<=-0.01)] #294 cases, with q 312
    # # # #wet_filter
    # # msgopen = msgopen[(msgopen['SMmean0']>=0.3) & (msgopen['SMmean-1']>=-0.01)] #295 cases, with q 318


    #basic filter
    msgopen = msgopen[(msgopen['lat']>9.5) & (msgopen['lat']<20.5) & (msgopen['topo']<=450) & (msgopen['dtime']<=2)]
    #propagation filter
    msgopen = msgopen[(msgopen['xdiff']>=100) | ((msgopen['initTime'] <10 )&(msgopen['SMmean0'] >=-12)&(msgopen['SMmean0'] <=-7 ) & (msgopen['SMmean-1']<=-1) & (msgopen['SMdry0']==1))]# | ((msgopen['initTime'] <11 )&(msgopen['initTime'] >2 ))]
    #msgopen = msgopen[(msgopen['xdiff'] >= 120) | ((msgopen['initTime'] < 10))]# & (msgopen['initTime'] > 2))]# & (msgopen['SMmean0'] >=0.1)  & (msgopen['SMmean-1'] >=0))]  # | ((msgopen['initTime'] <11 )&(msgopen['initTime'] >2 ))]

    #lsta filter
    #msgopen = msgopen[msgopen['LSTAslotfrac']>=0.05]
    #wetness_filter
    msgopen = msgopen[np.isfinite(msgopen['SMmean0'])]# & np.isfinite(msgopen['SMmean-1'])]
    #eraq_filter
    #msgopen = msgopen[(msgopen['ERAqmean'] >= 14)] #14
    # #dry_filter
    msgopen = msgopen[(msgopen['SMmean0']<=-1.8)]# -1.8 & (msgopen['SMmean-1']<=-0.5)]# & (msgopen['SMmean-1']<=-1)]#&(msgopen['SMmean0']<=-7.19)] #294 cases, with q 312
    # # #wet_filter
    #msgopen = msgopen[(msgopen['SMmean0']>=0.31) ]#& (msgopen['SMmean0']>=1.8)] #295 cases, with q 318, 0.16-> 317, noMCS filter



    msgin = msgopen

    print('Number of cores', len(msgin))
    #ipdb.set_trace()

    # calculate the chunk size as an integer
    #'chunk_size = int(msg.shape[0] / pnumber)
    msgin.sort_values(by='date')
    #ipdb.set_trace()

    msgy = msgin

    chunk, chunk_ind, chunk_count = np.unique(msgy.date, return_index=True, return_counts=True)

    chunks = [msgy.loc[msgy.index[ci:ci + cc]] for ci, cc in zip(chunk_ind, chunk_count)] # daily chunks

    # res = []
    # for m in chunks[115:120]:
    #     out = file_loop(m)
    #     res.append(out)
    #
    # return
    #ipdb.set_trace()
    pool = multiprocessing.Pool(processes=3)

    res = pool.map(file_loop, chunks)
    pool.close()

    print('Returned from parallel')


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

    outpath = cnst.network_data + '/figs/LSTA/corrected_LSTA/new/wavelet_coefficients/'
    pkl.dump(dic, open(outpath+"coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_WAVELET_INIT_SCALES_" + key + ".p", "wb"))
    print('Save file written!')



def cut_kernel(xpos, ypos, arr, date, lat, lon, rotate=False):

    dist = 200

    kernel = u_arrays.cut_kernel_3d(arr,xpos, ypos,dist)

    # if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
    #     return

    if kernel.shape != (arr.shape[0], 2*dist+1, 2*dist+1):
        print('Wavelet kernel wrong shape')
        return

    if rotate:
        kernel = u_met.era_wind_rotate3d(kernel, date, lat, lon, level=850, ref_angle=270)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        wav_ns = np.nanmean(kernel[:,:, 90:111], axis=2) #kernel[:,:,50] #

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        wav_we = np.nanmean(kernel[:,90:111,:], axis=1) #kernel[:,50,:] #


    # filler = kernel[[1,2,4,6,10], :, :]
    # slices = [(1,3), (3,5), (5,7), (7,9), (9,11)]
    # for ids, sl in enumerate(slices):
    #     filler[ids,]
    #

    return wav_ns, wav_we, kernel[[2,4,6,8], :, :]  #[1,4,6,10]


def cut_kernel_lsta(xpos, ypos, arr):

    dist = 200

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    kernel = kernel - np.mean(kernel)

    # if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
    #     print('Not enough valid in kernel')
    #     return

    if kernel.shape != (2*dist+1, 2*dist+1):
        print('Kernel wrong shape')
        return

    # plt.figure()
    # plt.imshow(kernel, origin='lower')

    return kernel


def get_previous_hours_msg(lsta_date):

    edate = lsta_date.replace(hour=12)  # make 12 reference hour for MCS filter

    t1 = edate - pd.Timedelta('1 hours')
    t2 = edate + pd.Timedelta('1 hours')

    file = cnst.MCS_15K# MCS_15K #_POINTS_DOM
    msg = xr.open_dataarray(file)
    try:
        msg = msg.sel(time=slice(t1.strftime("%Y-%m-%dT%H"), t2.strftime("%Y-%m-%dT%H")))
    except OverflowError:
        return None

    #print(prev_time.strftime("%Y-%m-%dT%H"), date.strftime("%Y-%m-%dT%H"))
    pos = np.where((msg.values <= -40) ) #(msg.values >= 5) & (msg.values < 65)) # #

    out = np.zeros_like(msg)
    out[pos] = 1
    out = np.sum(out, axis=0)
    out[out>0]=1

    msg = msg.sum(axis=0)*0

    xout = msg.copy()
    del msg
    xout.name = 'probs'
    xout.values = out

    return xout


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

    lsta = xr.open_dataset(cnst.AMSRE_ANO_DAY + 'sma_' + fdate + '.nc')
    lsta = lsta.sel(time=str(daybefore.year)+'-'+str(daybefore.month)+'-'+str(daybefore.day))
    lsta = lsta.sel(lon=slice(-13, 13), lat=slice(7, 23))
    print('Doing '+ 'AMSR_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    lsta_da = lsta['SM'].squeeze()

    topo = xr.open_dataset(cnst.WA_TOPO_3KM)
    topo = topo.sel(lon=slice(-12, 12), lat=slice(7.5, 22.5))

    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])


    try:
        lsta_da = topo.salem.transform(lsta_da)
    except RuntimeError:
        print('lsta_da on LSTA interpolation problem')
        return None


    lsta_da.values[ttopo.values >= 450] = np.nan
    lsta_da.values[gradsum > 30] = np.nan

    wav_input = lsta_da.values
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

    wavpos = wutil.applyHat_pure(wav_input, dataset='METSRFC_LS')

    wavarr = wavpos['coeffs'].copy()
    scales = wavpos['scales'].copy()

    for inds, wava in enumerate(wavarr):
        wava[inter1] = np.nan
        wavarr[inds, :, :] = (wava) / (np.nanstd(wava))

        #wavarr[inds, inter1[0], inter1[1]] = np.nan
    del wavpos
    wav_input[inter1] = np.nan

    # probs_msg = get_previous_hours_msg(daybefore)
    # probsm_on_lsta = topo.salem.transform(probs_msg, interp='nearest')
    #del probs_msg
    del topo


    rnspos = []
    rwepos = []
    rkernel = []
    lsta = []


    ###############################Blob loop

    nspos = []
    wepos = []
    vkernel = []

    for dids, dit in df.iterrows():

        lat = dit.lat
        lon = dit.lon

        try:
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.04)
        except KeyError:
            print('Nearest point finding error')
            continue
        #ipdb.set_trace()
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])+10
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            lsta_kernel = cut_kernel_lsta(xpos, ypos, wav_input)
        except TypeError:
            print('LSTA kernel error')
            continue

        # dist=100
        # if np.nansum(msg_kernel[dist-30:dist+30,dist-30:dist+67])>=2:   # filter out cases with MCSs at 12 [dist-50:dist+50,dist-30:dist+100]
        #     print('Meteosat MCS continue')
        #     continue


        try:
            wavpos_ns, wavpos_we, vvkernel = cut_kernel(xpos, ypos, wavarr, daybefore, plon, plat, rotate=False)
        except TypeError:
            print('Wavelt kernel error')
            continue


        nspos.append(wavpos_ns)  # north-south wavelet
        wepos.append(wavpos_we)  # west-east wavelet
        vkernel.append(vvkernel)
        lsta.append(lsta_kernel)

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

        for ry, rx in zip(randy, randx+10):

            if ry < 0:
                continue
            if ry > lsta_da.shape[0] - 1:
                continue

            if rx < 0:
                continue
            if rx > lsta_da.shape[1] - 1:
                continue

            try:
                lat = lsta_da['lat'][ry]
            except IndexError:
                ipdb.set_trace()
            try:
                lon = lsta_da['lon'][rx]
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



def plot(hour):
    name = 'UTC_15000_WAVELET_INIT2_NEWTRACKING'#'UTC_15000_AMSLST_3slot_noOverplotFilter'#'UTC_15000_-60_AMSRE'
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+name+".p", "rb"))

    scales = dic['scales']
    nbcores = dic['nbcores']
    nbrcores = dic['nbrcores']
    del dic['scales']
    del dic['nbcores']
    del dic['nbrcores']
    del dic['kernel']


    keys = list(dic.keys())
    cnt = (dic['SN-pos'][0]).shape[0]


    l=0

    snblob = (dic[keys[l]])[0] #/ (dic[keys[l]])[3]).T#-(dic[keys[l+2]])[0]
    snrandom = (dic[keys[l]])[1] #/ (dic[keys[l]])[4]).T#-(dic[keys[l+2]])[1]
    snmask = (dic[keys[l]])[2]#-(dic[keys[l+2]])[2]
    # snblob_std = (dic[keys[l]])[3]
    # snrandom_std = (dic[keys[l]])[4]

    weblob = (dic[keys[l+1]])[0] #/ (dic[keys[l+1]])[3]).T#-(dic[keys[l+3]])[0]
    werandom = (dic[keys[l+1]])[1] #/ (dic[keys[l+1]])[4]).T#-(dic[keys[l+3]])[1]
    wemask = (dic[keys[l+1]])[2]#-(dic[keys[l+3]])[2]
    # weblob_std = (dic[keys[l+1]])[3]
    # werandom_std = (dic[keys[l+1]])[4]

    l=2
    dist=200

    f = plt.figure(figsize=(9, 9))
    ax = f.add_subplot(221)
    levels= np.linspace(-0.5,0.5,20)#[-0.14, -0.12,-0.1, -0.08,-0.06, -0.04, -0.02,0.02,0.04,0.06,0.08,0.1, 0.12, 0.14]
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales ,  (snrandom) , cmap='RdBu_r', levels = levels, extend='both')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels = [0.5,1])

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(222)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3,scales,   (werandom) , cmap='RdBu_r', levels =levels, extend='both')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, wemask, colors='none', hatches='.', levels = [0.5,1])

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)


    ax = f.add_subplot(223)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales,  snblob-snrandom,  cmap='RdBu_r', levels =levels, extend='both') #, vmin = -0.1, vmax=0.1)
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, snmask, colors='none', hatches='.', levels=[0.5, 1],
                 linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(nbrcores) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(224)

    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, weblob-werandom    , cmap='RdBu_r', levels =levels, extend='both')
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales,wemask, colors='none', hatches='.', levels=[0.5, 1],
                 linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('West-East scales', fontsize=10)

    plt.tight_layout()
    #plt.savefig(path + '/amsreVSlsta/'+name+'_cross_AMSRE_INIT2_WAVELET_'+str(hour)+'.png')
    plt.show()
    #plt.close()



def plot_single(hour):
    name = 'UTC_15000_WAVELET_INIT_SCALES_NEWTRACKING'#'UTC_15000_AMSLST_3slot_noOverplotFilter'#'UTC_15000_-60_AMSRE'
    path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
    dic = pkl.load(open(path+"/coeffs_nans_stdkernel_USE_"+str(hour)+name+".p", "rb"))

    scales = dic['scales']
    nbcores = dic['nbcores']
    nbrcores = dic['nbrcores']
    del dic['scales']
    del dic['nbcores']
    del dic['nbrcores']
    del dic['kernel']


    keys = list(dic.keys())
    cnt = (dic['SN-pos'][0]).shape[0]
    levels = np.linspace(-0.4, 0.4, 20)

    l=0

    snblob = (dic[keys[l]])[0] #/ (dic[keys[l]])[3]).T#-(dic[keys[l+2]])[0]
    snrandom = (dic[keys[l]])[1] #/ (dic[keys[l]])[4]).T#-(dic[keys[l+2]])[1]
    snmask = (dic[keys[l]])[2]#-(dic[keys[l+2]])[2]
    # snblob_std = (dic[keys[l]])[3]
    # snrandom_std = (dic[keys[l]])[4]

    weblob = (dic[keys[l+1]])[0] #/ (dic[keys[l+1]])[3]).T#-(dic[keys[l+3]])[0]
    werandom = (dic[keys[l+1]])[1] #/ (dic[keys[l+1]])[4]).T#-(dic[keys[l+3]])[1]
    wemask = (dic[keys[l+1]])[2]#-(dic[keys[l+3]])[2]
    # weblob_std = (dic[keys[l+1]])[3]
    # werandom_std = (dic[keys[l+1]])[4]

    weblob = ndimage.gaussian_filter(weblob, 3, mode='nearest')

    l=2
    dist=200

    f = plt.figure(figsize=(7, 5), dpi=200)
    ax = f.add_subplot(111)

    cs = ax.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales, weblob , cmap='RdBu_r', levels =levels, extend='both')
    plt.colorbar(cs, label='Wavelet coefficient (-)')
    plt.contourf((np.arange(0, 2*dist+1) - dist) * 3, scales,wemask, colors='none', hatches='.', levels=[0.5, 1],
                 linewidth=0.25)

    ax.set_xlabel('West-east distance from core (km)')
    ax.set_ylabel('SMA scales')
    ax.set_yscale('log')
    ax.set_yticks([30,50,100,150,200,250,300])
    ax.set_yticklabels([30,50,100,150,200,250,300])

    plt.title(str(hour).zfill(2) + '00 UTC', fontsize=10)

    plt.tight_layout()
    plt.savefig(path + '/2hOverlap/'+name+'_cross_AMSRE_INIT2_WAVELET_'+str(hour)+'.png')
    plt.show()
    #plt.close()

def plot_trio():
    f = plt.figure(figsize=(7, 5))
    for ids,h in enumerate([15,17,19]):

        name = 'UTC_15000_WAVELET_INIT_SCALES_NEWTRACKING'  # 'UTC_15000_AMSLST_3slot_noOverplotFilter'#'UTC_15000_-60_AMSRE'
        path = cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients'
        dic = pkl.load(open(path + "/coeffs_nans_stdkernel_USE_" + str(h) + name + ".p", "rb"))

        scales = dic['scales']
        nbcores = dic['nbcores']
        nbrcores = dic['nbrcores']
        del dic['scales']
        del dic['nbcores']
        del dic['nbrcores']
        del dic['kernel']

        keys = list(dic.keys())
        cnt = (dic['SN-pos'][0]).shape[0]
        levels = np.linspace(-0.4, 0.4, 20)

        l = 0

        snblob = (dic[keys[l]])[0]  # / (dic[keys[l]])[3]).T#-(dic[keys[l+2]])[0]
        snrandom = (dic[keys[l]])[1]  # / (dic[keys[l]])[4]).T#-(dic[keys[l+2]])[1]
        snmask = (dic[keys[l]])[2]  # -(dic[keys[l+2]])[2]
        # snblob_std = (dic[keys[l]])[3]
        # snrandom_std = (dic[keys[l]])[4]

        weblob = (dic[keys[l + 1]])[0]  # / (dic[keys[l+1]])[3]).T#-(dic[keys[l+3]])[0]
        werandom = (dic[keys[l + 1]])[1]  # / (dic[keys[l+1]])[4]).T#-(dic[keys[l+3]])[1]
        wemask = (dic[keys[l + 1]])[2]  # -(dic[keys[l+3]])[2]
        # weblob_std = (dic[keys[l+1]])[3]
        # werandom_std = (dic[keys[l+1]])[4]

        l = 2
        dist = 200


        ax = f.add_subplot(1,3,ids+1)

        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, weblob, cmap='RdBu_r', levels=levels,
                     extend='both')
        plt.colorbar(label='Power difference (Blob-random)')
        plt.contourf((np.arange(0, 2 * dist + 1) - dist) * 3, scales, wemask, colors='none', hatches='.',
                     levels=[0.5, 1],
                     linewidth=0.25)

        ax.set_xlabel('West-east distance from core (km)')
        ax.set_ylabel('SMA scales')

        plt.title(str(h).zfill(2)+'00 UTC', fontsize=10)

    plt.tight_layout()
    #plt.savefig(path + '/amsreVSlsta/'+name+'_cross_AMSRE_INIT2_WAVELET_TRIO.png')
    plt.show()