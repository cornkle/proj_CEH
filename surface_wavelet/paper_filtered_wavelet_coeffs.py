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
    # msgopen = pd.read_csv(
    #     cnst.network_data + 'figs/LSTA/corrected_LSTA/new/wavelet_coefficients/core_txt/cores_gt15000km2_table_1640_580_' + str(
    #         hour) + '.csv')

    msg = pd.DataFrame.from_dict(msgopen)

    msg['date'] = pd.to_datetime(msg[['year','month','day']])
    print('Start core number ', len(msg))

    msg = msg[(msg['LSTAslotfrac']>=0.025) & (msg['dtime']<=2) & (np.isfinite(msg['SMmean0'])) ] # & (msg['SMmean0']>1) & (np.isfinite(msg['SMmean0']))
    #msg = msg[(msg['LSTAslotfrac'] >= 0.5) & (msg['dtime'] <= 2) & (np.isfinite(msg['SMmean0']))]
    msgin = msg[(msg['lat']>8.5) & (msg['lat']<20.5) & (msg['topo']<450)]#[msg['initTime']<=3]#[msg['SMwet']==2]

    print('Number of cores', len(msgin))
    ipdb.set_trace()

    # calculate the chunk size as an integer
    #'chunk_size = int(msg.shape[0] / pnumber)
    msgin.sort_values(by='date')

    msgy = msgin#[msgin['year']==year]

    chunk, chunk_ind, chunk_count = np.unique(msgy.date, return_index=True, return_counts=True)

    chunks = [msgy.loc[msgy.index[ci:ci + cc]] for ci, cc in zip(chunk_ind, chunk_count)] # daily chunks

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

    slot_list = []

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

        slot_list.append(r[11])

    dic = collections.OrderedDict([('SN-pos' , [snpos_list, rsnpos_list]),
                                   ('WE-pos' , [wepos_list, rwepos_list]),
                                   ('kernel' , [vkernel_list, rkernel_list]),
                                   ('lsta' , [lsta_list, lsta_cnt]),
                                   ('cnt' , [vkernel_cnt, rkernel_cnt]),
                                   ('scales' , scales),
                                  ('slots', [slot_list])])

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
        if (l == 'scales') | (l in ['lsta', 'slots']):
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
    pkl.dump(dic, open(outpath+"coeffs_nans_stdkernel_USE_"+str(hour)+"UTC_15000_LSTAA_1330_" + key + ".p", "wb"))
    print('Save file written!')



def cut_kernel(xpos, ypos, arr, date, lat, lon, rotate=False):

    dist = 100

    kernel = u_arrays.cut_kernel_3d(arr,xpos, ypos,dist)

    # if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
    #     return

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


def cut_kernel_lsta(xpos, ypos, arr, msg=None, nbslot=None):

    dist = 100

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    kernel = kernel

    # if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
    #     return

    if kernel.shape != (2*dist+1, 2*dist+1):
        return

    if msg is not False:
        msgkernel = u_arrays.cut_kernel(msg, xpos, ypos, dist)
    else:
        msgkernel = kernel*np.nan

    if nbslot is not False:
        nbkernel = u_arrays.cut_kernel(nbslot, xpos, ypos, dist)
    else:
        nbkernel = kernel*np.nan

    return kernel, msgkernel, nbkernel


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
    # if np.sum(out>1) != 0:
    #     'Stop!!!'
    #     pdb.set_trace()

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

    try:
        lsta = xr.open_dataset(cnst.LSTA_1330 + 'lsta_daily_' + fdate + '.nc')
    except OSError:
        return None
    print('Doing ' + 'lsta_daily_' + fdate + '.nc')

    topo = xr.open_dataset(cnst.LSTA_TOPO)
    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0]) + abs(grad[1])

    lsta_da = lsta['LSTA'].squeeze()
    slot_da = lsta['NbSlot'].squeeze().values


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

    wavpos = wutil.applyHat_pure(wav_input, dataset='METSRFC')

    wavarr = wavpos['coeffs'].copy()
    scales = wavpos['scales'].copy()

    for inds, wava in enumerate(wavarr):
        wava[inter1] = np.nan

        wavarr[inds, :, :] = (wava) / (np.nanstd(wava))

        wavarr[inds, inter1[0], inter1[1]] = np.nan
    del wavpos
    wav_input[inter1] = np.nan


    # probs_msg = get_previous_hours_msg(daybefore)
    # probsm_on_lsta = topo.salem.transform(probs_msg, interp='nearest')
    del topo

    rnspos = []
    rwepos = []
    rkernel = []
    lsta = []
    slot = []

    ###############################Blob loop

    nspos = []
    wepos = []
    vkernel = []

    for dids, dit in df.iterrows():

        lat = dit.lat
        lon = dit.lon

        #initiation filter:
        # initpath = cnst.network_data + 'data/OBS/MSG_WA30/track_back_cores_vn1_'+str(hour)+'Z.txt'
        # if os.path.isfile(initpath):
        #     dic = pd.read_table(initpath, delim_whitespace=True, header=None,
        #                         names=['year', 'mon', 'day', 'i_core', 'j_core', 'i_initiation', 'j_initiation',
        #                                'core_time', 'initiation_time'])
        #     ddic = dic[
        #         (dic['i_core'] == dit['xloc']) & (dic['j_core'] == dit['yloc']) & (dic['year'] == dit['year']) & (
        #                     dic['mon'] == dit['month']) & (dic['day'] == dit['day'])]

            #ipdb.set_trace()
            # if len(ddic) == 0:
            #     continue
            #ipdb.set_trace()
            # if (ddic['initiation_time'].values >3 ) | (ddic['initiation_time'].values < 0  ): #& (ddic['initiation_time'].values >= 12):
            #     continue

            # if (ddic['initiation_time'].values <=12 ) | (ddic['initiation_time'].values >= hour ): #& (ddic['initiation_time'].values >= 12):
            #     continue

            ## init location filter

            # min_init = lon + 3
            # is_init = mlon[ddic['i_core']]
            #
            # print('Init location', is_init, 'Min location', min_init, 'init_hour', ddic['initiation_time'])

            # if (dit['xdiff'] < 100) & (dit['initTime']!=2):
            #     print('Initiation point too close')
            #     continue

            # if np.abs((dit['xdiff']) > 33) | (dit['xinit'] < 0):
            #     print('Initiation point too far', dit['initTime'])
            #     continue



        try:
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.04)
        except KeyError:
            continue
        #ipdb.set_trace()
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_da['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            lsta_kernel, msg_kernel, nbkernel = cut_kernel_lsta(xpos, ypos, lsta_da.values,msg=False, nbslot=slot_da)
        except TypeError:
            print('Kernel error')
            continue
        dist=100

        # if np.nansum(msg_kernel[dist-30:dist+30,dist-30:dist+67])>=2:   # filter out cases with MCSs at 12 [dist-50:dist+50,dist-30:dist+100]
        #     print('Meteosat MCS continue')
        #     continue


        try:
            wavpos_ns, wavpos_we, vvkernel = cut_kernel(xpos, ypos, wavarr, daybefore, plon, plat, rotate=False)
        except TypeError:
            print('Kernel error')
            continue


        nspos.append(wavpos_ns)  # north-south wavelet
        wepos.append(wavpos_we)  # west-east wavelet
        vkernel.append(vvkernel)
        lsta.append(lsta_kernel)
        slot.append(nbkernel)

        ##### random

        rdist = 1.5
        y = lat
        x = lon
        randy50 = [y-rdist, y-rdist, y-rdist, y, y, y+rdist, y+rdist, y+rdist]
        randx50 = [x-rdist, x, x+rdist, x-rdist, x+rdist, x-rdist, x, x + rdist]
        randy50_100 = [y - rdist,  y - rdist, y, y, y + rdist, y + rdist]

        rdist = 3
        randx100 = [x-rdist,  x+rdist, x-rdist, x+rdist, x-rdist, x + rdist]

        rdist = 4
        randx150 = [x-rdist,  x+rdist, x-rdist, x+rdist, x-rdist, x + rdist]


        randy = np.array(randy50 + randy50_100 + randy50_100)
        randx = np.array(randx50 + randx100 + randx150)
        #ipdb.set_trace()

        for ry, rx in zip(randy, randx):

            try:
                point = lsta_da.sel(lat=ry, lon=rx, method='nearest', tolerance=0.04)
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
        slot_sum = np.nansum(np.stack(slot, axis=0), axis=0)[np.newaxis, ...]


        rnspos_sum = np.squeeze(np.stack(rnspos, axis=0))
        rwepos_sum = np.squeeze(np.stack(rwepos, axis=0))
        rkernel_sum = np.nansum(np.stack(rkernel, axis=0), axis=0)[np.newaxis,...]
        rkernel_cnt = np.sum(np.isfinite(np.stack(rkernel, axis=0)), axis=0)[np.newaxis,...]


    print('Returning with kernel, success!!')

    return (nspos_sum, wepos_sum,
            rnspos_sum, rwepos_sum, vkernel_sum, rkernel_sum, scales, vkernel_cnt, rkernel_cnt, lsta_sum, lsta_cnt, slot_sum)
# rcns_sum, rcwe_sum, cns_sum, cwe_sum,


if __name__ == "__main__":
    run_hours()
