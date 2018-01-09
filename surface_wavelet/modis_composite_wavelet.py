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
from wavelet import util
from utils import u_plot
from scipy.stats import ttest_ind as ttest
from scipy.interpolate import griddata
import pickle as pkl
import collections

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
matplotlib.rcParams['hatch.linewidth'] = 0.1

def run_hours():

    l = [17,18,19,0,1,2,3]
    for ll in l:
        composite(ll)


def composite(hour):
    pool = multiprocessing.Pool(processes=8)

    nightp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_0-3UTC_centrePoint.nc'
    dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'
    path = '/users/global/cornkle/figs/LSTA-bullshit/scales/picks'
    # nightp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_0-3UTC.nc'
    # dayp = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/blob_map_35km_-70_15-18UTC.nc'

    hour = hour

    if hour > 6:
        file = dayp
    else:
        file = nightp

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 8) ]

    msg = msg.sel(lat=slice(11,18), lon=slice(-9, 9))

    res = pool.map(file_loop, msg)
    pool.close()

    res = [x for x in res if x is not None]

    snpos_list = []
    wepos_list = []
    snneg_list = []
    weneg_list = []
    snboth_list = []
    weboth_list = []

    rsnpos_list = []
    rwepos_list = []
    rsnneg_list = []
    rweneg_list = []
    rsnboth_list = []
    rweboth_list = []

    for r in res:
        snpos_list.append(np.squeeze(r[0]))
        wepos_list.append(np.squeeze(r[1]))
        snneg_list.append(np.squeeze(r[2]))
        weneg_list.append(np.squeeze(r[3]))
        snboth_list.append(np.squeeze(r[4]))
        weboth_list.append(np.squeeze(r[5]))

        rsnpos_list.append(np.squeeze(r[6]))
        rwepos_list.append(np.squeeze(r[7]))
        rsnneg_list.append(np.squeeze(r[8]))
        rweneg_list.append(np.squeeze(r[9]))
        rsnboth_list.append(np.squeeze(r[10]))
        rweboth_list.append(np.squeeze(r[11]))

        scales = r[12]


        dic = collections.OrderedDict([('SN-pos' , [snpos_list, rsnpos_list]),
                                       ('WE-pos' , [wepos_list, rwepos_list]),
                                       ('SN-neg' , [snneg_list, rsnneg_list]),
                                       ('WE-neg' , [weneg_list, rweneg_list]),
                                       ('SN-both', [snboth_list, snboth_list]),
                                       ('WE-both', [weboth_list, weboth_list]),
                                       ('scales' , scales)])

    keys = list(dic.keys())
    for l in keys:
        (dic[l])[0] = np.squeeze(np.vstack((dic[l])[0]))
        (dic[l])[1] = np.squeeze(np.vstack((dic[l])[1]))

    for l in keys:
        nsstat, nspvalue = ttest(dic[l][0], dic[l][1], axis=0, equal_var=False, nan_policy='omit')
        mask = nspvalue < 0.05
        dic[l].append(mask)

    for l in keys:
        (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
        (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)

    pkl.dump(dic, open(path+"/composite_pos_neg_both_"+str(hour)+"UTC.p", "wb"))



def plot(hour):

    path = '/users/global/cornkle/figs/LSTA-bullshit/scales/picks'
    dic = pkl.load(open(path+"/composite_pos_neg_both_"+str(hour)+"UTC.p", "rb"))
    scales = dic['scales']
    del dic['scales']

    keys = list(dic.keys())
    cnt = (dic['SN-pos'][0]).shape[0]

    for l in [0,2,4]:
        f = plt.figure(figsize=(9, 4))
        ax = f.add_subplot(121)

        snblob = (dic[keys[l]])[0]
        snrandom = (dic[keys[l]])[1]
        snmask = (dic[keys[l]])[2]

        weblob = (dic[keys[l+1]])[0]
        werandom = (dic[keys[l+1]])[1]
        wemask = (dic[keys[l+1]])[2]

        plt.contourf((np.arange(0, 101) - 50) * 3, scales ,snblob - snrandom  , cmap='RdBu_r', vmax=0.05, vmin=-0.05)
        plt.colorbar(label='Power difference (Blob-random)')
        plt.contourf((np.arange(0, 101) - 50) * 3, scales, snmask, colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

        ax.set_xlabel('km')
        ax.set_ylabel('Scales')

        plt.title('South-North scales, Nb cores: ' + str(cnt) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
                  fontsize=10)

        ax = f.add_subplot(122)

        plt.contourf((np.arange(0, 101) - 50) * 3,scales,  weblob- werandom  , cmap='RdBu_r', vmax=0.05, vmin=-0.05)
        plt.colorbar(label='Power difference (Blob-random)')
        plt.contourf((np.arange(0, 101) - 50) * 3, scales, wemask, colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

        ax.set_xlabel('km')
        ax.set_ylabel('km')

        plt.title('West-East scales', fontsize=10)

        plt.tight_layout()
        plt.savefig(path + '/lsta_hours_'+keys[l]+'_'+str(hour)+'.png')


def cut_kernel(xpos, ypos, wl):

    dist = 50
    #pdb.set_trace()
    try:
        kernel = wl[:, ypos - dist: ypos + dist + 1,  xpos - dist: xpos + dist + 1]

    except ValueError:
        # print('Point at boundary, pass')
        return

    if kernel.shape != (np.shape(wl)[0], 101, 101):
        return

    wav_ns = kernel[:,:,50]
    wav_we = kernel[:,50,:]

    return wav_ns, wav_we


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

    pos = np.where((fi.values >= 1) & (fi.values <= 30))
    if np.sum(pos) == 0:
        print('No blobs found')
        return

    wav_input = lsta_day2-lsta_day2.mean()
    wav_input = wav_input.values
    points = np.where(np.isfinite(wav_input))
    inter = np.where(np.isnan(wav_input))
    # interpolate over sea from land points
    try:
         wav_input[inter] = griddata(points, np.ravel(wav_input[points]), inter, method='nearest')
    except ValueError:
        wav_input[inter]=0
    wavpos = util.waveletLSTA(np.squeeze(wav_input), 3, method='pos')
    wavneg = util.waveletLSTA(np.squeeze(wav_input), 3, method='neg')
    wavboth = util.waveletLSTA(np.squeeze(wav_input), 3, method=None)

    wlpos = wavpos['power']
    wlneg = wavneg['power']
    wlboth = wavboth['power']

    wlpos[:,inter[0], inter[1]]=np.nan
    wlneg[:,inter[0], inter[1]]=np.nan
    wlboth[:,inter[0], inter[1]]=np.nan

    xfi = fi.shape[1]
    yfi = fi.shape[0]
    randx = np.random.randint(0, xfi, 250)
    randy = np.random.randint(0, yfi, 250)
    posr = (randy, randx)

##############################Random loop
    rnspos_list = []
    rwepos_list = []
    rnsneg_list = []
    rweneg_list = []
    rnsboth_list = []
    rweboth_list = []

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
            wavpos_ns, wavpos_we = cut_kernel(xpos, ypos, wlpos)
            wavneg_ns, wavneg_we = cut_kernel(xpos, ypos, wlneg)
            wavboth_ns, wavboth_we = cut_kernel(xpos, ypos, wlboth)
        except TypeError:
            continue

        rnspos_list.append(wavpos_ns)  # north-south wavelet
        rwepos_list.append(wavpos_we)  # west-east wavelet
        rnsneg_list.append(wavneg_ns)  # north-south wavelet
        rweneg_list.append(wavneg_we)  # west-east wavelet
        rnsboth_list.append(wavboth_ns)  # north-south wavelet
        rweboth_list.append(wavboth_we)  # west-east wavelet

###############################Blob loop
    nspos_list = []
    wepos_list = []
    nsneg_list = []
    weneg_list = []
    nsboth_list = []
    weboth_list = []

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
            wavpos_ns, wavpos_we = cut_kernel(xpos, ypos, wlpos)
            wavneg_ns, wavneg_we = cut_kernel(xpos, ypos, wlneg)
            wavboth_ns, wavboth_we = cut_kernel(xpos, ypos, wlboth)
        except TypeError:
            continue

        nspos_list.append(wavpos_ns)  # north-south wavelet
        wepos_list.append(wavpos_we)  # west-east wavelet
        nsneg_list.append(wavneg_ns)  # north-south wavelet
        weneg_list.append(wavneg_we)  # west-east wavelet
        nsboth_list.append(wavboth_ns)  # north-south wavelet
        weboth_list.append(wavboth_we)  # west-east wavelet

    if nspos_list == []:
        return None

    if len(nspos_list) == 1:
      return None
    else:
        nspos_sum = np.squeeze(np.stack(nspos_list, axis=0))
        wepos_sum = np.squeeze(np.stack(wepos_list, axis=0))
        nsneg_sum = np.squeeze(np.stack(nsneg_list, axis=0))
        weneg_sum = np.squeeze(np.stack(weneg_list, axis=0))
        nsboth_sum = np.squeeze(np.stack(nsboth_list, axis=0))
        weboth_sum = np.squeeze(np.stack(weboth_list, axis=0))

        rnspos_sum = np.squeeze(np.stack(rnspos_list, axis=0))
        rwepos_sum = np.squeeze(np.stack(rwepos_list, axis=0))
        rnsneg_sum = np.squeeze(np.stack(rnsneg_list, axis=0))
        rweneg_sum = np.squeeze(np.stack(rweneg_list, axis=0))
        rnsboth_sum = np.squeeze(np.stack(rnsboth_list, axis=0))
        rweboth_sum = np.squeeze(np.stack(rweboth_list, axis=0))

    scales = wavpos['scales']

    print('Returning')

    return (nspos_sum, wepos_sum,  nsneg_sum, weneg_sum,  nsboth_sum, weboth_sum,
            rnspos_sum, rwepos_sum, rnsneg_sum, rweneg_sum, rnsboth_sum, rweboth_sum,
            scales)
# rcns_sum, rcwe_sum, cns_sum, cwe_sum,


if __name__ == "__main__":
    run_hours()
