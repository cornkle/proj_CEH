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

    file = '/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'
    #file = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'
    path = '/users/global/cornkle/figs/LSTA-bullshit/scales/new'

    hour = hour

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] ==hour) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 8) ]

    msg = msg.sel(lat=slice(10,18), lon=slice(-10, 10))

    res = pool.map(file_loop, msg)
    pool.close()

    # for m in msg[2:5]:
    #     file_loop(m)
    # return
    res = [x for x in res if x is not None]

    snpos_list = []
    wepos_list = []

    rsnpos_list = []
    rwepos_list = []


    for r in res:
        snpos_list.append(np.squeeze(r[0]))
        wepos_list.append(np.squeeze(r[1]))

        rsnpos_list.append(np.squeeze(r[2]))
        rwepos_list.append(np.squeeze(r[3]))


        scales = r[4]


        dic = collections.OrderedDict([('SN-pos' , [snpos_list, rsnpos_list]),
                                       ('WE-pos' , [wepos_list, rwepos_list]),
                                       ('scales' , scales)])

    keys = list(dic.keys())

    for l in keys:
        if l == 'scales':
            continue
        (dic[l])[0] = np.squeeze(np.vstack((dic[l])[0]))
        (dic[l])[1] = np.squeeze(np.vstack((dic[l])[1]))

    for l in keys:
        if l == 'scales':
            continue
        nsstat, nspvalue = ttest(dic[l][0], dic[l][1], axis=0, equal_var=False, nan_policy='omit')
        mask = nspvalue < 0.05
        dic[l].append(mask)

    for l in keys:
        if l == 'scales':
            continue

        (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
        (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)


    pkl.dump(dic, open(path+"/composite_wet_"+str(hour)+"UTC.p", "wb"))



def plot(hour):

    path = '/users/global/cornkle/figs/LSTA-bullshit/scales/new'
    dic = pkl.load(open(path+"/composite_wet_"+str(hour)+"UTC.p", "rb"))


    scales = dic['scales']
    del dic['scales']

    keys = list(dic.keys())
    cnt = (dic['SN-pos'][0]).shape[0]

    # for l in keys:
    #     if keys == 'scales':
    #         continue
    #     (dic[l])[0] = np.nanmean((dic[l])[0], axis=0)
    #     (dic[l])[1] = np.nanmean((dic[l])[1], axis=0)

    l=0
    f = plt.figure(figsize=(9, 4))
    ax = f.add_subplot(121)

    snblob = (dic[keys[l]])[0]
    snrandom = (dic[keys[l]])[1]
    snmask = (dic[keys[l]])[2]

    weblob = (dic[keys[l+1]])[0]
    werandom = (dic[keys[l+1]])[1]
    wemask = (dic[keys[l+1]])[2]

    plt.contourf((np.arange(0, 101) - 50) * 3, scales ,snblob - snrandom   , cmap='RdBu_r', vmin=-0.05, vmax=0.05)
    plt.colorbar(label='Power difference (Blob-random)')
    plt.contourf((np.arange(0, 101) - 50) * 3, scales, snmask, colors='none', hatches='.', levels = [0.5,1], linewidth=0.25)

    ax.set_xlabel('km')
    ax.set_ylabel('Scales')

    plt.title('South-North scales, Nb cores: ' + str(cnt) + '| ' + str(hour).zfill(2) + '00UTC, Jul-Sep',
              fontsize=10)

    ax = f.add_subplot(122)

    plt.contourf((np.arange(0, 101) - 50) * 3,scales,  weblob - werandom  , cmap='RdBu_r', vmin=-0.05, vmax=0.05)
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
        print('Point at boundary, pass')
        return

    if kernel.shape != (np.shape(wl)[0], 101, 101):

        kernel =  np.zeros([np.shape(wl)[0], 101, 101])*np.nan


        if xpos - dist >= 0:
            xmin = 0
            xmindist = dist
        else:
            xmin = (xpos - dist)*-1
            xmindist = dist + (xpos - dist)

        if ypos - dist >= 0:
            ymin = 0
            ymindist = dist
        else:
            ymin = (ypos - dist)*-1
            ymindist = dist + (ypos - dist)

        if xpos + dist < wl.shape[2]:
            xmax = kernel.shape[2]
            xmaxdist = dist + 1
        else:
            xmax = dist -  (xpos - wl.shape[2])
            xmaxdist = dist - ( xpos + dist - wl.shape[2])

        if ypos + dist < wl.shape[1]:
            ymax = kernel.shape[1]
            ymaxdist = dist +1
        else:
            ymax = dist - (ypos - wl.shape[1])
            ymaxdist = dist - (ypos + dist - wl.shape[1])

        cutk = wl[:, ypos - ymindist: ypos + ymaxdist  , xpos - xmindist: xpos + xmaxdist]

        kernel[:, ymin : ymax,  xmin:xmax] = cutk
        #print(ypos, xpos, wl.shape)
        # f = plt.figure()
        # plt.imshow(kernel[0,:,:])
        # plt.plot(51,51, 'ro')

        # print('Point at boundary, pass')
        # return

    wav_ns = np.nanmean(kernel[:,:,49:51], axis=2) #kernel[:,:,50] #
    wav_we = np.nanmean(kernel[:,49:51,:], axis=1) #kernel[:,50,:] #

    return wav_ns, wav_we


def file_loop(fi):


    date = pd.to_datetime(str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))
    dayd = pd.Timedelta('1 days')
    if fi['time.hour'].values.size != 1:
        'hour array too big, problem!'

    if (fi['time.hour'].values) <= 7:
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

    pos = np.where((fi.values >=1) & (fi.values <=30)) #< -40 )) (fi.values >=1) & (fi.values <=20)
    if np.sum(pos) == 0:
        print('No blobs found')
        return

    wav_input = np.squeeze(lsta_day2.values)
    points = np.where(np.isfinite(wav_input))
    inter1 = np.where(np.isnan(wav_input))
    # interpolate over sea from land points
    #wav_input[inter1] = 0

    try:
         wav_input[inter1] = griddata(points, np.ravel(wav_input[points]), inter1, method='linear')
    except ValueError:
        pass

    inter = np.where(np.isnan(wav_input))
    try:
         wav_input[inter] = griddata(points, np.ravel(wav_input[points]), inter, method='nearest')
    except ValueError:
        wav_input[inter]=0

    wavpos = util.waveletLSTA(wav_input, 3, method='wet')

    wlpos = wavpos['power']

    wlpos[:,inter1[0], inter1[1]]=np.nan
    #wlpos[np.abs(wlpos)<3] = np.nan

    xfi = fi.shape[1]
    yfi = fi.shape[0]
    randx = np.random.randint(0, xfi, 250)
    randy = np.random.randint(0, yfi, 250)
    posr = (randy, randx)

##############################Random loop
    rnspos_list = []
    rwepos_list = []

    for y, x in zip(posr[0], posr[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        point = lsta_day2.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_day2['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_day2['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            wavpos_ns, wavpos_we = cut_kernel(xpos, ypos, wlpos)
        except TypeError:
            print('Kernel random error')
            continue

        rnspos_list.append(wavpos_ns)  # north-south wavelet
        rwepos_list.append(wavpos_we)  # west-east wavelet

###############################Blob loop
    nspos_list = []
    wepos_list = []


    for y, x in zip(pos[0], pos[1]):


        lat = fi['lat'][y]
        lon = fi['lon'][x]

        point = lsta_day2.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta_day2['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta_day2['lat'].values == plat)
        ypos = int(ypos[0])


        try:
            wavpos_ns, wavpos_we = cut_kernel(xpos, ypos, wlpos)
        except TypeError:
            print('Kernel error')
            continue

        nspos_list.append(wavpos_ns)  # north-south wavelet
        wepos_list.append(wavpos_we)  # west-east wavelet


    if nspos_list == []:
        print('NorthSouth empty')
        return None

    if len(nspos_list) == 1:
        print('NorthSouth 1')
        return None
    else:
        nspos_sum = np.squeeze(np.stack(nspos_list, axis=0))
        wepos_sum = np.squeeze(np.stack(wepos_list, axis=0))


        rnspos_sum = np.squeeze(np.stack(rnspos_list, axis=0))
        rwepos_sum = np.squeeze(np.stack(rwepos_list, axis=0))

    scales = wavpos['scales']

    print('Returning')

    return (nspos_sum, wepos_sum,
            rnspos_sum, rwepos_sum,
            scales)
# rcns_sum, rcwe_sum, cns_sum, cwe_sum,


if __name__ == "__main__":
    run_hours()
