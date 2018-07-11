# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import pdb
import pandas as pd
import salem
from utils import u_met, u_parallelise, u_gis, u_arrays, constants, u_grid
from scipy.interpolate import griddata

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def diurnal_loop():

    for l in np.arange(17,23):
        print('Doing '+str(l))
        composite(l)

def composite(h):
    #pool = multiprocessing.Pool(processes=8)


    file = constants.MCS_POINTS_DOM

    hour = h

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == hour ) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >=6) ]

    msg = msg.sel(lat=slice(10.9,19), lon=slice(-9.8,9.8))

    dic = u_parallelise.run_arrays(7,file_loop,msg,['ano', 'regional', 'cnt',  'prob', 'cntp']) #'rano', 'rregional', 'rcnt',

    for k in dic.keys():
       dic[k] = np.nansum(dic[k], axis=0)


    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/cm_prevday_-70/composite_backtrack_CMORPH_"+str(hour).zfill(2)+".p", "wb"))


def cut_kernel(xpos, ypos, arr, date, lon, lat, t, parallax=False, rotate=False, probs=False):

    if parallax:
        km, coords = u_gis.call_parallax_era(date.month, t, lon, lat, 0, 0)
        lx, ly = km

        lx = int(np.round(lx / 3.))
        ly = int(np.round(ly / 3.))  # km into pixels
        xpos = xpos - lx
        ypos = ypos - ly

    dist = 100

    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    if np.sum(probs) > 0:
        prob = u_arrays.cut_kernel(probs,xpos, ypos,dist)
        cnt2 = np.zeros_like(kernel)
        cnt2[np.isfinite(prob)] = 1
    else:
        prob = 0
        cnt2 = 0

    if rotate:
        kernel = u_met.era_wind_rotate(kernel,date,lat,lon,level=700, ref_angle=90)

    if (np.sum(np.isfinite(kernel)) < 0.01 * kernel.size):
        return

    kernel3 = kernel - np.nanmean(kernel)

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1



    if kernel.shape != (201, 201):
        pdb.set_trace()



    return kernel, kernel3, cnt, prob, cnt2


def get_previous_hours(date):


    tdic = {18 : ('36 hours', '15 hours'),
            21: ('39 hours', '18 hours'),
            0: ('42 hours', '21 hours'),
            3: ('45 hours', '24 hours'),
            6: ('48 hours', '27 hours')}
    before = pd.Timedelta(tdic[date.hour][0])
    before2 = pd.Timedelta(tdic[date.hour][1])

    t1 = date - before
    t2 = date - before2

    file = constants.CMORPH
    try:
        cmm = xr.open_dataarray(file + 'CMORPH_WA_onLSTA' + str(date.year) + '.nc')
    except:
        return None
    cmm = cmm.sel(lat=slice(10.9,19), lon=slice(-9.8,9.8), time=slice(t1, t2)).sum(dim='time') #
    # cmm = cmm.sel(lat=slice(10.9, 19), lon=slice(-9.8, 9.8))
    # posi = np.where(pd.to_datetime(cmm.time.values) == date)
    #
    # cmm = cmm.isel(time=slice(posi[0][0]-2,posi[0][0]))
    cm = cmm
    pos = np.where(cm.values>=10) #(msg.values >= 5) & (msg.values < 65)) # #

    out = np.zeros_like(cm)
    out[pos] = 1
    #out = np.sum(out, axis=0) / out.shape[0]



    xout = cm.copy()
    xout.name = 'probs'
    xout.values = out

    return xout



def file_loop(fi):
    print('Doing day: ', fi.time.values)

    date = pd.Timestamp(fi.time.values)

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
    if fdate == '20060708':
        return None
    try:
        lsta = xr.open_dataset(constants.LSTA_NEW + 'lsta_daily_' + fdate + '.nc')
    except OSError:
        return None
    print('Doing '+ 'lsta_daily_' + fdate + '.nc')

    topo = xr.open_dataset(constants.LSTA_TOPO)
    ttopo = topo['h']
    grad = np.gradient(ttopo.values)
    gradsum = abs(grad[0])+abs(grad[1])


    lsta_da = lsta['LSTA'].squeeze()
    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.10:
        print('Not enough valid')
        return None

    lsta_da.values[ttopo.values>=450] = np.nan
    lsta_da.values[gradsum>30] = np.nan
    pos = np.where((fi.values >= 5) & (fi.values < 65) )  #(fi.values >= 5) & (fi.values < 65)

    if (np.sum(pos) == 0) | (len(pos[0]) < 3):
        print('No blobs found')
        return None

    kernel2_list = []
    kernel3_list = []
    cnt_list = []
    cntp_list = []
    prkernel_list = []

    probs = get_previous_hours(date)

    try:
        probs_on_lsta = lsta.salem.transform(probs)
    except RuntimeError:
        return None

    f = plt.figure()
    ax = f.add_subplot(111)

    plt.contourf(probs_on_lsta, cmap='viridis',extend='both')
    plt.colorbar(label='K')
    plt.savefig('/users/global/cornkle/test/' + fdate + '_single.png')
    plt.close()



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
            kernel2, kernel3, cnt, prkernel, cntp = cut_kernel(xpos, ypos, lsta_da, daybefore.month, plon, plat, -40, parallax=False, rotate=False, probs=probs_on_lsta)
        except TypeError:
            continue


        kernel2_list.append(kernel2)
        kernel3_list.append(kernel3)
        prkernel_list.append(prkernel)
        cnt_list.append(cnt)
        cntp_list.append(cnt)

    if kernel2_list == []:
        return None

    if len(kernel2_list) == 1:
      return None
    else:

        kernel2_sum = np.nansum(np.stack(kernel2_list, axis=0), axis=0)
        kernel3_sum = np.nansum(np.stack(kernel3_list, axis=0), axis=0)
        cnt_sum = np.nansum(np.stack(cnt_list, axis=0), axis=0)
        cntp_sum = np.nansum(np.stack(cntp_list, axis=0), axis=0)
        pr_sum = np.nansum(np.stack(prkernel_list, axis=0), axis=0)


    print('Returning')

    return (kernel2_sum, kernel3_sum, cnt_sum, pr_sum, cntp_sum)


def plot_gewex(h):
    hour=h
    dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/cm_prevday_-70/composite_backtrack_CMORPH_"+str(hour).zfill(2)+".p", "rb"))

    extent = (dic['ano'].shape[1]-1)/2

    f = plt.figure(figsize=(7, 5))
    ax = f.add_subplot(111)

    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu_r',  levels=[ -0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5], extend='both') #-(rkernel2_sum / rcnt_sum)
    #plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    #pdb.set_trace()

    contours = plt.contourf((dic['prob']/ dic['cntp'])*100, extend='both', levels=np.arange(10,41,2)) # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.0f')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title(str(hour).zfill(2)+'00 UTC | '+str(np.max(dic['cnt']))+' cores', fontsize=17)


    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/cm_prevday_-70/'+str(hour).zfill(2)+'_single.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_all():

    hours = [16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]

    for h in hours:
        plot_gewex(h)

