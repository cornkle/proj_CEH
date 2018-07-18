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
from collections import OrderedDict
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


    file = constants.MCS_CENTRE70_SMALL

    hour = h

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == 17 ) & ((msg['time.minute'] == 0) & (
        msg['time.year'] >= 2008) & (msg['time.year'] <= 2010) & (msg['time.month'] >=6)) ]

    msg = msg.sel(lat=slice(10.9,19), lon=slice(-9.8,9.8))

    dic = u_parallelise.run_arrays(5,file_loop,msg[0:40], ['lsta', 'cnt', 'cntp', 't', 'u', 'v', 'shear']) #'rano', 'rregional', 'rcnt',

    # for k in dic.keys():
    #    dic[k] = np.nansum(dic[k], axis=0)


    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/ERA/composite_backtrack_ERA"+str(hour).zfill(2)+".p", "wb"))


def cut_kernel(xpos, ypos, arr, date, lon, lat, t, dist, probs=False):


    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)


    t = u_arrays.cut_kernel(probs['t'].values,xpos, ypos,dist)
    t = t - np.nanmean(t)

    u = u_arrays.cut_kernel(probs['u950'].values,xpos, ypos,dist)
    v = u_arrays.cut_kernel(probs['v950'].values, xpos, ypos, dist)
    shear = u_arrays.cut_kernel(probs['shear'].values, xpos, ypos, dist)

    cnt2 = np.zeros_like(u)
    cnt2[np.isfinite(u)] = 1


    if (np.sum(np.isfinite(kernel)) < 2):
        return

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1



    if kernel.shape != (dist*2+1, dist*2+1):
        pdb.set_trace()


    return kernel,  cnt, cnt2, t, u, v, shear


def get_previous_hours(date):

    ehour = 12

    if (date.hour) <= 16:
        print('Nighttime')
        edate = date - pd.Timedelta('1 days')

        edate = edate.replace(hour=ehour)

    else:
        edate = date
        edate = edate.replace(hour=ehour)


    t1 = edate
    t2 = edate + pd.Timedelta('3 hours')

    file = constants.ERA5

    try:
        cmm = xr.open_dataset(file + 'ERA5_'+str(date.year)+'_pl.nc')
    except:
        return None
    cmm = cmm.sel(time=t1)

    cm = cmm['t'].sel(level=950).squeeze() -273.15 #* 1000

    cm = cm.to_dataset()

    shear =  (cmm['u'].sel(level=600).squeeze() - cmm['u'].sel(level=925).squeeze() ) #

    vwind_srfc = cmm['v'].sel(level=950).squeeze()
    uwind_srfc = cmm['u'].sel(level=950).squeeze()

    cm['shear'] = shear
    cm['u950'] = uwind_srfc
    cm['v950'] = vwind_srfc

    pdb.set_trace()
    return cm



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
    if (np.sum(np.isfinite(lsta_da)) / lsta_da.size) < 0.05:
        print('Not enough valid')
        return None

    lsta_da.values[ttopo.values>=450] = np.nan
    lsta_da.values[gradsum>30] = np.nan
    pos = np.where((fi.values == 2) ) #(fi.values >= 5) & (fi.values < 65)

    if (np.sum(pos) == 0):
        print('No blobs found')
        return None


    dist = 200

    kernel2_list = []
    cnt_list = []
    cntp_list = []
    shear_list = []
    u_list = []
    v_list = []
    t_list = []

    probs = get_previous_hours(date)

    try:
        probs_on_lsta = lsta.salem.transform(probs)
    except RuntimeError:
        return None

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
            kernel2, cnt, cntp, t, u ,v, shear = cut_kernel(xpos, ypos, lsta_da, daybefore.month, plon, plat, -40, dist, probs=probs_on_lsta)
        except TypeError:
            continue

        kernel2_list.append(kernel2)
        cnt_list.append(cnt)
        cntp_list.append(cntp)
        t_list.append(t)
        u_list.append(u)
        v_list.append(v)
        shear_list.append(shear)

        if kernel2_list == []:
            return None

        if len(kernel2_list) == 1:
            return None
        else:

            lsta_sum = np.nansum(np.stack(kernel2_list, axis=0), axis=0)
            cnt_sum = np.nansum(np.stack(cnt_list, axis=0), axis=0)
            cntp_sum = np.nansum(np.stack(cntp_list, axis=0), axis=0)
            t_sum = np.nansum(np.stack(t_list, axis=0), axis=0)
            u_sum = np.nansum(np.stack(u_list, axis=0), axis=0)
            v_sum = np.nansum(np.stack(v_list, axis=0), axis=0)
            shear_sum = np.nansum(np.stack(shear_list, axis=0), axis=0)

        print('Returning')

        return (lsta_sum, cnt_sum, cntp_sum, t_sum, u_sum, v_sum, shear_sum)


def plot_gewex(h):
    hour=h
    dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/ERA/composite_backtrack_ERA"+str(hour).zfill(2)+".p", "rb"))
    dic2 = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/ERA/composite_backtrack_CMORPH_"+str(hour).zfill(2)+".p", "rb"))
    extent = (dic['lsta'].shape[1]-1)/2


    f = plt.figure(figsize=(12, 5))
    ax = f.add_subplot(121)

    plt.contourf((dic['lsta'] / dic['cnt']), cmap='RdBu_r',  levels=[ -0.7,-0.6,-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6, 0.7], extend='both') #-(rkernel2_sum / rcnt_sum)
    #plt.plot(extent, extent, 'bo')
    plt.colorbar(label='K')
    #pdb.set_trace()

    contours = plt.contour((dic2['prob']/ dic2['cntp'])*100, extend='both', levels=np.arange(10,70,10), cmap='viridis') # #, levels=np.arange(1,5, 0.5)
    plt.clabel(contours, inline=True, fontsize=11, fmt='%1.0f')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('23-01UTC | '+str(np.max(dic['cnt']))+' cores, LSTA & 06-06UTC antecedent rain', fontsize=9)


    ax = f.add_subplot(122)
    plt.contourf((dic['t']/ dic['cntp']), extend='both',  cmap='RdBu') # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='s-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-15,-10,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')

    ax.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')

    plt.title('ERA5 950hPa divergence (shading) & 700-950hPa wind shear', fontsize=9)


    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/ERA/'+str(hour).zfill(2)+'_single.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_all():

    hours = [16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]

    for h in hours:
        plot_gewex(h)

