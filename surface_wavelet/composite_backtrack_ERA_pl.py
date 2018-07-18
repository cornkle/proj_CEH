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

    dic = u_parallelise.era_run_arrays(5,file_loop,msg) #'rano', 'rregional', 'rcnt',

    # for k in dic.keys():
    #    dic[k] = np.nansum(dic[k], axis=0)


    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/ERA/composite_backtrack_ERA"+str(hour).zfill(2)+".p", "wb"))


def cut_kernel(xpos, ypos, arr, dist, probs=False):


    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    vdic = {}

    for d in probs.data_vars:

        var = u_arrays.cut_kernel(probs[d].values,xpos, ypos,dist)
        vdic[d] = var

    cnt2 = np.zeros_like(kernel)
    cnt2[np.isfinite(vdic[list(vdic.keys())[0]])] = 1


    if (np.sum(np.isfinite(kernel)) < 2):
        return

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1


    if kernel.shape != (dist*2+1, dist*2+1):
        pdb.set_trace()


    return kernel,  cnt, cnt2, vdic


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

    try:
        css = xr.open_dataset(file + 'ERA5_'+str(date.year)+'_srfc.nc')
    except:
        return None


    pl_clim = xr.open_dataset(file + 'CLIM/ERA5_2008-2010_CLIM_'+str(edate.month)+'-'+str(ehour)+'_pl.nc')
    srfc_clim = xr.open_dataset(file + 'CLIM/ERA5_2008-2010_CLIM_'+str(edate.month)+'-'+str(ehour)+'_srfc.nc')

    cmm = cmm.sel(time=t1)
    css = css.sel(time=t1)
    # pl_clim = pl_clim.sel(month=cmm['time.month'].values)
    # srfc_clim = srfc_clim.sel(month=cmm['time.month'].values)

    cm = cmm['t'].sel(level=950).squeeze() - pl_clim['t'].sel(level=950).squeeze() #* 1000

    cm = cm.to_dataset()

    shear =  (cmm['u'].sel(level=600).squeeze() - cmm['u'].sel(level=925).squeeze() ) #


    vwind_srfc = cmm['v'].sel(level=950).squeeze() - pl_clim['v'].sel(level=950).squeeze()
    uwind_srfc = cmm['u'].sel(level=950).squeeze() - pl_clim['u'].sel(level=950).squeeze()
    div = cmm['d'].sel(level=950).squeeze()

    cape = css['cape'].squeeze() - srfc_clim['cape'].squeeze()
    surface_pressure = css['sp'].squeeze() / 100 - srfc_clim['sp'].squeeze() / 100
    sl_pressure = css['msl'].squeeze() / 100 - srfc_clim['msl'].squeeze()/ 100
    sh = css['ishf'].squeeze() - srfc_clim['ishf'].squeeze()
    t2 = css['t2m'].squeeze() - srfc_clim['t2m'].squeeze()
    q = cmm['q'].sel(level=950).squeeze() - pl_clim['q'].sel(level=950).squeeze()

    cm['shear'] = shear
    cm['u950'] = uwind_srfc
    cm['v950'] = vwind_srfc
    cm['cape'] = cape
    cm['sf'] = surface_pressure
    cm['slp'] = sl_pressure
    cm['sh'] = sh
    cm['t2'] = t2
    cm['div'] = div *1000
    cm['q'] = q
    srfc_clim.close()
    pl_clim.close()
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

    kernel2_sum = np.zeros((dist*2+1, dist*2+1))
    cnt_sum = np.zeros((dist*2+1, dist*2+1))
    cntp_sum = np.zeros((dist*2+1, dist*2+1))
    edic = {}

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
            kernel2, cnt, cntp, vdic = cut_kernel(xpos, ypos, lsta_da, dist, probs=probs_on_lsta)
        except TypeError:
            continue


        kernel2_sum = np.nansum(np.stack([kernel2_sum, kernel2]), axis=0)
        cntp_sum = np.nansum(np.stack([cntp_sum, cntp]), axis=0)
        cnt_sum = np.nansum(np.stack([cnt_sum, cnt]), axis=0)

        for ks in vdic.keys():
            if ks in edic:
                edic[ks] = np.nansum(np.stack([edic[ks], vdic[ks]]), axis=0)
            else:
                edic[ks] = vdic[ks]

    outlist = [kernel2_sum, cnt_sum, cntp_sum]
    outnames = ['lsta',  'cnt', 'cntp']
    for ek in edic.keys():
        outnames.append(ek)
        outlist.append(edic[ek])

    print('Returning')

    return outlist, outnames


def plot_gewex(h):
    hour=h
    dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/ERA/composite_backtrack_ERA"+str(hour).zfill(2)+".p", "rb"))
    dic2 = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/ERA/composite_backtrack_CMORPH_"+str(hour).zfill(2)+".p", "rb"))
    extent = (dic['lsta'].shape[1]-1)/2
    xlen = dic['lsta'].shape[1]
    ylen = dic['lsta'].shape[0]

    xv, yv = np.meshgrid(np.arange(ylen), np.arange(xlen))
    st=30
    xquiv = xv[4::st, 4::st]
    yquiv = yv[4::st, 4::st]

    u = (dic['u950']/ dic['cntp'])[4::st, 4::st]
    v = (dic['v950']/ dic['cntp'])[4::st, 4::st]


    f = plt.figure(figsize=(15, 8))
    ax = f.add_subplot(231)

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


    ax1 = f.add_subplot(232)
    plt.contourf(((dic['lsta'])/ dic['cnt']), extend='both',  cmap='RdBu_r',levels=[ -0.7,-0.6,-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6, 0.7]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='s-1')
    contours = plt.contour((dic['t2'] / dic['cntp']), extend='both',levels=[ -0.7,-0.6,-0.5,-0.4,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6, 0.7], cmap='RdBu') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.1f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')

    ax1 = f.add_subplot(233)
    plt.contourf(((dic['cape'])/ dic['cntp']), extend='both',  cmap='RdBu') # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='J kg-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    #qu = ax1.quiver(xquiv, yquiv, u, v, scale=50)

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')

    ax1 = f.add_subplot(234)
    plt.contourf(((dic['t'])/ dic['cntp']), extend='both',  cmap='RdBu') # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='s-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=30)

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')

    ax1 = f.add_subplot(235)
    plt.contourf(((dic['t2'])/ dic['cntp']), extend='both',  cmap='RdBu', levels=np.arange(-1,1.1,0.2)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='m s-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=30)

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')

    ax1 = f.add_subplot(236)
    plt.contourf(((dic['q'])/ dic['cntp']), extend='both',  cmap='RdBu') # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='m s-1')
    contours = plt.contour((dic['shear'] / dic['cntp']), extend='both',levels=np.arange(-17,-12,0.5), cmap='viridis') #np.arange(-15,-10,0.5)
    plt.clabel(contours, inline=True, fontsize=9, fmt='%1.2f')
    qu = ax1.quiver(xquiv, yquiv, u, v, scale=30)

    ax1.set_xticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_yticklabels(np.array((np.linspace(0, extent*2, 9) - extent) * 3, dtype=int))
    ax1.set_xlabel('km')
    ax1.set_ylabel('km')

    plt.title('ERA5 950hPa divergence (shading) & 700-950hPa wind shear', fontsize=9)


    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/ERA/'+str(hour).zfill(2)+'_single.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_all():

    hours = [16,17,18,19,20,21,22,23,0,1,2,3,4,5,6,7]

    for h in hours:
        plot_gewex(h)

