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


    file = constants.MCS_CENTRE70

    hour = h

    msg = xr.open_dataarray(file)
    msg = msg[((msg['time.hour'] >= 23 ) | (msg['time.hour'] <= 1 )) & ((msg['time.minute'] == 0) & (
        msg['time.year'] >= 2008) & (msg['time.year'] <= 2010) & (msg['time.month'] >=6)) ]

    msg = msg.sel(lat=slice(10.9,19), lon=slice(-9.8,9.8))

    dic = u_parallelise.era_run_arrays(5,file_loop,msg) #'rano', 'rregional', 'rcnt',

    # for k in dic.keys():
    #    dic[k] = np.nansum(dic[k], axis=0)


    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/composite_backtrack_ERA_pl_"+str(hour).zfill(2)+".p", "wb"))


def cut_kernel(xpos, ypos, expos, eypos, arr, dist, probs=False):


    kernel = u_arrays.cut_kernel(arr,xpos, ypos,dist)

    edist = int(dist/10)

    vdic = {}
    cnt2 = 0

    for d in probs.data_vars:
        #print(d)
        var = u_arrays.cut_kernel_3d(probs[d].values,expos, eypos,edist)
        var = np.mean(var[:,:, edist-1:edist+2], axis=2)
        #var = np.mean(var[:, edist - 1:edist + 2, :], axis=1)
        vdic[d] = var


    cnt2 = np.zeros_like(vdic[list(vdic.keys())[0]])
    cnt2[np.isfinite(vdic[list(vdic.keys())[0]])] = 1

    if (np.sum(np.isfinite(kernel)) < 2):
        return

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1


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
        cmm = xr.open_dataset(file + 'ERA5_'+str(date.year)+'_pls.nc')
    except:
        return None

    pl_clim = xr.open_dataset(file + 'CLIM/ERA5_2008-2010_CLIM_'+str(edate.month)+'-'+str(ehour)+'_pls.nc')


    cmm = cmm.sel(time=t1)

    cm = cmm['t'] - pl_clim['t']

    cm = cm.to_dataset()

    vwind_srfc = cmm['v'] - pl_clim['v']
    uwind_srfc = cmm['u'] - pl_clim['u']
    div = cmm['d'] - pl_clim['d']

    q = cmm['q'] - pl_clim['q']

    cm['u950'] = uwind_srfc
    cm['v950'] = vwind_srfc

    cm['div'] = div *1000
    cm['q'] = q

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



    probs = get_previous_hours(date)


    kernel2_sum = np.zeros((dist*2+1, dist*2+1))
    cnt_sum = np.zeros((dist*2+1, dist*2+1))
    cntp_sum = np.zeros((len(probs.level), int(dist/10)*2+1))
    edic = {}

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


        epoint = probs.sel(latitude=lat, longitude=lon, method='nearest')
        elat = epoint['latitude'].values
        elon = epoint['longitude'].values

        expos = np.where(probs['longitude'].values == elon)
        expos = int(expos[0])
        eypos = np.where(probs['latitude'].values == elat)
        eypos = int(eypos[0])
        try:
            kernel2, cnt, cntp, vdic = cut_kernel(xpos, ypos, expos, eypos, lsta_da, dist, probs=probs)
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

    outlist = [cntp_sum]
    outnames = ['cntp']
    for ek in edic.keys():
        outnames.append(ek)
        outlist.append(edic[ek])

    print('Returning')

    return outlist, outnames



def plot_doug(h):
    hour=h
    chour=17
    dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/composite_backtrack_ERA_pl_"+str(hour).zfill(2)+"NS.p", "rb"))

    levels= [400,450,500,550,600,650,700,750,825,850,875,900,925,950,975]


    f = plt.figure(figsize=(10,8))
    ax = f.add_subplot(221)
    plt.contourf(np.arange(-20,21)*30, levels, (dic['u950'] / dic['cntp']), cmap='RdBu_r', levels=[-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.title('u wind', fontsize=9)
    plt.xlabel('South-North extent')
    plt.ylabel('Pressure level (hPa)')


    ax1 = f.add_subplot(222)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['v950'])/ dic['cntp']), extend='both',  cmap='RdBu_r',levels=[-1.5,-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1.5]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.title('v wind', fontsize=9)
    plt.xlabel('South-North extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(223)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.title('specific humidity', fontsize=9)
    plt.xlabel('South-North extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(224)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['t'])/ dic['cntp']), extend='both',  cmap='RdBu_r', vmin=-0.5, vmax=0.5) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.title('Divergence', fontsize=9)
    plt.xlabel('South-North extent')
    plt.ylabel('Pressure level (hPa)')


    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/'+str(hour).zfill(2)+'_NS_pl.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


def plot_doug_big(h):
    hour=h
    chour=17
    dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/composite_backtrack_ERA_pl_"+str(hour).zfill(2)+"NS.p", "rb"))

    levels= [400,450,500,550,600,650,700,750,825,850,875,900,925,950,975]


    f = plt.figure(figsize=(15,8))
    ax = f.add_subplot(231)
    plt.contourf(np.arange(-20,21)*30, levels, (dic['u950'] / dic['cntp']), cmap='RdBu_r', levels=[-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.title('u wind', fontsize=9)
    plt.xlabel('South-North extent')
    plt.ylabel('Pressure level (hPa)')


    ax1 = f.add_subplot(232)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['v950'])/ dic['cntp']), extend='both',  cmap='RdBu_r',levels=[-1.5,-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1.5]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.title('v wind', fontsize=9)
    plt.xlabel('South-North extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(233)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.title('specific humidity', fontsize=9)
    plt.xlabel('South-North extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(234)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['div'])/ dic['cntp']), extend='both',  cmap='RdBu_r', levels=np.arange(-0.005, 0.0051, 0.001)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='s-1')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.title('Divergence', fontsize=9)
    plt.xlabel('South-North extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(235)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['t'])/ dic['cntp']), extend='both',  cmap='RdBu_r', levels=np.arange(-0.5,0.6,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.title('temperature', fontsize=9)
    plt.xlabel('South-North extent')
    plt.ylabel('Pressure level (hPa)')


    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/'+str(hour).zfill(2)+'_NS_pl_big.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()



#np.arange(-20,21)*30, levels,

def plot_doug_WE(h):
    hour=h
    chour=17
    dic = pkl.load(open("/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/composite_backtrack_ERA_pl_"+str(hour).zfill(2)+"WE.p", "rb"))

    levels= [400,450,500,550,600,650,700,750,825,850,875,900,925,950,975]


    f = plt.figure(figsize=(10,8))
    ax = f.add_subplot(221)
    plt.contourf(np.arange(-20,21)*30, levels, (dic['u950'] / dic['cntp']), cmap='RdBu_r', levels=[-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1], extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    #plt.gca().invert_xaxis()
    plt.title('u wind', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')


    ax1 = f.add_subplot(222)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['v950'])/ dic['cntp']), extend='both',  cmap='RdBu_r',levels=[-1.5,-1, -0.8, -0.6,-0.4,-0.2,0.2,0.4,0.6,  0.8, 1.5]) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='m s-1')
    plt.gca().invert_yaxis()
    #plt.gca().invert_xaxis()
    plt.title('v wind', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(223)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['q'])*1000/ dic['cntp']), extend='both',  cmap='RdBu',levels=np.arange(-1,1.1,0.1)) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='g kg-1')
    plt.gca().invert_yaxis()
    #plt.gca().invert_xaxis()
    plt.title('specific humidity', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')

    ax1 = f.add_subplot(224)
    plt.contourf(np.arange(-20,21)*30, levels,((dic['div'])/ dic['cntp']), extend='both',  cmap='RdBu_r', vmin=-0.004, vmax=0.004) # #, levels=np.arange(1,5, 0.5), levels=np.arange(10,70,5)
    plt.colorbar(label='K')
    plt.gca().invert_yaxis()
    #plt.gca().invert_xaxis()
    plt.title('Divergence', fontsize=9)
    plt.xlabel('West-East extent')
    plt.ylabel('Pressure level (hPa)')


    plt.tight_layout()
    plt.savefig('/users/global/cornkle/figs/LSTA-bullshit/corrected_LSTA/system_scale/doug/'+str(hour).zfill(2)+'_WE_pl.png')#str(hour).zfill(2)+'00UTC_lsta_fulldomain_dominant<60.png')
    plt.close()


#np.arange(-20,21)*30, levels,