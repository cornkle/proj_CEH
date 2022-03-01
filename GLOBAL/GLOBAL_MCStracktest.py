# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import ipdb
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
import matplotlib.pylab as pylab
import glob
import os
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_gis, u_arrays as ua, constants as cnst, u_grid, u_darrays
from scipy.interpolate import griddata
import multiprocessing
import metpy
from metpy import calc
from metpy.units import units

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)



MREGIONS = {'WAf' : [[-15,15,10,25], 'spac', 0], # last is hourly offset to UCT # 12    # [-18,25,4,25]
 'SAf' : [[20,35, -35,-15], 'spac', 2], # 10
 'india' : [[70,90, 5,30], 'asia', 5], # 7
 'china' : [[105,115,25,40], 'asia', 8 ], # 4
 'australia' : [[120,140,-23, -11], 'asia', 9], # 3
 'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4] , # 16
 'trop_SA' : [[-75, -50, -20, -5], 'spac', -5], # 17
 'GPlains' : [[-100,-90,32,47], 'nam', -6] # # 18

}

REGION = 'WAf'

def composite(h):
    #pool = multiprocessing.Pool(processes=8)

    h = h - (MREGIONS[REGION])[2]

    # if h >= 24:
    #     shift = (12 + h) * (-1)
    # else:
    #     shift = 12 - h

    print('Hour: ', h)

    print(REGION)
    path = cnst.network_data + 'data/GLOBAL_MCS/save_composites/'

    for y in np.arange(2012, 2020):
        msg = pkl.load(
            open('/home/ck/DIR/cornkle/data/GLOBAL_MCS/save_files/'+REGION+'_noinit__mcs_tracks_extc_'+str(y)+'0101_'+str(y)+'1231.p', "rb"))
        hour = h

        # file = cnst.MCS_POINTS_DOM
        # msg = xr.open_dataarray(file)
        # msg = msg[(msg['time.hour'] == h) & (msg['time.minute'] == 0) & (
        #         msg['time.year'] == y) & (msg['time.month'] >= 6)]
        domain = (MREGIONS[REGION])[0]
        m1 = 6
        m2 = 9
        lontag = 'meanlon'
        lattag = 'meanlat'

        for k in msg.keys():

            msg[k] = np.array(msg[k])

        # inmask = (msg[lontag]>=domain[0])&(msg[lontag]<=domain[1])&(msg[lattag]>=domain[2])\
        #        &(msg[lattag]<=domain[3]) & (msg['hour'] == h) & (msg['mcs_status'] >= 0) & (msg['pf_landfrac'] > 0.8)  & (msg['month']>=m1) & (msg['month']<=m2)
        # if (np.sum(inmask) == 0) & (REGION in ['GPlains']):
        #     inmask = (msg[lontag]>=domain[0])&(msg[lontag]<=domain[1])&(msg[lattag]>=domain[2])\
        #        &(msg[lattag]<=domain[3]) & (msg['hour'] == h) & (msg['mcs_status'] >= 0)


        inmask = (msg[lontag] >= domain[0]) & (msg[lontag] <= domain[1]) & (msg[lattag] >= domain[2]) \
                 & (msg[lattag] <= domain[3]) & (msg['hour'] >= h-1) & (msg['hour'] <= h+1) & (msg['pf_mcsstatus'] > 0) & (
                             msg['pf_landfrac'] > 0.9)  & (msg['month']>=m1) & (msg['month']<=m2)

        if (np.sum(inmask) == 0) & (REGION in ['GPlains']):
            inmask = (msg[lontag] >= domain[0]) & (msg[lontag] <= domain[1]) & (msg[lattag] >= domain[2]) \
                     & (msg[lattag] <= domain[3]) & (msg['hour'] >= h - 1) & (msg['hour'] <= h + 1) & (
                                 msg['pf_mcsstatus'] > 0) & (msg['month']>=m1) & (msg['month']<=m2)

        #msc_status > 0 : MCS  = (-32C over 40000km2)
        # pf_landfrac > 0.8: 80% of rain fields over land

        mask = np.where(inmask)
        #ipdb.set_trace()
        for k in msg.keys():
            msg[k] = (msg[k])[mask]

        msg['date'] = []
        msg['day'] = []
        for yi, k in enumerate(msg['base_time']):
            bbt = pd.to_datetime(k)

            msg['date'].append(bbt.replace(hour=0, minute=0))
            msg['day'].append(bbt.day)
        msg['date'] = np.array(msg['date'])
        msg['day'] = np.array(msg['day'])
        #ipdb.set_trace()
        # msg.attrs['refhour'] = h
        # msg.attrs['eh'] = eh
        #
        ubt = np.unique(msg['date'])

        chunks = []
        for ut in ubt:
            daydir = {}

            pos = np.where(msg['date'] == ut)

            for k in msg.keys():
                 daydir[k] = (msg[k])[pos]
            chunks.append(daydir)

        dic = u_parallelise.run_arrays(4, file_loop, chunks, ['ano', 'regional', 'cnt', 'allcnt'])

        # res = []
        # #ipdb.set_trace()
        # for m in chunks:
        #     out = file_loop(m)
        #     res.append(out)
        # #ipdb.set_trace()
        # return

        for k in dic.keys():
           dic[k] = np.nansum(dic[k], axis=0)
        print('File written')
        #pkl.dump(dic, open(path + "/"+REGION+"_AMSR2-global_2012-2019_SM_MCSTRACK_BOX-ANNUAL_PFbased_"+str(y)+'_h'+str(hour).zfill(2)+".p", "wb"))
        pkl.dump(dic, open(path + "/"+REGION+"_AMSR2-global_2012-2019_SM_MCSTRACK_BOX-ANNUAL_PF_"+str(y)+'_h'+str(hour).zfill(2)+".p", "wb"))




def cut_kernel(xpos, ypos, arr, wd = None):

    dist = 11
    #dist=100


    kernel = ua.cut_kernel(arr,xpos, ypos,dist)

    if wd != None:
        kernel = ua.rotate(kernel, wd, ref_angle=90)
    #ipdb.set_trace()

    kernel3 = kernel - np.nanmean(kernel)

    cnt = np.zeros_like(kernel)
    cnt[np.isfinite(kernel)] = 1

    if kernel.shape != (dist*2+1, dist*2+1):
        return None

    return kernel, kernel3, cnt


def file_loop(fi):


    date = fi['date'][0]

    hour = fi['hour']
    print('Doing day: ', date)


    box = (MREGIONS[REGION])[0]

    dayd = pd.Timedelta('1 days')

    if (date.hour) <= 16:
        print('Nighttime')
        daybefore = date #- dayd
    else:
        print('Daytime')
        daybefore = date #- dayd

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)
    # alt_path = cnst.AMSRE_ANO_DAY
    # alt_path = cnst.AMSRE_ANO_DAY_CORR
    # try:
    #     lsta = xr.open_dataset(alt_path + 'sma_' + fdate + '.nc')
    # except:
    #     print('Surface file not found, return')
    #     return None

    alt_path = '/media/ck/Elements/global/AMSR2/daily/25km/day_anom/'
    try:
        lsta = xr.open_dataset(alt_path + 'amsr2_25km_anom_' + fdate + '.nc')
    except:
        print('Surface file not found, return')
        return None

    lsta = lsta.sel(time=str(daybefore.year) + '-' + str(daybefore.month) + '-' + str(daybefore.day))
    lsta = lsta.sel(lon=slice(box[0]-3, box[1]+3), lat=slice(box[2]-3, box[3]+3)) # define SM region box to consider

    print('Doing ' + 'AMSR2_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')


    lsta_da = lsta['soil_moisture_c1'].squeeze()  # soil_moisture_c1

    kernel2_list = []
    kernel3_list = []
    cnt_list = []
    all_cnt_list = []
    #ipdb.set_trace()
    cid = -1

############################
    for id in range(len(fi['direction'])):

        direction = fi['direction'][id]
        if np.isnan(direction):
            continue


        ##################
        # permcs = 0
        # for lat, lon in zip(fi['pf_lat'][id].squeeze().flatten(), fi['pf_lon'][id].squeeze().flatten()):  #zip(fi['meanlat'], fi['meanlon'])
        #     cid += 1
        #     if ((fi['pf_maxrainrate'][id]).squeeze().flatten()[permcs] < 1) | np.isnan((fi['pf_maxrainrate'][id]).squeeze().flatten()[permcs]):
        #         print('Pf rainrate below 5mm, continue')
        #         permcs +=1
        #         continue
        #     permcs += 1
        #####################
        maxrrpos = np.argmax(fi['pf_maxrainrate'][id].squeeze().flatten())
        maxrr = np.max(fi['pf_maxrainrate'][id].squeeze().flatten())
        cid += 1
        if (maxrr < 1) | np.isnan(maxrr):
            continue
        lat = fi['pf_lat'][id].squeeze().flatten()[maxrrpos]
        lon = fi['pf_lon'][id].squeeze().flatten()[maxrrpos]
       # ipdb.set_trace()
    ############################
    # for lat, lon in zip(fi['meanlat'], fi['meanlon']):
    #     cid += 1
    #     #direction = (fi[direction])[cid-1]
    ####################
        try:
            point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.2)
        except KeyError:
            print('point index error')
            #ipdb.set_trace()
            continue

        plat = point['lat'].values
        plon = point['lon'].values
        #ipdb.set_trace()
        xpos = np.where(lsta_da['lon'].values == plon)
        try:
            xpos = int(xpos[0])
        except:
            continue
        ypos = np.where(lsta_da['lat'].values == plat)
        ypos = int(ypos[0])

        try:
            kernel2, kernel3, cnt = cut_kernel(xpos, ypos, lsta_da, wd=direction)
        except TypeError:
            continue

        # f = plt.figure()
        # ax = f.add_subplot(121)
        # cb = ax.contourf(kernel2, cmap='RdBu')
        # plt.colorbar(cb)
        # ax = f.add_subplot(122)
        # cb = ax.contourf(cnt, cmap='viridis')
        # plt.colorbar(cb)
        # f.show()

        cnt_all = np.zeros_like(cnt)+1
        kernel2_list.append(kernel2)
        kernel3_list.append(kernel3)
        cnt_list.append(cnt)
        all_cnt_list.append(cnt_all)

    if kernel2_list == []:
        return None

    if len(kernel2_list) == 1:
        return None
    else:

        kernel2_sum = np.nansum(np.stack(kernel2_list, axis=0), axis=0)
        kernel3_sum = np.nansum(np.stack(kernel3_list, axis=0), axis=0)
        cnt_sum = np.nansum(np.stack(cnt_list, axis=0), axis=0)
        allcnt_sum = np.nansum(np.stack(all_cnt_list, axis=0), axis=0)

    print('Returning')
    del lsta
    return (kernel2_sum, kernel3_sum, cnt_sum, allcnt_sum)


def plot(h):
    h = h - (MREGIONS[REGION])[2]
    print('Hour: ', h)
    pin = cnst.network_data + 'data/GLOBAL_MCS/save_composites/' + "/" + REGION + "_AMSR2-global_2012-2019_SM_MCSTRACK_BOX-ANNUAL_PF_"
    #pin = cnst.network_data + 'data/GLOBAL_MCS/save_composites/' + "/"+REGION+"_AMSR2-global_2012-2019_SM_MCSTRACK_BOX-ANNUAL_PFbased_"
    y1 = 2012
    y2 = 2020

    dic = pkl.load(open(
            pin+str(y1)+'_h'+str(h).zfill(2)+".p", "rb"))

    def coll(dic, h, year):
        print(h)
        core = pkl.load(open(
            pin+str(year)+'_h'+str(h).zfill(2)+".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]

    for y in range(y1+1, y2):
        print('Coll', y)
        coll(dic, h, y)

    extent = 11

    f = plt.figure(figsize=(15, 4))
    ax = f.add_subplot(131)

    thresh = 2.5

    plt.contourf(dic['regional'] / dic['cnt'], cmap='RdBu',  levels=np.linspace(thresh*-1,thresh,10), extend='both')
    plt.plot(extent, extent, 'bo')
    ax.axvline(extent, linestyle='dashed', color='k')
    ax.axhline(extent, linestyle='dashed', color='k')
    print('Extent', extent)
    ax.set_xticks((np.linspace(0, 2 * extent, 9)))
    ax.set_xticklabels(((np.linspace(0, (2*extent), 9)-extent)*30).astype(int))
    #ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
    ax.set_yticks((np.linspace(0, 2 * extent, 9)))
    ax.set_yticklabels(((np.linspace(0, (2*extent), 9)-extent)*30).astype(int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    plt.colorbar(label='%')
    plt.title(REGION+' '+str(y1)+'-'+str(y2)+' box anom., Nb: ' + str(np.max(dic['allcnt'])) + '| ' + str(h).zfill(2) + '00UTC',
              fontsize=10)


    ax = f.add_subplot(132)

    plt.contourf((dic['ano'] / dic['cnt']), cmap='RdBu',  levels=np.linspace(thresh*-1,thresh,10), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticks((np.linspace(0, 2 * extent, 9)))
    ax.set_xticklabels(((np.linspace(0, (2*extent), 9)-extent)*30).astype(int))
    ax.set_yticks((np.linspace(0, 2 * extent, 9)))
    ax.set_yticklabels(((np.linspace(0, (2*extent), 9)-extent)*30).astype(int))
    ax.set_xlabel('km')
    ax.set_ylabel('km')
    ax.axvline(extent, linestyle='dashed', color='k')
    ax.axhline(extent, linestyle='dashed', color='k')
    plt.colorbar(label='%')
    plt.title('Seasonal anomaly',
              fontsize=10)


    ax = f.add_subplot(133)

    plt.contourf(dic['cnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticks((np.linspace(0, 2 * extent, 9)))
    ax.set_xticklabels(((np.linspace(0, (2*extent), 9)-extent)*30).astype(int))
    ax.set_yticks((np.linspace(0, 2 * extent, 9)))
    ax.set_yticklabels(((np.linspace(0, (2*extent), 9)-extent)*30).astype(int))
    ax.set_xlabel('Degrees')
    ax.set_ylabel('Degrees')
    ax.axvline(extent, linestyle='dashed', color='k')
    ax.axhline(extent, linestyle='dashed', color='k')
    plt.colorbar(label='n')
    plt.title('Valid count',
              fontsize=10)


    plt.tight_layout()
    f.savefig('/home/ck/DIR/cornkle/figs/GLOBAL_MCS/'+REGION+'_'+str(y1)+'-'+str(y2-1)+'_SM_composite_AMSR2-global_BOX-ANNUAL_PF.jpg')