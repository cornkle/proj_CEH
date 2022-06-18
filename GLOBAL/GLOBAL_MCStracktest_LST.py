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
from utils import u_met, u_parallelise, u_gis, u_arrays as ua, constants as cnst, u_grid, u_darrays as uda
from scipy.interpolate import griddata
import multiprocessing
import metpy
from metpy import calc
from metpy.units import units

import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)



MREGIONS = {'WAf' : [[-18,25,4,25], 'spac', 0, (1,7), (8,12), (1,12)], # last is hourly offset to UCT # 12    # [-18,25,4,25]
 'SAf' : [[20,35, -35,-15], 'spac', 2, (9,12), (1,5), (1,12)], # 10
 'india' : [[70,90, 5,30], 'asia', 5, (1,7), (8,12), (1,12)], # 7
 'china' : [[105,115,25,40], 'asia', 8 , (1,7), (8,12), (1,12)], # 4
 'australia' : [[120,140,-23, -11], 'asia', 9, (10,12), (1,5), (1,12)], # 3
 'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4, (9,12), (1,5), (1,12)] , # 16
 'trop_SA' : [[-75, -50, -20, -5], 'spac', -5, (1,12), (1,12), (1,12)], # 17
 'GPlains' : [[-100,-90,32,47], 'nam', -6, (1,7), (8,12), (1,12)] # # 18

}

REGIONS = ['GPlains', 'sub_SA', 'WAf', 'china', 'india', 'australia']
SENSOR = 'terra'
SENSOP = 10 # (LT overpass)
extag = '_day0init' #
INIT_DISTANCE = 0
AREA = 0 #1000
TRACKTIME = 1

def composite(rawhour):

    h = rawhour - (MREGIONS[REGION])[2]
    def h_checker(h):
        if h >= 24:
            h = h-24
        if h == 24:
            h = 0
        if h < 0:
            h = h+24
        return h
    h = h_checker(h)
    h2 = h_checker(h+1)
    h1 = h_checker(h-1)

    print('Hour: ', h, h1, h2)

    print(REGION)
    path = cnst.network_data + 'data/GLOBAL_MCS/save_composites/'

    for y in np.arange(2003, 2019):

        msg = pd.read_csv(cnst.lmcs_drive +'/save_files/'+REGION+'_initTime__mcs_tracks_extc_'+str(y)+'0101_'+str(y)+'1231.csv')
        msg = msg.to_dict(orient='list')

        domain = (MREGIONS[REGION])[0]

        m1 = MONTHS[0]
        m2 = MONTHS[1]

        msg['date'] = []
        msg['utc_date'] = []
        msg['day'] = []
        msg['lt_hour'] =[]
        for yi, k in enumerate(msg['base_time']):

            bbt = pd.to_datetime(k) #- pd.Timedelta('1 days')
            hourchange = (MREGIONS[REGION])[2]

            if hourchange < 0:
                bbl = bbt - pd.Timedelta(str(np.abs(hourchange))+' hours')
            elif hourchange > 0:
                    bbl = bbt + pd.Timedelta(str(np.abs(hourchange)) + ' hours')
            else:
                bbl = bbt

            msg['utc_date'].append(bbt.replace(hour=0, minute=0))
            msg['lt_hour'].append(bbl.hour)
            msg['date'].append(bbl.replace(hour=0, minute=0))
            msg['day'].append(bbt.day)

        for k in msg.keys():

            msg[k] = np.array(msg[k])

        inmask =   ((msg['lt_hour'] >= SENSOP+2) & (msg['tracktime'] >= TRACKTIME) & (msg['ccs_area'] >= AREA) & \
                   ((msg['hour']==h)  | (msg['hour']==h1)  | (msg['hour']==h2)) & \
                    (msg['pf_landfrac'] > 0.99)  & (msg['month']>=m1) & (msg['month']<=m2) & ((np.abs(msg['londiff_loc-init'])>=INIT_DISTANCE)|(np.abs(msg['latdiff_loc-init'])>=INIT_DISTANCE)))


        if (np.sum(inmask) == 0) & (REGION in ['GPlains']):
            inmask = ((msg['lt_hour'] >= SENSOP+2)  & (msg['tracktime'] >= TRACKTIME) & (msg['ccs_area'] >= AREA) & \
                      ((msg['hour'] == h) | (msg['hour'] == h1) | (msg['hour'] == h2)) & \
                       (msg['month'] >= m1) & (msg['month'] <= m2) & (
                                  (np.abs(msg['londiff_loc-init']) >= INIT_DISTANCE) | (np.abs(msg['latdiff_loc-init']) >= INIT_DISTANCE)))

        #msc_status > 0 : MCS  = (-32C over 40000km2)
        # pf_landfrac > 0.8: 80% of rain fields over land

        mask = np.where(inmask)
        for k in msg.keys():
            msg[k] = (msg[k])[mask]

        ubt = np.unique(msg['date'])
        msg['date'] = np.array(msg['date'])

        chunks = []
        for ut in ubt:
            daydir = {}

            pos = np.where(msg['date'] == ut)[0]

            for k in msg.keys():
                 daydir[k] = np.array(msg[k])[pos]
            chunks.append(daydir)

        try:
            dic = u_parallelise.run_arrays(4, file_loop, chunks, ['ano', 'regional', 'cnt', 'allcnt', 'init'])
        except:
            continue
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
        pkl.dump(dic, open(path + "/"+REGION+"_LSTA_"+SENSOR+"_SWA_2012-2019_MCSTRACK_BOX-ANNUAL_PF_DAY_"+str(m1).zfill(2)+'-'+str(m2).zfill(2)+'_'+str(y)+'_h'+str(rawhour).zfill(2)+extag+".p", "wb"))


def cut_kernel(xpos, ypos, arr, inits, wd = None):

    dist = 85
    #dist=100

    res = 0.05

    kernel = ua.cut_kernel(arr,xpos, ypos,dist)
    if np.sum(np.isfinite(kernel)) < 0.1 * kernel.size:
        print('Not enough data')
        return None

    ilon = np.round(inits[0] / res).astype(int)
    ilat = np.round(inits[1] / res).astype(int)

    if (np.abs(ilon)<29) & (np.abs(ilat)<29):
        return None

    #print('lonlat', ilon, ilat)

    cnt = np.zeros_like(kernel)
    init = np.zeros_like(kernel)

    try:
        init[dist+1-ilat, dist+1-ilon] = 1
    except:
        pass

    if kernel.shape != (dist*2+1, dist*2+1):
        return None

    if wd != None:
        kernel = ua.rotate(kernel, wd, ref_angle=90)
        init = ua.rotate(init, wd, ref_angle=90)

    cnt[np.isfinite(kernel)] = 1
    kernel3 = kernel - np.nanmean(kernel)

    cut = 25
    print('shape', kernel[cut:-cut, cut:-cut].shape)
    return kernel[cut:-cut, cut:-cut], kernel3[cut:-cut,cut:-cut], cnt[cut:-cut,cut:-cut], init[cut:-cut,cut:-cut]


def file_loop(fi):


    date = fi['date'][0]

    hour = fi['hour']
    print('Doing day: ', date)



    # ulist = []
    # vlist = []
    # for direc in fi['direction']:
    #     ulist.append(np.sin(np.deg2rad(direc)))
    #     vlist.append(np.cos(np.deg2rad(direc)))
    #
    # umean = np.sum(ulist)
    # vmean = np.sum(vlist)
    #
    # avg_wd = np.rad2deg(np.arctan2(umean,vmean))
    # if avg_wd<0:
    #     avg_wd = avg_wd+360
    # #print('mean WD', avg_wd, fi['direction'])
    # direct = avg_wd


    box = (MREGIONS[REGION])[0]

    dayd = pd.Timedelta('1 days')

    # if (date.hour) <= 16:
    #     print('Nighttime')
    #     daybefore = date #- dayd
    # else:
    #     print('Daytime')
    daybefore = date #- dayd

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)

    alt_path = cnst.lmcs_drive+'MODIS_LST/'+SENSOR+'/anom_v61/'
    try:
        lsta = xr.open_dataset(alt_path +SENSOR+'_05deg_anom_' + fdate + '.nc')

    except:
        print('Surface file not found, return')
        return None

    lsta = lsta.sel(time=str(daybefore.year) + '-' + str(daybefore.month) + '-' + str(daybefore.day))
    lsta = uda.flip_lat(lsta)
    lsta = lsta.sel(lon=slice(box[0]-3, box[1]+3), lat=slice(box[2]-3, box[3]+3)) # define SM region box to consider

    print('Doing ' + SENSOR+'_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc', REGION)


    lsta_da = lsta['LST_Day'].squeeze()
    lsta_da = lsta_da.where((np.isfinite(lsta_da)) & (lsta_da > -8000)) / 100

    #lsta_da = lsta_da - np.nanmean(lsta_da)

    kernel2_list = []
    kernel3_list = []
    cnt_list = []
    init_list = []
    all_cnt_list = []

    cid = 0

############################
    # for id in range(len(fi['direction'])):
    #
    #     direct = fi['direction'][id]
    #     if np.isnan(direct):
    #         continue
    #
    #
    #         #################
    #     permcs = 0
    #     for lat, lon in zip(fi['pf_lat'][id].squeeze().flatten(), fi['pf_lon'][id].squeeze().flatten()):  #zip(fi['meanlat'], fi['meanlon'])
    #         cid += 1
    #         # if ((fi['pf_maxrainrate'][id]).squeeze().flatten()[permcs] < 1) | np.isnan((fi['pf_maxrainrate'][id]).squeeze().flatten()[permcs]):
    #         #     print('Pf rainrate below 5mm, continue')
    #         #     permcs +=1
    #         #     continue
    #         permcs += 1
            #####################
            #takes just strongest rain feature per MCS
        # maxrrpos = np.argmax(fi['pf_maxrainrate'][id].squeeze().flatten())
        # maxrr = np.max(fi['pf_maxrainrate'][id].squeeze().flatten())
        # cid += 1
        # if (maxrr < 1) | np.isnan(maxrr):
        #     continue
        # lat = fi['pf_lat'][id].squeeze().flatten()[maxrrpos]
        # lon = fi['pf_lon'][id].squeeze().flatten()[maxrrpos]
       # ipdb.set_trace()
    ############################
    # if len(fi['pf_lat']) == 0:
    #     return None
    for id in range(3): # 3 precip entries
        idss = 0
        for lat, lon in zip(fi['pf_lat'+str(id+1)], fi['pf_lon'+str(id+1)]):

            ####################
            ulist = []
            vlist = []
            for direc in [fi['direction1'][idss], fi['direction0'][idss],fi['direction-1'][idss] , fi['direction-2'][idss]]: #
                ulist.append(np.sin(np.deg2rad(direc)))
                vlist.append(np.cos(np.deg2rad(direc)))

            umean = np.sum(ulist)
            vmean = np.sum(vlist)

            avg_wd = np.rad2deg(np.arctan2(umean, vmean))
            if avg_wd < 0:
                avg_wd = avg_wd + 360
            # print('mean WD', avg_wd, fi['direction'])
            direct = avg_wd

            init_lon = fi['londiff_loc-init'][idss]
            init_lat = fi['latdiff_loc-init'][idss]

            idss = idss+1

            try:
                point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.05)
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
                kernel2, kernel3, cnt, init = cut_kernel(xpos, ypos, lsta_da, [init_lon,init_lat], wd=direct)
            except TypeError:
                continue

            cnt_all = np.zeros_like(cnt)+1
            kernel2_list.append(kernel2)
            kernel3_list.append(kernel3)
            cnt_list.append(cnt)
            init_list.append(init)
            all_cnt_list.append(cnt_all)
            cid += 1
    #####################################
    if kernel2_list == []:
        return None

    if len(kernel2_list) == 1:
        return None
    else:

        kernel2_sum = np.nansum(np.stack(kernel2_list, axis=0), axis=0)
        kernel3_sum = np.nansum(np.stack(kernel3_list, axis=0), axis=0)
        cnt_sum = np.nansum(np.stack(cnt_list, axis=0), axis=0)
        init_sum = np.nansum(np.stack(init_list, axis=0), axis=0)
        allcnt_sum = np.nansum(np.stack(all_cnt_list, axis=0), axis=0)


    print('Returning')
    del lsta
    return (kernel2_sum, kernel3_sum, cnt_sum, allcnt_sum, init_sum)



def plot(rawhour):
    h = rawhour - (MREGIONS[REGION])[2]
    if h >= 24:
        h = h-24
    if h == 24:
        h = 0
    if h < 0:
        h = h+24
    print('Hour: ', h)
    pin = cnst.network_data + 'data/GLOBAL_MCS/save_composites/' + "/" +REGION+"_LSTA_"+SENSOR+"_SWA_2012-2019_MCSTRACK_BOX-ANNUAL_PF_DAY_"+str(MONTHS[0]).zfill(2)+'-'+str(MONTHS[1]).zfill(2)+'_'

    y1 = 2003
    y2 = 2020

    dic = pkl.load(open(
            pin+str(y1)+'_h'+str(rawhour).zfill(2)+extag+".p", "rb"))

    def coll(dic, h, year):
        print(h)

        core = pkl.load(open(
                pin+str(year)+'_h'+str(rawhour).zfill(2)+extag+".p", "rb"))
        for id, k in enumerate(core.keys()):
            try:
                dic[k] = dic[k] + core[k]
            except KeyError:
                dic[k] = core[k]

    for y in range(y1+1, y2):
        print('Coll', y)
        try:
            coll(dic, h, y)
        except:
            print(h, y, 'missing')
            continue

    extent = 60

    f = plt.figure(figsize=(15, 4))
    ax = f.add_subplot(131)

    #thresh = np.percentile(dic['init'], 99)


    plt.contourf(dic['init'], cmap='RdBu_r',  levels=np.arange(0,11), extend='both')
    plt.plot(extent, extent, 'bo')
    ax.axvline(extent, linestyle='dashed', color='k')
    ax.axhline(extent, linestyle='dashed', color='k')
    print('Extent', extent)
    ax.set_xticks((np.linspace(0, 2 * extent, 9)))
    ax.set_xticklabels(((np.linspace(0, (2*extent), 9)-extent)/20).round(1).astype(float))
    #ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
    ax.set_yticks((np.linspace(0, 2 * extent, 9)))
    ax.set_yticklabels(((np.linspace(0, (2*extent), 9)-extent)/20).round(1).astype(float))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.colorbar(label='Count')
    plt.title(REGION+' '+str(y1)+'-'+str(y2)+'| Initiations, Nb: ' + str(np.max(dic['allcnt'])) + '| ' + str(h).zfill(2) + '00UTC',
              fontsize=10)

    ano = (dic['ano'] / dic['cnt']) - np.nanmean((dic['ano'] / dic['cnt']))  # / dic['cnt']
    thresh = np.abs(np.percentile(ano, 95)) #np.linspace(thresh*-1,thresh,20)
    div = 1/thresh #np.array([-1,-0.75,-0.5,-0.25,-0.1,0.1,0.25,0.5,0.75,1])/div
    ax = f.add_subplot(132)

    plt.contourf(ano, cmap='RdBu_r',  levels=np.linspace(thresh*-1,thresh,20), extend='both') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticks((np.linspace(0, 2 * extent, 9)))
    ax.set_xticklabels(((np.linspace(0, (2*extent), 9)-extent)/20).round(1).astype(float))
    ax.set_yticks((np.linspace(0, 2 * extent, 9)))
    ax.set_yticklabels(((np.linspace(0, (2*extent), 9)-extent)/20).round(1).astype(float))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.axvline(extent, linestyle='dashed', color='k')
    ax.axhline(extent, linestyle='dashed', color='k')
    plt.colorbar(label='K')
    plt.title('Seasonal anomaly',
              fontsize=10)


    ax = f.add_subplot(133)

    plt.contourf(dic['cnt'], cmap='viridis') #-(rkernel2_sum / rcnt_sum)
    plt.plot(extent, extent, 'bo')
    ax.set_xticks((np.linspace(0, 2 * extent, 9)))
    ax.set_xticklabels(((np.linspace(0, (2*extent), 9)-extent)/20).round(1).astype(float))
    ax.set_yticks((np.linspace(0, 2 * extent, 9)))
    ax.set_yticklabels(((np.linspace(0, (2*extent), 9)-extent)/20).round(1).astype(float))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.axvline(extent, linestyle='dashed', color='k')
    ax.axhline(extent, linestyle='dashed', color='k')
    plt.colorbar(label='n')
    plt.title('Valid count',
              fontsize=10)


    plt.tight_layout()
    #f.savefig('/home/ck/DIR/cornkle/figs/GLOBAL_MCS/'+REGION+'_'+str(y1)+'-'+str(y2-1)+'_SM_composite_AMSR2-global_BOX-ANNUAL_PF.jpg')
    plt.savefig(cnst.elements_drive + '/'+str(rawhour)+'/' + REGION + '_'+SENSOR+'_LSTA_'+str(MONTHS[0]).zfill(2)+'-'+str(MONTHS[1]).zfill(2)+'_2003-2019_allPF_'+str(rawhour).zfill(2)+'h'+extag+'.jpg')

    plt.close('all')


for regs in REGIONS:
    REGION = regs
    MONTHS = (MREGIONS[REGION])[5]
    composite(19)
    plot(19)