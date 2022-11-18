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


MREGIONS = {'WAf' : [[-18,25,4,25], 'spac', 0, (1,7), (8,12), (1,12)], # last is hourly offset to UCT # 12    # [-18,25,4,25]
 'SAf' : [[20,35, -35,-15], 'spac', 2, (9,12), (1,5), (1,12)], # 10
 'india' : [[70,90, 5,30], 'asia', 5, (1,7), (8,12), (1,12)], # 7
 'china' : [[105,115,25,40], 'asia', 8 , (1,7), (8,12), (1,12)], # 4
 'australia' : [[120,140,-23, -11], 'asia', 9, (10,12), (1,5), (1,12)], # 3
 'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4, (9,12), (1,5), (1,12)] , # 16
 'trop_SA' : [[-75, -50, -20, -5], 'spac', -5, (1,12), (1,12), (1,12)], # 17
 'GPlains' : [[-100,-90,32,47], 'nam', -6, (1,7), (8,12), (1,12)] # # 18

}

OUT = '/home/ck/DIR/cornkle/figs/GLOBAL_MCS/'

REGIONS = ['SAf','GPlains', 'sub_SA', 'WAf', 'china', 'india', 'australia']
SENSOR = 'AMSR2'
SENSOP = 13 # (LT overpass sensor)
extag = 'propagate_AMSR_Smfiltercheck'
INIT_DISTANCE = 2 # degrees
AREA = 1000 #1000
TRACKTIME = 1 # tracktime = 0 is initiation

TOPO = xr.open_dataarray('/home/ck/DIR/cornkle/data/ancils_python/gtopo_1min.nc').sel(lon=slice(-110,125), lat=slice(-50,55))


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

    for y in np.arange(2012, 2020):

        m1 = MONTHS[0]
        m2 = MONTHS[1]
        outfile = path + "/"+REGION+"_SM_"+SENSOR+"_SWA_2012-2019_MCSTRACK_BOX-ANNUAL_PF_DAY_"+str(m1).zfill(2)+'-'+\
                  str(m2).zfill(2)+'_'+str(y)+'_h'+str(rawhour).zfill(2)+'_'+extag+".p"

        msg = pd.read_csv(cnst.lmcs_drive +'/save_files/'+REGION+'_initTime__mcs_tracks_extc_'+str(y)+'0101_'+str(y)+'1231.csv')
        msg = msg.to_dict(orient='list')

        domain = (MREGIONS[REGION])[0]

        msg['date'] = []
        msg['utc_date'] = []
        msg['day'] = []
        msg['lt_hour'] =[]
        msg['lt_init'] = []
        for yi, k in enumerate(msg['base_time']):

            bbt = pd.to_datetime(k) #- pd.Timedelta('1 days')
            lt_init = bbt.replace(hour=msg['init_hour'][yi], minute=0) # real date for lt_init may well be previous UTC day, lt_init just used for hour calculation to LT init hour
            hourchange = (MREGIONS[REGION])[2]

            if hourchange < 0:
                bbl = bbt - pd.Timedelta(str(np.abs(hourchange))+' hours')
                binit = lt_init - pd.Timedelta(str(np.abs(hourchange))+' hours')
            elif hourchange > 0:
                    bbl = bbt + pd.Timedelta(str(np.abs(hourchange)) + ' hours')
                    binit = lt_init + pd.Timedelta(str(np.abs(hourchange)) + ' hours')
            else:
                bbl = bbt
                binit = lt_init

            msg['utc_date'].append(bbt.replace(hour=0, minute=0))
            msg['lt_hour'].append(bbl.hour)
            msg['lt_init'].append(binit.hour)
            msg['date'].append(bbl.replace(hour=0, minute=0))
            msg['day'].append(bbt.day)
        #ipdb.set_trace()
        for k in msg.keys():

            msg[k] = np.array(msg[k])

        inmask =   ((msg['lt_init'] >= SENSOP+2) & (msg['tracktime'] >= TRACKTIME)  & (msg['ccs_area'] >= AREA) & \
                   ((msg['hour']==h)  | (msg['hour']==h1)  | (msg['hour']==h2)) & \
                    (msg['pf_landfrac'] > 0.99)  & (msg['month']>=m1) & (msg['month']<=m2) & \
                    ((np.abs(msg['londiff_loc-init'])>=INIT_DISTANCE)|(np.abs(msg['latdiff_loc-init'])>=INIT_DISTANCE)))


        if (np.sum(inmask) == 0) & (REGION in ['GPlains']):
            inmask = ((msg['lt_init'] >= SENSOP+2)  & (msg['tracktime'] >= TRACKTIME) & (msg['ccs_area'] >= AREA) & \
                      ((msg['hour'] == h) | (msg['hour'] == h1) | (msg['hour'] == h2)) & \
                       (msg['month'] >= m1) & (msg['month'] <= m2) & \
                      ((np.abs(msg['londiff_loc-init']) >=INIT_DISTANCE) | (np.abs(msg['latdiff_loc-init']) >=INIT_DISTANCE)))

        #(msg['lt_init'] >= SENSOP+2) & & (msg['pf_mcsstatus'] > 0)

        # inmask =   ( (msg['tracktime'] >= TRACKTIME) & \
        #            ((msg['hour']==h)  | (msg['hour']==h1)  | (msg['hour']==h2)) & \
        #             (msg['pf_landfrac'] > 0.99)  & np.isfinite(msg['pf_lat1']) & \
        #             ((np.abs(msg['londiff_loc-init'])>=INIT_DISTANCE)|(np.abs(msg['latdiff_loc-init'])>=INIT_DISTANCE)))
        #
        #
        #
        # if (np.sum(inmask) == 0) & (REGION in ['GPlains']):
        #     inmask = ((msg['tracktime'] >= TRACKTIME) & \
        #               ((msg['hour'] == h) | (msg['hour'] == h1) | (msg['hour'] == h2)) & \
        #                np.isfinite(msg['pf_lat1']) & \
        #             ((np.abs(msg['londiff_loc-init'])>=INIT_DISTANCE)|(np.abs(msg['latdiff_loc-init'])>=INIT_DISTANCE)))


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

            pos = np.where(msg['date'] == ut)#[0]).astype(int)
            #ipdb.set_trace()
            for k in msg.keys():
                 daydir[k] = np.array(msg[k])[pos]
            chunks.append(daydir)

        # try:
        #     os.remove(outfile)
        #     print('Removed file')
        # except:
        #     pass

        if os.path.isfile(outfile):
            print('File exists, continue')
            continue


        print('Doing', y)

        try:
            dic = u_parallelise.run_arrays(4, file_loop, chunks, ['ano', 'regional', 'cnt', 'allcnt', 'ano_ts', 'regional_ts', 'cnt_ts', 'init'])
        except:
            print('##########NO DATA WAS SAVED, RETURN#####', y)
            continue
        # res = []
        #
        # for m in chunks:
        #     out = file_loop(m)
        #     res.append(out)
        # ipdb.set_trace()
        # return

        for k in dic.keys():
           dic[k] = np.nansum(dic[k], axis=0)
        print('File written', y)
        #pkl.dump(dic, open(path + "/"+REGION+"_AMSR2-global_2012-2019_SM_MCSTRACK_BOX-ANNUAL_PFbased_"+str(y)+'_h'+str(hour).zfill(2)+".p", "wb"))
        pkl.dump(dic, open(outfile, "wb"))





def cut_kernel(xpos, ypos, arr, inits, wd = None):

    dist = 16
    #dist = 35
    #dist=100

    res = 0.25
    #res = 0.1

    kernel = ua.cut_kernel(arr,xpos, ypos,dist)
    if np.sum(np.isfinite(kernel)) < 0.05 * kernel.size:
        print('Not enough data')
        return None

    ilon = np.round(inits[0] / res).astype(int)
    ilat = np.round(inits[1] / res).astype(int)

    # if (np.abs(ilon)<4) & (np.abs(ilat)<4):
    #     return None

    #print('lonlat', ilon, ilat)

    cnt = np.zeros_like(kernel)
    init = np.zeros_like(kernel)

    try:
        init[dist+1-ilat, dist+1-ilon] = 1
    except:
        #print('Init index error, pass')
        pass

    if kernel.shape != (dist*2+1, dist*2+1):
        print('Kernel shape error, pass')
        return None

    if wd != None:
        try:
            kernel = ua.rotate(kernel, wd, ref_angle=90)
        except:
            ipdb.set_trace()
        init = ua.rotate(init, wd, ref_angle=90)

    cnt[np.isfinite(kernel)] = 1
    kernel3 = kernel - np.nanmean(kernel)

    cut = 5
    #print('shape', kernel[cut:-cut, cut:-cut].shape)
    return kernel[cut:-cut, cut:-cut], kernel3[cut:-cut,cut:-cut], cnt[cut:-cut,cut:-cut], init[cut:-cut,cut:-cut]


def file_loop(fi):

    print('Entering func')
    date = fi['date'][0]

    hour = fi['hour']

    print('Doing day: ', date, hour[0])

    box = (MREGIONS[REGION])[0]

    daybefore = date

    print(daybefore)

    fdate = str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2)
    edate = str(daybefore.year) +'_' + str(daybefore.month).zfill(2) + '_' + str(daybefore.day).zfill(2)

    alt_path = cnst.lmcs_drive + 'AMSR2/daily/25km/day_anom/'
    try:
        lsta = xr.open_dataset(alt_path + 'amsr2_25km_anom_' + fdate + '.nc')
    except:
        print('Surface file not found, return')
        return None

    lsta = lsta.sel(time=str(daybefore.year) + '-' + str(daybefore.month) + '-' + str(daybefore.day))
    lsta = lsta.sel(lon=slice(box[0]-3, box[1]+3), lat=slice(box[2]-3, box[3]+3)) # define SM region box to consider

    print('Doing ' + 'AMSR2_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc', date)


    lsta_da = lsta['soil_moisture_c1'].squeeze()  # soil_moisture_c1
    ts_da = lsta['ts'].squeeze()

    mask = (lsta_da<-1) & (ts_da<-0.5)
    lsta_da.values[mask] = np.nan
    ts_da.values[mask] = np.nan
    #ipdb.set_trace()
    # try:
    #  era_pl = xr.open_dataset('/media/ck/LStorage/ERA5/hourly/pressure_levels/'+REGION+'/ERA5_'+edate+'_'+REGION+'_'+'pl.nc')
    # except:
    #     print('ERA5 file not found, return')
    #     return None
    # ttopo = lsta.salem.lookup_transform(TOPO, lut=ILUT)
    #
    # grad = np.gradient(ttopo.values)
    # gradsum = abs(grad[0]) + abs(grad[1])
    # lsta_da.values[ttopo.values >= 500] = np.nan
    # lsta_da.values[gradsum > 30] = np.nan

    kernel2_list = []
    kernel3_list = []
    ts2_list = []
    ts3_list = []
    cnt_list = []
    cntts_list = []
    init_list = []
    all_cnt_list = []

    cid = 0

############################

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
            #ipdb.set_trace()


            init_lon = fi['londiff_loc-init'][idss]
            init_lat = fi['latdiff_loc-init'][idss]

            try:
                point = lsta_da.sel(lat=lat, lon=lon, method='nearest', tolerance=0.2)
            except KeyError:
                print('point index error')
                idss = idss + 1
                continue

            # try:
            #     umean = era_pl['u'].sel(time=(era_pl['time.hour']==hour[idss]), level=650).sel(longitude=lon, latitude=lat, method='nearest', tolerance=0.2).values
            # except:
            #     idss = idss + 1
            #     continue
            # vmean = era_pl['v'].sel(time=(era_pl['time.hour']==hour[idss]), level=650).sel(longitude=lon, latitude=lat, method='nearest', tolerance=0.2).values

            avg_wd = np.rad2deg(np.arctan2(umean, vmean))
            if avg_wd < 0:
                avg_wd = avg_wd + 360
            # print('mean WD', avg_wd, fi['direction'])
            direct = avg_wd#[0]

            idss = idss+1

            plat = point['lat'].values
            plon = point['lon'].values

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

            try:
                tskernel2, tskernel3, cntts, init = cut_kernel(xpos, ypos, ts_da, [init_lon,init_lat], wd=direct)
            except TypeError:
                continue

            cnt_all = np.zeros_like(cnt)+1
            kernel2_list.append(kernel2)
            kernel3_list.append(kernel3)
            ts2_list.append(tskernel2)
            ts3_list.append(tskernel3)
            cnt_list.append(cnt)
            cntts_list.append(cntts)
            init_list.append(init)
            all_cnt_list.append(cnt_all)

            cid += 1

    if kernel2_list == []:
        return None

    if len(kernel2_list) == 1:
        return None
    else:

        kernel2_sum = np.nansum(np.stack(kernel2_list, axis=0), axis=0)
        kernel3_sum = np.nansum(np.stack(kernel3_list, axis=0), axis=0)
        ts2_sum = np.nansum(np.stack(ts2_list, axis=0), axis=0)
        ts3_sum = np.nansum(np.stack(ts3_list, axis=0), axis=0)
        cnt_sum = np.nansum(np.stack(cnt_list, axis=0), axis=0)
        cntts_sum = np.nansum(np.stack(cntts_list, axis=0), axis=0)
        init_sum = np.nansum(np.stack(init_list, axis=0), axis=0)
        allcnt_sum = np.nansum(np.stack(all_cnt_list, axis=0), axis=0)

    print('Returning')
    del lsta

    return (kernel2_sum, kernel3_sum, cnt_sum, allcnt_sum, ts2_sum, ts3_sum, cntts_sum, init_sum)


def plot_sm_ts(rawhour):
        h = rawhour - (MREGIONS[REGION])[2]
        print('Hour: ', h)
        pin = cnst.network_data + 'data/GLOBAL_MCS/save_composites/' + "/" + REGION + "_SM_" + SENSOR + "_SWA_2012-2019_MCSTRACK_BOX-ANNUAL_PF_DAY_" + str(
        MONTHS[0]).zfill(2) + '-' + str(MONTHS[1]).zfill(2) + '_'
        y1 = 2013
        y2 = 2020

        dic = pkl.load(open(
            pin + str(y1) + '_h' + str(rawhour).zfill(2) + extag+".p", "rb"))

        def coll(dic, h, year):
            print(h)
            core = pkl.load(open(
                pin + str(year) + '_h' + str(rawhour).zfill(2) + extag+ ".p", "rb"))
            for id, k in enumerate(core.keys()):
                try:
                    dic[k] = dic[k] + core[k]
                except KeyError:
                    dic[k] = core[k]

        for y in range(y1 + 1, y2):
            print('Coll', y)
            coll(dic, h, y)

        extent = 11

        f = plt.figure(figsize=(15, 8))
        ax = f.add_subplot(231)

        ano = (dic['ano']/ dic['cnt'] ) -np.mean((dic['ano'])/ dic['cnt'])  #/ dic['cnt']

        thresh = np.max(np.abs(np.percentile(ano, [5, 95])))

        plt.contourf(dic['init'], cmap='RdBu', levels=np.arange(0,30,2),
                     extend='both')
        plt.plot(extent, extent, 'bo')
        ax.axvline(extent, linestyle='dashed', color='k')
        ax.axhline(extent, linestyle='dashed', color='k')
        print('Extent', extent)
        ax.set_xticks((np.linspace(0, 2 * extent, 9)))
        ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        # ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_yticks((np.linspace(0, 2 * extent, 9)))
        ax.set_yticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_xlabel('km')
        ax.set_ylabel('km')
        plt.colorbar(label='%')
        plt.title(REGION + ' ' + str(y1) + '-' + str(y2) + ' initiations, Nb: ' + str(np.max(dic['allcnt'])) + '| ' + str(
            h).zfill(2) + '00UTC',
                  fontsize=10)

        ax = f.add_subplot(232)

        plt.contourf(ano, cmap='RdBu', levels=np.linspace(thresh * -1, thresh, 10),
                     extend='both')  # -(rkernel2_sum / rcnt_sum)
        plt.plot(extent, extent, 'bo')
        ax.set_xticks((np.linspace(0, 2 * extent, 9)))
        ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_yticks((np.linspace(0, 2 * extent, 9)))
        ax.set_yticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_xlabel('km')
        ax.set_ylabel('km')
        ax.axvline(extent, linestyle='dashed', color='k')
        ax.axhline(extent, linestyle='dashed', color='k')
        plt.colorbar(label='%')
        plt.title('Spatial anomaly',
                  fontsize=10)

        ax = f.add_subplot(233)

        plt.contourf(dic['cnt'], cmap='viridis')  # -(rkernel2_sum / rcnt_sum)
        plt.plot(extent, extent, 'bo')
        ax.set_xticks((np.linspace(0, 2 * extent, 9)))
        ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_yticks((np.linspace(0, 2 * extent, 9)))
        ax.set_yticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_xlabel('Degrees')
        ax.set_ylabel('Degrees')
        ax.axvline(extent, linestyle='dashed', color='k')
        ax.axhline(extent, linestyle='dashed', color='k')
        plt.colorbar(label='n')
        plt.title('Valid count',
                  fontsize=10)


        ###########################

        ax = f.add_subplot(234)

        ano_ts = (dic['ano_ts'] / dic['cnt']) - np.mean((dic['ano_ts'] / dic['cnt']))
        thresh = np.max(np.abs(np.percentile(ano_ts, [5, 95])))



        plt.contourf(dic['init'], cmap='RdBu',levels=np.arange(0,30,2),
                     extend='both')
        plt.plot(extent, extent, 'bo')
        ax.axvline(extent, linestyle='dashed', color='k')
        ax.axhline(extent, linestyle='dashed', color='k')
        print('Extent', extent)
        ax.set_xticks((np.linspace(0, 2 * extent, 9)))
        ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        # ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_yticks((np.linspace(0, 2 * extent, 9)))
        ax.set_yticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_xlabel('km')
        ax.set_ylabel('km')
        plt.colorbar(label='%')
        plt.title(REGION + ' ' + str(y1) + '-' + str(y2) + ' initiations, Nb: ' + str(np.max(dic['allcnt'])) + '| ' + str(
            h).zfill(2) + '00UTC',
                  fontsize=10)

        ax = f.add_subplot(235)

        plt.contourf(ano_ts, cmap='RdBu_r', levels=np.linspace(thresh * -1, thresh, 10),
                     extend='both')  # -(rkernel2_sum / rcnt_sum)
        plt.plot(extent, extent, 'bo')
        ax.set_xticks((np.linspace(0, 2 * extent, 9)))
        ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_yticks((np.linspace(0, 2 * extent, 9)))
        ax.set_yticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_xlabel('km')
        ax.set_ylabel('km')
        ax.axvline(extent, linestyle='dashed', color='k')
        ax.axhline(extent, linestyle='dashed', color='k')
        plt.colorbar(label='%')
        plt.title('Spatial anomaly',
                  fontsize=10)

        ax = f.add_subplot(236)

        plt.contourf(dic['cnt_ts'], cmap='viridis')  # -(rkernel2_sum / rcnt_sum)
        plt.plot(extent, extent, 'bo')
        ax.set_xticks((np.linspace(0, 2 * extent, 9)))
        ax.set_xticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_yticks((np.linspace(0, 2 * extent, 9)))
        ax.set_yticklabels(((np.linspace(0, (2 * extent), 9) - extent) * 30).astype(int))
        ax.set_xlabel('Degrees')
        ax.set_ylabel('Degrees')
        ax.axvline(extent, linestyle='dashed', color='k')
        ax.axhline(extent, linestyle='dashed', color='k')
        plt.colorbar(label='n')
        plt.title('Valid count',
                  fontsize=10)

        plt.tight_layout()


        plt.savefig(OUT+ '/'+str(rawhour) + '/' + REGION + '_' + SENSOR + '_SM_TS_' + str(MONTHS[0]).zfill(2) + '-' + str(
            MONTHS[1]).zfill(2) + '_2012-2019_allPF_' + str(rawhour).zfill(2) + 'h'+ extag+'.jpg')
        plt.close('all')



for regs in REGIONS:
    REGION = regs
    MONTHS = (MREGIONS[REGION])[5]

    alt_path = cnst.lmcs_drive + 'AMSR2/daily/25km/day_anom/'
    bbox = (MREGIONS[REGION])[0]

    #
    # lsta = xr.open_dataset(glob.glob(alt_path + 'amsr2_25km_anom_*.nc')[0]).isel(time=0)
    # lsta = lsta.sel(lon=slice(bbox[0] - 3, bbox[1] + 3), lat=slice(bbox[2] - 3, bbox[3] + 3))
    #
    # #ilut_file ='/home/ck/DIR/cornkle/data/GLOBAL_MCS/save_topo/'+REGION+'_topo.p'
    # # if os.path.isfile(ilut_file):
    # #     ilut = pkl.load(open(ilut_file, 'rb'))
    # #     try:
    # #         ttopo = lsta.salem.lookup_transform(TOPO, lut=ilut)
    # #     except:
    # #         print('Topo file fail, writing new file')
    # #         ttopo, ilut = lsta.salem.lookup_transform(TOPO, return_lut=True)
    # #         pkl.dump(ilut, open(ilut_file, 'wb'))
    # # else:
    # ttopo, ilut = lsta.salem.lookup_transform(TOPO, return_lut=True)
    # #pkl.dump(ilut, open(ilut_file, 'wb'))
    #
    # ILUT = ilut

    composite(19)
    #plot_sm_ts(19)