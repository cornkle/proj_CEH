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
from GLOBAL import glob_util
from utils import u_met
#import metpy
#from metpy import calc
#from metpy.units import units


import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


MREGIONS = glob_util.MREGIONS


REGIONS = ['sub_SA', 'WAf', 'china', 'india', 'australia','SAf', 'GPlains']
SENSOR = 'ERA5'
SENSOP = 12 # (LT overpass)
extag = 'FullYear' #
INIT_DISTANCE = 0
AREA = 5000 #1000
TRACKTIME = 0
OUTPATH = cnst.lmcs_drive + '/ERA5_MCS_saveFiles/trackv2/'
INPATH = cnst.lmcs_drive + '/save_files_v2/'


#TOPO = xr.open_dataarray('/home/ck/DIR/cornkle/data/ancils_python/gtopo_1min.nc').sel(longitude=slice(-110,125), latitude=slice(-50,55))

def dictionary():

    dic = {}
    vars = ['hour', 'month', 'year', 'day', 'date', 'core_area', 'ccs_area', 'pf_area',
            'lon', 'lat', 'clon', 'clat', 'direction',
            'init_lon', 'init_lat', 'tlon', 'tlat',
            'tmin', 'tmean_core', 'tmean_ccs', 'tcwv', 'tgrad2m', 'tgrad925',
            'pmax', 'pmean', 'ptot',
            'q925', 'q650', 'q850', 'era_precip', 'sm'
            'u925', 'u650',
            'v925', 'v650',
            'w925', 'w650',
            'rh925','rh650',
            't925', 't650',
            'div925','div850', 'div650',
            'pv925', 'pv650',
            'ushear925_650', 'ushear850_650', 'ushear925_850',
            'vshear925_650', 'vshear850_650', 'vshear925_850',
            'shear925_650', 'shear850_650', 'shear925_850',
            'cape', 't2m']

    for v in vars:
        dic[v] = []
    
    dummy = pd.read_csv(glob.glob(INPATH + '*.csv')[0])
    for v in dummy.keys():
        dic[v] = []
    
    return dic

def composite(lt_hour):

    h = glob_util.LT_to_UTC_hour(lt_hour, REGION)
    h2 = glob_util.LT_to_UTC_hour(lt_hour+1, REGION)
    h1 = glob_util.LT_to_UTC_hour(lt_hour-1, REGION)

    print('Hour: ', h, h1, h2)

    print(REGION)

    for y in np.arange(2000, 2021):

        m1 = MONTHS[0]
        m2 = MONTHS[1]
        outfile = OUTPATH + REGION + "_" + SENSOR + "_2000-2019_MCSTRACK_localBulk_" + str(m1).zfill(
            2) + '-' + \
                  str(m2).zfill(2) + '_' + str(y) + '_h' + str(lt_hour).zfill(2) + '_' + extag + ".csv"

        msg = pd.read_csv(
            glob.glob(INPATH + REGION + '_mcs_tracks_final_' + str(y) + '*.0000_' + str(y+1) + '0101.0000.csv')[0])
        msg = msg.to_dict(orient='list')

        for k in msg.keys():
            msg[k] = np.array(msg[k])
        inmask =   ((msg['lt_init_hour'] >= SENSOP+2) & (msg['tracktime'] >= TRACKTIME) &  (msg['utc_month'] >= m1) & (msg['utc_month'] <= m2) &\
                   ((msg['lt_hour']==h)  | (msg['lt_hour']==h1)  | (msg['lt_hour']==h2)) & \
                    (msg['pf_landfrac'] > 0.95)  & (msg['pf_area1'] > AREA)) & np.isfinite(msg['pf_lat1'])


        mask = np.where(inmask)
        for k in msg.keys():
            msg[k] = (msg[k])[mask]

        ubt = np.unique(msg['utc_date'])
        msg['utc_date'] = np.array(msg['utc_date'])

        chunks = []
        for ut in ubt:
            daydir = {}

            pos = np.where(msg['utc_date'] == ut)  # [0]).astype(int)
            # ipdb.set_trace()
            for k in msg.keys():
                daydir[k] = np.array(msg[k])[pos]
            chunks.append(daydir)

        try:
            os.remove(outfile)
            print('Removed file')
        except:
            pass

        print('Doing', y)

        #pool = multiprocessing.Pool(processes=3)
        #mdic = dictionary()  # defaultdict(list)
        #res = pool.map(file_loop, chunks)
        #pool.close()

        res = []
        for f in chunks[0:10]:
        
             out = file_loop(f)
             res.append(out)
        
        ipdb.set_trace()

        print('Back from multiproc')
        keys = mdic.keys()
        for vv in res:
            if vv == None:
                continue
            for k in keys:
                try:
                    mdic[k].extend(vv[k])
                except TypeError:
                    ipdb.set_trace()
                    continue

        for kk in mdic.keys():
            print(kk,len(mdic[kk]))

        df = pd.DataFrame.from_dict(mdic)
        df.to_csv(outfile, index=False)
        print('Saved ', outfile)


def file_loop(fi):

    print('Entered loop')
    date = pd.to_datetime(fi['utc_time'][0])

    hour = fi['utc_hour']
    print('Doing day: ', date, hour[0])

    box = (MREGIONS[REGION])[0]
    daybefore = fi['lt_date'][0]
    edate = pd.to_datetime(daybefore)
    edate = edate.replace(hour=MHOUR, minute=0)
    hourchange = (MREGIONS[REGION])[2]


    if hourchange < 0:
        edate = edate + pd.Timedelta(str(np.abs(hourchange)) + ' hours') # translated back to UTC
    elif hourchange > 0:
        edate = edate - pd.Timedelta(str(np.abs(hourchange)) + ' hours')

    print(daybefore)

    try:
        era_pl = xr.open_dataset(
            cnst.lmcs_drive + 'ERA5/hourly/pressure_levels/'+REGION+'/ERA5_' + str(edate.year) + '_' + str(edate.month).zfill(
                2) + '_'+ str(edate.day).zfill(2) + '_'+REGION+'_pl.nc')
    except:
        print('ERA5 missing')
        return
    try:
        era_srfc = xr.open_dataset(
            cnst.lmcs_drive + 'ERA5/hourly/surface/'+REGION+'/ERA5_' + str(edate.year) + '_' + str(edate.month).zfill(
                2) + '_'+ str(edate.day).zfill(2) + '_'+REGION+ '_srfc.nc')
    except:
        print('ERA5 srfc missing')
        return

    try:
        mfile = glob.glob(cnst.lmcs_drive+'/MCS_Feng/global_v2/2d_fields/'+str(date.year)+'*/mcstrack_'+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'*'+'_'+str(date.hour).zfill(2)+'30.nc')[0]

        mcs_2d = xr.open_dataset(mfile)
    except:
        print('MCS 2d field file missing')
        return

    era_pl = u_darrays.flip_lat(era_pl)
    era_srfc = u_darrays.flip_lat(era_srfc)

    edate = edate.replace(hour=12, minute=0).floor('S')
    era_pl_day = era_pl.sel(time=edate)
    era_srfc_day = era_srfc.sel(time=edate)

    # if isinstance(fi['pf_maxrainrate1'], int):
    #     pr = np.array([fi['pf_maxrainrate1']])
    # else:
    pr = fi['pf_maxrainrate1']

    print('Nb per day', pr.size)
    out = dictionary()
    for ids in range(pr.size):

        cloudnumber = fi['cloudnumber'][ids]
        single_mcs = mcs_2d['tb'].where((mcs_2d['cloudnumber']==cloudnumber)).squeeze()
        pos = np.unravel_index(single_mcs.argmin().values, single_mcs.shape)
        print('Minimum MCS temp', single_mcs.values[pos], 'from track', fi['corecold_mintb'])
        tminlat = single_mcs.lat.values[pos[0]]
        tminlon = single_mcs.lon.values[pos[1]]
        ipdb.set_trace() 
        out['tlon'].append(tminlon)
        out['tlat'].append(tminlat)

        elon = fi['pf_lon1'][ids]
        elat = fi['pf_lat1'][ids]

        if (elon < box[0]) | (elon > box[1]) | (elat < box[2]) | (elat > box[3]):
            continue

        ####################

        umean = fi['movement_distance_x']
        vmean = fi['movement_distance_y']

        ws, direct = u_met.u_v_to_ws_wd(umean,vmean)


        ############
        for k in fi.keys():
            out[k] = fi[k][ids]

#        out['init_lon'].append(fi['londiff_loc-init'][ids])
#        out['init_lat'].append(fi['latdiff_loc-init'][ids])
        out['direction'].append(direct)

#        out['hour'].append(fi['lt_hour'][ids])
#        out['month'].append(fi['lt_month'][ids])
#        out['year'].append(fi['lt_year'][ids])
#        out['day'].append(fi['day'][ids])
#        out['date'].append(fi['date'][ids])

#        out['clat'].append(fi['meanlat'][ids])
#        out['clon'].append(fi['meanlon'][ids])

        # era_day = era_pl_day.sel(latitude=elat, longitude=elon, method='nearest')  # take point of minimum T
        # era_day_srfc = era_srfc_day.sel(latitude=elat, longitude=elon, method='nearest')  # take point of minimum T

        era_day = era_pl_day.sel(latitude=slice(elat-0.5, elat+0.5), longitude=slice(elon-0.5,elon+0.5)).mean(['latitude','longitude'])
        era_day_srfc = era_srfc_day.sel(latitude=slice(elat - 0.5, elat + 0.5), longitude=slice(elon - 0.5, elon + 0.5)).mean(['latitude','longitude'])



        tgrad = era_srfc_day['t2m'].sel(latitude=slice(elat-1.5, elat+1.5), longitude=slice(elon-0.5,elon+0.5)).mean('longitude').squeeze.polyfit(dim='latitude', deg=1)
        print('tgradlen', era_srfc_day['t2m'].sel(latitude=slice(elat-1.5, elat+1.5), longitude=slice(elon-0.5,elon+0.5)).mean('longitude').squeeze.size)

        tgrad_pl = era_pl_day['t'].sel(latitude=slice(elat-1.5, elat+1.5), longitude=slice(elon-0.5,elon+0.5)).mean('longitude').squeeze.polyfit(dim='latitude', deg=1)

        e925 = era_day.sel(level=925).mean()
        e650 = era_day.sel(level=650).mean()
        e850 = era_day.sel(level=850).mean()
        srfc = era_day_srfc.mean()

        del era_day
        del era_day_srfc

#        out['ccs_area'].append(fi['ccs_area'][ids])
#        out['core_area'].append(fi['core_area'][ids])
#        out['pf_area'].append(fi['pf_area1'][ids])

       out['tmin'].append(fi['mintb'][ids])
       out['tmean_core'].append(fi['core_meantb'][ids])
       out['tmean_ccs'].append(fi['meantb'][ids])

        out['pmax'].append(fi['pf_maxrainrate1'][ids])
        out['pmean'].append(fi['pf_rainrate1'][ids])
        out['ptot'].append(fi['pf_accumrain1'][ids])
#        out['lon'].append(elon)
#        out['lat'].append(elat)

        try:
            out['q925'].append(float(e925['q']))
        except TypeError:
            return

        out['q650'].append(float(e650['q']))
        out['v925'].append(float(e925['v']))
        out['v650'].append(float(e925['v']))
        out['u925'].append(float(e925['u']))
        out['u650'].append(float(e650['u']))
        out['w925'].append(float(e925['w']))
        out['w650'].append(float(e650['w']))
        out['rh925'].append(float(e925['r']))
        out['rh650'].append(float(e650['r']))
        out['t925'].append(float(e925['t']))
        out['t650'].append(float(e650['t']))
        out['pv925'].append(float(e925['pv']))
        out['pv650'].append(float(e650['pv']))
        out['div925'].append(float(e925['d']))
        out['div650'].append(float(e650['d']))
        out['div850'].append(float(e850['d']))
        out['q850'].append(float(e850['q']))
        out['tcwv'].append(float(srfc['tcwv']))
        out['cape'].append(float(srfc['cape']))
        out['t2m'].append(float(srfc['t2m']))
        out['era_precip'].append(float(srfc['mtpr'])*3600)
        out['sm'].append(float(srfc['swvl1']))

        out['tgrad2m'].append(float(tgrad))
        out['tgrad925'].append(float(tgrad_pl))

        u9256 = float(e650['u'] - e925['u'])
        u8506 =float(e650['u'] - e850['u'])
        u8509 = float(e850['u'] - e925['u'])

        v9256 = float(e650['v'] - e925['v'])
        v8506 =float(e650['v'] - e850['v'])
        v8509 = float(e850['v'] - e925['v'])

        out['ushear925_650'].append(u9256)
        out['ushear850_650'].append(u8506)
        out['ushear925_850'].append(u8509)

        out['vshear925_650'].append(v9256)
        out['vshear850_650'].append(v8506)
        out['vshear925_850'].append(v8509)

        out['shear925_650'].append(np.sqrt(u9256**2+v9256**2))
        out['shear850_650'].append(np.sqrt(u8506**2+v8506**2))
        out['shear925_850'].append(np.sqrt(u8509**2+v8509**2))

    return out



for regs in REGIONS:
    REGION = regs
    MONTHS = (MREGIONS[REGION])[4]
    MHOUR = 17
    composite(MHOUR)
