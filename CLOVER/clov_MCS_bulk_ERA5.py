import numpy as np
import xarray as xr
from utils import u_arrays as ua, u_darrays as uda
import matplotlib.pyplot as plt
import multiprocessing
import pickle as pkl
from collections import defaultdict
from utils import constants as cnst


import pdb
import glob
import ipdb
import pandas as pd

def dictionary():

    dic = {}
    vars = ['hour', 'month', 'year', 'day', 'date', 'area',
            'lon', 'lat', 'clon', 'clat', 'minlon', 'minlat',
            'tmin', 'tmean', 'tcwv', #'tgrad', 'tbox',
            'pmax', 'pmean',
            'q925', 'q650',
            'u925', 'u650',
            'v925', 'v650',
            'w925', 'w650',
            'rh925','rh650',
            't925', 't650',
            'div925','div650',
            'pv925', 'pv650',
            'shear',
            'pgt30', 'pgt01'
            'isvalid',
             't', 'p' ,
            'cape', 't2m']

    for v in vars:
        dic[v] = []
    return dic

def perSys():

    pool = multiprocessing.Pool(processes=3)
    tthresh = '-50'
    files = glob.glob(cnst.network_data + 'MCSfiles/WA5000_5-25N_17W-15E_-50_afternoon_GPM/*.nc')

    #ipdb.set_trace()

    print('Nb files', len(files))
    for y in range(2016,2019):

        yfiles = []
        for f in files:
            if str(y) in f:
                yfiles.append(f)

        pool = multiprocessing.Pool(processes=3)
        mdic = dictionary() #defaultdict(list)
        res = pool.map(file_loop, yfiles)
        pool.close()

        # res = []
        # for f in files[500:1000]:
        #     out = file_loop(f)
        #     res.append(out)
        #
        #
        # ipdb.set_trace()

        print('Back from multiproc')
        keys = mdic.keys()
        for v in res:
            for k in keys:
                try:
                    mdic[k].append(v[k])
                except TypeError:
                    continue


        pkl.dump(mdic, open(cnst.network_data + 'data/CLOVER/saves/bulk_'+tthresh+'_5000km2_GPM_5-26N_16W17E_p15_ERA0.7_TCWV_hourly_wDate_'+str(y)+'.p',
                               'wb'))


def file_loop(f):
    print('Doing file: ' + f)

    dic = xr.open_dataset(f)
    edate = pd.Timestamp(dic.time.values)

    out = dictionary()
    res = []
    outt = dic['tc_lag0'].values
    outp = dic['p'].values

    out['lon'] = dic['lon'].values
    out['lat'] = dic['lat'].values
    out['hour'] = dic['time.hour'].item()
    out['month'] = dic['time.month'].item()
    out['year'] = dic['time.year'].item()
    out['day'] = dic['time.day'].item()
    out['date'] = dic['time'].values

    if np.nanmin(dic['tc_lag0'].values) > -53:
        return
    #ipdb.set_trace()
    out['clat'] = np.min(out['lat'])+((np.max(out['lat'])-np.min(out['lat']))*0.5)
    out['clon'] = np.min(out['lon']) + ((np.max(out['lon']) - np.min(out['lon'])) * 0.5)

    if (out['clat']<4) | (out['clon']<-18) | (out['clon']>25):
        print('MCS out of box')
        return


    # if edate.hour < 17:
    #     return

    try:
        era_pl = xr.open_dataset(cnst.ERA5_HOURLY_PL+'ERA5_'+str(dic['time.year'].values)+'_'+str(dic['time.month'].values).zfill(2)+'_pl.nc')
    except:
        print('ERA5 missing')
        return
    try:
        era_srfc = xr.open_dataset(cnst.ERA5_HOURLY_SRFC+'ERA5_'+str(dic['time.year'].values)+'_'+str(dic['time.month'].values).zfill(2)+'_srfc.nc')
    except:
        print('ERA5 srfc missing')
        return
    era_pl = uda.flip_lat(era_pl)
    era_srfc = uda.flip_lat(era_srfc)

    edate = edate.replace(hour=12, minute=0)

    era_pl_day = era_pl.sel(time=edate, longitude=slice(-17,20), latitude=slice(3,28))
    era_srfc_day = era_srfc.sel(time=edate, longitude=slice(-17,20), latitude=slice(3,28))


    tminpos = np.where(dic['tc_lag0'].values == np.nanmin(dic['tc_lag0'].values)) # era position close to min temp
    if len(tminpos[0])>1:
        ptmax = np.nanmax((dic['p'].values)[tminpos])
        if ptmax > 0:
            prpos = np.where((dic['p'].values)[tminpos] == ptmax)
            tminpos = ((tminpos[0])[prpos], (tminpos[1])[prpos] )
        else:
            tminpos = ((tminpos[0])[0], (tminpos[1])[0])

    elon = dic['lon'].values[tminpos]
    elat = dic['lat'].values[tminpos]

    out['minlon'] = elon
    out['minlat'] = elat

    era_day = era_pl_day.sel(latitude=elat, longitude=elon , method='nearest') # take point of minimum T
    era_day_srfc = era_srfc_day.sel(latitude=elat, longitude=elon , method='nearest') # take point of minimum T

    del era_srfc_day

    e925 = era_day.sel(level=925).mean()

    e850 = era_pl_day['t'].sel(level=850)
    elow = era_day.sel(level=slice(925,850)).mean('level').mean()
    e650 = era_day.sel(level=650).mean()
    emid = era_day.sel(level=slice(600,700)).mean('level').mean()
    srfc = era_day_srfc.mean()


    t_thresh = -50  # -40C ~ 167 W m-2
    mask = np.isfinite(outp) & (outt<=t_thresh) & np.isfinite(outt)
    mask_area = (outt<=t_thresh) & np.isfinite(outt)
    mask70 = (outt<=-70) & np.isfinite(outt)

    if np.sum(mask) < 3:
        return

    print(np.nanmax(outt[mask]))   # can be bigger than cutout threshold because of interpolation to 5km grid after cutout

    out['area'] = np.sum(mask_area)
    out['area70'] = np.sum(mask70)

    out['tmin'] = np.min(outt[mask])
    out['tmean'] = np.mean(outt[mask])

    maxpos = np.unravel_index(np.nanargmax(outp), outp.shape)
    out['pmax'] = np.nanmean(ua.cut_kernel(outp,maxpos[1], maxpos[0],1)) #np.max(outp[mask])
    out['pmean'] = np.mean(outp[mask])

    dbox = e850.copy(deep=True)
    minlon = era_pl_day.sel(latitude=8, longitude=np.min(out['lon']), method='nearest')
    maxlon = era_pl_day.sel(latitude=8, longitude=np.max(out['lon']), method='nearest')

    del era_pl_day

    tgrad = dbox.sel(longitude=slice(minlon.longitude.values, maxlon.longitude.values)).mean('longitude')

    tmin = np.nanargmin(tgrad.values)
    tmax = np.nanargmax(tgrad.values)
    tgrad = tgrad.isel(latitude=slice(tmin, tmax))

    lingress = uda.linear_trend_lingress(tgrad)

    out['tgrad'] = lingress['slope'].values

    tgrad2 = dbox.sel(longitude=slice(np.min(out['lon']), np.max(out['lon'])), latitude=slice(10, 20)).mean(
        ['longitude', 'latitude']) - \
             dbox.sel(longitude=slice(np.min(out['lon']), np.max(out['lon'])), latitude=slice(5, 7)).mean(['longitude', 'latitude'])
    out['tbox'] = tgrad2.values

    try:
        out['q925'] =float(e925['q'])
    except TypeError:
        return

    out['q650'] = float(e650['q'])
    out['v925'] = float(e925['v'])
    out['v650'] = float(e925['v'])
    out['u925'] = float(e925['u'])
    out['u650'] = float(e650['u'])
    out['w925'] = float(e925['w'])
    out['w650'] = float(e650['w'])
    out['rh925'] = float(e925['r'])
    out['rh650'] = float(e650['r'])
    out['t925'] = float(e925['t'])
    out['t650'] = float(e650['t'])
    out['pv925'] = float(e925['pv'])
    out['pv650'] = float(e650['pv'])
    out['div925'] = float(e925['d'])
    out['div650'] = float(e650['d'])
    out['q925'] = float(e925['q'])
    out['q650'] = float(e650['q'])
    out['tcwv'] = float(srfc['tcwv'])
    out['cape'] = float(srfc['cape'])
    out['t2m'] = float(srfc['t2m'])

    out['shear'] = float(e650['u']-e925['u'])

    # theta_down = u_met.theta_e(925,e925['t']-273.15, e925['q'])
    # theta_up = u_met.theta_e(650,e650['t']-273.15, e650['q'])
    #
    # out['dtheta'] =  (theta_down-theta_up).values
    # out['thetaup'] = theta_up.values
    # out['thetadown'] = theta_down.values

    out['pgt30'] = np.sum(outp[mask]>=30)
    out['isvalid'] = np.sum(mask)
    out['pgt01'] = np.sum(outp[mask]>=0.1)
    #
    out['p'] = outp[mask]
    out['t'] = outt[mask]
    #ipdb.set_trace()
    dic.close()

    return out
