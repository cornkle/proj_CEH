# -*- coding: utf-8 -*-


import numpy as np
import datetime as dt
import xarray as xr
import os
import ipdb
import matplotlib.pyplot as plt
from utils import u_grid, u_arrays as ua, u_darrays as uda, u_arrays as ua
import ipdb
import pandas as pd
from utils import constants as cnst
import multiprocessing
import glob
import pickle as pkl


def dictionary_storm():

    dic = {}
    vars = ['q_col','u_col','r_col','v_col', 'd_col', 't_col', 'q_col_s','u_col_s','r_col_s','v_col_s','d_col_s', 't_col_s',
            'u925_s','u650_s',
           'q925_s','q700_s','t2_s', 'divMoist_s','slp_s', 'd925_s','tcwv_s','CAPE_s','tcwv','CAPE','dates','tmin','tmean','t10','area',
            'area70','lat','lon','u925','u650','q925','q700', 't2','divMoist', 'slp', 'd925']

    for v in vars:
        dic[v] = []
    return dic

def dictionary():

    dic = {}
    vars = ['q_col','u_col','r_col','v_col', 'd_col',
           'tcwv','CAPE','dates','tmin','tmean','t10','area',
            'area70','lat','lon','u925','u650','q925','q700']

    for v in vars:
        dic[v] = []
    return dic


def perSys(clim=False):

    pool = multiprocessing.Pool(processes=5)

    pdf = pkl.load(open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA.p', 'rb'))
    pdf_all = pdf.where((pdf.clat >= 5) & (pdf.clat <= 8) & (pdf.clon >= -12) & (pdf.clon <= 10) &
                        (pdf.year >= 2000) & (pdf.year <= 2014))
    pdf_all = pdf_all.dropna()

    #ipdb.set_trace()

    dates = np.unique(np.array(pdf_all.date))
    print(dates)

    inputs = []
    for d in dates:
        pos = np.where(pdf_all.date==d)

        inputs.append((pd.Timestamp(d),pdf_all,pos[0],clim)) #,era_pl,era_srfc,

    #
    res = pool.map(get_ERA5, inputs)
    pool.close()
    # res = []
    # for ins in inputs:
    #     res.append(get_ERA5(ins))

    merged = ua.merge_dicts(res, merge_lists=True)

    #test = pd.DataFrame.from_dict(merged, orient='index')
    if clim:
        pkl.dump(merged, open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA_ERA5_allmonth_5-8N_2000-2014_18UTC_front_CLIM.p',
                           'wb'))
    else:
        pkl.dump(merged, open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA_ERA5_allmonth_5-8N_2000-2014_18UTC_front.p',
                           'wb'))

def get_ERA5(inputs):

    date = inputs[0]
    indic = inputs[1]
    # era_pl = inputs[2]
    # era_srfc = inputs[3]
    ids = inputs[2]
    clim = inputs[3]

    dic = dictionary_storm()

    if clim:
        time = str(date.month).zfill(2) + '-12'
        stormtime = str(date.month).zfill(2) + '-18'

        pl_str = 'ERA5/monthly/synop_selfmade/CLIM_2000-2014/pressure_levels/ERA5_2000-2014_CLIM_'
        srfc_str = 'ERA5/monthly/synop_selfmade/CLIM_2000-2014/surface/ERA5_2000-2014_CLIM_'

        try:
            print('Open '+ cnst.local_data + pl_str + time + '_pl.nc')
            era_day_pl = xr.open_dataset(
                cnst.local_data + pl_str + time + '_pl.nc')
        except (TypeError, IndexError, KeyError):
            print('Era missing:', date)
            #             for k in dic.keys():
            #                 dic[k].append(np.nan)
            return
        dic['level'] = era_day_pl.level.values
        era_day_sf = xr.open_dataset(
            cnst.local_data + srfc_str + time + '_srfc.nc')
        try:
            print('Open ' + cnst.local_data + srfc_str + stormtime + '_srfc.nc')
            era_day_sft = xr.open_dataset(
                cnst.local_data + srfc_str + stormtime + '_srfc.nc')
        except IndexError:
            return
        era_day_plt = xr.open_dataset(
            cnst.local_data + pl_str + stormtime + '_pl.nc')

    else:

        era_pl = xr.open_dataset(cnst.local_data + 'ERA5/hourly/pressure_levels/ERA5_' +str(date.year) + '_' + str(date.month).zfill(2) + '_pl.nc')
        era_srfc = xr.open_dataset(
            cnst.local_data + 'ERA5/hourly/surface/ERA5_' + str(date.year) + '_' + str(date.month).zfill(2) + '_srfc.nc')

        era_pl = era_pl.rename({'longitude' : 'lon', 'latitude' : 'lat'})
        era_srfc = era_srfc.rename({'longitude' : 'lon', 'latitude' : 'lat'})

        time = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + 'T12'
        stormtime = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + 'T18'

        try:
            era_day_pl = era_pl.sel(time=time).isel(time=0)
        except (TypeError, IndexError, KeyError):
            print('Era missing:', date)
            #             for k in dic.keys():
            #                 dic[k].append(np.nan)
            return
        dic['level'] = era_pl.level.values
        era_day_sf = era_srfc.sel(time=time).isel(time=0)
        try:
            era_day_sft = era_srfc.sel(time=stormtime).isel(time=0)
        except IndexError:
            return
        era_day_plt = era_pl.sel(time=stormtime).isel(time=0)


    for id in ids:

        print('Doing', date)


        # elat = indic.clat[id]
        # elon = indic.clon[id]

        elat = indic.clat[id]
        elon = indic.minlon[id]

        dic['dates'].append(date)
        dic['lat'].append(elat)
        dic['lon'].append(elon)
        # ipdb.set_trace()
        point = era_day_pl.sel(lat=elat, lon=elon, method='nearest')

        posx = int(np.where(era_day_sf.lon == point.lon)[0])
        posy = int(np.where(era_day_sf.lat == point.lat)[0])

        posxx = int(np.where(era_day_pl.lon == point.lon)[0])
        posyy = int(np.where(era_day_pl.lat == point.lat)[0])

        #         xxx = 0
        #         yy = slice(posy-xx,posy+xx)
        #         xx = slice(posx-xx,posx+xx)
        yy = posy
        xx = posx

        xxx = posxx
        yyy = posyy

        try:
           # ipdb.set_trace()
            dic['u925'].append(np.asscalar((era_day_pl['u'].isel(lat=yyy, lon=xxx).sel(level=925).mean().values)))
        except ValueError:
            print('Value problem')
            continue
        # ipdb.set_trace()
        dic['u650'].append(np.asscalar((era_day_pl['u'].isel(lat=yyy, lon=xxx).sel(level=650).mean().values)))
        dic['q925'].append(np.asscalar((era_day_pl['q'].isel(lat=yyy, lon=xxx).sel(level=925).mean().values)))
        dic['q700'].append(np.asscalar((era_day_pl['q'].isel(lat=yyy, lon=xxx).sel(level=700).mean().values)))
        dic['d925'].append(np.asscalar((era_day_pl['d'].isel(lat=yyy, lon=xxx).sel(level=925).mean().values)))
        dic['CAPE'].append(np.asscalar((era_day_sf['cape'].isel(lat=yy, lon=xx).mean().values)))
        dic['tcwv'].append(np.asscalar((era_day_sf['tcwv'].isel(lat=yy, lon=xx).mean().values)))
        dic['t2'].append(np.asscalar((era_day_sf['t2m'].isel(lat=yy, lon=xx).mean().values)))
        dic['divMoist'].append(np.asscalar((era_day_sf['p84.162'].isel(lat=yy, lon=xx).mean().values)))
        dic['slp'].append(np.asscalar((era_day_sf['msl'].isel(lat=yy, lon=xx).mean().values)))

        dic['q_col'].append((era_day_pl['q'].isel(lat=yyy, lon=xxx).values))
        dic['u_col'].append((era_day_pl['u'].isel(lat=yyy, lon=xxx).values))
        dic['r_col'].append((era_day_pl['r'].isel(lat=yyy, lon=xxx).values))
        dic['v_col'].append((era_day_pl['v'].isel(lat=yyy, lon=xxx).values))
        dic['t_col'].append((era_day_pl['t'].isel(lat=yyy, lon=xxx).values))
        dic['d_col'].append((era_day_pl['d'].isel(lat=yyy, lon=xxx).values))


        dic['u925_s'].append(
            np.asscalar((era_day_plt['u'].isel(lat=yyy, lon=xxx).sel(level=925).mean().values)))
        dic['u650_s'].append(
            np.asscalar((era_day_plt['u'].isel(lat=yyy, lon=xxx).sel(level=650).mean().values)))
        dic['q925_s'].append(
            np.asscalar((era_day_plt['q'].isel(lat=yyy, lon=xxx).sel(level=925).mean().values)))
        dic['q700_s'].append(
            np.asscalar((era_day_plt['q'].isel(lat=yyy, lon=xxx).sel(level=700).mean().values)))
        dic['d925_s'].append(
            np.asscalar((era_day_plt['d'].isel(lat=yyy, lon=xxx).sel(level=925).mean().values)))
        dic['CAPE_s'].append(np.asscalar((era_day_sft['cape'].isel(lat=yy, lon=xx).mean().values)))
        dic['tcwv_s'].append(np.asscalar((era_day_sft['tcwv'].isel(lat=yy, lon=xx).mean().values)))
        dic['t2_s'].append(np.asscalar((era_day_sft['t2m'].isel(lat=yy, lon=xx).mean().values)))
        dic['divMoist_s'].append(np.asscalar((era_day_sft['p84.162'].isel(lat=yy, lon=xx).mean().values)))
        dic['slp_s'].append(np.asscalar((era_day_sft['msl'].isel(lat=yy, lon=xx).mean().values)))


        dic['q_col_s'].append((era_day_plt['q'].isel(lat=yyy, lon=xxx).values))
        dic['u_col_s'].append((era_day_plt['u'].isel(lat=yyy, lon=xxx).values))
        dic['r_col_s'].append((era_day_plt['r'].isel(lat=yyy, lon=xxx).values))
        dic['v_col_s'].append((era_day_plt['v'].isel(lat=yyy, lon=xxx).values))
        dic['d_col_s'].append((era_day_plt['d'].isel(lat=yyy, lon=xxx).values))
        dic['t_col_s'].append((era_day_plt['t'].isel(lat=yyy, lon=xxx).values))

        dic['tmin'].append(indic.tmin[id])
        dic['tmean'].append(indic.tmean[id])
        dic['t10'].append(indic.t10[id])
        dic['area'].append(indic.area[id])
        dic['area70'].append((indic['70area'])[id])

        print('DID', date)

    return dic

