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
    vars = ['q_col','u_col','r_col','v_col','q_col_s','u_col_s','r_col_s','v_col_s','u925_s','u650_s',
           'q925_s','q700_s','tcwv_s','CAPE_s','tcwv','CAPE','dates','tmin','tmean','t10','area',
            'area70','lat','lon','u925','u650','q925','q700']

    for v in vars:
        dic[v] = []
    return dic

def dictionary():

    dic = {}
    vars = ['q_col','u_col','r_col','v_col',
           'tcwv','CAPE','dates','tmin','tmean','t10','area',
            'area70','lat','lon','u925','u650','q925','q700']

    for v in vars:
        dic[v] = []
    return dic


def perSys(clim=False):

    pool = multiprocessing.Pool(processes=4)

    pdf = pkl.load(open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA.p', 'rb'))
    pdf_all = pdf.where((pdf.clat >= 5) & (pdf.clat <= 9) & (pdf.clon >= -12) & (pdf.clon <= 12) &
                        (pdf.year >= 2000) & (pdf.year <= 2014))
    pdf_all = pdf_all.dropna()

    era_pl = xr.open_mfdataset('/home/ck/DIR/mymachine/ERA5/pressure_levels/*.nc')
    era_srfc = xr.open_mfdataset('/home/ck/DIR/mymachine/ERA5/surface/*.nc')

    era_pl = uda.flip_lat(era_pl)
    era_srfc = uda.flip_lat(era_srfc)

    #ipdb.set_trace()

    dates = np.unique(np.array(pdf_all.date))

    inputs = []
    for d in dates:
        pos = np.where(pdf_all.date==d)

        inputs.append((pd.Timestamp(d),pdf_all,era_pl,era_srfc,pos[0],clim))


    res = pool.map(get_ERA5, inputs)
    pool.close()
    # res = []
    # for ins in inputs:
    #     res.append(get_ERA5(ins))

    merged = ua.merge_dicts(res, merge_lists=True)

    #test = pd.DataFrame.from_dict(merged, orient='index')

    pkl.dump(merged, open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA_ERA5_allmonth_2000-2014.p',
                           'wb'))


def get_ERA5(inputs):

    date = inputs[0]
    indic = inputs[1]
    era_pl = inputs[2]
    era_srfc = inputs[3]
    ids = inputs[4]
    clim = inputs[5]

    dic = dictionary_storm()

    for id in ids:

        print('Doing', date)

        if clim:

            time = str(date.year) + '-' + str(date.month) + '-' + '12'
            stormtime = str(date.year) + '-' + str(date.month) + '-' + '18'

        else:
            time = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + 'T12'
            stormtime = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + 'T18'

        try:
            era_day_pl = era_pl.sel(time=time).isel(time=0)
        except (TypeError, IndexError, KeyError):
            print('Era missing:', date)
            #             for k in dic.keys():
            #                 dic[k].append(np.nan)
            continue

        era_day_sf = era_srfc.sel(time=time).isel(time=0)
        era_day_sft = era_srfc.sel(time=stormtime).isel(time=0)
        era_day_plt = era_pl.sel(time=stormtime).isel(time=0)

        elat = indic.clat[id]
        elon = indic.clon[id]
        dic['dates'].append(date)
        dic['lat'].append(elat)
        dic['lon'].append(elon)

        point = era_day_sf.sel(latitude=elat, longitude=elon, method='nearest')

        posx = int(np.where(era_day_sf.longitude == point.longitude)[0])
        posy = int(np.where(era_day_sf.latitude == point.latitude)[0])
        
        xx = 1

        try:
            dic['u925'].append(np.asscalar((era_day_pl['u'].isel(latitude=slice(posy - xx, posy + xx +1),
                                                                 longitude=slice(posx - xx, posx + xx+1)).sel(
                level=925).mean().values)))
        except ValueError:
            continue
            

        dic['u650'].append(np.asscalar((era_day_pl['u'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                             longitude=slice(posx - xx, posx + xx+1)).sel(
            level=650).mean().values)))
        dic['q925'].append(np.asscalar((era_day_pl['q'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                             longitude=slice(posx - xx, posx + xx+1)).sel(
            level=925).mean().values)))
        dic['q700'].append(np.asscalar((era_day_pl['q'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                             longitude=slice(posx - xx, posx + xx+1)).sel(
            level=700).mean().values)))
        dic['CAPE'].append(np.asscalar((era_day_sf['cape'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                                longitude=slice(posx - xx, posx + xx+1)).mean().values)))
        dic['tcwv'].append(np.asscalar((era_day_sf['tcwv'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                                longitude=slice(posx - xx, posx + xx+1)).mean().values)))

        dic['q_col'].append((era_day_pl['q'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                  longitude=slice(posx - xx, posx + xx+1)).mean(
            ['latitude', 'longitude']).values))
        dic['u_col'].append((era_day_pl['u'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                  longitude=slice(posx - xx, posx + xx+1)).mean(
            ['latitude', 'longitude']).values))
        dic['r_col'].append((era_day_pl['r'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                  longitude=slice(posx - xx, posx + xx+1)).mean(
            ['latitude', 'longitude']).values))
        dic['v_col'].append((era_day_pl['v'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                  longitude=slice(posx - xx, posx + xx+1)).mean(
            ['latitude', 'longitude']).values))

        dic['u650_s'].append(np.asscalar((era_day_plt['u'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                                longitude=slice(posx - xx, posx + xx+1)).sel(
            level=650).mean().values)))
        dic['q925_s'].append(np.asscalar((era_day_plt['q'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                                longitude=slice(posx - xx, posx + xx)).sel(
            level=925).mean().values)))
        dic['q700_s'].append(np.asscalar((era_day_plt['q'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                                longitude=slice(posx - xx, posx + xx+1)).sel(
            level=700).mean().values)))
        dic['CAPE_s'].append(np.asscalar((era_day_sft['cape'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                                   longitude=slice(posx - xx, posx + xx+1)).mean().values)))
        dic['tcwv_s'].append(np.asscalar((era_day_sft['tcwv'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                                   longitude=slice(posx - xx, posx + xx+1)).mean().values)))

        dic['q_col_s'].append((era_day_plt['q'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                     longitude=slice(posx - xx, posx + xx+1)).mean(
            ['latitude', 'longitude']).values))
        dic['u_col_s'].append((era_day_plt['u'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                     longitude=slice(posx - xx, posx + xx+1)).mean(
            ['latitude', 'longitude']).values))
        dic['r_col_s'].append((era_day_plt['r'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                     longitude=slice(posx - xx, posx + xx+1)).mean(
            ['latitude', 'longitude']).values))
        dic['v_col_s'].append((era_day_plt['v'].isel(latitude=slice(posy - xx, posy + xx+1),
                                                     longitude=slice(posx - xx, posx + xx+1)).mean(
            ['latitude', 'longitude']).values))

        dic['tmin'].append(indic.tmin[id])
        dic['tmean'].append(indic.tmean[id])
        dic['t10'].append(indic.t10[id])
        dic['area'].append(indic.area[id])
        dic['area70'].append((indic['70area'])[id])
    dic['level'] = era_pl.level.values

    return dic

