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
    vars = ['q','u','r','v','q_s','u_s','r_s','v_s','tcwv_s','CAPE_s','tcwv','CAPE','dates','tmin','tmean','t10','area',
            'area70','clon','lon']

    for v in vars:
        dic[v] = []
    return dic

def dictionary():

    dic = {}
    vars = ['q','u','r','v',
           'tcwv','CAPE','dates','tmin','tmean','t10','area',
            'area70','lat','lon']

    for v in vars:
        dic[v] = []
    return dic


def perSys(clim=False):

    pool = multiprocessing.Pool(processes=5)

    pdf = pkl.load(open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA.p', 'rb'))
    pdf_all = pdf.where((pdf.clat >= 5) & (pdf.clat <= 9) & (pdf.clon >= -12) & (pdf.clon <= 12) &
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

    pkl.dump(merged, open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA_ERA5_allmonth_2000-2014_18UTC_cross_CLIM.p',
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

        pl_str = 'ERA5/CLIM_2000-2014/pressure_levels/ERA5_2000-2014_CLIM_'
        srfc_str = 'ERA5/CLIM_2000-2014/surface/ERA5_2000-2014_CLIM_'

        try:
            print('Open '+ cnst.local_data + pl_str + time + '_pl.nc')
            era_day_pl = xr.open_dataset(
                cnst.local_data + pl_str + time + '_pl.nc')
            era_day_pl = era_day_pl.isel(time=0)
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

        era_pl = xr.open_dataset(cnst.local_data + 'ERA5/pressure_levels/ERA5_' +str(date.year) + '_' + str(date.month).zfill(2) + '_pl.nc')
        era_srfc = xr.open_dataset(
            cnst.local_data + 'ERA5/surface/ERA5_' + str(date.year) + '_' + str(date.month).zfill(2) + '_srfc.nc')

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
        clon = indic.clon[id]

        dic['dates'].append(date)
        dic['lat'].append(elat)
        dic['lon'].append(elon)
        dic['clon'].append(clon)
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
           dic['CAPE'].append(np.array((era_day_sf['cape'].isel(lon=xx).values)))
        except ValueError:
            print('Value problem')
            continue
        dic['tcwv'].append(np.array((era_day_sf['tcwv'].isel(lon=xx).values)))

        dic['q'].append((era_day_pl['q'].isel(lon=xxx).values))
        dic['u'].append((era_day_pl['u'].isel(lon=xxx).values))
        dic['r'].append((era_day_pl['r'].isel(lon=xxx).values))
        dic['v'].append((era_day_pl['v'].isel(lon=xxx).values))


        dic['CAPE_s'].append(np.array((era_day_sft['cape'].isel(lon=xx).values)))
        dic['tcwv_s'].append(np.array((era_day_sft['tcwv'].isel(lon=xx).values)))

        dic['q_s'].append((era_day_plt['q'].isel(lon=xxx).values))
        dic['u_s'].append((era_day_plt['u'].isel(lon=xxx).values))
        dic['r_s'].append((era_day_plt['r'].isel(lon=xxx).values))
        dic['v_s'].append((era_day_plt['v'].isel(lon=xxx).values))

        dic['tmin'].append(indic.tmin[id])
        dic['tmean'].append(indic.tmean[id])
        dic['t10'].append(indic.t10[id])
        dic['area'].append(indic.area[id])
        dic['area70'].append((indic['70area'])[id])

        print('DID', date)

    return dic

