# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import multiprocessing
import pdb
import pandas as pd
from scipy import ndimage
from CLOVER import era_geop_t3d as era_geop
from utils import u_gis
import pickle as pkl
import statsmodels.stats.proportion as prop

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def composite(h):
    pool = multiprocessing.Pool(processes=8)


    file = '/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == 17) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2009) & (msg['time.month'] >= 6) ]

    msg = msg.sel(lat=slice(10.5,17.5), lon=slice(-9.5,9.5))

    res = pool.map(file_loop, msg)
    pool.close()

    # res = []
    # for m in msg[0:50]:
    #     r = file_loop(m)
    #     res.append(r)



    res = [x for x in res if x is not None]

    scales = res

    scales = [item for sublist in scales for item in sublist]  # flatten list of lists

    scales = np.concatenate(scales)


    return scales



def cut_kernel(xpos, ypos, array):


    # kernel = lsta_day2.isel(lon=slice(xpos+3 , xpos + 4 + 1),
    #                         lat=slice(ypos-1 , ypos + 1 + 1)).values
    # #
    # #
    # kernel =  np.mean(kernel)

    kernel = array.isel(lon=xpos+5,
                             lat=ypos).values

    # if (kernel == np.nan):
    #     return

    return kernel


def file_loop(fi):
    print('Doing day: ', fi)

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))
    dayd = pd.Timedelta('1 days')

    if (fi['time.hour'].values) <= 15:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date

    print('Opening ', '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/power_maps/' \
                           'lsta_daily_powerMaxscale_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2) + '.nc' )
    try:
        lsta = xr.open_dataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/power_maps/' \
                           'lsta_daily_powerMaxscale_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2) + '.nc')
    except OSError:
        print('Cant find file')
        return

    print('Doing '+ 'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    temp_da =  xr.open_dataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/scale_maps/' \
                           'lsta_daily_scale_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2) + '.nc')


    lsta = lsta['LSTA']
    temp_da = temp_da['LSTA']
    lsta = lsta.where(temp_da > -800)
    temp_da = temp_da.where(temp_da > -800)

    lsta = lsta.where(temp_da <= 20)


    # lsta = lsta.where(np.abs(temp) > 0.2)
    # temp_da = temp_da.where(np.abs(temp) > 0.2)

    if np.sum(lsta) == np.nan:
        return

    pos = np.where( (fi.values >= 5) & (fi.values <= 30) )

    if np.sum(pos) == 0:
        print('No blobs found')
        return

    scale_list = []

    for y, x in zip(pos[0], pos[1]):

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        t = fi.sel(lat=lat, lon=lon)

        point = lsta.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            ano = cut_kernel(xpos, ypos, lsta)
        except TypeError:
            continue

        if ano == np.nan:
            continue


        scale_list.append(ano)


    if scale_list == None:
        return None

    # if np.isnan(scale_list).all():
    #     return None

    print('Returning')

    return scale_list

if __name__ == "__main__":
    composite()

