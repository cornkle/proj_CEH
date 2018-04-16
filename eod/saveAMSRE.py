# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import pdb
import os
import xarray as xr
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

def saveNetcdf():

    sm_folder = '/users/global/cornkle/data/OBS/AMSRE/day_aqua/raw_day'
    pool = multiprocessing.Pool(processes=7)
    files = glob.glob(sm_folder+'/AMSR*.gra')

    # for f in files:
    #     ds = rewrite_data.rewrite_AMSRE(f, day=True)

    res = pool.map(rewrite_data.rewrite_AMSRE, files)

def saveMonthly():

    bla = xr.open_mfdataset('/users/global/cornkle/data/OBS/AMSRE/day_aqua/nc/*.nc')
    monthly = bla.resample('m', dim='time', how='mean')
    monthly.to_netcdf('/users/global/cornkle/data/OBS/AMSRE/day_aqua/amsre_day_monthly.nc')

def saveAnomaly():
    mf = xr.open_mfdataset('/users/global/cornkle/data/OBS/AMSRE/day_aqua/nc_day/AMSR*.nc', concat_dim='time')
    mf = mf.sel(lon=slice(-11,11), lat=slice(9,21))

    mf = mf['SM'][(mf['time.month'] >= 6) & (mf['time.month'] <= 9)]

    mf['ymonth'] = ('time', [str(y)+'-'+str(m) for (y,m) in zip(mf['time.year'].values,mf['time.month'].values)])
    # minus =  mf.groupby('ymonth').mean(dim='time')
    # dso = mf.groupby('ymonth') - minus
    # dso = dso.drop('ymonth')
    #
    # dso.to_netcdf('/users/global/cornkle/data/OBS/AMSRE/day_aqua/amsre_monthly_anomaly.nc')

    grouped='ymonth'

    valid_days = mf.groupby(grouped).count(dim='time') # number of valid days per month

    minus =  mf.groupby(grouped).mean(dim='time')
    arr = minus.values

    arr[valid_days.values<10] = np.nan
    minus.values = arr

    dso = mf.groupby(grouped) - minus

    for d in dso['time']:
        try:
            arr = dso.sel(time=d.values).drop('ymonth')
        except ValueError:
            arr = dso.sel(time=d.values)
        day = arr['time.day'].values
        month = arr['time.month'].values
        year = arr['time.year'].values

        date = [pd.datetime(year, month, day, 0, 0)]
        da = xr.DataArray(arr.values[None, ...],
                          coords={'time': date, 'lat': arr.lat, 'lon': arr.lon},
                          dims=['time', 'lat', 'lon'])  # [np.newaxis, :]
        ds = xr.Dataset({'SM': da})

        date = str(arr['time.year'].values)+str(arr['time.month'].values).zfill(2)+str(arr['time.day'].values).zfill(2)
        ds.to_netcdf('/users/global/cornkle/data/OBS/AMSRE/day_aqua/sma_nc_day/smaD_'+date+'.nc')
