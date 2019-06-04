# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import pdb
import os
import xarray as xr
from utils import constants as cnst
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import ipdb

def saveNetcdf(day=True):

    if day:
        daystring = '_A_'
        dstring = 'day'
    else:
        daystring = '_D_'
        dstring = 'night'


    sm_folder = cnst.network_data + 'data/OBS/AMSRE/aqua/raw_'+dstring
    pool = multiprocessing.Pool(processes=5)
    files = glob.glob(sm_folder+'/AMSR*.gra')
    print('start loop')
    # ipdb.set_trace()
    # for f in files:
    #     ds = rewrite_data.rewrite_AMSRE(f, day=True)

    res = pool.map(rewrite_data.rewrite_AMSRE, files)

def saveMonthly():

    bla = xr.open_mfdataset(cnst.network_data + 'data/OBS/AMSRE/aqua/nc/*.nc')
    monthly = bla.resample('m', dim='time', how='mean')
    monthly.to_netcdf(cnst.network_data + 'data/OBS/AMSRE/aqua/amsre_day_monthly.nc')

def saveAnomaly():
    mf = xr.open_mfdataset(cnst.network_data + 'data/OBS/AMSRE/aqua/nc_night/AMSR*.nc', concat_dim='time')
    #mf = mf.sel(lon=slice(-11,11), lat=slice(9,21))
    mf = mf['SM'][(mf['time.month'] >= 3) & (mf['time.month'] <= 11)]

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
        print('Doing ', d.values)
        day = arr['time.day'].values
        month = arr['time.month'].values
        year = arr['time.year'].values

        date = [pd.datetime(year, month, day, 0, 0)]
        da = xr.DataArray(arr.values[None, ...],
                          coords={'time': date, 'lat': arr.lat, 'lon': arr.lon},
                          dims=['time', 'lat', 'lon'])  # [np.newaxis, :]
        ds = xr.Dataset({'SM': da})

        date = str(arr['time.year'].values)+str(arr['time.month'].values).zfill(2)+str(arr['time.day'].values).zfill(2)

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(path=cnst.network_data + 'data/OBS/AMSRE/aqua/sma_nc_night/sma_'+date+'.nc', mode='w', encoding=encoding, format='NETCDF4')
