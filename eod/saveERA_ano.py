# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import ipdb
import os
import xarray as xr
import numpy as np
import pandas as pd
from utils import constants as cnst, u_met
from scipy.interpolate import griddata


def saveAnomaly():
    mf = xr.open_dataset('/localscratch/wllf030/cornkle/ERA5/ERA5_2010_12UTCpl.nc' )

    u = mf['u'].values
    v = mf['v'].values

    ws, wd = u_met.u_v_to_ws_wd(u,v)
    mf['ws'] = (('time', 'level', 'latitude', 'longitude'), ws)
    mf['wd'] = (('time', 'level', 'latitude', 'longitude'), wd)

    mf['ymonth'] = ('time', [str(y)+'-'+str(m) for (y,m) in zip(mf['time.year'].values,mf['time.month'].values)])
    # minus =  mf.groupby('ymonth').mean(dim='time')
    # dso = mf.groupby('ymonth') - minus
    # dso = dso.drop('ymonth')
    #
    # dso.to_netcdf('/users/global/cornkle/data/OBS/AMSRE/day_aqua/amsre_monthly_anomaly.nc')

    grouped='ymonth'

    valid_days = mf.groupby(grouped).count(dim='time') # number of valid days per month

    minus =  mf.groupby(grouped).mean(dim='time')
    # arr = minus.values
    #
    # arr[valid_days.values<10] = np.nan
    # minus.values = arr

    dso = mf.groupby(grouped) - minus
    dso = dso.drop('ymonth')

    dso.to_netcdf('/localscratch/wllf030/cornkle/ERA5/ERA5_2010_12UTCpl_anomaly.nc')

def rewrite_ERA_latflip(file):

    with xr.open_dataset(file) as era:
        era = era.sel(latitude=slice(None, None, -1))
        era_write = era.load().copy()
    era_write.to_netcdf(file)


def saveClimatology():
    flist = []
    for y in range(2000,2015):
        fs = glob.glob(cnst.local_data + '/ERA5/pressure_levels/ERA5_*'+str(y)+'*_pl.nc')
        flist.extend(fs)

    mf = xr.open_mfdataset(flist, concat_dim='time')
    #mf['ymonth'] = ('time', [str(y) + '-' + str(m) for (y, m) in zip(mf['time.year'].values, mf['time.month'].values)])

    mf['monthHour'] = (
    'time', [str(y).zfill(2) + '-' + str(m).zfill(2) for (y, m) in zip(mf['time.month'].values, mf['time.hour'].values)])
    clim = mf.groupby('monthHour').mean(dim='time')
    for _, sl in clim.groupby('monthHour'):

        ds = xr.Dataset()

        date = [pd.datetime(2014, np.array(str(sl.monthHour.values)[0:2], dtype=int), 15, np.array(str(sl.monthHour.values)[3:5], dtype=int), 0)]
        for svar in sl.data_vars:

            da = xr.DataArray(sl[svar].values[None, ...],
                              coords={'time': date, 'level' : clim.level.values, 'lat': sl.latitude.values, 'lon': sl.longitude.values},
                              dims=['time', 'level', 'lat', 'lon'])  # [np.newaxis, :]
            da.attrs = sl[svar].attrs
            ds[svar] = da
            ds.attrs = clim.attrs
        # comp = dict(zlib=True, complevel=5)
        # encoding = {var: comp for var in ds.data_vars}
        #ipdb.set_trace()
        ds.to_netcdf(cnst.local_data + '/ERA5/CLIM_test/ERA5_2000-2014_CLIM_'+str(ds['time.month'].values[0]).zfill(2)+'-'+str(ds['time.hour'].values[0]).zfill(2)+'_pl.nc')#, mode='w', encoding=encoding, format='NETCDF4')
