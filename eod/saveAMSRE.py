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


def saveAnomalyDay():

    mf = xr.open_mfdataset(cnst.network_data + 'data/OBS/AMSRE/aqua/nc_day/AMSR_*.nc', concat_dim='time', combine='by_coords')
    #mf = mf.sel(lon=slice(-11,11), lat=slice(9,21))
    mf = mf['SM'][(mf['time.month'] >= 3) & (mf['time.month'] <= 11)]

    mf['monthDay'] = (
    'time', [ str(m).zfill(2) + '-' + str(d).zfill(2) for (m,d) in zip(mf['time.month'].values,  mf['time.day'].values)])

    clim = mf.groupby('monthDay').sum(dim='time')
    valid_days = mf.groupby('monthDay').count(dim='time')  # number of valid days per month

    times=[]
    for tt in clim['monthDay']:
        hh = str(tt.values)

        mmonth = int(hh[0:2])
        dday = int(hh[3:5])
        date = pd.datetime(2008, mmonth, dday, 0, 0)
        times.append(date)
    tseries = pd.to_datetime(times)


    climda = xr.DataArray(clim.values,
                      coords={'time': tseries, 'lat': clim.lat.values,
                              'lon': clim.lon.values},
                      dims=['time', 'lat', 'lon'])

    countda = xr.DataArray(valid_days.values,
                      coords={'time': tseries, 'lat': clim.lat.values,
                              'lon': clim.lon.values},
                      dims=['time', 'lat', 'lon'])


    for tstep in climda.time:

        dt = pd.to_datetime([tstep.values])


        window1 = dt - pd.Timedelta('5days')
        window2 = dt + pd.Timedelta('5days')

        climsliced = climda.sel(time=slice(window1[0],window2[0]))
        countsliced = countda.sel(time=slice(window1[0], window2[0]))

        outclim = climsliced.sum('time')
        outcount = countsliced.sum('time')

        allclim = outclim / outcount
        allclim.values[outcount.values < 10] = np.nan

        date = [pd.datetime(2008, dt.month[0], dt.day[0], 0, 0)]

        da = xr.DataArray(allclim.values[None, ...],
                              coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
                              dims=['time', 'lat', 'lon'])  # [np.newaxis, :]
        da.attrs = climda.attrs

        ds = xr.Dataset({'SM': da})
        odate = date[0]

        outdate = str(odate.year)+str(odate.month).zfill(2)+str(odate.day).zfill(2)

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(path=cnst.network_data + 'data/OBS/AMSRE/aqua/sma_clim_day/sma_'+outdate+'.nc', mode='w', encoding=encoding, format='NETCDF4')



def writeAnomaly():

    tag = 'night'

    files = glob.glob(cnst.network_data + 'data/OBS/AMSRE/aqua/nc_'+tag+'/AMSR_*.nc')

    for f in files:

        basename = os.path.basename(f)

        day = xr.open_dataset(f)
        climpath = cnst.network_data + 'data/OBS/AMSRE/aqua/sma_clim_'+tag+'/'

        try:
            clim = xr.open_dataset(climpath + 'sma_2008'+basename[-7:-3]+'.nc')
        except FileNotFoundError:
            continue

        out = day.copy()

        out['SM'].values = day['SM'].values - clim['SM'].values

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in out.data_vars}
        out.to_netcdf(path=cnst.network_data + 'data/OBS/AMSRE/aqua/sma_nc_'+tag+'_new/'+ 'sma_'+basename[-11:-3]+'.nc', mode='w', encoding=encoding, format='NETCDF4')

