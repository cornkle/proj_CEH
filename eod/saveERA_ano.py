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


def saveClimatology_pl():
    flist = []
    for y in range(2006,2014):
        fs = glob.glob(cnst.local_data + '/ERA5/hourly/pressure_levels/ERA5_*'+str(y)+'*_pl.nc')
        flist.extend(fs)

    mf = xr.open_mfdataset(flist, concat_dim='time')
    #mf['ymonth'] = ('time', [str(y) + '-' + str(m) for (y, m) in zip(mf['time.year'].values, mf['time.month'].values)])

    mf['monthHour'] = (
    'time', [str(y).zfill(2) + '-' + str(m).zfill(2) for (y, m) in zip(mf['time.month'].values, mf['time.hour'].values)])
    clim = mf.groupby('monthHour').mean(dim='time')
    for _, sl in clim.groupby('monthHour'):

        ds = xr.Dataset()

        date = [pd.datetime(2010, np.array(str(sl.monthHour.values)[0:2], dtype=int), 15, np.array(str(sl.monthHour.values)[3:5], dtype=int), 0)]
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
        ds.to_netcdf(cnst.local_data + '/ERA5/monthly/synop_selfmade/CLIM_2006-2010/ERA5_2006-2010_CLIM_'+str(ds['time.month'].values[0]).zfill(2)+'-'+str(ds['time.hour'].values[0]).zfill(2)+'_pl.nc')#, mode='w', encoding=encoding, format='NETCDF4')


def saveClimatology_srfc():
    flist = []
    for y in range(2006,2011):
        fs = glob.glob(cnst.local_data + '/ERA5/hourly/surface/ERA5_*'+str(y)+'*_srfc.nc')
        flist.extend(fs)

    mf = xr.open_mfdataset(flist, concat_dim='time')
    #mf['ymonth'] = ('time', [str(y) + '-' + str(m) for (y, m) in zip(mf['time.year'].values, mf['time.month'].values)])

    mf['monthHour'] = (
    'time', [str(y).zfill(2) + '-' + str(m).zfill(2) for (y, m) in zip(mf['time.month'].values, mf['time.hour'].values)])
    clim = mf.groupby('monthHour').mean(dim='time')
    for _, sl in clim.groupby('monthHour'):

        ds = xr.Dataset()

        date = [pd.datetime(2010, np.array(str(sl.monthHour.values)[0:2], dtype=int), 15, np.array(str(sl.monthHour.values)[3:5], dtype=int), 0)]
        for svar in sl.data_vars:

            da = xr.DataArray(sl[svar].values[None, ...],
                              coords={'time': date, 'lat': sl.latitude.values, 'lon': sl.longitude.values},
                              dims=['time', 'lat', 'lon'])  # [np.newaxis, :]
            da.attrs = sl[svar].attrs
            ds[svar] = da
            ds.attrs = clim.attrs
        # comp = dict(zlib=True, complevel=5)
        # encoding = {var: comp for var in ds.data_vars}
        #ipdb.set_trace()
        ds.to_netcdf(cnst.local_data + '/ERA5/monthly/synop_selfmade/CLIM_2006-2010/ERA5_2006-2010_CLIM_'+str(ds['time.month'].values[0]).zfill(2)+'-'+str(ds['time.hour'].values[0]).zfill(2)+'_srfc.nc')#, mode='w', encoding=encoding, format='NETCDF4')



def saveClimatology_pl_day():
    flist = []
    for y in range(2006,2011):
        fs = glob.glob(cnst.local_data + '/ERA5/hourly/pressure_levels/ERA5_'+str(y)+'*_pl.nc')
        flist.extend(fs)

    mf = xr.open_mfdataset(flist, concat_dim='time', combine='by_coords')  #
    mf['monthHourDay'] = (
    'time', [str(y).zfill(2) + '-' + str(m).zfill(2) + '-' + str(d).zfill(2) for (y, m,d) in zip(mf['time.month'].values, mf['time.hour'].values, mf['time.day'].values)])

    clim = mf.groupby('monthHourDay').mean(dim='time')
    times=[]
    for tt in clim['monthHourDay']:
        hh = str(tt.values)

        mmonth = int(hh[0:2])
        dday = int(hh[6:8])
        hhour = int(hh[3:5])
        date = pd.datetime(2008, mmonth, dday, hhour, 0)
        times.append(date)
    tseries = pd.to_datetime(times)



    for hh in np.unique(tseries.hour):


            pos = np.where(tseries.hour==hh)[0]

            hdata = xr.Dataset()

            for svar in clim.data_vars:
                da = xr.DataArray(clim[svar].values[pos,:,:],
                                  coords={'time': tseries[pos], 'level': clim.level.values, 'lat': clim.latitude.values,
                                          'lon': clim.longitude.values},
                                  dims=['time', 'level', 'lat', 'lon'])

                hdata[svar] = da


            for tstep in hdata.time:

                dt = pd.to_datetime([tstep.values])

                if (dt.month < 6) | (dt.month > 9):
                    continue


                window1 = dt - pd.Timedelta('5days')
                window2 = dt + pd.Timedelta('5days')

                sliced = hdata.sel(time=slice(window1[0],window2[0]))

                out = sliced.mean('time')

                ds = xr.Dataset()
                #ipdb.set_trace()
                date = [pd.datetime(2008, dt.month[0], dt.day[0], dt.hour[0], 0)]

                for svar in out.data_vars:

                    da = xr.DataArray(out[svar].values[None, ...],
                                      coords={'time': date, 'level' : clim.level.values, 'lat': out.lat.values, 'lon': out.lon.values},
                                      dims=['time', 'level', 'lat', 'lon'])  # [np.newaxis, :]
                    da.attrs = out[svar].attrs
                    ds[svar] = da
                ds.attrs = clim.attrs
                # comp = dict(zlib=True, complevel=5)
                # encoding = {var: comp for var in ds.data_vars}


                ds.to_netcdf(cnst.local_data + '/ERA5/monthly/synop_selfmade/CLIM_2006-2010/ERA5_2006-2010_CLIM_'+str(ds['time.month'].values[0]).zfill(2)+'-'+str(ds['time.day'].values[0]).zfill(2)+'-'+str(ds['time.hour'].values[0]).zfill(2)+'_pl.nc')#, mode='w', encoding=encoding, format='NETCDF4')
                print('Wrote '+'ERA5_2006-2010_CLIM_'+str(ds['time.month'].values[0]).zfill(2)+'-'+str(ds['time.day'].values[0]).zfill(2)+'-'+str(ds['time.hour'].values[0]).zfill(2)+'_pl.nc')
                del sliced


def saveStddev_pl_day():
    flist = []
    for y in range(2006,2011):
        fs = glob.glob(cnst.local_data + '/ERA5/hourly/pressure_levels/ERA5_'+str(y)+'*_pl.nc')
        flist.extend(fs)

    mf = xr.open_mfdataset(flist, concat_dim='time', combine='by_coords')  #
    mf['monthHourDay'] = (
    'time', [str(y).zfill(2) + '-' + str(m).zfill(2) + '-' + str(d).zfill(2) for (y, m,d) in zip(mf['time.month'].values, mf['time.hour'].values, mf['time.day'].values)])

    clim = mf.groupby('monthHourDay').std(dim='time')
    times=[]
    for tt in clim['monthHourDay']:
        hh = str(tt.values)

        mmonth = int(hh[0:2])
        dday = int(hh[6:8])
        hhour = int(hh[3:5])
        date = pd.datetime(2008, mmonth, dday, hhour, 0)
        times.append(date)
    tseries = pd.to_datetime(times)



    for hh in np.unique(tseries.hour):


            pos = np.where(tseries.hour==hh)[0]

            hdata = xr.Dataset()

            for svar in clim.data_vars:
                da = xr.DataArray(clim[svar].values[pos,:,:],
                                  coords={'time': tseries[pos], 'level': clim.level.values, 'lat': clim.latitude.values,
                                          'lon': clim.longitude.values},
                                  dims=['time', 'level', 'lat', 'lon'])

                hdata[svar] = da


            for tstep in hdata.time:

                dt = pd.to_datetime([tstep.values])

                if (dt.month < 6) | (dt.month > 9):
                    continue


                window1 = dt - pd.Timedelta('5days')
                window2 = dt + pd.Timedelta('5days')

                sliced = hdata.sel(time=slice(window1[0],window2[0]))

                out = sliced.mean('time')

                ds = xr.Dataset()
                #ipdb.set_trace()
                date = [pd.datetime(2008, dt.month[0], dt.day[0], dt.hour[0], 0)]

                for svar in out.data_vars:

                    da = xr.DataArray(out[svar].values[None, ...],
                                      coords={'time': date, 'level' : clim.level.values, 'lat': out.lat.values, 'lon': out.lon.values},
                                      dims=['time', 'level', 'lat', 'lon'])  # [np.newaxis, :]
                    da.attrs = out[svar].attrs
                    ds[svar] = da
                ds.attrs = clim.attrs
                # comp = dict(zlib=True, complevel=5)
                # encoding = {var: comp for var in ds.data_vars}


                ds.to_netcdf(cnst.local_data + '/ERA5/monthly/synop_selfmade/CLIM_2006-2010/ERA5_2006-2010_STD_'+str(ds['time.month'].values[0]).zfill(2)+'-'+str(ds['time.day'].values[0]).zfill(2)+'-'+str(ds['time.hour'].values[0]).zfill(2)+'_pl.nc')#, mode='w', encoding=encoding, format='NETCDF4')
                print('Wrote '+'ERA5_2006-2010_STD_'+str(ds['time.month'].values[0]).zfill(2)+'-'+str(ds['time.day'].values[0]).zfill(2)+'-'+str(ds['time.hour'].values[0]).zfill(2)+'_pl.nc')
                del sliced



def saveClimatology_srfc_day():
    flist = []
    for y in range(2006,2011):
        fs = glob.glob(cnst.local_data + '/ERA5/hourly/surface/ERA5_'+str(y)+'*_srfc.nc')
        flist.extend(fs)

    mf = xr.open_mfdataset(flist, concat_dim='time', combine='by_coords')  #
    mf['monthHourDay'] = (
    'time', [str(y).zfill(2) + '-' + str(m).zfill(2) + '-' + str(d).zfill(2) for (y, m,d) in zip(mf['time.month'].values, mf['time.hour'].values, mf['time.day'].values)])

    clim = mf.groupby('monthHourDay').mean(dim='time')
    times=[]
    for tt in clim['monthHourDay']:
        hh = str(tt.values)

        mmonth = int(hh[0:2])
        dday = int(hh[6:8])
        hhour = int(hh[3:5])
        date = pd.datetime(2008, mmonth, dday, hhour, 0)
        times.append(date)
    tseries = pd.to_datetime(times)



    for hh in np.unique(tseries.hour):


            pos = np.where(tseries.hour==hh)[0]

            hdata = xr.Dataset()

            for svar in clim.data_vars:
                da = xr.DataArray(clim[svar].values[pos,:,:],
                                  coords={'time': tseries[pos],'lat': clim.latitude.values,
                                          'lon': clim.longitude.values},
                                  dims=['time', 'lat', 'lon'])

                hdata[svar] = da


            for tstep in hdata.time:

                dt = pd.to_datetime([tstep.values])

                if (dt.month < 6) | (dt.month > 9):
                    continue


                window1 = dt - pd.Timedelta('5days')
                window2 = dt + pd.Timedelta('5days')

                sliced = hdata.sel(time=slice(window1[0],window2[0]))

                out = sliced.mean('time')

                ds = xr.Dataset()
                #ipdb.set_trace()
                date = [pd.datetime(2008, dt.month[0], dt.day[0], dt.hour[0], 0)]

                for svar in out.data_vars:

                    da = xr.DataArray(out[svar].values[None, ...],
                                      coords={'time': date, 'lat': out.lat.values, 'lon': out.lon.values},
                                      dims=['time','lat', 'lon'])  # [np.newaxis, :]
                    da.attrs = out[svar].attrs
                    ds[svar] = da
                ds.attrs = clim.attrs
                # comp = dict(zlib=True, complevel=5)
                # encoding = {var: comp for var in ds.data_vars}


                ds.to_netcdf(cnst.local_data + '/ERA5/monthly/synop_selfmade/CLIM_2006-2010/ERA5_2006-2010_CLIM_'+str(ds['time.month'].values[0]).zfill(2)+'-'+str(ds['time.day'].values[0]).zfill(2)+'-'+str(ds['time.hour'].values[0]).zfill(2)+'_srfc.nc')#, mode='w', encoding=encoding, format='NETCDF4')
                print('Wrote '+'ERA5_2006-2010_CLIM_'+str(ds['time.month'].values[0]).zfill(2)+'-'+str(ds['time.day'].values[0]).zfill(2)+'-'+str(ds['time.hour'].values[0]).zfill(2)+'_srfc.nc')
                del sliced