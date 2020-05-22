# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import ipdb
import os
import xarray as xr
import numpy as np
from utils import u_arrays as ua
from utils import constants as cnst
import pandas as pd

def saveYearly():

    out = cnst.GRIDSAT_PERU #cnst.local_data + 'GRIDSAT/MCS18_peru/'
    infolder = '/media/ck/Elements/SouthAmerica/GRIDSAT/3hourly/www.ncei.noaa.gov/data/geostationary-ir-channel-brightness-temperature-gridsat-b1/access/'

    years = np.arange(1984,2018)  # list(next(os.walk(msg_folder))[1])

    for y in years:
        filename = 'gridsat_WA_-40_1000km2' + str(y) + '.nc'
        da = None
        if os.path.isfile(out + filename):
            continue

        files = glob.glob(infolder + str(y) + '/GRIDSAT-AFRICA_CP*.nc')
        files.sort()
        for f in files:
            print('Doing ' + f)

            df = xr.open_dataset(f)
            if (df['time.hour']<15) | (df['time.hour']>21):
                continue

            df.rename({'irwin_cdr':'tir'}, inplace=True)
            df['tir'].values = df['tir'].values-273.15
            labels, goodinds = ua.blob_define(df['tir'].values, -40, minmax_area=[16, 25000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
            df['tir'].values[labels == 0] = 0
            df['tir'].values[df['tir'].values < -110] = 0
            df['tir'].values = (np.round(df['tir'].values, decimals=2)*100).astype(np.int16)
            try:
                da = xr.concat([da, df ], dim='time')
            except TypeError:
                da = df.copy()

        enc = {'tir': {'complevel': 5, 'shuffle': True, 'zlib': True}}
        da.to_netcdf(out + filename, encoding=enc)
        da.close()



def saveYearly_parallel():


    years = np.arange(2018, 2020)  # list(next(os.walk(msg_folder))[1])

    pool = multiprocessing.Pool(processes=4)

    res = pool.map(loop, years)
    pool.close()


def loop(y):

    out = cnst.GRIDSAT_PERU
    infolder = cnst.elements_drive + 'SouthAmerica/GRIDSAT/3hourly/www.ncei.noaa.gov/data/geostationary-ir-channel-brightness-temperature-gridsat-b1/access/'
    filename = 'gridsat_WA_-40_5000km2_13-19UTC' + str(y) + '.nc'
    da = None
    #ipdb.set_trace()
    if os.path.isfile(out + filename):
        return

    files = glob.glob(infolder + str(y) + '/GRIDSAT-SouthAmerica_CP*.nc')
    #ipdb.set_trace()
    if files == []:
        return
    files.sort()
    for f in files:
        print('Doing ' + f)

        df = xr.open_dataset(f)

        if df['time.hour'] not in [18,21,0]: # LT UTC-5 [13,16,19]
            continue

        df = df.rename({'irwin_cdr': 'tir'}) #, inplace=True
        df['tir'].values = df['tir'].values-273.15
        labels, goodinds = ua.blob_define(df['tir'].values, -40, minmax_area=[83, 25000],
                                          max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
        df['tir'].values[labels == 0] = 0
        df['tir'].values[df['tir'].values < -110] = 0
        df['tir'].values = (np.round(df['tir'].values, decimals=2)*100).astype(np.int16)
        try:
            da = xr.concat([da, df], dim='time')
        except TypeError:
            da = df.copy()

    enc = {'tir': {'complevel': 5, 'shuffle': True, 'zlib': True}}
    da.to_netcdf(out + filename, encoding=enc)
    da.close()

def rewrite(file):

    ds = xr.open_dataset(file).load()
    # ds = ds.where(ds['time.day'] == 1, drop=True)
    #ds = ds.rename({'lat': 'latitude', 'lon': 'longitude'})
    out = file.replace('.nc', '_compressed.nc')
    ipdb.set_trace()

    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(out, mode='w', format='NETCDF4')

def rewrite_day():

    for y in np.arange(1994,1995):#(1985,2017):

        path = cnst.GRIDSAT_PERU +'daily_LT/gridsat_WA_-40Min_5000km2_13-19UTCperDay_'
        try:
            dataset = xr.open_mfdataset(cnst.GRIDSAT_PERU+'*'+str(y)+'*.nc')
        except:
            continue

        ds = None

        for date in dataset['time'].sel(time=((dataset['time.hour']==18)&(dataset['time.month']==4))):


            date_today = pd.to_datetime(date.values)
            print('Doing', date_today)
            date_today = date_today.replace(hour=18)

            date_end = date_today + pd.Timedelta('6 hours')

            dat = dataset['tir'].sel(time=slice(date_today,date_end)).min('time')
            data = np.reshape(dat.values, (1,dat.values.shape[0], dat.values.shape[1]))
            #ipdb.set_trace()
            da = xr.DataArray(data, coords=[[date_today],  dat['lat'], dat['lon']], dims=['time', 'lat', 'lon'])

            da.name = 'tir'

            try:
                ds = xr.concat([ds, da], dim='time')
            except TypeError:
                ds = da.copy()

            #if date_end.month > date_today.month:
        out = path + str(date_today.year)+'-'+str(date_today.month).zfill(2)+'.nc'
        comp = dict(zlib=True, complevel=5)
        encoding = {'tir': comp}
        ds.to_netcdf(out, mode='w', format='NETCDF4', encoding=encoding)

        ds = None


def rewrite_UTCday():

    for y in np.arange(2018,2020): #(1984,2018)

        print('Doing', y)
        inpath = cnst.elements_drive + 'SouthAmerica/GRIDSAT/3hourly/www.ncei.noaa.gov/data/geostationary-ir-channel-brightness-temperature-gridsat-b1/access/'+str(y)+'/'
        #file = glob.glob(cnst.GRIDSAT_PERU + '*'+str(y)+'*.nc')

        path = cnst.GRIDSAT_PERU +'daily_ALLkm2_UTC_DAY/gridsat_WA_-40Min_ALLkm2_UTCDay_'
        try:
            dataset = xr.open_mfdataset(inpath + '*.nc', combine='nested', concat_dim='time')
        except:
            print('FAIL')
            continue

        #ipdb.set_trace()

        datout = dataset.resample(time='1D').min()
        ipdb.set_trace()
        for m in np.unique(datout['time.month']):

            monthly = datout['irwin_cdr'].sel(time=datout['time.month']==m).load()-273.15
            monthly.values[(monthly.values>-40) | (monthly.values<-109)] = 0

            monthly.values = (np.round(monthly.values, decimals=2) * 100).astype(np.int16)
            monthly.name = 'tir'

            comp = dict(zlib=True, complevel=5)
            #ipdb.set_trace()
            encoding = {'tir': comp}
            out = path + str(y) + '-' + str(m).zfill(2) + '.nc'
            monthly.to_netcdf(out, mode='w', format='NETCDF4', encoding=encoding)
        del datout
        del dataset


def rewrite_UTCday_afternoon():

    for y in np.arange(1994,1995):#(1984,2018):

        print('Doing', y)
        #inpath = cnst.elements_drive + 'SouthAmerica/GRIDSAT/3hourly/www.ncei.noaa.gov/data/geostationary-ir-channel-brightness-temperature-gridsat-b1/access/'+str(y)+'/'
        inpath = cnst.GRIDSAT_PERU + '*'+str(y)+'*.nc'

        path = cnst.GRIDSAT_PERU +'daily_5000km2_UTC_afternoon/gridsat_WA_-40Min_5000km2_UTCDay_'
        try:
            dataset = xr.open_mfdataset(inpath).chunk({'time':365})
        except:
            continue

        #ipdb.set_trace()

        datout = dataset.resample(time='1D').min()

        for m in np.unique(datout['time.month']):

            monthly = datout.sel(time=datout['time.month']==m).load()

            comp = dict(zlib=True, complevel=5)
            #ipdb.set_trace()
            encoding = {'tir': comp}
            out = path + str(y) + '-' + str(m).zfill(2) + '.nc'
            monthly.to_netcdf(out, mode='w', format='NETCDF4', encoding=encoding)

        del datout
        del dataset


def rewrite_UTCday_fullDay():  ######### FULL DAY MCS NEEDED

    for y in np.arange(2018,2020):

        print('Doing', y)
        #inpath = cnst.elements_drive + 'SouthAmerica/GRIDSAT/3hourly/www.ncei.noaa.gov/data/geostationary-ir-channel-brightness-temperature-gridsat-b1/access/'+str(y)+'/'
        inpath = cnst.GRIDSAT_PERU + '*'+str(y)+'*.nc'

        path = cnst.GRIDSAT_PERU +'daily_5000km2_UTC_DAY/gridsat_WA_-40Min_5000km2_UTCDay_'
        try:
            dataset = xr.open_mfdataset(inpath)#.chunk({'time':365})
        except:
            ipdb.set_trace()
            continue

        #ipdb.set_trace()

        datout = dataset.resample(time='1D').min()

        for m in np.unique(datout['time.month']):

            monthly = datout['irwin_cdr'].sel(time=datout['time.month']==m).load()

            comp = dict(zlib=True, complevel=5)
            #ipdb.set_trace()
            encoding = {'tir': comp}
            out = path + str(y) + '-' + str(m).zfill(2) + '.nc'
            monthly.to_netcdf(out, mode='w', format='NETCDF4', encoding=encoding)
            print('Saving', out)

        del datout
        del dataset


