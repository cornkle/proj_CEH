# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import ipdb
import os
import xarray as xr
import numpy as np
from utils import u_darrays as uda
from utils import constants as cnst
import pandas as pd


def saveCHIRPS():
    chirpsbox = [-81,-68,-18.5,0]  # peru daily

    chirpsall = xr.open_dataset(cnst.elements_drive + 'SouthAmerica/CHIRPS/chirps-v2.0.daily.peru.nc').chunk({'time':365})
    #chirpsm = xr.open_dataset(cnst.elements_drive + 'SouthAmerica/CHIRPS/chirps-v2.0.monthly.nc')
    chirpsall = chirpsall['precip']

    date = '2016-01-13'
    dt = pd.to_datetime(date)
    era5 = xr.open_dataset(cnst.ERA5_HOURLY_PL_HU+'/ERA5_'+str(dt.year)+'_'+str(dt.month).zfill(2)+'_'+str(dt.day).zfill(2)+'_pl.nc')

    u200 = era5['u'].sel(longitude=slice(chirpsbox[0], chirpsbox[1]), latitude=slice(chirpsbox[3], chirpsbox[2]), level=200, time=(era5['time.hour']==15)).squeeze().load()
    u200 = uda.flip_lat(u200)

    chirps = chirpsall.sel(time=date).load().squeeze()
    ch_on_e, lut = u200.salem.lookup_transform(chirps, return_lut=True)

    ch_on_e_all = u200.salem.lookup_transform(chirpsall, lut=lut)

    comp = dict(zlib=True, complevel=5)
    encoding = {'precip': comp}
    ch_on_e_all.to_netcdf('/media/ck/Elements/SouthAmerica/CHIRPS/CHIRPS_peru_onERA5.nc', mode='w', encoding=encoding, format='NETCDF4')


def saveGRIDSAT():
    chirpsbox = [-81,-68,-18.5,0]  # peru daily

    chirpsall_list = glob.glob(cnst.GRIDSAT_PERU + '/daily_-15ALLkm2_UTC_DAY/*.nc')

    for ids, ch in enumerate(chirpsall_list):

        fname = os.path.basename(ch)
        fname = fname.replace('UTCDay', 'UTCDay_onERA')
        outpath = cnst.GRIDSAT_PERU + '/daily_-15ALLkm2_UTC_DAY_onBIGERA/' + fname + '.nc'

        if os.path.isfile(outpath):
            print('File exists, continue')
            continue

        chirpsall = xr.open_dataset(ch)

        if ids ==0:
            date = '2016-01-13'
            dt = pd.to_datetime(date)
            era5 = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/uv_15UTC/uv_15UTC_1985_peru.nc')

            u200 = era5['u'].sel( level=250, time=(era5['time.hour']==15)).squeeze().load()
            u200 = uda.flip_lat(u200)

            chirps = chirpsall.isel(time=0).squeeze()
            #ipdb.set_trace()
            ch_on_e, lut = u200.salem.lookup_transform(chirps, return_lut=True)
        else:
            ch_on_e = u200.salem.lookup_transform(chirpsall, lut=lut)
        #ch_on_e.name = 'tir'
        comp = dict(zlib=True, complevel=5)
        encoding = {'tir': comp}
        fname = os.path.basename(ch)
        fname = fname.replace('UTCDay', 'UTCDay_onERA')
        #ipdb.set_trace()
        ch_on_e.to_netcdf(outpath, mode='w', encoding=encoding, format='NETCDF4')
        #ipdb.set_trace()

def saveGPM():
    chirpsbox = [-81,-68,-18.5,0]  # peru daily

    chirpsall_list = glob.glob('/media/ck/Elements/SouthAmerica/GPM/daily/*.nc4')
    for ids, ch in enumerate(chirpsall_list):

        chirpsall = xr.open_dataset(ch)
        #ipdb.set_trace()

        if ids ==0:
            date = '2016-01-13'
            dt = pd.to_datetime(date)
            era5 = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/uv_15UTC/uv_15UTC_1985_peru.nc')

            u200 = era5['u'].isel(time=0).sel(level=250).squeeze().load()
            u200 = uda.flip_lat(u200)

            chirps = chirpsall['HQprecipitation'].T.squeeze()
            #ipdb.set_trace()
            ch_on_e, lut = u200.salem.lookup_transform(chirps, return_lut=True)

        else:
            ch_on_e = u200.salem.lookup_transform(chirpsall['HQprecipitation'].T.squeeze(), lut=lut)
        ch_on_e.name = 'precip'
        comp = dict(zlib=True, complevel=5)
        encoding = {'precip': comp}
        fname = os.path.basename(ch)
        fname = fname.replace('.V06.nc4.SUB.nc4', '_onERA.nc')
        #ipdb.set_trace()
        ch_on_e.to_netcdf('/media/ck/Elements/SouthAmerica/GPM/daily_onERA/'+fname+'.nc', mode='w', encoding=encoding, format='NETCDF4')
        #ipdb.set_trace()


def saveERA5():
    #chirpsbox = [-81, -68, -18.5, 0]  # peru daily

    u200orig = xr.open_mfdataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/uv_15UTC/q*.nc')
    u200orig = u200orig['q'].sel(level=250).load()
    u200orig.name = 'q'
    comp = dict(zlib=True, complevel=5)
    encoding = {'q': comp}
    u200orig.to_netcdf('/media/ck/Elements/SouthAmerica/ERA5/hourly/q_15UTC_1981-2019_peru_big.nc', mode='w',
                  encoding=encoding, format='NETCDF4')




