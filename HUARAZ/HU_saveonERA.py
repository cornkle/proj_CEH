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
import time
import datetime


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


def saveCHIRPS_SA():
    chirpsbox = [-81,-65,-25,0]  # peru daily
    closer = [-90, -46, -35, 10]


    #chirpsm = xr.open_dataset(cnst.elements_drive + 'SouthAmerica/CHIRPS/chirps-v2.0.monthly.nc')
    path = cnst.elements_drive + 'SouthAmerica/CHIRPS/Global_Daily/'

    date = '2016-01-13'
    dt = pd.to_datetime(date)
    era5 = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/v850_15UTC_1981-2019_peru_big.nc')

    u200 = era5['v'].isel(time=0).squeeze().load()
    u200 = uda.flip_lat(u200)

    chirpsall = xr.open_dataset(
        glob.glob(cnst.elements_drive + 'SouthAmerica/CHIRPS/Global_Daily/chirps-*.nc')[0])
    chirpsall = chirpsall['precip'].sel(longitude=slice(closer[0], closer[1]), latitude=slice(closer[2], closer[3])).isel(time=0).squeeze()
    ch_on_e_all, lut = u200.salem.lookup_transform(chirpsall, return_lut=True)

    del chirpsall
    del ch_on_e_all

    for ids, y in enumerate(range(1981, 2020)):
        chirpsall = xr.open_dataset(glob.glob(cnst.elements_drive + 'SouthAmerica/CHIRPS/Global_Daily/chirps-*'+str(y)+'*.nc')[0]).chunk({'time':5})
        chirpsall = chirpsall['precip'].sel(longitude=slice(closer[0], closer[1]), latitude=slice(closer[2], closer[3]))
        #ipdb.set_trace()
        print('Doing year', y)

        ch_on_e_all = u200.salem.lookup_transform(chirpsall, lut=lut)

        comp = dict(zlib=True, complevel=5)
        encoding = {'precip': comp}
        ch_on_e_all.to_netcdf('/media/ck/Elements/SouthAmerica/CHIRPS/SA_daily_onERA/CHIRPS_daily_onERA_'+str(y)+'.nc', mode='w', encoding=encoding, format='NETCDF4')

        del ch_on_e_all


def saveGRIDSAT():
    chirpsbox = [-81,-65,-25,0]  # peru daily

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

    chirpsall_list = glob.glob('/media/ck/Elements/SouthAmerica/GPM/daily_precipCal/*.nc4')
    for ids, ch in enumerate(chirpsall_list):

        chirpsall = xr.open_dataset(ch)
        #ipdb.set_trace()

        if ids ==0:
            date = '2016-01-13'
            dt = pd.to_datetime(date)
            era5 = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/uv_15UTC/uv_15UTC_1985_peru.nc')

            u200 = era5['u'].isel(time=0).sel(level=250).squeeze().load()
            u200 = uda.flip_lat(u200)

            chirps = chirpsall['precipitationCal'].T.squeeze()
            #ipdb.set_trace()
            ch_on_e, lut = u200.salem.lookup_transform(chirps, return_lut=True)

        else:
            ch_on_e = u200.salem.lookup_transform(chirpsall['precipitationCal'].T.squeeze(), lut=lut)
        ch_on_e.name = 'precip'
        comp = dict(zlib=True, complevel=5)
        encoding = {'precip': comp}
        fname = os.path.basename(ch)
        fname = fname.replace('.V06.nc4.SUB.nc4', '_onERA.nc')
        #ipdb.set_trace()
        ch_on_e.to_netcdf('/media/ck/Elements/SouthAmerica/GPM/daily_precipCal_onERA/'+fname+'.nc', mode='w', encoding=encoding, format='NETCDF4')



def saveSM():
    chirpsbox = [-81,-65,-25,0]  # peru daily

    chirpsall_list = glob.glob('/media/ck/Elements/global/SM_CSI_ESA/daily/year_files_v6.1_combined_GLOBAL/*.nc')
    for ids, ch in enumerate(chirpsall_list):

        chirpsall = xr.open_dataset(ch)['sm']
        #chirpsall = uda.flip_lat(chirpsall)
        #chirpsall = chirpsall.sel(lon=slice(chirpsbox[0], chirpsbox[1]), lat=slice(chirpsbox[2], chirpsbox[3]))


        if ids ==0:
            date = '2016-01-13'
            dt = pd.to_datetime(date)
            era5 = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/qr_15UTC/uv_15UTC_1985_peru.nc')

            u200 = era5['u'].isel(time=0).sel(level=250).squeeze().load()
            u200 = uda.flip_lat(u200)

            ch_on_e = u200.salem.transform(chirpsall, interp='nearest')

        else:
            ch_on_e = u200.salem.transform(chirpsall, interp='nearest')

        ch_on_e.name = 'sm'
        comp = dict(zlib=True, complevel=5)
        encoding = {'sm': comp}
        fname = os.path.basename(ch)
        fname = fname.replace('.nc', '_onERA.nc')
        #ipdb.set_trace()
        ch_on_e.to_netcdf('/media/ck/Elements/SouthAmerica/SM_CSI_ESA/onERA/'+fname, mode='w', encoding=encoding, format='NETCDF4')
        print('Written ', ch)



def saveVOD():
    chirpsbox = [-81,-65,-25,0]  # peru daily

    chirpsall_list = glob.glob('/media/ck/Elements/global/VOD/*.nc')
    for ids, ch in enumerate(chirpsall_list):

        chirpsall = xr.open_dataset(ch)['vod']
        #chirpsall = uda.flip_lat(chirpsall)
        #chirpsall = chirpsall.sel(lon=slice(chirpsbox[0], chirpsbox[1]), lat=slice(chirpsbox[2], chirpsbox[3]))


        if ids ==0:
            date = '2016-01-13'
            dt = pd.to_datetime(date)
            era5 = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/uv_15UTC/uv_15UTC_1985_peru.nc')

            u200 = era5['u'].isel(time=0).sel(level=250).squeeze().load()
            u200 = uda.flip_lat(u200)

            ch_on_e = u200.salem.transform(chirpsall, interp='nearest')

        else:
            ch_on_e = u200.salem.transform(chirpsall, interp='nearest')

        ch_on_e.name = 'vod'
        comp = dict(zlib=True, complevel=5)
        encoding = {'vod': comp}
        fname = os.path.basename(ch)
        fname = fname.replace('_filtered', 'filtered_onERA_')
        #ipdb.set_trace()
        ch_on_e.to_netcdf('/media/ck/Elements/SouthAmerica/VOD/onERA/'+fname, mode='w', encoding=encoding, format='NETCDF4')


def saveVOD_native():
    chirpsbox = [-81,-65,-25,0]  # peru daily

    chirpsall_list = glob.glob('/media/ck/Elements/global/VOD/*.nc')
    for ids, ch in enumerate(chirpsall_list):

        chirpsall = xr.open_dataset(ch)['vod']
        #chirpsall = uda.flip_lat(chirpsall)

        ch_on_e = chirpsall.sel(lon=slice(chirpsbox[0], chirpsbox[1]), lat=slice(chirpsbox[2], chirpsbox[3]))
        ch_on_e.name = 'vod'
        comp = dict(zlib=True, complevel=5)
        encoding = {'vod': comp}
        fname = os.path.basename(ch)
        fname = fname.replace('_filtered', 'filtered_onERA_')
        #ipdb.set_trace()
        ch_on_e.to_netcdf('/media/ck/Elements/SouthAmerica/VOD/cutout/'+fname, mode='w', encoding=encoding, format='NETCDF4')

def saveNDVI():
    chirpsbox = [-81,-65,-25,0]  # peru daily

    chirpsall_list = glob.glob('/media/ck/Elements/global/VOD/*.nc')
    for ids, ch in enumerate(chirpsall_list):

        chirpsall = xr.open_dataarray(ch)
        chirpsall = uda.flip_lat(chirpsall)
        #chirpsall = chirpsall.sel(lon=slice(chirpsbox[0], chirpsbox[1]), lat=slice(chirpsbox[2], chirpsbox[3]))


        if ids ==0:
            date = '2016-01-13'
            dt = pd.to_datetime(date)
            era5 = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/uv_15UTC/uv_15UTC_1985_peru.nc')

            u200 = era5['u'].isel(time=0).sel(level=250).squeeze().load()
            u200 = uda.flip_lat(u200)

            ch_on_e, lut = u200.salem.lookup_transform(chirpsall, return_lut=True)

        else:
            ch_on_e = u200.salem.lookup_transform(chirpsall, lut=lut)

        ch_on_e.name = 'precip'
        comp = dict(zlib=True, complevel=5)
        encoding = {'precip': comp}
        fname = os.path.basename(ch)
        fname = fname.replace('MOD13C2', 'MOD13C2_onERA_')
        #ipdb.set_trace()
        ch_on_e.to_netcdf('/media/ck/Elements/SouthAmerica/NDVI/onERA/'+fname, mode='w', encoding=encoding, format='NETCDF4')

def saveNDVI_onCHIRPSbox():
    chirpsbox = [-81,-65,-25,0]  # peru daily


    chirpsall_list = glob.glob('/media/ck/Elements/global/NDVI_monthly/*.nc')
    for ids, ch in enumerate(chirpsall_list):

        chirpsall = xr.open_dataarray(ch)
        chirpsall = uda.flip_lat(chirpsall)
        chirpsall = chirpsall.sel(lon=slice(chirpsbox[0], chirpsbox[1]), lat=slice(chirpsbox[2], chirpsbox[3]))

        if ids == 0:
            date = '2016-01-13'
            dt = pd.to_datetime(date)
            era5 = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/uv_15UTC/uv_15UTC_1985_peru.nc')

            u200 = era5['u'].isel(time=0).sel(longitude=slice(chirpsbox[0], chirpsbox[1]),
                                 latitude=slice(chirpsbox[3], chirpsbox[2]), level=250).squeeze().load()
            u200 = uda.flip_lat(u200)

            chirps = chirpsall.squeeze()
            ch_on_e, lut = u200.salem.lookup_transform(chirpsall, return_lut=True)
        else:
            ch_on_e = u200.salem.lookup_transform(chirpsall, lut=lut)

        ch_on_e.name = 'precip'
        comp = dict(zlib=True, complevel=5)
        encoding = {'precip': comp}
        fname = os.path.basename(ch)
        fname = fname.replace('MOD13C2', 'MOD13C2_onCHIRPSbox_')
        #ipdb.set_trace()
        ch_on_e.to_netcdf('/media/ck/Elements/SouthAmerica/NDVI/onCHIRPSbox/'+fname+'.nc', mode='w', encoding=encoding, format='NETCDF4')



def saveERA5():
    #chirpsbox = [-81, -68, -18.5, 0]  # peru daily

    #u200orig = xr.open_mfdataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/uv_15UTC/q*.nc')
    var = 'tp'
    u200orig = xr.open_mfdataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/qr_15UTC/rain*peru2.nc', combine='nested', concat_dim='time')
    time.sleep(2)
    #u200orig = xr.open_mfdataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/peru/surface/ERA5*.nc', lock=False, combine='nested', concat_dim='time')
    time.sleep(2)
    print('Dataset loaded')

    u200orig = u200orig[var].isel(time=u200orig['time.hour']==15).load() #.sel(level=550).load()#.sel(level=850)#.load() .sel(level=650)
    u200orig.name = var
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp}
    u200orig.to_netcdf('/media/ck/Elements/SouthAmerica/ERA5/hourly/rain_15UTC_1981-2019_peru_big.nc', mode='w',
                  encoding=encoding, format='NETCDF4')


def saveERA5_multi():

    def get_ERA5(file):
        var = 'tcwv'
        u200orig = xr.open_dataset(file)
        try:
            u200orig = u200orig[var].isel(time=u200orig['time.hour']==15,expver=0).load()#.sel(level=200).isel(time=u200orig['time.hour']==15).load()
        except:
            u200orig = u200orig[var].isel(time=u200orig['time.hour'] == 15).load()
        u200orig.name = var

        return u200orig

       # ipdb.set_trace()

    var = 'tcwv'
    #inputs = glob.glob('/media/ck/Elements/SouthAmerica/ERA5/hourly/peru/pressure_levels/ERA5*.nc')
    inputs = glob.glob('/media/ck/Elements/SouthAmerica/ERA5/hourly/peru/surface/ERA5*.nc')
   #  pool = multiprocessing.Pool(processes=5)
   # # ipdb.set_trace()
   #  res = pool.map(get_ERA5, inputs[0:5])
   #  pool.close()

    res = []
    for file in inputs:
        print('Doing ', file)
        out = get_ERA5(file)
        res.append(out)

    #ipdb.set_trace()

    ds = xr.concat(res, dim='time')
    try:
        ds = ds.drop('expver')
    except:
        pass
    ds = ds.load()

    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp}
    ds.to_netcdf('/media/ck/Elements/SouthAmerica/ERA5/hourly/tcwv_15UTC_1981-2019_peru_big.nc', mode='w',
                       encoding=encoding, format='NETCDF4')


def saveERA5_multifix():

    def get_ERA5(file):
        var = 'v'
        u200orig = xr.open_dataset(file)
        try:
            u200orig = u200orig[var].sel(level=400).isel(time=u200orig['time.hour']==15,expver=0).load()#.sel(level=200).isel(time=u200orig['time.hour']==15).load()
        except:
            u200orig = u200orig[var].sel(level=400).isel(time=u200orig['time.hour'] == 15).load()
        u200orig.name = var

        return u200orig

       # ipdb.set_trace()

    var = 'v'
    #inputs = glob.glob('/media/ck/Elements/SouthAmerica/ERA5/hourly/peru/pressure_levels/ERA5*.nc')
    inputs = np.sort(glob.glob('/media/ck/Elements/SouthAmerica/ERA5/hourly/peru/pressure_levels/ERA5*.nc'))

   #  pool = multiprocessing.Pool(processes=5)
   # # ipdb.set_trace()
   #  res = pool.map(get_ERA5, inputs[0:5])
   #  pool.close()

    res = []
    dummy = get_ERA5(inputs[0])
    dummy_new = dummy.copy(deep=True)
    #ipdb.set_trace()
    del dummy
    for file in inputs:
        print('Doing ', file)
        out = get_ERA5(file)

        duout = dummy_new.copy(deep=True)
        duout.values = out.values

        # y =  int(file[-18:-14])
        # m = int(file[-13:-11])
        # d = int(file[-10:-8])
       # ipdb.set_trace()
        #date = datetime.datetime(y,m,d,15)

        duout['time'] = ('time', out.time) #[date]
        print(out.time)
        #
        res.append(duout)
    #ipdb.set_trace()
    ds = xr.concat(res, dim='time')

    ds = ds.load()

    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp}
    ds.to_netcdf('/media/ck/Elements/SouthAmerica/ERA5/hourly/v400_15UTC_1981-2019_peru.nc', mode='w',
                       encoding=encoding, format='NETCDF4')


def saveERA5_multifix_srfc():

    def get_ERA5(file):
        var = 'cape'
        u200orig = xr.open_dataset(file)
        try:
            u200orig = u200orig[var].isel(time=u200orig['time.hour']==15,expver=0).load()#.sel(level=200).isel(time=u200orig['time.hour']==15).load()
        except:
            u200orig = u200orig[var].isel(time=u200orig['time.hour'] == 15).load()
        u200orig.name = var

        return u200orig

       # ipdb.set_trace()

    var = 'cape'
    #inputs = glob.glob('/media/ck/Elements/SouthAmerica/ERA5/hourly/peru/pressure_levels/ERA5*.nc')
    inputs = np.sort(glob.glob('/media/ck/Elements/SouthAmerica/ERA5/hourly/peru/surface/ERA5*.nc'))

   #  pool = multiprocessing.Pool(processes=5)
   # # ipdb.set_trace()
   #  res = pool.map(get_ERA5, inputs[0:5])
   #  pool.close()

    res = []
    dummy = get_ERA5(inputs[0])
    dummy_new = dummy.copy(deep=True)
    #ipdb.set_trace()
    del dummy
    for file in inputs:
        print('Doing ', file)
        out = get_ERA5(file)

        duout = dummy_new.copy(deep=True)
        duout.values = out.values

        # y =  int(file[-18:-14])
        # m = int(file[-13:-11])
        # d = int(file[-10:-8])
       # ipdb.set_trace()
        #date = datetime.datetime(y,m,d,15)

        duout['time'] = ('time', out.time) #[date]
        print(out.time)
        #
        res.append(duout)
    #ipdb.set_trace()
    ds = xr.concat(res, dim='time')

    ds = ds.load()

    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp}
    ds.to_netcdf('/media/ck/Elements/SouthAmerica/ERA5/hourly/cape_15UTC_1981-2019_peru_big.nc', mode='w',
                       encoding=encoding, format='NETCDF4')





