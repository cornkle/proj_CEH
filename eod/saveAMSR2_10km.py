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
import datetime
from scipy.interpolate import griddata
import ipdb

def saveNetcdf(day=True):

    if day:
        daystring = '_A_'
        dstring = 'day'
    else:
        daystring = '_D_'
        dstring = 'night'


    sm_folder = cnst.elements_drive + 'global/AMSR2/daily/10km/'+dstring+'/raw'
    pool = multiprocessing.Pool(processes=3)
    files = glob.glob(sm_folder+'/LPRM*.nc4')
    print('start loop')
    #ipdb.set_trace()
    # ipdb.set_trace()
    # for f in files:
    #     ds = rewrite_data.rewrite_AMSRE(f, day=True)

    res = pool.map(rewrite_data.rewrite_AMSR2_10km, files)

def saveMonthly():

    bla = xr.open_mfdataset(cnst.network_data + 'data/OBS/AMSRE/aqua/nc/*.nc')
    monthly = bla.resample('m', dim='time', how='mean')
    monthly.to_netcdf(cnst.network_data + 'data/OBS/AMSRE/aqua/amsre_day_monthly.nc')

def saveAnomaly():
    mf = xr.open_mfdataset(cnst.elements_drive + 'global/AMSR2/daily/10km/day/LPRM*.nc4', concat_dim='time')
    #mf = mf.sel(lon=slice(-11,11), lat=slice(9,21))
    mf = mf['SM'][(mf['time.month'] >= 6) & (mf['time.month'] <= 9)].load()

    mf.values[mf.values<1] = np.nan

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

        date = [datetime.datetime(year, month, day, 0, 0)]
        da = xr.DataArray(arr.values[None, ...],
                          coords={'time': date, 'lat': arr.lat, 'lon': arr.lon},
                          dims=['time', 'lat', 'lon'])  # [np.newaxis, :]
        ds = xr.Dataset({'SM': da})

        date = str(arr['time.year'].values)+str(arr['time.month'].values).zfill(2)+str(arr['time.day'].values).zfill(2)

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(path=cnst.elements_drive + 'global/AMSR2/daily/10km/day/sma_10km_'+date+'.nc', mode='w', encoding=encoding, format='NETCDF4')

from utils import u_arrays
def saveAnomalyDay():

    #allfiles = glob.glob(cnst.elements_drive + 'global/AMSR2/daily/10km/day/nc/AMSR2*.nc')
    allfiles = u_arrays.locate('.nc', cnst.elements_drive + 'global/AMSR2/daily/10km/day/nc/')
    #ipdb.set_trace()
    mstruc = {}
    for m in range(1,13):
        mstruc[m]=[]

    only_months=list()
    for af in allfiles:
        ismonth = af[-7:-5]
        only_months.append(int(ismonth))

    for m in range(1,13):

        goodpos = []
        pos = np.where(np.array(only_months) == m)
        for pp in pos[0]:
            #ipdb.set_trace()
            gpos = np.arange(pp-5,pp+6)
            gpos[gpos<0]=0
            gpos[gpos>len(only_months)-1]=len(only_months)-1
            goodpos.extend(gpos)
        goodpos = np.unique(goodpos)
       # ipdb.set_trace()
        goodfiles = list(np.array(allfiles)[goodpos])

        mstruc[m]=goodfiles


    for m in range(1,13):

        mf = xr.open_mfdataset(mstruc[m], concat_dim='time') #, combine='by_coords'

        # mf['monthDay'] = (
        # 'time', [ str(m).zfill(2) + '-' + str(d).zfill(2) for (m,d) in zip(mf['time.month'].values,  mf['time.day'].values)])


        clim = mf['soil_moisture_c1'].groupby('time.dayofyear').sum(dim='time')
        clim2 = mf['soil_moisture_c2'].groupby('time.dayofyear').sum(dim='time')
        ts = mf['ts'].groupby('time.dayofyear').sum(dim='time')

        valid_days = mf['soil_moisture_c1'].groupby('time.dayofyear').count(dim='time')  # number of valid days per month

        times=[]
        #ipdb.set_trace()
        for tt in clim['dayofyear']:
            hh = tt.values-1

            # mmonth = int(hh[0:2])
            # dday = int(hh[3:5])
            #ipdb.set_trace()
            if (m == 1) & (hh >=360):
                    date = datetime.datetime(2007, 1, 1, 0, 0) + pd.Timedelta(str(hh)+' days')
            elif (m == 12) & (hh <=7):
                    date = datetime.datetime(2009, 1, 1, 0, 0) + pd.Timedelta(str(hh) + ' days')
            else:
                date = datetime.datetime(2008, 1, 1, 0, 0) + pd.Timedelta(str(hh) + ' days')

            times.append(date)
        tseries = pd.to_datetime(times)


        climda = xr.DataArray(clim.values,
                          coords={'time': tseries, 'lat': clim.lat.values,
                                  'lon': clim.lon.values},
                          dims=['time', 'lat', 'lon'])

        climda = climda.sortby(climda.time)

        climda2 = xr.DataArray(clim2.values,
                          coords={'time': tseries, 'lat': clim.lat.values,
                                  'lon': clim.lon.values},
                          dims=['time', 'lat', 'lon'])

        climda2 = climda2.sortby(climda2.time)

        temp = xr.DataArray(ts.values,
                          coords={'time': tseries, 'lat': clim.lat.values,
                                  'lon': clim.lon.values},
                          dims=['time', 'lat', 'lon'])
        temp = temp.sortby(temp.time)

        countda = xr.DataArray(valid_days.values,
                          coords={'time': tseries, 'lat': clim.lat.values,
                                  'lon': clim.lon.values},
                          dims=['time', 'lat', 'lon'])
        countda = countda.sortby(countda.time)


        for tstep in climda.time[5:-5]:

            print('Doing', tstep)

            dt = pd.to_datetime([tstep.values])

            try:
                window1 = dt - pd.Timedelta('5days')
            except:
                window1= dt
            try:
                window2 = dt + pd.Timedelta('5days')
            except:
                window2 = dt

            #ipdb.set_trace()
            climsliced = climda.sel(time=slice(window1[0],window2[0]))
            climsliced2 = climda2.sel(time=slice(window1[0], window2[0]))
            tsliced = temp.sel(time=slice(window1[0],window2[0]))
            countsliced = countda.sel(time=slice(window1[0], window2[0]))

            outclim = climsliced.sum('time')
            outclim2 = climsliced2.sum('time')
            outtemp = tsliced.sum('time')
            outcount = countsliced.sum('time')

            allclim = outclim / outcount
            alltclim = outtemp / outcount
            allclim2 = outclim2/outcount
            allclim.values[outcount.values < 10] = np.nan
            alltclim.values[outcount.values < 10] = np.nan
            allclim2.values[outcount.values < 10] = np.nan


            date = [pd.datetime(2008, dt.month[0], dt.day[0], 0, 0)]

            da = xr.DataArray(allclim.values[None, ...],
                                  coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
                                  dims=['time', 'lat', 'lon'])  # [np.newaxis, :]

            da2 = xr.DataArray(allclim.values[None, ...],
                                  coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
                                  dims=['time', 'lat', 'lon'])

            dtt = xr.DataArray(alltclim.values[None, ...],
                                  coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
                                  dims=['time', 'lat', 'lon'])

            da.attrs = climda.attrs
            dtt.attrs = temp.attrs
            da2.attrs = climda2.attrs

            ds = xr.Dataset({'SM': da})
            ds['LST'] = dtt
            ds['SM2'] = da2
            odate = date[0]

            outdate = str(odate.year)+str(odate.month).zfill(2)+str(odate.day).zfill(2)

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in ds.data_vars}
            ds.to_netcdf(path=cnst.elements_drive + 'global/AMSR2/daily/10km/day_clim/amsr2_10km_clim_'+outdate+'.nc', mode='w', encoding=encoding, format='NETCDF4')

            del ds, da, da2, dtt
        del climda, climda2, temp, countda

def writeAnomaly():

    tag = 'day'

    files = glob.glob(cnst.elements_drive + 'global/AMSR2/daily/10km/'+tag+'/nc/*.nc')

    for f in files:

        basename = os.path.basename(f)

        day = xr.open_dataset(f)
        climpath = cnst.elements_drive + 'global/AMSR2/daily/10km/'+tag+'_clim/'

        try:
            clim = xr.open_dataset(climpath + 'amsr2_10km_clim_2008'+basename[-7:-3]+'.nc')
        except FileNotFoundError:
            continue

        out = day.copy()

        var_orig = ['ts', 'soil_moisture_c1', 'soil_moisture_c2']
        var_new = ['LST', 'SM', 'SM2']
        for vo, vn in zip(var_orig,var_new):

            out[vo].values = day[vo].values - clim[vn].values

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in out.data_vars}
        #outf = basename.replace(".nc", "_anom.nc")
        out.to_netcdf(path=cnst.elements_drive + 'global/AMSR2/daily/10km/day_anom/amsr2_10km_anom_'+basename[-11:-3]+'.nc', mode='w', encoding=encoding, format='NETCDF4')

