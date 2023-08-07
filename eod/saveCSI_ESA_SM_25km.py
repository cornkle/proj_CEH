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



from utils import u_arrays

def saveAnomaly():

    mf = xr.open_mfdataset(cnst.elements_drive + 'global/SM_CSI_ESA/daily/year_files_v6.1_combined_GLOBAL/20*.nc', concat_dim='time')
    #mf = mf.sel(lon=slice(-11,11), lat=slice(9,21))
    mf = mf['sm'][(mf['time.month'] >= 6) & (mf['time.month'] <= 9)]#.load()
    ipdb.set_trace()
    mf.values[mf.values<1] = np.nan

    mf['ymonth'] = ('time', [str(y)+'-'+str(m) for (y,m) in zip(mf['time.year'].values,mf['time.month'].values)])
    # minus =  mf.groupby('ymonth').mean(dim='time')
    # dso = mf.groupby('ymonth') - minus
    # dso = dso.drop('ymonth')
    #
    # dso.to_netcdf('/users/global/cornkle/data/OBS/AMSRE/day_aqua/amsre_monthly_anomaly.nc')

    grouped='ymonth'
    ipdb.set_trace()
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
        ds.to_netcdf(path=cnst.elements_drive + 'global/SM_CSI_ESA/daily/clim/volSM_25km_clim_'+date+'.nc', mode='w', encoding=encoding, format='NETCDF4')




def saveAnomalyDay():

    #allfiles = glob.glob(cnst.elements_drive + 'global/AMSR2/daily/10km/day/nc/AMSR2*.nc')
    allfiles = glob.glob(cnst.elements_drive + 'global/SM_CSI_ESA/daily/year_files_v6.1_combined_GLOBAL/20*.nc')



    for m in range(1,13):

        mf = xr.open_mfdataset(allfiles, concat_dim='time') #, combine='by_coords'

        # mf['monthDay'] = (
        # 'time', [ str(m).zfill(2) + '-' + str(d).zfill(2) for (m,d) in zip(mf['time.month'].values,  mf['time.day'].values)])


        clim = mf['sm'].groupby('time.dayofyear').sum(dim='time')
        valid_days = mf['sm'].groupby('time.dayofyear').count(dim='time')  # number of valid days per month

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
            ds.to_netcdf(path=cnst.elements_drive + 'global/SM_CSI_ESA/daily/clim/volSM_25km_clim_'+outdate+'.nc', mode='w', encoding=encoding, format='NETCDF4')

            del ds, da
        del climda, countda


def writeAnomaly():

    tag = 'day'

    files = glob.glob(cnst.elements_drive + 'global/AMSR2/daily/25km/'+tag+'/nc/*.nc')

    for f in files:
        print('Doing', f)

        basename = os.path.basename(f)

        day = xr.open_dataset(f)
        climpath = cnst.elements_drive + 'global/AMSR2/daily/25km/'+tag+'_clim/'

        try:
            clim = xr.open_dataset(climpath + 'amsr2_25km_clim_2008'+basename[-7:-3]+'.nc')
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
        out.to_netcdf(path=cnst.elements_drive + 'global/AMSR2/daily/25km/day_anom/amsr2_25km_anom_'+basename[-11:-3]+'.nc', mode='w', encoding=encoding, format='NETCDF4')

