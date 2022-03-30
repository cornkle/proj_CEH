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
#from scipy.interpolate import griddata
import ipdb

def saveNetcdf(day=True):

    if day:
        daystring = '_A_'
        dstring = 'day'
    else:
        daystring = '_D_'
        dstring = 'night'


    sm_folder = cnst.lmcs_drive + '/AMSR2/daily/25km/'+dstring+'/raw'
    pool = multiprocessing.Pool(processes=3)
    files = glob.glob(sm_folder+'/LPRM*.nc4')

    print('start loop')
    #ipdb.set_trace()
    # ipdb.set_trace()
    # for f in files:
    #     ds = rewrite_data.rewrite_AMSR2_10km(f)

    res = pool.map(rewrite_data.rewrite_AMSR2_10km, files)





from utils import u_arrays
def saveAnomalyDay():

    dtag = 'night'

    #allfiles = glob.glob(cnst.elements_drive + 'global/AMSR2/daily/10km/day/nc/AMSR2*.nc')
    allfiles = u_arrays.locate('.nc', cnst.lmcs_drive + '/AMSR2/daily/25km/'+dtag+'/nc/')
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
            gpos = np.arange(pp-15,pp+16)
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


        clim = mf['soil_moisture_c1'].where(mf['soil_moisture_c1']>=10).groupby('time.dayofyear').sum(dim='time')
        #clim2 = mf['soil_moisture_c2'].groupby('time.dayofyear').sum(dim='time')
        ts = mf['ts'].where(mf['soil_moisture_c1']>=10).groupby('time.dayofyear').sum(dim='time')

        valid_days = mf['soil_moisture_c1'].where(mf['soil_moisture_c1']>=10).groupby('time.dayofyear').count(dim='time')  # number of valid days per month

        times=[]
        #ipdb.set_trace()
        for tt in clim['dayofyear']:
            hh = tt.values-1

            # mmonth = int(hh[0:2])
            # dday = int(hh[3:5])
            #ipdb.set_trace()
            if (m == 1) & (hh >=350):
                    date = datetime.datetime(2007, 1, 1, 0, 0) + pd.Timedelta(str(hh)+' days')
            elif (m == 12) & (hh <=17):
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

        # climda2 = xr.DataArray(clim2.values,
        #                   coords={'time': tseries, 'lat': clim.lat.values,
        #                           'lon': clim.lon.values},
        #                   dims=['time', 'lat', 'lon'])

        # climda2 = climda2.sortby(climda2.time)

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


        for tstep in climda.time[15:-15]:

            print('Doing', tstep)

            dt = pd.to_datetime([tstep.values])

            try:
                window1 = dt - pd.Timedelta('15days')
            except:
                window1= dt
            try:
                window2 = dt + pd.Timedelta('15days')
            except:
                window2 = dt

            #ipdb.set_trace()
            climsliced = climda.sel(time=slice(window1[0],window2[0]))
            # climsliced2 = climda2.sel(time=slice(window1[0], window2[0]))
            tsliced = temp.sel(time=slice(window1[0],window2[0]))
            countsliced = countda.sel(time=slice(window1[0], window2[0]))

            outclim = climsliced.sum('time')
            # outclim2 = climsliced2.sum('time')
            outtemp = tsliced.sum('time')
            outcount = countsliced.sum('time')

            allclim = outclim / outcount
            alltclim = outtemp / outcount
            # allclim2 = outclim2/outcount
            allclim.values[outcount.values < 10] = np.nan
            alltclim.values[outcount.values < 10] = np.nan
            # allclim2.values[outcount.values < 10] = np.nan


            date = [pd.datetime(2008, dt.month[0], dt.day[0], 0, 0)]

            da = xr.DataArray(allclim.values[None, ...],
                                  coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
                                  dims=['time', 'lat', 'lon'])  # [np.newaxis, :]

            # da2 = xr.DataArray(allclim2.values[None, ...],
            #                       coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
            #                       dims=['time', 'lat', 'lon'])

            dtt = xr.DataArray(alltclim.values[None, ...],
                                  coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
                                  dims=['time', 'lat', 'lon'])

            da.attrs = climda.attrs
            dtt.attrs = temp.attrs
            # da2.attrs = climda2.attrs

            ds = xr.Dataset({'SM': da})
            ds['LST'] = dtt
            #ds['SM2'] = da2
            odate = date[0]

            outdate = str(odate.year)+str(odate.month).zfill(2)+str(odate.day).zfill(2)

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in ds.data_vars}
            ds.to_netcdf(path=cnst.lmcs_drive + 'AMSR2/daily/25km/'+dtag+'_clim_v2/amsr2_25km_clim_'+outdate+'.nc', mode='w', encoding=encoding, format='NETCDF4')
            print('Written')
            del ds, da,  dtt #da2,
        del climda,  temp, countda #climda2,




def saveStddevDay():
    tag = 'night'
    #allfiles = glob.glob(cnst.elements_drive + 'global/AMSR2/daily/10km/day/nc/AMSR2*.nc')
    allfiles = u_arrays.locate('.nc', cnst.lmcs_drive + '/AMSR2/daily/25km/'+tag+'/nc/')
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


        clim = mf['soil_moisture_c1'].groupby('time.dayofyear').mean(dim='time')
        #clim2 = mf['soil_moisture_c2'].groupby('time.dayofyear').mean(dim='time')
        ts = mf['ts'].groupby('time.dayofyear').mean(dim='time')

        climstd = mf['soil_moisture_c1'].groupby('time.dayofyear').std(dim='time')
        #clim2std = mf['soil_moisture_c2'].groupby('time.dayofyear').std(dim='time')
        tsstd = mf['ts'].groupby('time.dayofyear').std(dim='time')

        times=[]
        #ipdb.set_trace()
        for tt in clim['dayofyear']:
            hh = tt.values-1

            # mmonth = int(hh[0:2])
            # dday = int(hh[3:5])
            #ipdb.set_trace()
            if (m == 1) & (hh >=350):
                    date = datetime.datetime(2007, 1, 1, 0, 0) + pd.Timedelta(str(hh)+' days')
            elif (m == 12) & (hh <=17):
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

        # climda2 = xr.DataArray(clim2.values,
        #                   coords={'time': tseries, 'lat': clim.lat.values,
        #                           'lon': clim.lon.values},
        #                   dims=['time', 'lat', 'lon'])

        # climda2 = climda2.sortby(climda2.time)

        temp = xr.DataArray(ts.values,
                          coords={'time': tseries, 'lat': clim.lat.values,
                                  'lon': clim.lon.values},
                          dims=['time', 'lat', 'lon'])
        temp = temp.sortby(temp.time)

 #############################
        climda_std = xr.DataArray(climstd.values,
                          coords={'time': tseries, 'lat': clim.lat.values,
                                  'lon': clim.lon.values},
                          dims=['time', 'lat', 'lon'])

        climda_std = climda_std.sortby(climda_std.time)

        # climda2_std = xr.DataArray(clim2std.values,
        #                   coords={'time': tseries, 'lat': clim.lat.values,
        #                           'lon': clim.lon.values},
        #                   dims=['time', 'lat', 'lon'])

        # climda2_std = climda2_std.sortby(climda2_std.time)

        temp_std = xr.DataArray(tsstd.values,
                          coords={'time': tseries, 'lat': clim.lat.values,
                                  'lon': clim.lon.values},
                          dims=['time', 'lat', 'lon'])
        temp_std = temp_std.sortby(temp_std.time)




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
            # climsliced2 = climda2.sel(time=slice(window1[0], window2[0]))
            tsliced = temp.sel(time=slice(window1[0],window2[0]))

            outclim = climsliced.std('time')
            # outclim2 = climsliced2.std('time')
            outtemp = tsliced.std('time')

            allclim = (outclim + climda_std) /2
            alltclim = (outtemp + temp_std) / 2
            # allclim2 = (outclim2 + climda2_std)/ 2


            date = [pd.datetime(2008, dt.month[0], dt.day[0], 0, 0)]

            da = xr.DataArray(allclim.values[None, ...],
                                  coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
                                  dims=['time', 'lat', 'lon'])  # [np.newaxis, :]

            # da2 = xr.DataArray(allclim2.values[None, ...],
            #                       coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
            #                       dims=['time', 'lat', 'lon'])

            dtt = xr.DataArray(alltclim.values[None, ...],
                                  coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
                                  dims=['time', 'lat', 'lon'])

            da.attrs = climda.attrs
            dtt.attrs = temp.attrs
            # da2.attrs = climda2.attrs

            ds = xr.Dataset({'SM': da})
            ds['LST'] = dtt
            # ds['SM2'] = da2
            odate = date[0]

            outdate = str(odate.year)+str(odate.month).zfill(2)+str(odate.day).zfill(2)

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in ds.data_vars}
            ds.to_netcdf(path=cnst.lmcs_drive + 'AMSR2/daily/25km/'+tag+'_clim/amsr2_25km_stddev_'+outdate+'.nc', mode='w', encoding=encoding, format='NETCDF4')

            del ds, da, dtt #da2,
        del climda, temp #climda2,


def writeAnomaly():

    tag = 'night'

    files = glob.glob(cnst.lmcs_drive + 'AMSR2/daily/25km/'+tag+'/nc/*.nc')

    for f in files:
        print('Doing', f)

        basename = os.path.basename(f)

        day = xr.open_dataset(f)
        climpath = cnst.lmcs_drive + 'AMSR2/daily/25km/'+tag+'_clim_v2/'

        try:
            clim = xr.open_dataset(climpath + 'amsr2_25km_clim_2008'+basename[-7:-3]+'.nc')
        except FileNotFoundError:
            continue

        out = day.copy()

        var_orig = ['ts', 'soil_moisture_c1'] # 'soil_moisture_c2'
        var_new = ['LST', 'SM'] #'SM2'
        for vo, vn in zip(var_orig,var_new):

            out[vo].values = day[vo].where(day['soil_moisture_c1']>=10).values - clim[vn].values

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in out.data_vars}
        #outf = basename.replace(".nc", "_anom.nc")
        out.to_netcdf(path=cnst.lmcs_drive + 'AMSR2/daily/25km/'+tag+'_anom_v2/amsr2_25km_anom_'+basename[-11:-3]+'.nc', mode='w', encoding=encoding, format='NETCDF4')

