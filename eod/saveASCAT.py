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
import datetime

def saveNetcdf():

    path = cnst.lmcs_drive + 'ASCAT/raw/'
    outpath = cnst.lmcs_drive + 'ASCAT/nc/'
    for y in range(2007,2013):
        if not os.path.isdir(outpath + str(y)):
            os.mkdir(outpath+str(y))

        for ff in glob.glob(path+str(y)+'/*.gra'):
                rewrite_data.rewrite_ASCAT(ff)


def saveAnomalyDay():
    dtag = 'am'

    allfiles = sorted(glob.glob(cnst.lmcs_drive + 'ASCAT/nc/*/*_' + dtag + '.nc'))

    #ipdb.set_trace()
    mstruc = {}
    for m in range(1, 13):
        mstruc[m] = []

    only_months = list()
    for af in allfiles:
        ismonth = af[-10:-8]
        #ipdb.set_trace()
        only_months.append(int(ismonth))

    for m in range(1, 13):

        goodpos = []
        pos = np.where(np.array(only_months) == m)
        for pp in pos[0]:
            gpos = np.arange(pp - 15, pp + 16)
            gpos[gpos < 0] = 0
            gpos[gpos > len(only_months) - 1] = len(only_months) - 1
            goodpos.extend(gpos)
        goodpos = np.unique(goodpos)

        goodfiles = list(np.array(allfiles)[goodpos])

        mstruc[m] = goodfiles

    for m in range(1, 13):

        mf = xr.open_mfdataset(mstruc[m], concat_dim='time')  # , combine='by_coords'

        clim = mf['SM'].groupby('time.dayofyear').sum(dim='time')
        valid_days = mf['SM'].groupby('time.dayofyear').count(dim='time')  # number of valid days per month

        times = []
        # ipdb.set_trace()
        for tt in clim['dayofyear']:
            hh = tt.values - 1

            if (m == 1) & (hh >= 350):
                date = datetime.datetime(2007, 1, 1, 0, 0) + pd.Timedelta(str(hh) + ' days')
            elif (m == 12) & (hh <= 17):
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

        for tstep in climda.time[15:-15]:

            print('Doing', tstep)

            dt = pd.to_datetime([tstep.values])
            date = [datetime.datetime(2008, dt.month[0], dt.day[0], 0, 0)]
            odate = date[0]
            outdate = str(odate.year) + str(odate.month).zfill(2) + str(odate.day).zfill(2)
            outpath = cnst.lmcs_drive + 'ASCAT/clim_' + dtag
            if os.path.isfile(outpath + '/ascat_' + outdate + '_' + dtag + '.nc'):
                print('File exists, continue')
                continue

            try:
                window1 = dt - pd.Timedelta('15days')
            except:
                window1 = dt
            try:
                window2 = dt + pd.Timedelta('15days')
            except:
                window2 = dt

            climsliced = climda.sel(time=slice(window1[0], window2[0]))
            countsliced = countda.sel(time=slice(window1[0], window2[0]))

            outclim = climsliced.sum('time')
            outcount = countsliced.sum('time')

            allclim = outclim / outcount
            allclim.values[outcount.values < 10] = np.nan



            da = xr.DataArray(allclim.values[None, ...],
                              coords={'time': date, 'lat': allclim.lat.values, 'lon': allclim.lon.values},
                              dims=['time', 'lat', 'lon'])

            da.attrs = climda.attrs

            ds = xr.Dataset({'SM': da})

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in ds.data_vars}

            if not os.path.isdir(outpath):
                os.mkdir(outpath)

            ds.to_netcdf(path=outpath + '/ascat_' + outdate + '_' + dtag + '.nc',
                         mode='w', encoding=encoding, format='NETCDF4')
            print('Written')
            del ds, da
        del climda, countda


def writeAnomaly():

    tag = 'am'

    files = glob.glob(cnst.lmcs_drive + 'ASCAT/nc/*/*_'+tag+'.nc')

    for f in files:
        print('Doing', f)

        basename = os.path.basename(f)
        outpath = cnst.lmcs_drive + 'ASCAT/anom_'+tag+'/'+basename

        if os.path.isfile(outpath):
            print('File exists, continue')
            continue
        #ipdb.set_trace()
        day = xr.open_dataset(f)
        climpath = cnst.lmcs_drive + 'ASCAT/clim_'+tag+'/'

        try:
            clim = xr.open_dataset(climpath + 'ascat_2008'+basename[-10:-6]+'_'+tag+'.nc')
        except FileNotFoundError:
            continue

        out = day.copy()

        out['SM'].values = day['SM'].values - clim['SM'].values

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in out.data_vars}

        out.to_netcdf(path=outpath, mode='w', encoding=encoding, format='NETCDF4')

