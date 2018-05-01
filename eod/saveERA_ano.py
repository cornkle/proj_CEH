# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import pdb
import os
import xarray as xr
import numpy as np
import pandas as pd
from utils import constants, u_met
from scipy.interpolate import griddata


def saveAnomaly():
    mf = xr.open_dataset('/local/localscratch/cornkle/ERA-I/daily_2006-2010_12UTCpl.nc')

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

    dso.to_netcdf('/local/localscratch/cornkle/ERA-I/daily_2006-2010_12UTCpl_anomaly.nc')

def rewrite_ERA_latflip(file):

    with xr.open_dataset(file) as era:
        era = era.sel(latitude=slice(None, None, -1))
        era_write = era.load().copy()
    era_write.to_netcdf(file)
