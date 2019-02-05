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

def saveYearly():

    out = cnst.local_data + 'GRIDSAT/MCS18/'
    infolder = cnst.local_data + 'GRIDSAT/www.ncei.noaa.gov/data/geostationary-ir-channel-brightness-temperature-gridsat-b1/access/'

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


    years = np.arange(1983, 2018)  # list(next(os.walk(msg_folder))[1])

    pool = multiprocessing.Pool(processes=4)

    res = pool.map(loop, years)
    pool.close()


def loop(y):

    out = cnst.local_data + 'GRIDSAT/MCS18/'
    infolder = cnst.local_data + 'GRIDSAT/www.ncei.noaa.gov/data/geostationary-ir-channel-brightness-temperature-gridsat-b1/access/'
    filename = 'gridsat_WA_-40_1000km2_15-21UTC' + str(y) + '.nc'
    da = None
    if os.path.isfile(out + filename):
        return

    files = glob.glob(infolder + str(y) + '/GRIDSAT-AFRICA_CP*.nc')
    files.sort()
    for f in files:
        print('Doing ' + f)

        df = xr.open_dataset(f)

        if (df['time.hour']<15) | (df['time.hour']>21):
            continue

        df.rename({'irwin_cdr': 'tir'}, inplace=True)
        df['tir'].values = df['tir'].values-273.15
        labels, goodinds = ua.blob_define(df['tir'].values, -40, minmax_area=[16, 25000],
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
