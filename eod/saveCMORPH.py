# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import pdb
import os
import xarray as xr
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

def saveNetcdf():

    sm_folder = '/users/global/cornkle/data/OBS/CMORPH/CMORPH_raw'
    pool = multiprocessing.Pool(processes=7)
    for y in range(2006,2011):
        files = glob.glob(sm_folder+'/'+ str(y) + '/' +'*.gra')

        for f in files:
            ds = rewrite_data.rewrite_CMORPH(f)

    #res = pool.map(rewrite_data.rewrite_AMSRE, files)


def mergeCMORPH():

    sm_folder = '/users/global/cornkle/data/OBS/CMORPH/CMORPH_nc/'

    for y in range(2006, 2011):
        files = sm_folder + str(y) + '/' + '*.nc'
        ds = xr.open_mfdataset(files)

        enc = {'pr': {'complevel': 5, 'zlib': True}}
        ds.to_netcdf(sm_folder + 'CMORPH_WA_' + str(y) + '.nc', encoding=enc, format='NETCDF4')

        print('Wrote ' + sm_folder + 'CMORPH_WA_' + str(y) + '.nc')



def to_datetime():

    sm_folder = '/users/global/cornkle/data/OBS/CMORPH/CMORPH_nc/'

    for y in range(2006, 2011):
        files = sm_folder + 'CMORPH_WA_' + str(y) + '.nc'
        dsa = xr.open_dataset(files)

        date = pd.to_datetime(dsa.time.values)

        da = xr.DataArray(dsa.pr, coords={'time': date,
                                      'lat': dsa.lat,
                                      'lon': dsa.lon},
                          dims=['time', 'lat', 'lon'])  # .isel(time=0)

        ds = xr.Dataset({'pr': da})

        enc = {'pr': {'complevel': 5, 'zlib': True}}
        ds.to_netcdf(sm_folder + 'CMORPH_WA_T' + str(y) + '.nc', encoding=enc, format='NETCDF4')

        print('Wrote ' + sm_folder + 'CMORPH_WA_T' + str(y) + '.nc')