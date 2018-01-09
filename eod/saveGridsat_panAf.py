# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import pdb
import os
import xarray as xr
import numpy as np

def saveNetcdf():
    #msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_raw_binary/'
    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_raw_18Z_Africa'
    pool = multiprocessing.Pool(processes=7)
    files = glob.glob(msg_folder+'/*/*/irwin_cdr_*18.gra')

    #res = pool.map(rewrite_data.rewriteGridsat_toNetcdf, files)
    res = pool.map(rewrite_data.rewriteGridsat_extract18Z, files)


def saveYearly():

    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_raw_18Z_Africa'

    years = np.arange(1982,2017)#list(next(os.walk(msg_folder))[1])

    for y in years:

        if os.path.isfile(msg_folder+'gridsat_WA_'+str(y)+'.nc'):
            continue

        files = glob.glob(msg_folder+str(y)+'/*/irwin_cdr_*.nc')
        files.sort()

        da = xr.open_dataset(files[0])

        for f in files[1::]:
            print('Doing '+f)

            da=xr.concat([da, xr.open_dataset(f)], dim='time')
        enc = {'t':{'complevel': 5, 'shuffle': True, 'zlib': True}}
        da.to_netcdf(msg_folder+'gridsat_WA_'+str(y)+'.nc', encoding=enc)

def saveMonthly18():
    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/z18_panAfrica/'

    da = xr.open_mfdataset(msg_folder+'gridsat_WA_*18UTC.nc')
    da = da.where((da<=-40) & (da>=-110))
    da = da.resample('m', dim='time', how='mean')
    da.to_netcdf(msg_folder+'gridsat_monthly_18UTC.nc')



def saveColdClimatology():
    years = range(1984,2016)
    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'


    for y in years:
        da = xr.open_dataset(msg_folder+'gridsat_WA_'+str(y)+'.nc')