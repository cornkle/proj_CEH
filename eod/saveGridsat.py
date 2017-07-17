# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import pdb
import os
import xarray as xr

def saveNetcdf():
    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_raw_binary/'
    pool = multiprocessing.Pool(processes=7)
    files = glob.glob(msg_folder+'/1990/*/irwin_cdr_*.gra')

    res = pool.map(rewrite_data.rewriteGridsat_toNetcdf, files)


def saveYearly():
    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'

    years = list(next(os.walk(msg_folder))[1])

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




def saveColdClimatology():
    years = range(1984,2016)
    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'


    for y in years:
        da = xr.open_dataset(msg_folder+'gridsat_WA_'+str(y)+'.nc')