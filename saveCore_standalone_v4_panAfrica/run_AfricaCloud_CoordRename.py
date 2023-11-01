# -*- coding: utf-8 -*-


import numpy as np
import xarray as xr
import os
import ccores.cores as cores
import datetime
import multiprocessing
from utils import constants as cnst
from eod import msg_panAfrica
import ipdb
import glob

### Function to loop over ch9 grads files.
def _loop(passit):

    print('Doing file: ' + passit)

    # Define outfile name
    outfile = passit.replace('ch9_wavelet', 'ch9_wavelet_coord')

    if os.path.isfile(outfile):
        print('File exists, continue')
        return

    # Calls package to read grads files, returns mdic object containing native MSG tir data and lat/lon coordinates.
    try:
        mdic = xr.open_dataset(passit)
    except FileNotFoundError:
        print('File not found')
        return

    newdic = mdic.rename({'lat':'y1d_dummy', 'lon' : 'x1d_dummy'})

    ## Netcdf compression step.
    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in mdic.data_vars}

    path = os.path.dirname(outfile)
    if not os.path.exists(path):
       os.makedirs(path)

    newdic.to_netcdf(path=outfile, mode='w', encoding=enc, format='NETCDF4')

    print('Saved ' + outfile)


############################
############################
# Loop initiaition for defined years, calling multiprocessing.

meteosat_folder = '/prj/Africa_cloud/wavelet/*.nc'

files = glob.glob(meteosat_folder)

pool = multiprocessing.Pool(processes=5)
res = pool.map(_loop, files)



