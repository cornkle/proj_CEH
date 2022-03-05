# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import ipdb
import os
import xarray as xr
import numpy as np
import pandas as pd
from utils import u_darrays as uda
from scipy.interpolate import griddata


VARS = {  'rlut' : 'toa_outgoing_longwave_flux'


}


def extr(var, dir, out):

    files = glob.glob(dir+os.sep+var+os.sep+'*.nc')

    for svar in files:

        ds = xr.open_dataset(svar)#[VARS[var]]

        ds = uda.shift_lons(ds, lon_dim='longitude')
        ds = ds.sel(longitude=slice(-18,25), latitude=slice(9,22))

        outfolder = out+os.sep+var
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)

        name = os.path.basename(svar)
        #ds.name = VARS[var]

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(path=outfolder+os.sep+name, mode='w', encoding=encoding, format='NETCDF4')