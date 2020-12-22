# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import ipdb
import os
import xarray as xr
import numpy as np
import pandas as pd
from utils import constants as cnst, u_met
from scipy.interpolate import griddata


def extr(var, dir, out):

    files = glob.glob(dir+'*.nc')

    for svar in files:

        ds = xr.open_dataset(svar).sel(longitude=slice(-18,25), latitude=slice(9,22))

        outfolder = out+os.sep+var
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(path=outfolder, mode='w', encoding=encoding, format='NETCDF4')