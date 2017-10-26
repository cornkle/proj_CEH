# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import pdb
import os
import xarray as xr
import numpy as np

def saveNetcdf():

    modis_folder = '/users/global/cornkle/data/OBS/modis_LST/modis_raw_binary'
    pool = multiprocessing.Pool(processes=7)
    files = glob.glob(modis_folder+'/lsta_daily_2*.gra')

    for f in files:
        rewrite_data.rewriteModis_toNetcdf(f)


def saveClimNetcdf():

    modis_folder = '/users/global/cornkle/data/OBS/modis_LST/modis_raw_binary'
    pool = multiprocessing.Pool(processes=7)
    files = glob.glob(modis_folder+'/seasonal*.gra')

    #res = pool.map(rewrite_data.rewriteGridsat_toNetcdf, files)
    #res = pool.map(rewrite_data.rewriteModis_toNetcdf, files)

    for f in files:
        rewrite_data.rewriteModisClim_toNetcdf(f)
