# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import multiprocessing
import pdb
import pandas as pd
from scipy import ndimage
####these functions create plots of the average of the LSTAs and all blobs over the time period
####AVERAGE PLOTS!


def modis():
    files = '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/'

    msg = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/*.nc')
    msg = msg['LSTA']
    msg = msg.sel(lat=slice(10, 20), lon=slice(-10, 10))
    #msg = msg[ (msg['time.month']>= 8) ]

    msg = msg.where(msg>-900)

    dat = msg.mean(dim='time')

    f = plt.figure()
    dat.plot.contourf(cmap='RdBu_r', vmin=-1, vmax=1)

    dat = dat.mean(dim='lon')
    f = plt.figure()
    dat.plot()


def blobs():
    dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_17-19UTC_centrePoint.nc'
    dayp = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_0-3UTC_centrePoint.nc'
    msg = xr.open_dataarray(dayp)
    msg = msg.sel(lat=slice(10, 20), lon=slice(-10, 5))
    msg = msg[ (msg['time.month'] >= 8 ) & (msg['time.hour'] == 3 ) ]
    msg = msg.where((msg.values >= 1) & (msg.values <= 35))
    msg.values[msg.values>=1] = 1
    msg = msg.sum(dim='time')

    f = plt.figure()
    msg.plot.contourf()

    msg = msg.sum(dim='lon')
    f = plt.figure()
    msg.plot()