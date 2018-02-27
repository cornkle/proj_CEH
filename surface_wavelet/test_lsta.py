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
import salem
from scipy import ndimage
####these functions create plots of the average of the LSTAs and all blobs over the time period
####AVERAGE PLOTS!



def modis():
    files = '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/'

    msg = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/*.nc')
    msg = msg['LSTA']
    msg = msg.sel(lat=slice(10, 20), lon=slice(-10, 10))
    msg = msg[ (msg['time.month']>= 7) ]

    msg = msg.where(msg>-900)

    dat = msg.mean(dim='time')


    f = plt.figure()
    dat.plot.contourf(cmap='RdBu_r', vmin=-1, vmax=1)

    dat = dat.mean(dim='lon')
    f = plt.figure()
    dat.plot()


def blobs():
    file = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_points.nc'
    fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
    msg = xr.open_dataarray(file)
    msg = msg.sel(lat=slice(10, 20), lon=slice(-10, 10))
    msg = msg[ (msg['time.month'] >= 7 )  ]
    msg = msg.where(msg > -900)
    msg.values[msg.values<=-40] = 1
    msg = msg.sum(dim='time')


    map = msg.salem.get_map(cmap='viridis')
    top = xr.open_dataarray(fpath)
    f = plt.figure()
    z = map.set_topography(top, relief_factor=1.4)
    map.set_contour(z, levels=(200,400,600,800), cmap='Reds' )

    map.set_data(msg)
    map.visualize(title='Blobs and topo')

    msg = msg.sum(dim='lon')
    f = plt.figure()
    msg.plot()