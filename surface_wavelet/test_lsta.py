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
from utils import constants
from scipy import ndimage
####these functions create plots of the average of the LSTAs and all blobs over the time period
####AVERAGE PLOTS!


def lsta():

    msgf = xr.open_mfdataset(constants.LSTA_NEW + '*.nc')
    msg = msgf#['LSTA']
    #msg = msg.sel(lat=slice(10, 20), lon=slice(-10, 10))
    #msg = msg[ (msg['time.month']==9) ]

    #msg = msg.where(msg>-900)

    dat = (msg['LSTA']).sum('time') / msg['NbSlot'].sum('time') #- msg.mean(dim=['lat', 'lon'])
    #dat = dat.mean('time')
    # t = msg.values
    # count = msgf['NbSlot'].values
    # mean = np.nansum(t * count, axis=0)/ np.nansum(count, axis=0)
    f = plt.figure()
    dat.plot.contourf(cmap='RdBu_r',  extend='both', vmin=-0.05, vmax=0.05) #vmin=-1, vmax=1,

    # f = plt.figure()
    # plt.contourf(mean, cmap='RdBu_r', vmin=-1, vmax=1)

    dat = dat.mean(dim='lon')
    f = plt.figure()
    dat.plot()
    plt.show()


def blobs():
    #file = '/users/global/cornkle/MCSfiles/blob_map_30km_-67_JJAS_points.nc'
    file = '/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'
    fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
    msg = xr.open_dataarray(file)
    msg = msg.sel(lat=slice(10, 20), lon=slice(-10, 10))
    msg = msg[ (msg['time.month'] >= 6 )  ]
    msg = msg.where(msg > 6)
    msg.values[msg.values>6] = 1
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
