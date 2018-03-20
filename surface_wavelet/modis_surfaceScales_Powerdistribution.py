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
from cold_cloud_trend import era_geop_t3d as era_geop
from utils import u_gis
import pickle as pkl

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def scaleVSpower():

    power = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/power_maps/' \
                           'lsta_daily_power*.nc')


    scale = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/scale_maps/' \
                           'lsta_daily_scale*.nc')


    scales = np.unique(scale['LSTA'].values[0:300,:,:])
    scales = scales[np.isfinite(scales)]

    power_arr = power['LSTA'][0:300]
    scale_arr = scale['LSTA'][0:300]

    mlist = []

    for s in scales:
        print('Doing '+str(s))
        mean = np.nanmean(power_arr.where(scale_arr.values == s).values)
        mlist.append(mean)


    f= plt.figure()

    plt.scatter(scales,mlist)

def powerDistribution():

    power = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/power_maps/' \
                           'lsta_daily_power*.nc')


    scale = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/scale_maps/' \
                           'lsta_daily_scale*.nc')


    scales = np.unique(np.abs(scale['LSTA'].values[0:300,:,:]))
    scales = scales[np.isfinite(scales)]

    power_arr = power['LSTA'][0:300]
    scale_arr = scale['LSTA'][0:300]

    mlist = []

    for s in scales:
        print('Doing '+str(s))
        all = power_arr.where(np.abs(scale_arr.values) == s).values
        all = all[np.isfinite(all)]
        f = plt.figure()
        plt.hist(all, bins=[-50, -25, -15, -8, -4, -2, 0,2,4,8,16,25,50])







