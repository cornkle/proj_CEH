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
from CLOVER import era_geop_t3d as era_geop
from utils import u_gis
import pickle as pkl

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def file_loop():

    lsta = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/scale_maps/' \
                           'lsta_daily_scale_*.nc')


    lsta_check = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/' \
                           'lsta_daily_*.nc')

    lsta_check = lsta_check.sel(lat=slice(lsta['lat'].values.min(),lsta['lat'].values.max()), lon=slice(lsta['lon'].values.min(),lsta['lon'].values.max()))


    lsta_checks = lsta_check['LSTA'].where(lsta_check['LSTA']>-800)
    lsta_checks = lsta_checks.where(lsta.time==lsta_checks.time)

    bins = np.arange(-20,20,2)
    f=plt.figure()
    plt.hist(lsta_checks.values[np.isfinite(lsta_checks.values)], bins=bins, edgecolor='k')

    bins = np.arange(-140, 141, 10)

    ll = []

    for i, b in enumerate(bins[0:-1]):

        b1 = bins[i+1]

        lmean = np.percentile(lsta_checks.where((lsta['LSTA'].values>=b) &  (lsta['LSTA'].values<b1)), 90)

        ll.append(lmean)

    pdb.set_trace()
    f = plt.figure()
    plt.scatter(bins[1::], ll)


