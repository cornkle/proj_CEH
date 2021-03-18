# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import pdb
import pandas as pd
from utils import u_met, u_parallelise, u_gis, u_arrays, constants as cnst
import ipdb
import pickle as pkl
import glob


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def monthly_mean():
    #pool = multiprocessing.Pool(processes=8)

    path = cnst.network_data + 'figs/NFLICS/LSTA_stats_study/'

    for y in range(2004,2016):


        for m in range(6,10):

            files = glob.glob(cnst.elements_drive + '/Africa/WestAfrica/cores_bigDomain/*_' + str(y) +  str(m).zfill(2) +'.nc')
            if len(files) >1:
                print('Found too many files, return')

            datout = xr.open_dataset(files[0])

            monthly = datout['small_scale'].sel(((datout['time.hour'] >= 12) & (datout['time.hour'] <= 23)) | (datout['time.hour'] == 0))

            monthly = monthly.where(monthly.values>10).resample(time='1D').sum()
            monthly.values[monthly.values>0] = 1
            monthly_means = monthly.mean('time')
            monthly_means.to_netcdf(cnst.elements_drive + '/Africa/WestAfrica/cores_bigDomain_dayFrequency12-0/coresPower_smallScaleMeans_' + str(y) +  str(m).zfill(2) +'.nc')



def clim_mean():

    for m in range(7,10):

        allmonths = []

        files = glob.glob(cnst.elements_drive + '/Africa/WestAfrica/cores_bigDomain/*_'+  str(m).zfill(2) +'.nc')

        for f in files:

            datout = xr.open_dataset(f)
            print('Doing ', f)

            monthly = datout['small_scale']#[((datout['time.hour'] >= 12) & (datout['time.hour'] <= 23)) | (datout['time.hour'] == 0)]
            monthly = monthly.resample(time='1D').max('time')
            #mask = monthly.values > 10
            mask = monthly.values > 10
            monthly.values = np.array(mask, dtype=int)
            allmonths.append(monthly)

            del datout

        means = xr.concat(allmonths, dim='time')

        monthly_means = means.mean('time')
        monthly_means.to_netcdf(
            cnst.elements_drive + '/Africa/WestAfrica/cores_climMean_dayFrequency12-0/coresPower_smallScaleMeans_' + str(
                m).zfill(2) + '.nc')

        del means