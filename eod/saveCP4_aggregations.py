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


def saveMonthly():
    tag = 'hist'
    data_path = '/media/ck/Elements/Africa/WestAfrica/CP4/CP4'+tag+'/'
    vars = [ 'lw_out_PBLtop', 'q_pl', 'lsRain'] #'tcwv', 't_pl', 'u_pl', 'v_pl'

    out = cnst.network_data + 'data/CP4/CLOVER/CP4_monthly_'+tag+'/'

    for v in vars:
        mf = xr.open_mfdataset(data_path+v+os.sep+"*.nc", combine='nested', concat_dim='time') #

        clim = mf.groupby('time.month').mean('time')

        for m in clim.month.values:
            #ipdb.set_trace()
            month = clim.sel(month=m).load()

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in month.data_vars}

            outname = out + v + os.sep
            if not os.path.isdir(outname):
                os.mkdir(outname)
            outfile = outname + v + '_monthlyClim_'+tag+'_4km_'+str(m).zfill(2)+'.nc'
            month.to_netcdf(outfile, mode='w', encoding=encoding, format='NETCDF4')


