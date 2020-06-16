# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import ipdb
import os
import xarray as xr
import numpy as np
from utils import u_darrays as uda
from utils import constants as cnst
import pandas as pd

def monthly():

    for y in np.arange(2020,2021): #(1984,2018)

        print('Doing', y)
        inpath = glob.glob('/media/ck/Elements/SouthAmerica/GPM/daily/*.'+str(y)+'*.nc4')
        #file = glob.glob(cnst.GRIDSAT_PERU + '*'+str(y)+'*.nc')

        path = '/media/ck/Elements/SouthAmerica/GPM/monthly/'

        try:
            datout = xr.open_mfdataset(inpath, combine='nested', concat_dim='time')
        except:
            print('FAIL')
            continue

        for m in np.unique(datout['time.month']):


            monthly = datout['HQprecipitation'].sel(time=datout['time.month']==m).T.squeeze().load()

            timeout = monthly.time[0]
            #ipdb.set_trace()
            monthly = monthly.mean('time')

            monthly.name = 'precipitation'
            monthly = monthly.expand_dims({'time' : [timeout.values]})

            comp = dict(zlib=True, complevel=5)
            #ipdb.set_trace()
            encoding = {'precipitation': comp}
            out = path + str(y) + '-' + str(m).zfill(2) + '.nc'
            monthly.to_netcdf(out, mode='w', format='NETCDF4', encoding=encoding)
        del datout



