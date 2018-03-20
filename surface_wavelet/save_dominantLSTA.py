import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from wavelet import util
import pdb
from scipy import ndimage
from utils import u_arrays as ua
import pandas as pd
import time
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
from utils import u_grid
import os
import itertools




def run_netcdf():

    for d, m, y in itertools.product(np.arange(1,32), np.arange(6,10), np.arange(2006,2010)):


        DATE = {'day': d,
                'month': m,
                'year': y}


        file = '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/lsta_daily_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'
        outfile = '/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/power_maps/lsta_daily_powerMaxscale_'+str(DATE['year'])+str(DATE['month']).zfill(2)+str(DATE['day']).zfill(2)+'.nc'

        ds = xr.open_dataset(file)
        ds = ds.sel(lon=slice(-10,10), lat=slice(10,20))
        lsta = ds['LSTA'][0,:,:]

        # f = plt.figure()
        # plt.imshow(lsta)

        lsta[lsta<-800] = np.nan


        lsta = lsta   - lsta.mean()

        points = np.where(np.isfinite(lsta.values))
        inter1 = np.where(np.isnan(lsta.values))
        try:
            lsta.values[inter1] = griddata(points, np.ravel(lsta.values[points]), inter1, method='linear')
        except ValueError:
            continue

        inter = np.where(np.isnan(lsta))
        try:
            lsta.values[inter] = griddata(points, np.ravel(lsta.values[points]), inter, method='nearest')
        except ValueError:
            continue
        #lsta[inter1]=0

        wav = util.waveletLSTA_dom(lsta.values,3)

        wl = wav['dominant']

        wl[inter[0], inter[1]] = np.nan
        wl[inter1[0], inter1[1]] = np.nan
        # f = plt.figure()
        # plt.imshow(wl, cmap='RdBu', vmin=9, vmax=120)
        scales = wav['scales']

        print(scales)

        ds['LSTA'].values = wl[None, ...]
        ds.attrs['scales'] = scales

        try:
            os.remove(outfile)
        except OSError:
            pass
        ds.to_netcdf(path=outfile, mode='w')


        print('Saved ' + outfile)


