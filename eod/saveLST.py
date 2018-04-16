# -*- coding: utf-8 -*-


import numpy as np
from wavelet import util
from eod import msg
import xarray as xr
import os
from utils import u_grid
from scipy.interpolate import griddata
from scipy import ndimage
from utils import u_arrays as ua
import multiprocessing
import datetime as dt
import matplotlib.pyplot as plt
import pdb
import pandas as pd
import glob
from scipy.ndimage.measurements import label
import cartopy
import cartopy.crs as ccrs



def run():
    #  (1174, 378)
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_WA30'
    pool = multiprocessing.Pool(processes=7)

    m = msg.ReadMsg(msg_folder, y1=2006, y2=2010)
    files  = m.fpath

    #files = files[1050:1057]
    mdic = m.read_data(files[0])
    # make salem grid
    grid = u_grid.make(mdic['lon'].values, mdic['lat'].values, 5000) #m.lon, m.lat, 5000)


    files_str = []

    for f in files:
        # if f[-8:-6] != 6:
        #     continue
        files_str.append(f[0:-8])  # we keep only daily file names and deal with hours and minutes in loop


    files_str = np.unique(files_str)

    passit = []
    for f in files_str:
        passit.append((grid,m, f))


    # res = pool.map(file_loop, passit)
    #
    # pool.close()
    for p in passit:
        file_loop(p)


def file_loop(passit):


    grid = passit[0]

    m = passit[1]
    files = passit[2]

    hfiles = glob.glob(files+'*[!_182].gra')

    for id, file in enumerate(hfiles):

        strr = file.split(os.sep)[-1]
        if ((np.int(strr[4:6]) > 9) | (np.int(strr[4:6])<6)):
            print('Skip month')
            continue

        if ((np.int(strr[8:10]) <7 ) | (np.int(strr[8:10]) > 13) ): #((np.int(strr[8:10]) > 3)): #not ((np.int(strr[8:10]) >= 16) & (np.int(strr[8:10]) <= 19) ): #& (np.int(strr[8:10]) < 18): #(np.int(strr[4:6]) != 6) & #(np.int(strr[8:10]) != 3) , (np.int(strr[8:10]) > 3)
        #     print('Skip hour')
             continue

        lon, lat = grid.ll_coordinates

        print('Doing file: ' + file)

        file2 = '/users/global/cornkle/data/OBS/MSG_LSTA/lsta_netcdf/lsta_daily_20060624.nc'
        ds = xr.open_dataset(file2)

        try:
            mdic = m.read_data(file, llbox=[ds.lon.min().values, ds.lon.max().values+0.1, ds.lat.min().values, ds.lat.max().values])
        except FileNotFoundError:
            print('File not found')
            return
        if not mdic:
            print('File missing')
            return

        mdic = mdic.sel(lon=slice(ds.lon.min().values, ds.lon.max().values), lat=slice(ds.lat.min().values, ds.lat.max().values))

       # outt = u_grid.quick_regrid(mdic['lon'].values, mdic['lat'].values,mdic['t'].values.flatten(), grid)
        outt = np.array(mdic['t'].values, dtype=float)

        outt[outt < 20] = np.nan

        try:
            astack = np.concatenate((astack, outt[None,...]), axis=0)
        except UnboundLocalError:
            astack = outt[None, ...]

    hour = 0
    minute = 0
    day = mdic['time.day'].values
    month = mdic['time.month'].values
    year = mdic['time.year'].values

    date = [pd.datetime(year, month, day, hour, minute)]
    figure = np.nanmean(astack, axis=0)
    isvalid = np.sum(np.isfinite(astack), axis=0)
    figure[isvalid<14] = np.nan
    pdb.set_trace()
    da = xr.DataArray(figure[None, ...], coords={'time': date, 'lat': mdic.lat.values, 'lon':mdic.lon.values},
                      dims=['time', 'lat', 'lon']) #[np.newaxis, :]
    ds = xr.Dataset({'LSTA': da})

    ds.to_netcdf('/users/global/cornkle/data/OBS/MSG_LSTA/lst_new/lst_daily_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'.nc')

    print('Did ', file)

