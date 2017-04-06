# -*- coding: utf-8 -*-


import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
from wavelet import util
from eod import msg
import xarray as xr
import os
import ipdb
import matplotlib.pyplot as plt
from utils import u_grid
from scipy.interpolate import griddata
from scipy import ndimage
from utils import u_arrays as ua
import multiprocessing
import datetime as dt

def run():
    #  (1174, 378)
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_WA30'
    pool = multiprocessing.Pool(processes=7)



    m = msg.ReadMsg(msg_folder)
    files  = m.fpath

    files = files

    # make salem grid
    grid = u_grid.make(m.lon, m.lat, 5000)

    files_str = []

    for f in files:
        files_str.append(f[0:-6])

    files_str = np.unique(files_str)


    passit = []
    for f in files_str:
        passit.append((grid,m, f))

    res = pool.map(file_loop, passit)


    #
    # for l in passit:
    #
    #     test = file_loop(l)

    pool.close()

    #return

    res = [x for x in res if x is not None]

    da = xr.concat(res, 'time')
    #da = da.sum(dim='time')

    savefile = '/users/global/cornkle/MCSfiles/blob_map_noscthresh.nc'

    try:
        os.remove(savefile)
    except OSError:
        pass
    da.to_netcdf(path=savefile, mode='w')
    print('Saved ' + savefile)



def file_loop(passit):



    grid = passit[0]

    m = passit[1]
    files = passit[2]

    min_list = ['00']#, '15','30', '45']

    strr = files.split(os.sep)[-1]

    if (np.int(strr[8:10]) != 18): #& (np.int(strr[8:10]) < 18): #(np.int(strr[4:6]) != 6) &
        print('Skip')
        return

    lon, lat = grid.ll_coordinates

    ds = xr.Dataset()


    file = files+'00'+'.gra'

    print('Doing file: ' + file)

    mdic = m.read_data(file)

    if not mdic:
        print('File missing')
        return

    # interpolate MSG to salem grid
    inter, mpoints = u_grid.griddata_input(mdic['lon'].values, mdic['lat'].values, grid)

    # Interpolate TRMM using delaunay triangularization
    dummyt = griddata(mpoints, mdic['t'].values.flatten(), inter, method='linear')
    outt = dummyt.reshape((grid.ny, grid.nx))

    outt[outt >= -40] = np.nan
    outt[np.isfinite(outt)]=1

    if np.sum(outt) < 10:
        return

    hour = mdic['time.hour']
    minute = mdic['time.minute' ]
    day = mdic['time.day']
    month = mdic['time.month']
    year = mdic['time.year']

    date = dt.datetime(year, month, day, hour, minute)

    da = xr.DataArray(outt, coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon']) #[np.newaxis, :]


    print('Did ', file)

    return (da)

