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
    mdic = m.read_data(files[0], llbox=[-11, 11, 9, 20])
    # make salem grid
    grid = u_grid.make(mdic['lon'].values, mdic['lat'].values, 5000) #m.lon, m.lat, 5000)

    files_str = []

    for f in files:
        files_str.append(f[0:-6])

    files_str = np.unique(files_str)

    passit = []
    for f in files_str:
        passit.append((grid,m, f))

    #pdb.set_trace()

    res = pool.map(file_loop, passit)

    # for l in passit:
    #
    #     test = file_loop(l)

    pool.close()


    res = [x for x in res if x is not None]

    da = xr.concat(res, 'time')
    #da = da.sum(dim='time')

    savefile = '/users/global/cornkle/MCSfiles/blob_map_MCSs_-50_JJAS_points_dominant_minT.nc'

    try:
        os.remove(savefile)
    except OSError:
        pass
    da.to_netcdf(path=savefile, mode='w')
    #
    # das = da.sum(dim='time')
    #
    # das.to_netcdf('/users/global/cornkle/MCSfiles/blob_map_35km_-67_JJAS_sum_17-19UTC.nc')

    print('Saved ' + savefile)



def file_loop(passit):


    grid = passit[0]

    m = passit[1]
    files = passit[2]

    min = '00'

    strr = files.split(os.sep)[-1]

    if ((np.int(strr[4:6]) > 9) | (np.int(strr[4:6])<6)):
        print('Skip month')
        return

    # if not ((np.int(strr[8:10]) >= 20)): #& (np.int(strr[8:10]) <= 19) ): #((np.int(strr[8:10]) > 3)): #not ((np.int(strr[8:10]) >= 16) & (np.int(strr[8:10]) <= 19) ): #& (np.int(strr[8:10]) < 18): #(np.int(strr[4:6]) != 6) & #(np.int(strr[8:10]) != 3) , (np.int(strr[8:10]) > 3)
    #     print('Skip hour')
    #     return

    lon, lat = grid.ll_coordinates


    file = files+min+'.gra'

    print('Doing file: ' + file)
    try:
        mdic = m.read_data(file, llbox=[-11, 11, 9, 20])
    except FileNotFoundError:
        print('File not found')
        return

    if not mdic:
        print('File missing')
        return

    outt = u_grid.quick_regrid(mdic['lon'].values, mdic['lat'].values,mdic['t'].values.flatten(), grid)
    figure = np.zeros_like(outt)

    outt[outt > -50] = 0

    labels, numL = label(outt)

    u, inv = np.unique(labels, return_inverse=True)
    n = np.bincount(inv)

    badinds = u[(n < 600)]  # all blobs with more than 36 pixels = 18 km x*y = 324 km2 (meteosat ca. 3km)
    goodinds = u[(n >= 600)]

    for bi in badinds:
        inds = np.where(labels == bi)
        outt[inds] = 0

    for gi in goodinds:
        if gi == 0:
            continue
        pos = np.where(labels == gi)
        ismin = np.argmin(outt[pos])
        if (outt[pos])[ismin] < -50:
            (figure[pos])[ismin] = 10

    hour = mdic['time.hour']
    minute = mdic['time.minute' ]
    day = mdic['time.day']
    month = mdic['time.month']
    year = mdic['time.year']

    date = dt.datetime(year, month, day, hour, minute)

    da = xr.DataArray(figure, coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon']) #[np.newaxis, :]

    print('Did ', file)

    return (da)

