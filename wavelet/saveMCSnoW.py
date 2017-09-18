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


def run():
    #  (1174, 378)
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_WA30'
    pool = multiprocessing.Pool(processes=7)

    m = msg.ReadMsg(msg_folder)

    files  = m.fpath

    #files = files[0:1]

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

    savefile = '/users/global/cornkle/MCSfiles/blob_map_MCS_-73_JJAS_16-19UTC.nc'

    try:
        os.remove(savefile)
    except OSError:
        pass
    da.to_netcdf(path=savefile, mode='w')

    das = da.sum(dim='time')

    das.to_netcdf('/users/global/cornkle/MCSfiles/blob_map_MCS_-73_JJAS_sum_16-19UTC.nc')

    print('Saved ' + savefile)



def file_loop(passit):

    grid = passit[0]

    m = passit[1]
    files = passit[2]

    min_list = ['00']#, '15','30', '45']

    strr = files.split(os.sep)[-1]

    if ((np.int(strr[4:6]) > 9) | (np.int(strr[4:6])<6)):
        print('Skip month')
        return

    if not ((np.int(strr[8:10]) >= 16) & (np.int(strr[8:10]) <= 19) ): #((np.int(strr[8:10]) > 3)): #not ((np.int(strr[8:10]) >= 16) & (np.int(strr[8:10]) <= 19) ): #& (np.int(strr[8:10]) < 18): #(np.int(strr[4:6]) != 6) & #(np.int(strr[8:10]) != 3) , (np.int(strr[8:10]) > 3)
        print('Skip hour')
        return

    lon, lat = grid.ll_coordinates

    ds = xr.Dataset()

    for min in min_list:

        file = files+min+'.gra'

        print('Doing file: ' + file)
        try:
            mdic = m.read_data(file)
        except FileNotFoundError:
            print('File not found')
            return

        if not mdic:
            print('File missing')
            return

        # interpolate MSG to salem grid
        inter, mpoints = u_grid.griddata_input(mdic['lon'].values, mdic['lat'].values, grid)

        # Interpolate TRMM using delaunay triangularization
        dummyt = griddata(mpoints, mdic['t'].values.flatten(), inter, method='linear')
        outt = dummyt.reshape((grid.ny, grid.nx))

        # hour = mdic['time.hour']
        # minute = mdic['time.minute']
        # day = mdic['time.day']
        # month = mdic['time.month']
        # year = mdic['time.year']
        #
        # date = dt.datetime(year, month, day, hour, minute)
        #
        # da = xr.DataArray(outt, coords={'time': date, 'lat': lat[:, 0], 'lon': lon[0, :]},
        #                   dims=['lat', 'lon'])  # [np.newaxis, :]
        #
        # da.to_netcdf('/users/global/cornkle/test.nc')
        # return

        figure = np.zeros_like(outt, dtype=np.int)

        figure[outt<=-73] += 1

        figure[np.isnan(outt)] = 0

       # # figure[figure == 0] = np.nan
       #  f = plt.figure()
       #  f.add_subplot(133)
       #  plt.imshow(outt, cmap='inferno')
       #  plt.imshow(figure, cmap='viridis')
       #  ax = f.add_subplot(132, projection=ccrs.PlateCarree())
       #  plt.contourf(lon, lat, figure, cmap='viridis', transform=ccrs.PlateCarree())
       #  ax.coastlines()
       #  ax.add_feature(cartopy.feature.BORDERS, linestyle='--');
       #
       #  plt.colorbar()
       #  f.add_subplot(131)
       #  plt.imshow(outt, cmap='inferno')
       #
       #  plt.plot(xxx, yyy, 'yo', markersize=3)
       #  plt.show()

        print(np.sum(figure))

        if np.sum(figure) < 10:
            return

        hour = mdic['time.hour']
        minute = mdic['time.minute' ]
        day = mdic['time.day']
        month = mdic['time.month']
        year = mdic['time.year']

    date = dt.datetime(year, month, day, hour, minute)

    da = xr.DataArray(figure, coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon']) #[np.newaxis, :]

    #da.to_netcdf('/users/global/cornkle/MCSfiles/blob_maps_0-4UTC_-65/'+str(date)+'.nc')

    print('Did ', file)

    return (da)

