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
import cartopy
import cartopy.crs as ccrs

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

    savefile = '/users/global/cornkle/MCSfiles/blob_map_90km_3UTC.nc'

    try:
        os.remove(savefile)
    except OSError:
        pass
    da.to_netcdf(path=savefile, mode='w')

    das = da.sum(dim='time')

    das.to_netcdf('/users/global/cornkle/MCSfiles/blob_map_90km_sum_3UTC.nc')

    print('Saved ' + savefile)



def file_loop(passit):



    grid = passit[0]

    m = passit[1]
    files = passit[2]

    min_list = ['00']#, '15','30', '45']

    strr = files.split(os.sep)[-1]

    if (np.int(strr[8:10]) != 3): #& (np.int(strr[8:10]) < 18): #(np.int(strr[4:6]) != 6) &
        print('Skip')
        return

    lon, lat = grid.ll_coordinates

    ds = xr.Dataset()

    for min in min_list:

        file = files+min+'.gra'

        print('Doing file: ' + file)

        mdic = m.read_data(file)

        if not mdic:
            print('File missing')
            continue

        # interpolate MSG to salem grid
        inter, mpoints = u_grid.griddata_input(mdic['lon'].values, mdic['lat'].values, grid)

        # Interpolate TRMM using delaunay triangularization
        dummyt = griddata(mpoints, mdic['t'].values.flatten(), inter, method='linear')
        outt = dummyt.reshape((grid.ny, grid.nx))

        out = np.zeros_like(outt, dtype=np.int)

        outt[outt >= -40] = 150
        outt[np.isnan(outt)] = 150

        grad = np.gradient(outt)
        outt[outt == 150] = np.nan

        nogood = np.isnan(outt)

        outt[np.isnan(outt)] = -55
        nok = np.where(abs(grad[0]) > 80)
        d = 2
        i = nok[0]
        j = nok[1]

        for ii, jj in zip(i, j):
            kern = outt[ii - d:ii + d + 1, jj - d:jj + d + 1]
            outt[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

        wav = util.waveletT(outt, 5)

        outt[nogood] = np.nan

        arr = np.array(wav['scales'], dtype=str)

        scale_ind = range(arr.size)

        figure = np.zeros_like(outt)

        wll = wav['t']  # [nb, :, :]

        # maxoutt = (
        #     wll == ndimage.maximum_filter(wll, (5, 5), mode='reflect',
        #                                   cval=np.amax(wll) + 1))  # (np.round(orig / 5))

        yyy = []
        xxx = []
        scal = []
        for nb in scale_ind[::-1]:

            orig = float(arr[nb])

            if orig <90:#> 30:
                 continue

            scale = int(np.round(orig))

            print(np.round(orig))

            wl = wll[nb, :, :]
           # maxout = maxoutt[nb, :, :]

            maxout = (
                wl == ndimage.maximum_filter(wl, (5,5), mode='constant', cval=np.amax(wl) + 1))  # (np.round(orig / 5))

            try:
                yy, xx = np.where((maxout == 1) & (outt <= -75) & ((wl >= np.percentile(wl[wl >= 0.5], 90)) & (wl > orig**.5) ))  # )& (wl > orig**.5) (wl >= np.percentile(wl[wl >= 0.1], 90)) )#(wl > orig**.5))#  & (wlperc > orig**.5))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80)))# & (wlperc > np.percentile(wlperc[wlperc>=0.1], 80) ))  # & (wl100 > 5)
            except IndexError:
                continue

            for y, x in zip(yy, xx):

                ss = orig
                iscale = (np.ceil(ss / 2. / 5.)).astype(int)

                ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, outt)

                figure[ycirc, xcirc] = 1
                xxx.append(x)
                yyy.append(y)
                scal.append(orig)

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


    print('Did ', file)

    return (da)

