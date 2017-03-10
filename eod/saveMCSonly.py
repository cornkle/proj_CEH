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

    #files = files[0:5]

    # make salem grid
    grid = u_grid.make(m.lon, m.lat, 5000)

    files_str = []

    for f in files:
        files_str.append(f[0:-6])

    files_str = np.unique(files_str)


    passit = []
    for f in files_str:
        passit.append((grid,m, f))

    #res = pool.map(file_loop, passit)



    for l in passit:

        test = file_loop(l)

    #pool.close()

    return

    res = [x for x in res if x is not None]

    da = xr.concat(res, 'time')

    savefile = '/users/global/cornkle/MCSfiles/blob_map_June.nc'

    try:
        os.remove(savefile)
    except OSError:
        pass
    da.to_netcdf(path=savefile, mode='w')
    print('Saved ' + savefile)



def file_loop(passit):



    grid = passit[0]

    lon, lat = grid.ll_coordinates

    m = passit[1]
    files = passit[2]

    min_list = ['00']#, '15','30', '45']

    strr = files.split(os.sep)[-1]

    if (np.int(strr[4:6]) != 6):  #(np.int(strr[8:10]) != 18)
        #print('Skip')
        return

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
        torig = outt.copy()

        outt[outt >= -40] = 150
        outt[np.isnan(outt)] = 150

        grad = np.gradient(outt)
        outt[outt == 150] = -55

        nok = np.where(abs(grad[0]) > 80)
        d = 2
        i = nok[0]
        j = nok[1]

        for ii, jj in zip(i, j):
            kern = outt[ii - d:ii + d + 1, jj - d:jj + d + 1]
            outt[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

        wav = util.waveletT(outt, 5)
        arr = np.array(wav['scales'], dtype=str)
        origs = np.array(wav['scales'], dtype=np.float)
        scales = np.array(np.round(origs), dtype = np.int)
        wl = wav['t']  # [nb, :, :]

        maxout = (
            wl == ndimage.maximum_filter(wl, (10, 10, len(arr)-1 ), mode='constant',
                                         cval=np.amax(wl) + 1))



        scale_ind = range(arr.size)

        for nb in scale_ind:

            orig = origs[nb]
            scale = scales[nb]

            maxoutt = maxout[nb, :, :]

            try:
                yy, xx = np.where((maxoutt == 1) & (torig <= -50))
            except IndexError:
                continue
            #print(scale, yy, xx)


            for y, x in zip(yy, xx):

                ss = orig
                #ss = 15  # just looking at a fixed radius surrounding points defined by wavelet
                iscale = (np.ceil(ss / 2. / 5.)).astype(int)

                ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, out)

                out[ycirc, xcirc] = scale

        plt.imshow(out)
        plt.show()

        ipdb.set_trace()
        hour = mdic['time.hour']
        minute = mdic['time.minute' ]
        day = mdic['time.day']
        month = mdic['time.month']
        year = mdic['time.year']

    date = dt.datetime(year, month, day, hour, minute)

    da = xr.DataArray(out, coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon']) #[np.newaxis, :]


    print('Did ', file)

    return (da)

