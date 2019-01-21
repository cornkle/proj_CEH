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
from utils import u_grid, u_interpolate as u_int
import datetime as dt
import matplotlib.pyplot as plt
import pdb
from scipy.ndimage.measurements import label
from utils import constants as cnst



def run():

    msg_folder = cnst.network_data + 'data/OBS/meteosat_WA30'
    #msg_folder

    for yy in range(2013,2014):   # (2004,2016)

        for mm in [6,9]:

            pool = multiprocessing.Pool(processes=6)

            m = msg.ReadMsg(msg_folder, y1=yy, y2=yy, months=[mm])
            files  = m.fpath

            mdic = m.read_data(files[0], llbox=[-14, 2.5, 4, 11.5])

            # make salem grid
            grid = u_grid.make(mdic['lon'].values, mdic['lat'].values, 5000)
            inds, weights, shape = u_int.interpolation_weights_grid(mdic['lon'].values, mdic['lat'].values, grid)
            gridd = (inds,weights,shape, grid)

            files_str = []

            for f in files:
                files_str.append(f[0:-6])

            files_str = np.unique(files_str)

            passit = []
            for f in files_str:
                passit.append((gridd,m, f))
            #
            res = pool.map(file_loop, passit)

            # for l in passit:
            #
            #     test = file_loop(l)

            pool.close()

            res = [x for x in res if x is not None]

            ds = xr.concat(res, 'time')
            path =  '/prj/vera/cores/' # cnst.network_data + 'MCSfiles/VERA_blobs/'
            savefile = path + 'test_size2_'+str(yy) + '_'+str(mm).zfill(2)+'.nc'#'blobMap_-40-700km2_-50-points_dominant_'+str(yy) + '_'+str(mm).zfill(2)+'.nc'

            try:
                os.remove(savefile)
            except OSError:
                pass
            #da.name = 'blob'
            #enc = {'blob': {'complevel': 5, 'zlib': True}}

            comp = dict(zlib=True, complevel=5)
            enc = {var: comp for var in ds.data_vars}

            ds.to_netcdf(path=savefile, mode='w', encoding=enc, format='NETCDF4')

            print('Saved ' + savefile)



def file_loop(passit):

    gridd = passit[0]
    inds = gridd[0]
    weights = gridd[1]
    shape = gridd[2]
    grid = gridd[3]

    m = passit[1]
    files = passit[2]

    min = '00'

    strr = files.split(os.sep)[-1]

    if ((np.int(strr[4:6]) != 9) & (np.int(strr[4:6]) != 6)):
        print('Skip month')
        return

    # if not ((np.int(strr[8:10]) >= 20)): #& (np.int(strr[8:10]) <= 19) ): #((np.int(strr[8:10]) > 3)): #not ((np.int(strr[8:10]) >= 16) & (np.int(strr[8:10]) <= 19) ): #& (np.int(strr[8:10]) < 18): #(np.int(strr[4:6]) != 6) & #(np.int(strr[8:10]) != 3) , (np.int(strr[8:10]) > 3)
    #     print('Skip hour')
    #     return

    lon, lat = grid.ll_coordinates

    file = files+min+'.gra'

    print('Doing file: ' + file)
    try:
        mdic = m.read_data(file, llbox=[-14, 2.5, 4, 11.5])
    except FileNotFoundError:
        print('File not found')
        return

    if not mdic:
        print('File missing')
        return
    outt = u_int.interpolate_data(mdic['t'].values, inds, weights, shape)
    savet = outt.copy()

    t_thresh_size = -40
    t_thresh_cut = -50

    core_min = -50

    outt[outt>=t_thresh_size] = 0
    outt[np.isnan(outt)] = 0

    labels, numL = label(outt)

    u, inv = np.unique(labels, return_inverse=True)
    n = np.bincount(inv)

    pix_nb = 28

    badinds = u[(n < pix_nb)]  # all blobs with more than 1000 pixels = 25,000km2 (meteosat regridded 5km), 200pix = 5000km2, 8pix = 200km2
    # scale 30km, radius 15km ca. 700km2 circular area equals 28 pix

    for bi in badinds:
        inds = np.where(labels == bi)
        outt[inds] = 0

    outt[outt >=t_thresh_cut] = 150

    grad = np.gradient(outt)
    outt[outt == 150] = np.nan

    nogood = np.isnan(outt) # filters edge maxima later, no maxima in -40 edge area by definition!

    tdiff = np.nanmax(outt)-np.nanmin(outt) # define background temperature for image
    if tdiff > 28:  # temp difference of 28 degrees
        xmin = 15
    else:
        xmin = 10

    outt[nogood] = t_thresh_cut-xmin
    nok = np.where(abs(grad[0]) > 80)
    d = 2
    i = nok[0]
    j = nok[1]
    # edge smoothing for wavelet application
    for ii, jj in zip(i, j):
        kern = outt[ii - d:ii + d + 1, jj - d:jj + d + 1]
        outt[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

    wav = util.waveletT(outt, dataset='METEOSAT5K')

    outt[nogood] = np.nan

    arr = np.array(wav['scales'], dtype=str)

    scale_ind = range(arr.size)

    figure = np.zeros_like(outt)

    wll = wav['t']

    maxoutt = (
        wll == ndimage.maximum_filter(wll, (5,4,4), mode='reflect',   # 5,4,4
                                      cval=np.amax(wll) + 1))  # (np.round(orig / 5))

    yyy = []
    xxx = []
    scal = []
    for nb in scale_ind[::-1]:

        orig = float(arr[nb])

        scale = int(np.round(orig))

        print(np.round(orig))

        wl = wll[nb, :, :]
        maxout = maxoutt[nb, :, :]

        try:
            yy, xx = np.where((maxout == 1) & (outt <= core_min) &  (wl > orig**.5)) #&  (wl > orig**.5)) #  #  &
        except IndexError:
            continue

        print(outt[yy,xx])

        for y, x in zip(yy, xx):

            ss = orig
            iscale = (np.ceil(ss / 2. / 5.)).astype(int)

            ycirc, xcirc = ua.draw_cut_circle(x, y, iscale, outt)

            figure[ycirc, xcirc] = scale  #outt
            figure[y,x] = scale * -1
            xxx.append(x)
            yyy.append(y)
            scal.append(orig)

    figure[np.isnan(outt)] = 0

    # f = plt.figure()
    # plt.contourf(outt)
    # plt.contour(figure)

    # figure[figure == 0] = np.nan
    # f = plt.figure()
    # f.add_subplot(111)
    # plt.imshow(outt, cmap='inferno')
    # plt.imshow(figure, cmap='viridis')
    # plt.colorbar()
    # plt.plot(xxx, yyy, 'yo', markersize=3)
    # ax = f.add_subplot(132, projection=ccrs.PlateCarree())
    # plt.contourf(lon, lat, figure, cmap='viridis', transform=ccrs.PlateCarree())
    # ax.coastlines()
    # ax.add_feature(cartopy.feature.BORDERS, linestyle='--');
    #
    # plt.colorbar()
    # f.add_subplot(131)
    # plt.imshow(outt, cmap='inferno')
    #
    #
    #plt.show()

    hour = mdic['time.hour']
    minute = mdic['time.minute' ]
    day = mdic['time.day']
    month = mdic['time.month']
    year = mdic['time.year']

    date = dt.datetime(year, month, day, hour, minute)

    ds = xr.Dataset()

    isnan = np.isnan(savet)
    savet[isnan]=0
    new_savet = (savet*100).astype(np.int16)

    ds['blobs'] = xr.DataArray(figure.astype(np.int16), coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon']) #[np.newaxis, :]
    ds['tir'] = xr.DataArray(new_savet, coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon'])
 #   ds['tir_w'] = xr.DataArray(outt, coords={'time': date, 'lat': lat[:, 0], 'lon': lon[0, :]}, dims=['lat', 'lon'])
    ds['power15-19km'] = xr.DataArray((np.mean(wav['t'][0:5,:,:], axis=0)*100).astype(np.uint16), coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon'])
    ds['power32-38km'] = xr.DataArray((np.mean(wav['t'][13:17, :, :], axis=0)*100).astype(np.uint16), coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon'])
    ds['power80-90km'] = xr.DataArray((np.mean(wav['t'][29:32, :, :], axis=0)*100).astype(np.uint16), coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon'])
    ds['power160-170km'] = xr.DataArray((np.mean(wav['t'][-5:-3, :, :], axis=0)*100).astype(np.uint16), coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon'])


    ds.attrs['radii']=(np.ceil(wav['scales'] / 2. / 5.)).astype(np.uint8)
    ds.attrs['scales_rounded'] = np.round(wav['scales']).astype(np.uint8)
    ds.attrs['scales_original'] = wav['scales']
    ds.attrs['cutout_T'] = t_thresh_size
    ds.attrs['core_minT'] = core_min
    ds.attrs['cutout_minPixelNb'] = pix_nb


    print('Did ', file)

    return (ds)
