# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import multiprocessing
import pdb
import pandas as pd
from scipy import ndimage
from cold_cloud_trend import era_geop_t3d as era_geop
from utils import u_gis
import pickle as pkl
import statsmodels.stats.proportion as prop

matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


def loop():
    dic = { 'scale': [],
            'blob' : [],
            'temp' : [],
}

    for l in np.arange(16,23):
        print('Doing '+str(l))
        blobs, scales, temp = composite(l)

        dic['scale'].append(np.array(scales))
        dic['blob'].append(np.array(blobs))
        dic['temp'].append(np.array(temp))



    pkl.dump(dic, open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/scalesVSblob.p", "wb"))


def plot():



    lsta_all = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/scale_maps/*.nc')

    temp_all = xr.open_mfdataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/lsta_daily_*.nc')
    temp_all = temp_all.sel(lat=slice(10.5,17.5), lon=slice(-9.5,9.5))
    lsta_all = temp_all.sel(lat=slice(10.5, 17.5), lon=slice(-9.5, 9.5))

    temp_all = temp_all.where(temp_all['time'] == lsta_all['time'])

    lsta_all = lsta_all.where(temp_all > -800)
    temp_all = temp_all.where(temp_all > -800)

    lsta_all = lsta_all.where(np.abs(temp_all['LSTA'].values) > 0.2)
    temp_all = temp_all.where(np.abs(temp_all['LSTA'].values) > 0.2)

    dic = pkl.load( open("/users/global/cornkle/figs/LSTA-bullshit/scales/new/scalesVSblob.p", "rb"))

    blob = np.squeeze(np.concatenate(dic['blob']))
    scale = np.squeeze(np.concatenate(dic['scale']))
    temp = np.squeeze(np.concatenate(dic['temp']))


    scalei = scale[np.isfinite(scale) & np.isfinite(temp)]
    blobi = blob[np.isfinite(scale) & np.isfinite(temp)]
    tempi = temp[np.isfinite(scale) & np.isfinite(temp)]


    H, xbins, ybins = np.histogram2d(tempi,np.abs(scalei) , bins = [ np.arange(-10,11,2), np.arange(0,151,15)])
    H = H.transpose() #/ np.sum(H)

    H2, xbins, ybins = np.histogram2d(temp_all['LSTA'].values.flatten(), np.abs(lsta_all['LSTA'].values.flatten()), bins=[np.arange(-10, 11, 2), np.arange(0, 151, 15)])
    #H2 = H2.transpose() / np.sum(H2)

    X,Y = np.meshgrid(xbins, ybins)

    f = plt.figure()

    plt.pcolormesh(X,Y,H, cmap='viridis')
    plt.colorbar()
    #plt.imshow(H)



def composite(h):
    pool = multiprocessing.Pool(processes=8)


    file = '/users/global/cornkle/MCSfiles/blob_map_allscales_-50_JJAS_points_dominant.nc'

    msg = xr.open_dataarray(file)
    msg = msg[(msg['time.hour'] == h) & (msg['time.minute'] == 0) & (
        msg['time.year'] >= 2006) & (msg['time.year'] <= 2010) & (msg['time.month'] >= 6) ]

    msg = msg.sel(lat=slice(10.5,17.5), lon=slice(-9.5,9.5))

    res = pool.map(file_loop, msg)
    pool.close()

    # for m in msg[0:50]:
    #     file_loop(m)
    #
    # return

    res = [x for x in res if x is not None]

    blobs = []
    scales = []
    temp = []

    for r in res:
        scales.append(r[0])
        temp.append(r[1])
        blobs.append(r[2])



    blobs = [item for sublist in blobs for item in sublist]  # flatten list of lists
    scales = [item for sublist in scales for item in sublist]  # flatten list of lists
    temp = [item for sublist in temp for item in sublist]


    return blobs, scales, temp



def cut_kernel(xpos, ypos, array):


    # kernel = lsta_day2.isel(lon=slice(xpos+3 , xpos + 4 + 1),
    #                         lat=slice(ypos-1 , ypos + 1 + 1)).values
    # #
    # #
    # kernel =  np.mean(kernel)

    kernel = array.isel(lon=xpos+5,
                             lat=ypos).values

    # if (kernel == np.nan):
    #     return

    return kernel


def file_loop(fi):
    print('Doing day: ', fi)

    date = pd.to_datetime(
        str(fi['time.year'].values) + str(fi['time.month'].values).zfill(2) + str(fi['time.day'].values).zfill(2))
    dayd = pd.Timedelta('1 days')

    if (fi['time.hour'].values) <= 15:
        print('Nighttime')
        daybefore = date - dayd
    else:
        print('Daytime')
        daybefore = date

    try:
        lsta = xr.open_dataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/scale_maps/' \
                           'lsta_daily_scale_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(daybefore.day).zfill(2) + '.nc')
    except OSError:
        return

    print('Doing '+ 'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')

    temp_da = xr.open_dataset('/users/global/cornkle/data/OBS/modis_LST/modis_netcdf/' \
                           'lsta_daily_' + str(daybefore.year) + str(daybefore.month).zfill(2) + str(
        daybefore.day).zfill(2) + '.nc')
    try:
        temp_da = temp_da.sel(lat=slice(lsta['lat'].values.min(),lsta['lat'].values.max()), lon=slice(lsta['lon'].values.min(),lsta['lon'].values.max()))
    except OSError:
        return

    lsta = lsta['LSTA']
    temp_da = temp_da['LSTA']
    lsta = lsta.where(temp_da > -800)
    temp_da = temp_da.where(temp_da > -800)


    temp = temp_da.values

    lsta = lsta.where(np.abs(temp) > 0.2)
    temp_da = temp_da.where(np.abs(temp) > 0.2)

    if np.sum(lsta) == np.nan:
        return

    pos = np.where( (fi.values >= 5) & (fi.values >= 75) )

    if np.sum(pos) == 0:
        print('No blobs found')
        return

    blob_list = []
    scale_list = []
    temp_list = []

    for y, x in zip(pos[0], pos[1]):

        blob_scale = fi.values[y,x]

        lat = fi['lat'][y]
        lon = fi['lon'][x]

        t = fi.sel(lat=lat, lon=lon)

        point = lsta.sel(lat=lat, lon=lon, method='nearest')
        plat = point['lat'].values
        plon = point['lon'].values

        xpos = np.where(lsta['lon'].values == plon)
        xpos = int(xpos[0])
        ypos = np.where(lsta['lat'].values == plat)
        ypos = int(ypos[0])
        try:
            ano = cut_kernel(xpos, ypos, lsta)
        except TypeError:
            continue

        if ano == np.nan:
            continue

        try:
            anot = cut_kernel(xpos, ypos, temp_da)
        except TypeError:
            continue

        scale_list.append(ano)
        temp_list.append(anot)
        blob_list.append(blob_scale)

    if scale_list == None:
        return None

    print('Returning')

    return (scale_list, temp_list, blob_list)

if __name__ == "__main__":
    composite()

