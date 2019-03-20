# -*- coding: utf-8 -*-


import numpy as np
from wavelet import util
import xarray as xr
import os
from scipy import ndimage
import multiprocessing
from scipy.ndimage.measurements import label
from utils import constants as cnst
import ipdb
import glob


def run():

    met_folder = cnst.network_data + '/data/vera_test/'

    pool = multiprocessing.Pool(processes=5)

    files = glob.glob(met_folder + 'cores_-40_700km2*.nc')

    #res = pool.map(file_loop, files)

    for f in files:
        file_loop(f)

    return


def filter_img(inarr):

        outt = inarr.copy()
        print('outmin', np.nanmin(outt), np.nanmax(outt))

        t_thresh_size = -40
        t_thresh_cut = -50

        outt[outt >= t_thresh_size] = 0
        outt[np.isnan(outt)] = 0

        labels, numL = label(outt)

        u, inv = np.unique(labels, return_inverse=True)
        n = np.bincount(inv)

        pix_nb = 28

        badinds = u[(
                    n < pix_nb)]  # all blobs with more than 1000 pixels = 25,000km2 (meteosat regridded 5km), 200pix = 5000km2, 8pix = 200km2
        # scale 30km, radius 15km ca. 700km2 circular area equals 28 pix

        for bi in badinds:
            inds = np.where(labels == bi)
            outt[inds] = 0

        outt[outt >= t_thresh_cut] = 150

        grad = np.gradient(outt)
        outt[outt == 150] = np.nan

        nogood = np.isnan(outt)  # filters edge maxima later, no maxima in -40 edge area by definition!

        # tdiff = np.nanmax(outt) - np.nanmin(outt)  # define background temperature for image
        # if tdiff > 28:  # temp difference of 28 degrees
        #     xmin = 15
        # else:
        #     xmin = 10

        xmin = 10
        outt[nogood] = t_thresh_cut - xmin
        nok = np.where(abs(grad[0]) > 80)
        d = 2
        i = nok[0]
        j = nok[1]
        # edge smoothing for wavelet application
        for ii, jj in zip(i, j):
            kern = outt[ii - d:ii + d + 1, jj - d:jj + d + 1]
            outt[ii - d:ii + d + 1, jj - d:jj + d + 1] = ndimage.gaussian_filter(kern, 3, mode='nearest')

        return outt, nogood, t_thresh_size, t_thresh_cut, pix_nb


def find_scales_dominant(wav, outt, core_min=None, no_good=None, dataset=None):
    outt[no_good] = np.nan

    wll = wav['t']

    power_img = np.sum(wll, axis=0)
    power_img[no_good] = 0

    if 'MFG' in dataset:
        smaller = -17
        thresh_p = np.sum((wav['scales'] + smaller) ** .5)
        power_img[(power_img < np.percentile(power_img[power_img > 1], 25)) | (power_img < (thresh_p))] = 0
    else:
        smaller = -8
        thresh_p = np.sum((wav['scales'] + smaller) ** .5)
        power_img[(power_img < np.percentile(power_img[power_img > 1], 25)) | (power_img < (thresh_p))] = 0


    labels, numL = label(power_img)
    u, inv = np.unique(labels, return_inverse=True)

    for inds in u:
        if inds == 0:
            continue

        arr = power_img.copy()
        arr[np.where(labels != inds)] = 0
        power_img.flat[np.argmax(arr)] = -999

    return power_img



def file_loop(passit):

    ds = xr.open_dataset(passit)

    ds['tir'].values = ds['tir'] / 100

    bloblist = []
    tirlist = []
    lat = ds.lat
    lon = ds.lon

    for day in ds['tir']:
        date = day.time
        day.values = day / 100
        if np.sum(day.values) == 0:
            continue
        img, nogood, t_thresh_size, t_thresh_cut, pix_nb = filter_img(day.values)

        power = util.waveletT(img, dataset='METEOSAT5K_vera')
        power_mfg = find_scales_dominant(power, img, no_good=nogood, core_min=-50, dataset=passit)

        bloblist.append(xr.DataArray(power_mfg.astype(np.int16), coords={'time': date, 'lat': lat, 'lon': lon},
                                     dims=['lat', 'lon']))  # [np.newaxis, :])
        tirlist.append(xr.DataArray(day.values, coords={'time': date, 'lat': lat, 'lon': lon}, dims=['lat', 'lon']))

    ds_mfg = xr.Dataset()
    ds_mfg['blobs'] = xr.concat(bloblist, 'time')
    ds_mfg['tir'] = xr.concat(tirlist, 'time')
    ds_mfg.sel(lat=slice(5, 12), lon=slice(-13, 13))


    savefile = passit.replace('core_', 'corePower_')

    try:
        os.remove(savefile)
    except OSError:
        pass

    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in ds.data_vars}

    ds_mfg.to_netcdf(path=savefile, mode='w', encoding=enc, format='NETCDF4')

    print('Saved ' + savefile)
