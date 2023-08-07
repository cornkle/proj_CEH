# -*- coding: utf-8 -*-


import numpy as np
from saveCore_standalone_NFLICS import util
import xarray as xr
from saveCore_standalone_NFLICS import powerBlob_utils
import datetime as dt


def wavelet_analysis(meteosat_data, longitudes, latitudes, date, savefile, data_resolution=5):


    outt, nogood, t_thresh_size, t_thresh_cut, pix_nb, area_img = powerBlob_utils.filter_img(meteosat_data, data_resolution)

    wav = util.waveletT(outt, dataset='METEOSAT5K_vera')
    meteosat_data[nogood] = np.nan

    power_msg = powerBlob_utils.find_scales_dominant(wav, nogood, area_img, dataset='MSG')

    if power_msg is None:  # if power calculation failed
        print('Power calc fail, continue')
        return

    #date = dt.datetime(year, month, day, hour, minute)

    isnan = np.isnan(meteosat_data)
    meteosat_data[isnan] = 0
    new_savet = (meteosat_data * 100).astype(np.int16)

    ds = xr.Dataset()

    if latitudes.ndim == 2:
        latitudes = latitudes[:,0]
    if longitudes.ndim == 2:
        longitudes = longitudes[0,:]

    blob = xr.DataArray(power_msg[np.newaxis, :], coords={'time': date, 'lat': latitudes, 'lon': longitudes},
                                dims=['time', 'lat', 'lon'])  # [np.newaxis, :])
    tir = xr.DataArray(new_savet[np.newaxis, :], coords={'time': date, 'lat': latitudes, 'lon': longitudes},
                               dims=['time','lat', 'lon'])

    ds['blobs'] = blob
    ds['tir'] = tir

    ds.attrs['radii']=(np.floor(wav['scales'] / 2. / np.float(data_resolution))).astype(np.uint8)
    ds.attrs['scales_rounded'] = np.round(wav['scales']).astype(np.uint8)
    ds.attrs['scales_original'] = wav['scales']
    ds.attrs['cutout_T'] = t_thresh_size
    ds.attrs['cutout_minPixelNb'] = pix_nb


    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in ds.data_vars}

    ds.to_netcdf(path=savefile, mode='w', encoding=enc, format='NETCDF4')
    print('Saved ' + savefile)


    return (ds)
