# -*- coding: utf-8 -*-


import numpy as np
from wavelet import util
import xarray as xr
import os
import multiprocessing
import glob
from wav_mixed import powerBlob_utils


def run():

    met_folder = '/prj/vera/cores/' #cnst.network_data + '/data/vera_test/'

    pool = multiprocessing.Pool(processes=6)

    files = glob.glob(met_folder + 'cores_*-40_700km2*.nc')

    res = pool.map(file_loop, files)

    #for f in files:
    #    file_loop(f)

    return



def file_loop(passit):

    ds = xr.open_dataset(passit)

    dataset = passit[6:9]

    ds['tir'].values = ds['tir']

    bloblist = []
    tirlist = []
    lat = ds.lat
    lon = ds.lon

    for ids, day in enumerate(ds['tir']):

        print('id', ids)

        date = day.time
        day.values = day / 100
        if np.sum(day.values) == 0:
            continue
        img, nogood, t_thresh_size, t_thresh_cut, pix_nb = powerBlob_utils.filter_img(day.values, 5)

        power = util.waveletT(img, dataset='METEOSAT5K_vera')
        power_out = powerBlob_utils.find_scales_dominant(power, nogood, dataset=dataset)

        if power_out is None:
            continue

        new_savet = (day.values*100).astype(np.int16)
        bloblist.append(xr.DataArray(power_out.astype(np.int16), coords={'time': date, 'lat': lat, 'lon': lon},
                                     dims=['lat', 'lon']))  # [np.newaxis, :])
        tirlist.append(xr.DataArray(new_savet, coords={'time': date, 'lat': lat, 'lon': lon}, dims=['lat', 'lon']))

    ds_mfg = xr.Dataset()
    ds_mfg['blobs'] = xr.concat(bloblist, 'time')
    ds_mfg['tir'] = xr.concat(tirlist, 'time')
    ds_mfg.sel(lat=slice(5, 12), lon=slice(-13, 13))


    savefile = passit.replace('cores_', 'coresPower_')
    try:
        os.remove(savefile)
    except OSError:
        pass

    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in ds_mfg.data_vars}

    ds_mfg.to_netcdf(path=savefile, mode='w', encoding=enc, format='NETCDF4')

    print('Saved ' + savefile)