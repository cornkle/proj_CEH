# -*- coding: utf-8 -*-


import numpy as np
from wavelet import util
from eod import msg, mfg
import xarray as xr
import os
from wav_mixed import powerBlob_utils
import multiprocessing
from utils import u_grid, u_interpolate as u_int
import datetime as dt
from utils import constants as cnst
import pickle as pkl
import ipdb
import glob


def run():

    met_folder = '/prj/PORCELAIN/Fengyun/CTT/'


    pool = multiprocessing.Pool(processes=5)

    for y in range(2015,2019):

        if (y >=2016) & (y<=2018):
            gstr = 'FY2G'
        if (y >= 2014) & (y <= 2015):
            gstr = 'FY2F'

        gfile = '/prj/PORCELAIN/Fengyun/CTT/CTT_grid_'+gstr+'.nc'
        usefile = '/prj/PORCELAIN/Fengyun/CTT/CTT_grid_FY2G.nc'
        mdic = xr.open_dataset(gfile)
        usedic = xr.open_dataset(usefile)

        # make salem grid
        grid = u_grid.make(usedic['longitude'].values, usedic['latitude'].values, 5000)
        inds, weights, shape = u_int.interpolation_weights_grid(mdic['longitude'].values, mdic['latitude'].values, grid)
        gridd = (inds, weights, shape, grid)

        for m in range(6,10):

            files = glob.glob(met_folder + str(y) +'/' + str(m).zfill(2) +'/' + gstr +'*.nc')
            passit = []
            for f in files:
                passit.append((gridd,f))

            res = pool.map(file_loop, passit)
            #
            # for f in passit:
            #   file_loop(f)

    return

def _timeLoop(timeslice, inds_inter, weights_inter, shape_inter, tag):

    try:

        outt = u_int.interpolate_data(timeslice.values, inds_inter, weights_inter, shape_inter)
        #print('Interpolated', id)
    except ValueError:
        print('Interpolation value error!!')
        return

    outt_out = outt.copy()
    outt_out[outt_out > -40] = 0

    if np.sum(outt) == 0:   # if temperature is all empty
        print('Temperature empty, continue')
        return

    outt, nogood, t_thresh_size, t_thresh_cut, pix_nb = powerBlob_utils.filter_img(outt, 5)

    wav = util.waveletT(outt, dataset='METEOSAT5K_vera')
    outt[nogood] = np.nan

    power_msg = powerBlob_utils.find_scales_dominant(wav, nogood, dataset=tag)
    if power_msg is None:  # if power calculation failed
        print('Power calc fail, continue')
        return


    isnan = np.isnan(outt)
    outt[isnan] = 0

    new_savet = (outt_out * 100).astype(np.int16)

    return power_msg.astype(np.int16), new_savet, outt, nogood, t_thresh_size, t_thresh_cut, pix_nb, wav




def file_loop(passit):


    gridd = passit[0]
    inds_inter = gridd[0]
    weights_inter = gridd[1]
    shape_inter = gridd[2]
    grid_inter = gridd[3]

    tag = 'neutral'

    file = passit[1]

    lon, lat = grid_inter.ll_coordinates

    print('Doing file: ' + file)

    mdic = xr.open_dataset(file)


    hour = np.array(file[-7:-5]).astype(int)
    minute = np.array(file[-5:-3]).astype(int)
    day = np.array(file[-10:-8]).astype(int)
    month = np.array(file[-12:-10]).astype(int)
    year = np.array(file[-16:-12]).astype(int)

    date = dt.datetime(year, month, day, hour, minute)


    ds = xr.Dataset()
    mdic['CTT'].values[mdic['CTT'].values==0] = np.nan#-273.15
    mdic['CTT'].values = mdic['CTT'].values-273.15
    mdic['CTT'].values[mdic['CTT'].values>=-40] = 0

    timeslice = mdic['CTT']
    try:
        power_msg, new_savet, outt, nogood, t_thresh_size, t_thresh_cut, pix_nb, wav \
            = _timeLoop(timeslice, inds_inter, weights_inter, shape_inter, tag)
    except TypeError:
        return

    bloblist = xr.DataArray(power_msg, coords={'time': date, 'lat': lat[:, 0], 'lon': lon[0, :]},
                            dims=['lat', 'lon'])  # [np.newaxis, :])
    tirlist = xr.DataArray(new_savet, coords={'time': date, 'lat': lat[:, 0], 'lon': lon[0, :]},
                           dims=['lat', 'lon'])

    ds['blobs'] = bloblist
    ds['CTT'] = tirlist

    ds.attrs['radii']=(np.ceil(wav['scales'] / 2. / 5.)).astype(np.uint8)
    ds.attrs['scales_rounded'] = np.round(wav['scales']).astype(np.uint8)
    ds.attrs['scales_original'] = wav['scales']
    ds.attrs['cutout_T'] = t_thresh_size
    ds.attrs['cutout_minPixelNb'] = pix_nb

    #ds = ds.sel(lat=slice(2,17), lon=slice(-18,13))     #[-14, 2.5, 4, 11.5] cutout to remove dodgy boundaries

    out = cnst.network_data + 'data/emma_test/'
    outdir = out + str(date.year) + '/' + str(date.month).zfill(2) + '/'
    fname = os.path.basename(file)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fnew = fname.replace('MLT', 'POWER')

    savefile = outdir + fnew  # 'blobMap_-40-700km2_-50-points_dominant_'+str(yy) + '_'+str(mm).zfill(2)+'.nc'

    try:
        os.remove(savefile)
    except OSError:
        pass
    # da.name = 'blob'
    # enc = {'blob': {'complevel': 5, 'zlib': True}}

    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in ds.data_vars}

    ds = ds.sel(lat=slice(20,44), lon=slice(79,122))     #[-14, 2.5, 4, 11.5] cutout to remove dodgy boundaries
    ds.to_netcdf(path=savefile, mode='w', encoding=enc, format='NETCDF4')
    print('Saved ' + savefile)

    print('Did ', file)

    return (ds)
