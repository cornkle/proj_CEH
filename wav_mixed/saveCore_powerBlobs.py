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

filepath = {

    'MFG_JJAS' : [cnst.network_data + '/data/OBS/MFG_JJAS/', [6,7,8,9], (1982,2005)],
    'MFG_MAMON' : [cnst.network_data +'/data/OBS/MFG_MAMON/', [3,4,5,10,11], (1982,2003)],
    'MSG_JJAS' : [cnst.network_data +'/data/OBS/MSG_WA30/', [6,7,8,9], (2004,2017)],
    'MSG_MAMON' : [cnst.network_data +'/data/OBS/MSG_tropWA/', [3,4,5,10,11], (2004,2015)],
}

def run(dataset):

    for yy in range((filepath[dataset])[2][0],((filepath[dataset])[2][1])+1):   # (2004,2016)

        for mm in (filepath[dataset])[1]:

            tag = dataset[0:3].upper()

            pool = multiprocessing.Pool(processes=3)
            print('Reading '+filepath[dataset][0])
            meteosat_folder = (filepath[dataset])[0]

            if tag == 'MFG':
                m = mfg.ReadMfg(meteosat_folder, y1=yy, y2=yy, months=[mm])
            if tag == 'MSG':
                m = msg.ReadMsg(meteosat_folder, y1=yy, y2=yy, months=[mm])

            files  = m.fpath

            gridll = pkl.load( open (cnst.network_data + 'data/OBS/saves/VERA_msg_latlon_18W12E_1N17N.p', 'rb'))

            mdic = m.read_data(files[0], llbox=[-25, 20, 2, 25])  #[-14, 2.5, 4, 11.5]

            # make salem grid
            grid = u_grid.make(gridll['lon'].values, gridll['lat'].values, 5000)
            inds, weights, shape = u_int.interpolation_weights_grid(mdic['lon'].values, mdic['lat'].values, grid)
            gridd = (inds,weights,shape, grid)

            files_str = []

            for f in files:
                if tag == 'MSG':
                    files_str.append(f[0:-4])
                if tag == 'MFG':
                    files_str.append(f)

            files_str = np.unique(files_str)

            passit = []
            for f in files_str:
                passit.append((gridd,m, f,tag))

            res = pool.map(file_loop, passit)

            # res=[]
            # for l in passit:
            #
            #     res.append(file_loop(l))

            pool.close()

            res = [x for x in res if x is not None]
            try:
                ds = xr.concat(res, 'time')
            except ValueError:
                return

            path =  '/prj/vera/cores/' # cnst.network_data + 'MCSfiles/VERA_blobs/'
            savefile = path + 'coresPower_'+tag.upper()+'_-40_700km2_-50points_dominant_'+str(yy) + '_'+str(mm).zfill(2)+'.nc'#'blobMap_-40-700km2_-50-points_dominant_'+str(yy) + '_'+str(mm).zfill(2)+'.nc'

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


def _timeLoop(timeslice, inds_inter, weights_inter, shape_inter, tag):

    try:

        outt = u_int.interpolate_data(timeslice.values, inds_inter, weights_inter, shape_inter)
        #print('Interpolated', id)
    except ValueError:
        print('Interpolation value error!!')
        return

    outt_out = outt.copy()

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

    hour = timeslice['time.hour']
    minute = timeslice['time.minute']
    day = timeslice['time.day']
    month = timeslice['time.month']
    year = timeslice['time.year']

    date = dt.datetime(year, month, day, hour, minute)

    isnan = np.isnan(outt)
    outt[isnan] = 0
    new_savet = (outt_out * 100).astype(np.int16)

    return date, power_msg.astype(np.int16), new_savet, outt, nogood, t_thresh_size, t_thresh_cut, pix_nb, wav


def file_loop(passit):

    gridd = passit[0]
    inds_inter = gridd[0]
    weights_inter = gridd[1]
    shape_inter = gridd[2]
    grid_inter = gridd[3]

    m = passit[1]
    file = passit[2]
    tag = passit[3]


    if tag == 'MFG':
            print(' month')
    if tag == 'MSG':
        strr = file.split(os.sep)[-1]

        if ((strr[-2::]) != '00') & ((strr[-2::]) != '30'):
            print('Skip minute')
            return

    # if not ((np.int(strr[8:10]) >= 20)): #& (np.int(strr[8:10]) <= 19) ): #((np.int(strr[8:10]) > 3)): #not ((np.int(strr[8:10]) >= 16) & (np.int(strr[8:10]) <= 19) ): #& (np.int(strr[8:10]) < 18): #(np.int(strr[4:6]) != 6) & #(np.int(strr[8:10]) != 3) , (np.int(strr[8:10]) > 3)
    #     print('Skip hour')
    #     return

    lon, lat = grid_inter.ll_coordinates

    if tag != 'MFG':
        file = file +'.gra'

    print('Doing file: ' + file)
    try:
        mdic = m.read_data(file, llbox=[-25, 20, 2, 25])
    except FileNotFoundError:
        print('File not found')
        return

    if not mdic:
        print('File missing')
        return

    bloblist = []
    tirlist = []
    ds = xr.Dataset()

    if tag == 'MFG':

        for id, timeslice in enumerate(mdic['t']):

            try:
                date, power_msg, new_savet, outt, nogood, t_thresh_size, t_thresh_cut, pix_nb, wav \
                    = _timeLoop(timeslice, inds_inter, weights_inter, shape_inter,tag)
            except TypeError:
                continue

            bloblist.append(xr.DataArray(power_msg, coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon'])) #[np.newaxis, :])
            tirlist.append(xr.DataArray(new_savet, coords={'time': date, 'lat': lat[:,0], 'lon':lon[0,:]}, dims=['lat', 'lon']))
        try:
            ds['blobs'] = xr.concat(bloblist, 'time')
            ds['tir'] = xr.concat(tirlist, 'time')
        except ValueError:
            return

    if tag == 'MSG':
        timeslice = mdic['t']
        try:
            date, power_msg, new_savet, outt, nogood, t_thresh_size, t_thresh_cut, pix_nb, wav \
                = _timeLoop(timeslice, inds_inter, weights_inter, shape_inter, tag)
        except TypeError:
            return

        bloblist = xr.DataArray(power_msg, coords={'time': date, 'lat': lat[:, 0], 'lon': lon[0, :]},
                                dims=['lat', 'lon'])  # [np.newaxis, :])
        tirlist = xr.DataArray(new_savet, coords={'time': date, 'lat': lat[:, 0], 'lon': lon[0, :]},
                               dims=['lat', 'lon'])

        ds['blobs'] = bloblist
        ds['tir'] = tirlist

    ds.attrs['radii']=(np.ceil(wav['scales'] / 2. / 5.)).astype(np.uint8)
    ds.attrs['scales_rounded'] = np.round(wav['scales']).astype(np.uint8)
    ds.attrs['scales_original'] = wav['scales']
    ds.attrs['cutout_T'] = t_thresh_size
    ds.attrs['cutout_minPixelNb'] = pix_nb

    ds = ds.sel(lat=slice(2,17), lon=slice(-18,13))     #[-14, 2.5, 4, 11.5] cutout to remove dodgy boundaries

    print('Did ', file)

    return (ds)
