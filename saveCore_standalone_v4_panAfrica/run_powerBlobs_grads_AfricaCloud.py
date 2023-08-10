# -*- coding: utf-8 -*-


import numpy as np
import xarray as xr
import ipdb
import os
import glob
import ccores.cores as cores
import datetime
import multiprocessing
from utils import constants as cnst
from eod import msg_panAfrica


filepath = {

    'panAfrica' : [cnst.other_drive +'nflics/SSA_data/', np.arange(1,13), (2013,2022)]  # 2004, 2022
}

dataset = 'panAfrica'



def _loop(passit):

    m = passit[0]
    file = passit[1]

    print('Doing file: ' + file)

    outfile = file.replace('ch9', 'ch9_wavelet')
    outfile = outfile.replace('.gra', '.nc')

    if os.path.isfile(outfile):
        print('File exists, continue')
        return

    try:
        mdic = m.read_data(file)
    except FileNotFoundError:
        print('File not found')
        return

    dat = mdic['t']
    lls = [-25, 55, -38, 26]
    print('raw data shape', dat.shape)
    dat = dat.where((dat.lon >= lls[0]) & (dat.lon <= lls[1]) & (dat.lat >= lls[2]) & (dat.lat <= lls[3]), drop=True)
    print('filtered data shape', dat.shape)

    hour = dat['time.hour']
    minute = dat['time.minute']
    day = dat['time.day']
    month = dat['time.month']
    year = dat['time.year']

    date = [datetime.datetime(int(year), int(month), int(day), int(hour), int(minute))]

    data = dat.squeeze().values

    wObj = cores.dataset('METEOSAT3K_veraLS')

    wObj.read_img(data, dat.lon.values, dat.lat.values, edge_smoothing=False)
    
    wObj.applyWavelet(normed='scale')
    dummy, max_location = wObj.scaleWeighting(wtype='nflics3k')

    nflics3k_da = wObj.to_dataarray(date=date, names=['cores', 'tir'])

    ## Add power maxima locations and associated storm area size (in pixels, not km2)
    nflics3k_da['PixelNb_-40C'] = xr.DataArray(max_location['area_pixels'], coords={'storm_idx': np.arange(len(max_location['area_pixels']))}, dims=['storm_idx'])
    nflics3k_da['max_lon'] = xr.DataArray(max_location['lon'], coords={'storm_idx': np.arange(len(max_location['area_pixels']))}, dims=['storm_idx'])
    nflics3k_da['max_lat'] = xr.DataArray(max_location['lat'], coords={'storm_idx': np.arange(len(max_location['area_pixels']))}, dims=['storm_idx'])

    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in nflics3k_da.data_vars}

    path = os.path.dirname(outfile)
    if not os.path.exists(path):
       os.makedirs(path)

    nflics3k_da.to_netcdf(path=outfile, mode='w', encoding=enc, format='NETCDF4')

    print('Saved ' + outfile)


############################
############################

for yy in range((filepath[dataset])[2][0],((filepath[dataset])[2][1])+1):

    meteosat_folder = (filepath[dataset])[0]
    data_folder = meteosat_folder + '/data/'
    print('Reading ' +meteosat_folder)

    for mm in (filepath[dataset])[1]:

            m = msg_panAfrica.ReadMsg(meteosat_folder, y1=yy, y2=yy, months=[mm])

            if len(m.fpath) == 0:
                print('No files for this year or month, continue', yy, mm)
                continue

            files  = m.fpath

            passit = []
            for f in files:
                passit.append((m, f))

            #res = []
            #for pas in passit:
            #    back = _loop(pas)
            #    res.append(back)

            pool = multiprocessing.Pool(processes=5)
            res = pool.map(_loop, passit)



