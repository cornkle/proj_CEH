# -*- coding: utf-8 -*-


import numpy as np
import xarray as xr
import os
import ccores.cores as cores
import datetime
import multiprocessing
from utils import constants as cnst
from eod import msg_panAfrica

#### Defines path to original ch9 grads files and years to consider.
filepath = {

    'panAfrica' : [cnst.other_drive +'nflics/SSA_data/', np.arange(1,13), (2016,2022)]  # 2004, 2022
}

dataset = 'panAfrica'


### Function to loop over ch9 grads files.
def _loop(passit):

    m = passit[0]
    file = passit[1]

    print('Doing file: ' + file)

    # Define outfile name
    outfile = file.replace('ch9', 'ch9_wavelet')
    outfile = outfile.replace('.gra', '.nc')

    if os.path.isfile(outfile):
        print('File exists, continue')
        return

    # Calls package to read grads files, returns mdic object containing native MSG tir data and lat/lon coordinates.
    try:
        mdic = m.read_data(file)  # Returns MSG dataset
    except FileNotFoundError:
        print('File not found')
        return

    dat = mdic['t']             # get raw MSG data from xarray dataset object
    lls = [-25, 55, -38, 26]    # define domain for wavelet code - note, this seems to be needed to avoid errors with original domain extent.
    print('raw data shape', dat.shape)
    dat = dat.where((dat.lon >= lls[0]) & (dat.lon <= lls[1]) & (dat.lat >= lls[2]) & (dat.lat <= lls[3]), drop=True) # cut out domain
    print('filtered data shape', dat.shape)

    hour = dat['time.hour']
    minute = dat['time.minute']
    day = dat['time.day']
    month = dat['time.month']
    year = dat['time.year']

    date = [datetime.datetime(int(year), int(month), int(day), int(hour), int(minute))]

    data = dat.squeeze().values # Get numpy array from xarray data array

    ############## Start of convective core package use
    wObj = cores.dataset('METEOSAT3K_veraLS')                                   # initialises the 3km scale decomposition and defines scale range

    wObj.read_img(data, dat.lon.values, dat.lat.values, edge_smoothing=False)   # Prepares data image for wavelets. Input here: Native MSG data and native lat/lon coordinates (irregular 2d!)
    
    wObj.applyWavelet(normed='scale')
    try:
        dummy, max_location = wObj.scaleWeighting(wtype='nflics3k')
    except:
        print('Date failed, return', date)
        return

    nflics3k_da = wObj.to_dataarray(date=date, names=['cores', 'tir'])          # Returns the calculated cores and original input image in an xarray dataset as saved in the wavelet object.
                                                                                # This output can be used to subsequently save the data.

    ## Add power maxima locations and associated storm area size (in pixels, not km2) as additional information.
    nflics3k_da['PixelNb_-40C'] = xr.DataArray(max_location['area_pixels'], coords={'storm_idx': np.arange(len(max_location['area_pixels']))}, dims=['storm_idx'])
    nflics3k_da['max_lon'] = xr.DataArray(max_location['lon'], coords={'storm_idx': np.arange(len(max_location['area_pixels']))}, dims=['storm_idx'])
    nflics3k_da['max_lat'] = xr.DataArray(max_location['lat'], coords={'storm_idx': np.arange(len(max_location['area_pixels']))}, dims=['storm_idx'])

    ## Netcdf compression step.
    comp = dict(zlib=True, complevel=5)
    enc = {var: comp for var in nflics3k_da.data_vars}

    path = os.path.dirname(outfile)
    if not os.path.exists(path):
       os.makedirs(path)

    nflics3k_da.to_netcdf(path=outfile, mode='w', encoding=enc, format='NETCDF4')

    print('Saved ' + outfile)


############################
############################
# Loop initiaition for defined years, calling multiprocessing.

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

            pool = multiprocessing.Pool(processes=5)
            res = pool.map(_loop, passit)



