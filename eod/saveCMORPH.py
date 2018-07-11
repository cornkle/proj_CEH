# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import pdb
import os
import xarray as xr
import salem
from utils import constants
import numpy as np

def saveNetcdf():

    sm_folder = '/users/global/cornkle/data/OBS/CMORPH/CMORPH_raw'
    pool = multiprocessing.Pool(processes=7)
    for y in range(2006,2011):
        outfolder = sm_folder+'/'+ str(y)
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        files = glob.glob( outfolder + '/' +'*.gra')

        for f in files:
            ds = rewrite_data.rewrite_CMORPH(f)

    #res = pool.map(rewrite_data.rewrite_AMSRE, files)


def mergeCMORPH():

    sm_folder = '/users/global/cornkle/data/OBS/CMORPH/CMORPH_nc/'

    for y in range(2006, 2011):
        files = sm_folder + str(y) + '/' + '*.nc'
        ds = xr.open_mfdataset(files)

        enc = {'pr': {'complevel': 5, 'zlib': True}}
        ds.to_netcdf(sm_folder + 'CMORPH_WA_' + str(y) + '.nc', encoding=enc, format='NETCDF4')

        print('Wrote ' + sm_folder + 'CMORPH_WA_' + str(y) + '.nc')



def regrid(cmorph):

    dummy = xr.open_dataset(constants.LSTA_TESTFILE)
    cm = xr.open_dataarray(cmorph)

    out = cmorph.replace('WA_', 'WA_onLSTA_')

    arrays = []
    for c in cm:

        c_on_lsta = dummy.salem.transform(c)
        arrays.append(c_on_lsta)

    astack = np.stack(arrays, axis=0)
    da = xr.DataArray(astack, coords={'time': cm.time,
                                  'lat': dummy.lat,
                                  'lon': dummy.lon},
                      dims=['time', 'lat', 'lon'])  # .isel(time=0)

    da.to_netcdf(out)

def regrid_simpler(cmorph):

    dummy = xr.open_dataset(constants.LSTA_TESTFILE)
    cm = xr.open_dataarray(cmorph)

    out = cmorph.replace('WA_', 'WA_onLSTA_')

    cm_on_lst = dummy.salem.transform(cm)
    enc = {'pr': {'complevel': 5, 'zlib': True}}

    cm_on_lst.to_netcdf(out, encoding=enc, format='NETCDF4')