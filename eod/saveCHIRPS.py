# -*- coding: utf-8 -*-


import glob
import xarray as xr
from shared.utils import constants as cnst
from salem import GeoTiff
import numpy as np
import pandas as pd

def saveNetcdf():

    folder = cnst.CHIRPS + 'monthly/*.tif'
    files = glob.glob(folder)
    for f in files:
        loop(f)


def loop(lst):

    g = GeoTiff(lst)

    ls = g.get_vardata()
    lon = g.grid.ll_coordinates[0]
    lat = g.grid.ll_coordinates[1]

    ls = np.array(ls, dtype=float)
    ls[ls==-9999] = np.nan
    year = lst[-11:-7]
    month = lst[-6:-4]

    date = pd.Timestamp(year+'-'+month)

    da = xr.DataArray(ls, coords={'time': date, 'lat': lat[:, 0], 'lon': lon[0, :]},
                            dims=['lat', 'lon'])  # [np.newaxis, :])

    ds = xr.Dataset()
    ds['precip'] = da

    enc = {'precip': {'complevel': 5, 'zlib': True}}
    out = lst.replace('.tif', '.nc')
    ds.to_netcdf(path=out, mode='w', encoding=enc, format='NETCDF4')
