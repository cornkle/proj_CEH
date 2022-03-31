# -*- coding: utf-8 -*-


import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pyhdf.SD import SD, SDC
import cartopy.crs as ccrs
import xarray as xr
import ipdb
import datetime
import glob
import os
import pandas as pd
import multiprocessing
from utils import constants as cnst



def run():
    files = glob.glob(cnst.lmcs_drive+ 'MODIS_LST/aqua/raw_v61/*.hdf')
    multi_list = []
    #ipdb.set_trace()
    for ff in files:
        multi_list.append((ff,['LST_Day','Day_view_angle','QC_Day'], True))
    pool = multiprocessing.Pool(processes=1)
    #ipdb.set_trace()
    res = pool.starmap(read_modis_monthly, multi_list)


def read_modis_monthly(FILE_NAME, dnames, save):  # dnames=None, save=False):

    def getMODISflag(flags):
        bla = flags.astype(int)
        mask = []

        for b, i in zip(np.nditer(bla), range(bla.size)):
            bb = '{0:008b}'.format(int(b))

            take = 0

            if len(bb) == 8:
                mandat = bb[0:2]
                quali = bb[2:4]
                emierr = bb[4:6]
                lsterr = bb[6:8]

            if (mandat in ['00']) | ((mandat in ['01']) & (lsterr in ['10', '11', '01']) & (quali in ['00','10'])): #, '01']) & (emierr in ['10', '11']) & (lsterr in ['10', '11', '01']):     #, '01']
                    take = 1
            # print(bb)
            mask.append(take)
        mask = np.reshape(np.array(mask), bla.shape)
        # print(np.sum(mask))
        return mask

    # Identify the data field.

    #     fields = ['CMG 0.05 Deg Monthly NDVI']
    #     DATAFIELD_NAME = 'CMG 0.05 Deg Monthly NDVI'
    print('Doing ', FILE_NAME)

    out = FILE_NAME.replace('.hdf', '.nc')
    out = out.replace('raw', 'nc')
    if (save) & os.path.isfile(out):
        print('File exists, return')
        return
    try:
        hdf = SD(FILE_NAME, SDC.READ)
    except:
        print('read error with file ', FILE_NAME, 'continue')
        return

    if dnames is None:
        dnames = hdf.datasets().keys()

    ds = xr.Dataset()

    # Normally we would use the grid metadata to reconstruct the grid, but
    # the grid metadata is incorrect in this case, specifically the upper left
    # and lower right coordinates of the grid.  We'll construct the grid
    # manually, taking into account the fact that we're going to subset the
    # data by a factor of 10 (the grid size is 3600 x 7200).
    lon = np.linspace(-180, 180, 7200)
    lat = np.linspace(90, -90, 3600)
    # lon, lat = np.meshgrid(x, y)

    split_date = FILE_NAME.split('.')[1]
    year = split_date[1:5]

    date = datetime.datetime(int(year), 1, 1) + pd.Timedelta(str(int(split_date[5:8]) - 1) + ' days')

    # ipdb.set_trace()

    for DATAFIELD_NAME in dnames:

        # Read dataset.
        data2D = hdf.select(DATAFIELD_NAME)
        data = data2D[:, :].astype(np.float)

        # Read attributes.
        attrs = data2D.attributes(full=1)
        lna = attrs["long_name"]
        long_name = lna[0]
        vra = attrs["valid_range"]
        valid_range = vra[0]
        try:
            aoa = attrs["_Offset"]
            add_offset = aoa[0]
        except:
            add_offset = 0
        fva = attrs["_FillValue"]
        _FillValue = fva[0]
        try:
            sfa = attrs["_Scale"]
            scale_factor = sfa[0]
        except:
            scale_factor = 1
        # ua=attrs["units"]
        # units = ua[0]

        # Handle fill value.
        if type(_FillValue) != str:
            invalid = (data == _FillValue)
            invalid = np.logical_or(invalid, data < valid_range[0])
            invalid = np.logical_or(invalid, data > valid_range[1])
            data[invalid] = np.nan

        # Apply scale factor and offset.
        data = data * scale_factor + add_offset
        if 'LST' in DATAFIELD_NAME:
            # ipdb.set_trace()
            data = ((data-273.15) * 100).astype(int)

        if ('angle' in DATAFIELD_NAME):
            # ipdb.set_trace()
            data = ((data <= 40) & (data >= -40)).astype(int)

        # ipdb.set_trace()
        # data = np.ma.masked_array(data, np.isnan(data))

        da = xr.DataArray(data[None, ...], coords={'time': np.array([date]),
                                                   'lat': lat,
                                                   'lon': lon},
                          dims=['time', 'lat', 'lon'])  # .isel(time=0)

        ds[DATAFIELD_NAME] = da

    #ql_mask = getMODISflag(ds['QC_Day'])

    ds['LST_Day'].values[(ds['Day_view_angle'].values != 1) ] = -9999   #| (ql_mask != 1)
    ds = ds.drop(['Day_view_angle'])
    ds = ds.drop(['QC_Day'])
    if save:
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}

        # ipdb.set_trace()
        ds.to_netcdf(path=out, mode='w', encoding=encoding, format='NETCDF4')
        return

    return ds

