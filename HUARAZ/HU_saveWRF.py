import xarray as xr
import geopandas as gpd
import pandas as pd
import numpy as np
import salem
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from netCDF4 import num2date
import multiprocessing
import glob
import os
from utils import constants as cnst
import salem


def run():

    in_path = '/media/ck/Elements/SouthAmerica/WRF/RAW_WRF/d02/*.nc'


    pool = multiprocessing.Pool(processes=2)
    files = glob.glob(in_path)
    print('start loop')

    res = pool.map(loop, files)




def loop(f):

    out_path = '/home/ck/rio_santa/'#'/media/ck/Elements/SouthAmerica/WRF/rio_santa/'
    fname = '/home/ck/DIR/cornkle/data/HUARAZ/shapes/riosan_sel_one.shp'
    sdf = salem.read_shapefile(fname)
    sdf = salem.transform_geopandas(sdf, to_crs=salem.wgs84)

    varlist = ['ALBEDO', 'LH', 'HFX', 'OLR', 'PBLH',
               'PSFC', 'Q2', 'QFX', 'PRCP', 'SNOWC',
               'LWDNB', 'SWDNB', 'SWUPB', 'LWUPB', 'T2C', 'TH2', 'TSK', 'U10', 'V10', 'RH2', 'col_int_QVAPOR',
               'col_int_QICE']

    print('Doing', f)
    fbase = os.path.basename(f)

    if os.path.isfile(out_path+fbase):
        print(out_path+fbase, 'File exists, next!')
        return
    ds = salem.open_wrf_dataset(f, decode_times=False)
    ds = ds[varlist]
    new_lat = ds['lat'].data[:, 0]
    new_lon = ds['lon'].data[0, :]
    south_north = ds['south_north']
    west_east = ds['west_east']

    ds = ds.rename({'south_north': 'latitude', 'west_east': 'longitude'})
    ds = ds.assign_coords({'latitude': new_lat, 'longitude': new_lon})
    ds = ds.assign_coords({'south_north': ('latitude', south_north), 'west_east': ('longitude', west_east)})

    ds = ds.sel(latitude=slice(-10.25, -8.6), longitude=slice(-78.02, -76.98))
    # ipdb.set_trace()
    dsout = ds.salem.roi(shape=sdf)
    #     new_time=[]
    #     for t in dsout.time:
    #         ipdb.set_trace()
    #         ot = pd.to_datetime(t.values, format="%Y-%m-%d_%H:%M:%S")
    #         new_time.append(ot)
    #     dsout = dsout.assign_coords({'time' : new_time})
    dsout = dsout.drop('xtime')
    # dsout = dsout.assign_coords({'time' : ttime.values})
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in dsout.data_vars}

    dsout.to_netcdf(out_path + fbase, mode='w', encoding=encoding,
                    format='NETCDF4')

    print('Saved', out_path+fbase)
    del ds
    del dsout








