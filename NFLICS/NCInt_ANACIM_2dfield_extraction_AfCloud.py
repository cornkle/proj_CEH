import numpy as np
import xarray as xr
import glob
from wavelet import util
import os
from utils import constants as cnst, u_arrays as ua, u_grid, u_interpolate as u_int
import pandas as pd
import ipdb

import matplotlib.pyplot as plt



#box = [-3.6, 1.7, 4, 11.4]
#lls = [-4,2,3,13]

box = [-18.8,-3,5,20.5]
lls = [-22.7,-2,4,21]

dummy = glob.glob('/prj/Africa_cloud/ch9_wavelet/2023/08/*.nc')[0]

infiles = '/prj/Africa_cloud/ch9_wavelet/'
geoloc_f = '/prj/Africa_cloud/geoloc/lat_lon_2268_2080.npz'

latlon = np.load(geoloc_f)
mlon = latlon['lon'][64:-64,57:-150]
mlat = latlon['lat'][64:-64,57:-150]
#outfiles = '/scratch/cornkle/ANACIM_stations/2d/'
outfiles = '/users/global/cornkle/shared/data/nflics/ANACIM_stations/2d/'
geoloc = xr.open_dataset(dummy)
#lls = [geoloc.lon.min()+4, geoloc.lons_mid.max()-4, geoloc.lats_mid.min(), geoloc.lats_mid.max()-4]

data_resolution = 3  # in km
# make salem grid
grid3k = u_grid.make(np.arange(lls[0], lls[1]), np.arange(lls[2], lls[3]), data_resolution * 1000)
dlon = mlon
dlat = mlat
lonN, latN = grid3k.ll_coordinates

inds3, weights3, shape3 = u_int.interpolation_weights(dlon, dlat, lonN, latN, irregular_1d=True)
for yy in range(2014,2024): #range(2018,2025):
   for mm in range(1,13):
    curpath = outfiles+str(yy)+'/'+str(mm).zfill(2)
    if not os.path.exists(curpath):
        os.makedirs(curpath)

    filepath = glob.glob(infiles+str(yy)+'/'+str(mm).zfill(2)+'/*.nc')

    for ff in filepath:
        
        if ff[-5:-3] != '30':
            print('Wrong hour', ff)
            continue
        ds = xr.Dataset()
        base = os.path.basename(ff)

        #y = base[-15:-11]
        #m = base[-11:-9]
        #d = base[-9:-7]
        #h = base[-7:-5]
        #mi = base[-5:-3]

        outbase = base.replace('.nc', '_Senegal.nc')
        outf = curpath+'/'+outbase
        if os.path.isfile(outf):
            print('Files exist, continue')
            continue

        try:
            dat = xr.open_dataset(ff).squeeze()
        except:
            print('Could not read ',ff,' continue!')
            continue
        print('Doing ',ff)
        tir = u_int.interpolate_data(dat.tir.values[64:-64,57:-150]
, inds3, weights3, shape3)
        cores = u_int.interpolate_data(dat.cores.values[64:-64,57:-150]
, inds3, weights3, shape3)
        tir_ds = xr.DataArray(np.round(tir,2).astype(np.int32), dims=("lat", "lon",), coords=[latN[:, 0], lonN[0, :]])
        cores_ds = xr.DataArray(cores, dims=("lat", "lon",), coords=[latN[:, 0], lonN[0, :]])

        ds['tir'] = tir_ds
        ds['cores'] = cores_ds

        ds = ds.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3]))
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}

        ds.to_netcdf(outf, mode='w', encoding=encoding, format='NETCDF4')
        del ds

