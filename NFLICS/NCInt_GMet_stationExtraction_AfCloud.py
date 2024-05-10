import numpy as np
import xarray as xr
import glob
from wavelet import util
import os
from utils import constants as cnst, u_arrays as ua, u_grid, u_interpolate as u_int
import pandas as pd
import ipdb
import datetime
import matplotlib.pyplot as plt


stations = pd.read_csv(cnst.DATA+'/nflics/core_gauge_comparison_GMet/GMet_stations_lon_lat.csv', index_col='SHORT_STN_NAME')
columns = stations.index.values

box = [-3.6, 1.7, 4, 11.4]
lls = [-4,2,3,13]
dummy = glob.glob('/prj/Africa_cloud/ch9_wavelet/2023/08/*.nc')[0]

infiles = '/prj/Africa_cloud/ch9_wavelet/'
geoloc_f = '/prj/Africa_cloud/geoloc/lat_lon_2268_2080.npz'

latlon = np.load(geoloc_f)
mlon = latlon['lon'][150:-150,150:-150]
mlat = latlon['lat'][150:-150,150:-150]


outfiles = '/scratch/cornkle/gmet_stations/tables/'
outfiles_2d = '/scratch/cornkle/gmet_stations/2d/'
geoloc = xr.open_dataset(dummy)

data_resolution = 3  # in km
# make salem grid
grid3k = u_grid.make(np.arange(lls[0], lls[1]), np.arange(lls[2], lls[3]), data_resolution * 1000)
dlon = mlon
dlat = mlat
lonN, latN = grid3k.ll_coordinates

inds3, weights3, shape3 = u_int.interpolation_weights(dlon, dlat, lonN, latN, irregular_1d=True)

for yy in range(2004,2025): #range(2018,2025):
   for mm in range(1,13):
    curpath = outfiles_2d+str(yy)+'/'+str(mm).zfill(2)
    if not os.path.exists(curpath):
        os.makedirs(curpath)

    struc_tir = {}
    struc_wav = {}
    for col in columns:
            struc_tir[col] = []
            struc_wav[col] = []
            time_index = []

            outpath_tir = outfiles+'/tir_stationLoc_' + str(yy) + str(mm).zfill(2) +'.csv'
            outpath_core = outfiles+'/core_stationLoc_' +  str(yy) + str(mm).zfill(2) +'.csv'

    if len(glob.glob(outfiles+'*' + str(yy) + str(mm).zfill(2) + "*")) != 0:
            print(glob.glob(outfiles+'*' + str(yy) + str(mm).zfill(2) + "*"),'Files exist, continue')
            continue
    filepath = glob.glob(infiles+str(yy)+'/'+str(mm).zfill(2)+'/*.nc')

    for ff in filepath:
        ds = xr.Dataset()
        base = os.path.basename(ff)

        outbase = base.replace('.nc', '_Ghana.nc')
        outf = curpath+'/'+outbase

        y = base[-15:-11]
        m = base[-11:-9]
        d = base[-9:-7]
        h = base[-7:-5]
        mi = base[-5:-3]

        try:
            dat = xr.open_dataset(ff).squeeze()
        except:
            print('Could not read ',ff,' continue!')
            continue
        print('Doing ',ff)

        time_index.append(datetime.datetime(int(y), int(m), int(d), int(h), int(mi)))

        tir = u_int.interpolate_data(dat.tir.values[150:-150,150:-150]
, inds3, weights3, shape3)
        cores = u_int.interpolate_data(dat.cores.values[150:-150,150:-150]
, inds3, weights3, shape3)
        tir_ds = xr.DataArray(tir.astype(np.int32), dims=("lat", "lon",), coords=[latN[:, 0], lonN[0, :]])
        cores_ds = xr.DataArray(cores, dims=("lat", "lon",), coords=[latN[:, 0], lonN[0, :]])

        ds['tir'] = tir_ds
        ds['cores'] = cores_ds

        for kk in struc_tir.keys():

                    slon = stations.loc[kk].LON
                    slat = stations.loc[kk].LAT

                    dat = ds.sel(lon=slon, lat=slat, method='nearest')

                    point = int(np.array(dat.cores.values > 0).astype(int))
                    temp = dat.tir.values

                    struc_tir[kk].append(temp)
                    struc_wav[kk].append(point)

        ds = ds.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3]))
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}

        ds.to_netcdf(outf, mode='w', encoding=encoding, format='NETCDF4')
        del ds

    if len(struc_tir['NAV']) > 0:

                df = pd.DataFrame(struc_tir, index=time_index)
                df.to_csv(outpath_tir)

                df2 = pd.DataFrame(struc_wav, index=time_index)
                df2.to_csv(outpath_core)
