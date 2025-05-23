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

def run():
    stations = pd.read_csv(cnst.DATA+'/nflics/stations/ANACIM/Position_MOD.csv', index_col='station number')

    geoloc = xr.open_dataset(cnst.DATA+'/nflics/geoloc/nxny1640_580_nxnyds164580_blobdx0.04491576_area4_n23_20_32.nc')
    columns = stations.index.values

    lls = [geoloc.lons_mid.min(), geoloc.lons_mid.max(), geoloc.lats_mid.min(), geoloc.lats_mid.max()]

    data_resolution = 3 # in km
    # make salem grid
    grid3k = u_grid.make(np.arange(lls[0],lls[1]), np.arange(lls[2],lls[3]), data_resolution*1000)
    dlon = geoloc.lons_mid.values
    dlat = geoloc.lats_mid.values
    inds3, weights3, shape3 = u_int.interpolation_weights_grid(dlon, dlat, grid3k)

    lonN, latN = grid3k.ll_coordinates
    ds = xr.Dataset()

    for yy in range(2004,2023):
         for mm in range(1,12):
    #for yy in range(2010, 2011):
    #    for mm in range(8, 9):
            print('Doing ',yy, mm)

            filepath = sorted(glob.glob(cnst.other_drive+'/nflics/hist_cores/'+str(yy)+'/'+str(mm).zfill(2)+'/*/*.nc'))

            struc_tir = {}
            struc_wav = {}
            for col in columns:
                struc_tir[col] = []
                struc_wav[col] = []
            time_index = []

            outpath_tir = cnst.network_data+'data/nflics/core_gauge_comparison_ANACIM/stations/extract/tir_stationLoc_' + str(yy) + str(mm).zfill(2) +'.csv'
            outpath_core = cnst.network_data+'data/nflics/core_gauge_comparison_ANACIM/stations/extract/core_stationLoc_' +  str(yy) + str(mm).zfill(2) +'.csv'

            if len(glob.glob(cnst.network_data + 'data/nflics/core_gauge_comparison_ANACIM/stations/extract/*' + str(yy) + str(mm).zfill(2) + "*")) != 0:
                print('Files exist, continue')
                continue

            for ff in filepath:

                base = os.path.basename(ff)

                y = base[-15:-11]
                m = base[-11:-9]
                d = base[-9:-7]
                h = base[-7:-5]
                mi = base[-5:-3]

                try:
                    dat = xr.open_dataset(ff)
                except:
                    print('Could not read ',ff,' continue!')
                    continue
                print('Doing ',ff)

                time_index.append(datetime.datetime(int(y), int(m), int(d), int(h), int(mi)))

                tir = u_int.interpolate_data(dat.msg_Tir.values, inds3, weights3, shape3)
                cores = u_int.interpolate_data(dat.msg_cores.values, inds3, weights3, shape3)

                tir_ds = xr.DataArray(tir, dims=("lat", "lon",), coords=[latN[:,0], lonN[0,:]])
                cores_ds = xr.DataArray(cores, dims=("lat", "lon",), coords=[latN[:,0], lonN[0,:]])

                ds['tir'] = tir_ds
                ds['cores'] = cores_ds


                for kk in struc_tir.keys():

                    slon = stations.loc[kk].longitude
                    slat = stations.loc[kk].latitude

                    dat = ds.sel(lon=slon, lat=slat, method='nearest')

                    point = int(np.array(dat.cores.values > 0).astype(int))
                    temp = dat.tir.values/10000

                    struc_tir[kk].append(temp)
                    struc_wav[kk].append(point)

            if len(struc_tir[1]) > 0:

                df = pd.DataFrame(struc_tir, index=time_index)
                df.to_csv(outpath_tir)

                df2 = pd.DataFrame(struc_wav, index=time_index)
                df2.to_csv(outpath_core)


