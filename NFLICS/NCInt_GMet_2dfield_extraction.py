import numpy as np
import xarray as xr
import glob
from wavelet import util
import os
from utils import constants as cnst, u_arrays as ua, u_grid, u_interpolate as u_int
import pandas as pd
import ipdb

import matplotlib.pyplot as plt


#stations = pd.read_csv(cnst.ext_drive+'/nflics/stations/GMet_stations_lon_lat.csv', index_col='SHORT_STN_NAME')
dates = ['20200909', '20201010', '20200924', '20210819', '20220521', '20220524', '20220618', '20220621', '20220610', '20220705']

box = [-3.6, 1.7, 4, 11.4]
geoloc = xr.open_dataset(cnst.other_drive+'/nflics/geoloc/nxny1640_580_nxnyds164580_blobdx0.04491576_area4_n23_20_32.nc')
lls = [geoloc.lons_mid.min()+4, geoloc.lons_mid.max()-4, geoloc.lats_mid.min(), geoloc.lats_mid.max()-4]

data_resolution = 3  # in km
# make salem grid
grid3k = u_grid.make(np.arange(lls[0], lls[1]), np.arange(lls[2], lls[3]), data_resolution * 1000)
dlon = geoloc.lons_mid.values
dlat = geoloc.lats_mid.values
inds3, weights3, shape3 = u_int.interpolation_weights_grid(dlon, dlat, grid3k)
lonN, latN = grid3k.ll_coordinates

for yy in dates:
    filepath = glob.glob(cnst.other_drive+'/nflics/hist_cores/'+yy[0:4]+'/'+yy[4:6]+'/'+yy[6:8]+'/Hist_cores_wa_'+yy+'*.nc')

    for ff in filepath:
        ds = xr.Dataset()
        base = os.path.basename(ff)

        y = base[-15:-11]
        m = base[-11:-9]
        d = base[-9:-7]
        h = base[-7:-5]
        mi = base[-5:-3]

        outbase = base.replace('Hist_cores', 'Hist_cores_Ghana')
        outf = cnst.other_drive+'/nflics/hist_cores_Ghana/'+outbase
        if os.path.isfile(outf):
            print('Files exist, continue')
            continue

        try:
            dat = xr.open_dataset(ff)
        except:
            print('Could not read ',ff,' continue!')
            continue
        print('Doing ',ff)

        tir = u_int.interpolate_data(dat.msg_Tir.values, inds3, weights3, shape3)
        cores = u_int.interpolate_data(dat.msg_cores.values, inds3, weights3, shape3)
        tir_ds = xr.DataArray(tir.astype(np.int32), dims=("lat", "lon",), coords=[latN[:, 0], lonN[0, :]])
        cores_ds = xr.DataArray(cores, dims=("lat", "lon",), coords=[latN[:, 0], lonN[0, :]])

        ds['tir'] = tir_ds
        ds['cores'] = cores_ds

        ds = ds.sel(lon=slice(box[0], box[1]), lat=slice(box[2], box[3]))
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}

        ds.to_netcdf(outf, mode='w', encoding=encoding, format='NETCDF4')
        del ds


