import numpy as np
import xarray as xr
import glob
from wavelet import util
import os
from utils import constants as cnst, u_arrays as ua
import pandas as pd
import ipdb

import matplotlib.pyplot as plt


stations = pd.read_csv(cnst.ext_drive+'/nflics/stations/GMet_stations_lon_lat.csv', index_col='SHORT_STN_NAME')

geoloc = xr.open_dataset(cnst.ext_drive+'/nflics/geoloc/nxny1640_580_nxnyds164580_blobdx0.04491576_area4_n23_20_32.nc')

time_index = []
columns = stations.index.values

for yy in range(2004,2022):
    for mm in range(1,12):

        filepath = glob.glob(cnst.ext_drive+'/nflics/hist_cores/'+str(yy)+'/'+str(mm).zfill(2)+'/*/*.nc')

        struc_tir = {}
        struc_wav = {}
        for col in columns:
            struc_tir[col] = []
            struc_wav[col] = []

        for ff in filepath:

            base = os.path.basename(ff)

            y = base[-15:-11]
            m = base[-11:-9]
            d = base[-9:-7]
            h = base[-7:-5]
            mi = base[-5:-3]

            if len(glob.glob(cnst.ext_drive+'/nflics/stations/*'+str(y)+str(m).zfill(2)+"*")) != 0:
                print('Files exist, continue')
                continue

            outpath = stations.replace('GMet_stations_lon_lat', 'core_tir_stationLoc_'+str(y)+str(m).zfill(2))

            try:
                dat = xr.open_dataset(ff)
            except:
                print('Could not read ',ff,' continue!')
                continue
            print('Doing ',ff)

            for kk in struc_tir.keys():

                slon = stations.loc[kk].LON
                slat = stations.loc[kk].LAT

                inds = ua.closest_point([slon,slat], np.array(list(zip(geoloc.lons_mid.values.flat,geoloc.lats_mid.values.flat))))

                uinds = np.unravel_index(inds, geoloc.lats_mid.values.shape)

                msg_lon = geoloc.lons_mid.values[uinds]
                msg_lat = geoloc.lats_mid.values[uinds]

                point = int(np.array(dat.msg_cores.values[uinds] > 0).astype(int))
                temp = dat.msg_Tir.values[uinds]/10000

                struc_tir[kk].append(temp)
                struc_wav[kk].append(point)

                print(temp, point)
