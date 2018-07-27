
import numpy as np
from scipy.ndimage.measurements import label
import xarray as xr
import os
import ipdb
import matplotlib.pyplot as plt
import glob
from utils import constants
import pandas as pd
import pdb

def file_save(p_path, olr_path):
    p_arr = xr.open_dataset(p_path)
    olr_arr = xr.open_dataset(olr_path)
    landsea = xr.open_dataset(constants.CP4_LANDSEA, decode_times=False)
    out = '/users/global/cornkle/data/CP4/CLOVER/MCS/'

    p_arr = p_arr['lsRain'].sel(longitude=slice(-13,13), latitude=slice(4,8))
    olr_arr = olr_arr['lw_out_PBLtop'].sel(longitude=slice(-13, 13), latitude=slice(4, 8))
    ls_arr = landsea['lsm'].sel(rlon=slice(-13+360, 13+360), rlat=slice(4, 8))

    pos = np.where(ls_arr[0,0,:,:]==0)
    p_arr.values[:,pos[0], pos[1]] = 0
    olr_arr.values[:,pos[0], pos[1]] = 0

    lons, lats = np.meshgrid(p_arr.longitude.values, p_arr.latitude.values)

    dates = p_arr.time.values

    cnt=0

    for date in dates:

        if (pd.Timestamp(date).hour != 18):
            continue

        olr = olr_arr[olr_arr.time == date].squeeze()
        precip = p_arr[p_arr.time == date].squeeze()

        olr.values[olr.values >= 167] = 0  # T threshold -40 maskout
        labels, numL = label(olr.values)
        #
        # plt.figure()
        # plt.imshow(labels)

        u, inv = np.unique(labels, return_inverse=True)
        n = np.bincount(inv)

        goodinds = u[n > 258]  # defines minimum MCS size e.g. 350 km2 = 39 pix at 3x3km res (258 pix at 4.4km is 5000km2)
        print(goodinds)
        if not sum(goodinds) > 0:
            return

        for gi in goodinds:
            if gi == 0:  # index 0 is always background, ignore!
                continue

            inds = np.where(labels == gi)
            mask = np.where(labels != gi)

            pdummy = precip.copy(deep=True)
            odummy = olr.copy(deep=True)

            pdummy.values[mask] = np.nan
            odummy.values[mask] = np.nan

            # cut a box for every single blob from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
            latmax, latmin = lats[inds].max(), lats[inds].min()
            lonmax, lonmin = lons[inds].max(), lons[inds].min()

            olr_box = odummy.sel(latitude=slice(latmin,latmax), longitude=slice(lonmin, lonmax))
            precip_box = pdummy.sel(latitude=slice(latmin, latmax), longitude=slice(lonmin, lonmax))


            olr_box = olr_box.to_dataset()
            olr_box['lsRain'] = precip_box

            savefile = out + pd.Timestamp(date).strftime('%Y-%m-%d_%H:%M:%S') + '_' + str(gi) + '.nc'
            try:
                os.remove(savefile)
            except OSError:
                pass

            olr_box.to_netcdf(path=savefile, mode='w')
            print('Saved ' + savefile)

            cnt = cnt + 1

        print('Saved ' + str(cnt) + ' MCSs as netcdf.')



def file_loop():

    path = constants.CP4_PATH
    rain = constants.CP4_RAIN
    olr = constants.CP4_OLR

    files = glob.glob(path+'CLOVER/CP4hist/lsRain/*.nc')

    for f in files:

        fname = os.path.basename(f)
        split = fname.split('_')[-1]
        datestring = split[-15:-7]
        year = split[-15:-11]
        # if (year != '2002'):
        #     continue

        vars = ['lw_out_PBLtop']
        vlist = []
        for v in vars:
            f2 = glob.glob(path+'CLOVER/CP4hist/'+v+'/*'+datestring+'*.nc')[0]
            vlist.append(f2)

            file_save(f, vlist[0])

    #file_save(rain, olr)