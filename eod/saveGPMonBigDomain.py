# -*- coding: utf-8 -*-


import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage.measurements import label
import datetime as dt
from eod import msg, trmm, tm_utils, trmm_clover
import xarray as xr
import os
import ipdb
import matplotlib.pyplot as plt
import glob
from utils import u_grid, constants as cnst
import ipdb
import multiprocessing
from utils import u_arrays as ua
from utils import u_grid, u_interpolate as u_int

HOD = range(24)  # hours of day
#YRANGE = range(2004, 2015)


def multi():
    pool = multiprocessing.Pool(processes=5)

    res = pool.map(saveMCS_WA15, np.arange(2005,2010))



def saveMCS_WA15(year):
    trmm_folder = cnst.network_data + 'data/OBS/IMERG_HQ_precip'
    msg_folder = cnst.network_data + 'data/OBS/MSG_WA30' #meteosat_WA30'

    msg_folder = 'prj/vera/cores_bigDomain'

    _y = year

    for ids, _m in enumerate(np.arange(6,10)):

        files = glob.glob(trmm_folder + '/'+str(_y) + '/'+str(_m).zfill(2) +'/*.nc4')# area=[-12, 12, 4, 9])   # [-15, 15, 4, 21], [-10, 10, 10, 20]
        msg = xr.open_dataset(glob.glob(msg_folder + '*_' + str(_y) + '_' + str(_m).zfill(2) + '.nc')[0])
        #ds = xr.Dataset()

        for tf in files:
            

            t = xr.open_dataset(tf)
            
            _h = t['time.hour'].values[0]

            _d = t['time.day'].values[0]
            _mi = t['time.minute'].values[0]

            # if (_h <15) | (_h>21):
            #     print('Wrong hour')
            #     continue

            if (_m<6) | (_m>9):
                print('Wrong month')
                continue

            da = t['HQprecipitation'].squeeze()
            da = da.T
            da = da.sel(lat=slice(3.5,18))   #[-12, 15, 5, 25]

            if np.sum(da.values) <= 0.01:
                continue

            date = dt.datetime(_y, _m, _d, _h, _mi)
            arr = np.array([15, 30, 45, 60, 0])

            # get closest minute
            dm = arr - _mi
            if (dm < 0).any():
                dm = dm[dm < 0]
           
            try:
                ind = (np.abs(dm)).argmin()
            except ValueError:
                continue

            # set zero shift time for msg

            dt0 = dm[ind]
            ndate = date + dt.timedelta(minutes=int(dt0))

            try:
                mmsg = msg.sel(time=ndate).copy(deep=True)
                #msg = xr.open_dataset(glob.glob(msg_folder+'*_'+str(ndate.year)+'_'+str(ndate.month)+'.nc')[0])

            except:
                dm = np.delete(dm, np.argmin(np.abs(dm)), axis=0)
                try:
                    dummy = np.min(np.abs(dm)) > 15
                except ValueError:
                    print('Date missing')
                    continue
                if dummy:
                    print('Date missing')
                    continue
                ind = (np.abs(dm)).argmin()
                dt0 = dm[ind]
                ndate = date + dt.timedelta(minutes=int(dt0))

                try:
                    mmsg = msg.sel(time=ndate).copy(deep=True)
                except:
                    print('Date missing')
                    continue

            print('TRMM:', date, 'MSG:', mmsg.time )

            lon1 = da['lon'].values
            lat1 = da['lat'].values

            if ids == 0:

                inds, weights, shape = u_int.interpolation_weights_grid(lon1, lat1, mmsg.salem.grid)

            orig = da['HQprecipitation'].values

            try:
                outorig = u_int.interpolate_data(orig, inds, weights, shape)
            except IndexError:
                print('Interpolation problem, continue')


            outorig = (np.round(outorig, decimals=2) * 100).astype(np.int16)


            mmsg['HQprecip'] = xr.DataArray(outorig[np.newaxis, :],
                                     coords={'time': date, 'lat': lat1, 'lon': lon1},
                                     dims=['time', 'lat', 'lon'])  # [np.newaxis, :])

            try:
                ds = xr.concat([ds, mmsg ], dim='time')
            except:
                ds = mmsg.copy()


        savefile = cnst.network_data + 'MCSfiles/GPM_on_Cores/GPM_Cores_'+str(_y)+'-'+str(_m).zfill(2)+'.nc'
        try:
            os.remove(savefile)
        except OSError:
            print('OSError, no dir?')
            pass

        ds.to_netcdf(path=savefile, mode='w')
        print('Saved ' + savefile)
        ds.close()
        del ds

