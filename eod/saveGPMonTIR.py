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
    msg_folder2 = cnst.network_data + 'data/OBS/MSG_MAMON'

    mJJAS = msg.ReadMsg(msg_folder)
    mMAMON = msg.ReadMsg(msg_folder2)
    cnt = 0
    _y = year

    for ids, _m in enumerate(np.arange(3,12)):

        files = glob.glob(trmm_folder + '/'+str(_y) + '/'+str(_m).zfill(2) +'/*.nc4')# area=[-12, 12, 4, 9])   # [-15, 15, 4, 21], [-10, 10, 10, 20]

        for tf in files:
            

            t = xr.open_dataset(tf)
            
            _h = t['time.hour'].values[0]

            _d = t['time.day'].values[0]
            _mi = t['time.minute'].values[0]

            # if (_h <15) | (_h>21):
            #     print('Wrong hour')
            #     continue

            if (_m<3) | (_m>11):
                print('Wrong month')
                continue

            da = t['HQprecipitation'].squeeze()
            da = da.T
            da = da.sel(lat=slice(3.5,17))   #[-12, 15, 5, 25]

            if np.sum(da.values) <= 0.01:
                continue

            if _m in [3,4,5,10,11]:
                m = mMAMON
            else:
                m = mJJAS

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
            m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)

            mdic = m.get_data(llbox=[da['lon'].values.min(), da['lon'].values.max(), da['lat'].values.min(),
                                     da['lat'].values.max()])
            #mdic = m.get_data() 
            # check whether date is completely missing or just 30mins interval exists
            # if str(date) == '2004-05-02 13:15:00':
            #     pdb.set_trace()
            if not mdic:
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
                m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
                mdic = m.get_data(llbox=[da['lon'].values.min(), da['lon'].values.max(), da['lat'].values.min(),
                                         da['lat'].values.max()])

                if not mdic:
                    #ipdb.set_trace()
                    print('Date missing')
                    continue

            print('TRMM:', date, 'MSG:', ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute )

            lon1 = mdic['lon'].values
            lat1 = mdic['lat'].values

            if ids == 0:

                inds, weights, shape = u_int.interpolation_weights_grid(lon1, lat1, da.salem.grid)

            orig = mdic['t'].values

            try:
                outorig = u_int.interpolate_data(orig, inds, weights, shape)
            except IndexError:
                print('Interpolation problem, continue')


            outt = outorig.copy()

            labels, goodinds = ua.blob_define(outt, -40, minmax_area=[278, 40000],
                                              max_area=None)  # 7.7x7.7km = 64km2 per pix in gridsat?
            outt[labels == 0] = 0
            outt[outt < -115] = 0
            outt = (np.round(outt, decimals=2) * 100).astype(np.int16)

            outorig = (np.round(outorig, decimals=2) * 100).astype(np.int16)

            rain = (np.round(da.values, decimals=2) * 100).astype(np.int16)


            #ipdb.set_trace() 
            ds = xr.Dataset({'rain': (['lat', 'lon'], rain),
                             'mcs': (['lat', 'lon'], outt),
                             'tir': (['lat', 'lon'], outorig),
                             },
                            coords={'lon': da.lon,
                                    'lat': da.lat,
                                    'time': da.time})
            #ipdb.set_trace()
            # try:
            #     blob = xr.DataArray(power_msg[np.newaxis, :],
            #                         coords={'time': date, 'lat': latitudes, 'lon': longitudes},
            #                         dims=['time', 'lat', 'lon'])  # [np.newaxis, :])
            # except ValueError:
            #     ipdb.set_trace()
            # tir = xr.DataArray(new_savet[np.newaxis, :], coords={'time': date, 'lat': latitudes, 'lon': longitudes},
            #                    dims=['time', 'lat', 'lon'])

            try:
                ds = xr.concat([ds, df ], dim='time')
            except:
                ds = df.copy()


        savefile = cnst.network_data + 'MCSfiles/TIR_on_GPM/GPM_MCS_'+str(_y)+'-'+str(_m).zfill(2)+'.nc'
        try:
            os.remove(savefile)
        except OSError:
            print('OSError, no dir?')
            pass

        ds.to_netcdf(path=savefile, mode='w')
        print('Saved ' + savefile)
        ds.close()
        del ds

