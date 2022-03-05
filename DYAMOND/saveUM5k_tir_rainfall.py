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
from scipy import ndimage
import pdb
import multiprocessing

HOD = range(24)  # hours of day
#YRANGE = range(2004, 2015)

def olr_to_bt(olr):
    sigma = 5.670373e-8
    return ((olr/sigma)**0.25)-273.15


def saveMCS():

    ffiles = '/home/ck/DIR/cornkle/data/DYAMOND/UM-5km/WA/'
    vars = ['rlut', 'pr']

    metum_res = 5

    fnames = glob.glob(ffiles+'/rlut/*.nc')

    cnt = 0
    for f in fnames:

        orl = (xr.open_dataset(f)['toa_outgoing_longwave_flux']).load()
        vals = orl.values
        nvals = olr_to_bt(vals)
        orl.values = nvals


        bname = os.path.basename(f)
        strdate = bname[-11:-3]

        pfile = glob.glob(ffiles+'/pr/*'+strdate+'.nc')
        pcp = (xr.open_dataset(pfile[0])['precipitation_flux'])*3600
        pfile = glob.glob(ffiles + '/prw/*' + strdate + '.nc')
        tcw = (xr.open_dataset(pfile[0])['atmosphere_water_vapor_content'])
        w = (tcw[tcw['time.hour'] == 10]).squeeze()
        lon, lat = np.meshgrid(orl.longitude, orl.latitude)

        for p,t in zip(pcp, orl):


            # ipdb.set_trace()
            #
            # _y = t['time.year'].values[0]
            # _m = t['time.month'].values[0]
            _h = int(t['time.hour'].values)
            #
            # _d = t['time.day'].values[0]
            _mi = int(t['time.minute'].values)
            #

            # date = dt.datetime(_y, _m, _d, _h, 0)
            t.values = ndimage.gaussian_filter(t.values, 2, mode='nearest')
            #ipdb.set_trace()
            t.values[t.values >= -40] = 0  # T threshold -10
            labels, numL = label(t.values)

            u, inv = np.unique(labels, return_inverse=True)
            n = np.bincount(inv)

            goodinds = u[n >= 200]  # defines minimum MCS size e.g. 350 km2 = 39 pix at 3x3km res , 5000km2 = 556 pixel
            print(goodinds)
            if not sum(goodinds) > 0:
                continue

            for gi in goodinds:
                if gi == 0:  # index 0 is always background, ignore!
                    continue

                inds = np.where(labels == gi)

                # plt.pcolormesh(labels, cmap='jet')
                # #plt.contour(t)
                # plt.colorbar()
                #
                # plt.show()

                #print(p.shape, t.shape)
                # cut a box for every single blob from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
                try:
                    latmax, latmin = np.max(lat[inds]), np.min(lat[inds])
                except:
                    ipdb.set_trace()
                lonmax, lonmin = np.max(lon[inds]), np.min(lon[inds])

                tout = t.sel(latitude=slice(latmin,latmax), longitude=slice(lonmin,lonmax))
                pout = p.sel(latitude=slice(latmin,latmax), longitude=slice(lonmin,lonmax))
                wout = w.sel(latitude=slice(latmin, latmax), longitude=slice(lonmin, lonmax))


                ds = xr.Dataset()
                ds['tir'] = tout
                ds['tcw'] = wout
                ds['prcp'] = pout

                #ipdb.set_trace()


                #
                # if (sum(mmask.flatten())*25 < 350) | (outt.max()>250):# or (sum(mmask.flatten())*25 > 1500000): #or (outt.max()<0.1)
                #     continue
                #
                # if sum(mask2.flatten()) < 5:  # sum(mmask.flatten())*0.3:
                #     print('Kickout: TRMM MSG overlap less than 3pix of cloud area')
                #     continue
                #
                # print('Hit:', gi)
                #
                # da = xr.Dataset({'p': (['x', 'y'], outt),
                #                  't_lag0': (['x', 'y'], dummy),
                #                  'tc_lag0': (['x', 'y'], outl),
                #                  },
                #                 coords={'lon': (['x', 'y'], lon),
                #                         'lat': (['x', 'y'], lat),
                #                         'time': date})

                ds.attrs['meanT'] = np.mean(t.values[inds])
                ds.attrs['minT'] = np.min(t.values[inds])
                ds.attrs['meanP'] = np.mean(p.values[inds])
                ds.attrs['maxP'] = np.max(p.values[inds])
                ds.attrs['area'] = inds[0].size

                savefile = '/home/ck/DIR/cornkle/data/DYAMOND/UM-5km/MCS/WA/UM-5km/' + strdate + '_'+str(_h).zfill(2)+str(_mi).zfill(2)+'_' + str(gi) + '.nc'
                # try:
                #     os.remove(savefile)
                # except OSError:
                #     print('OSError, no dir?')
                #     pass
                print('Hour', ds['time.hour'])
                ds.to_netcdf(path=savefile, mode='w')
                print('Saved ' + savefile)


        print('Saved ' + str(cnt) + ' MCSs as netcdf.')
