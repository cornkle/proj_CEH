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
from utils import u_grid
import pdb

HOD = range(24)  # hours of day
YRANGE = range(2004, 2015)


def saveMCS_WA15():
    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_tropWA' #meteosat_WA30'

    t = trmm_clover.ReadWA(trmm_folder, yrange=YRANGE, area=[-14, 12, 4, 8])   # [-15, 15, 4, 21], [-10, 10, 10, 20]
    m = msg.ReadMsg(msg_folder)

    cnt = 0

    # define the "0 lag" frist
    arr = np.array([15, 30, 45, 60, 0])
    #mon = [3,4,5] # months march april may only
    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered      
    for _y, _m, _d, _h, _mi in zip(t.dates.dt.year,  t.dates.dt.month, t.dates.dt.day, t.dates.dt.hour, t.dates.dt.minute):

        if (_h <10) | (_h>19):
            continue

        date = dt.datetime(_y, _m, _d, _h, _mi)

        tdic = t.get_ddata(date, cut=[4, 8])


        #get closest minute
        dm = arr - _mi
        dm = dm[dm<0]
        try:
            ind = (np.abs(dm)).argmin()
        except ValueError:
            continue

        # set zero shift time for msg


        dt0 = dm[ind]
        ndate = date + dt.timedelta(minutes=int(dt0))
        m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)

        mdic = m.get_data(llbox=[tdic['lon'].values.min(),  tdic['lon'].values.max(), tdic['lat'].values.min(),tdic['lat'].values.max()])

        # check whether date is completely missing or just 30mins interval exists
        # if str(date) == '2004-05-02 13:15:00':
        #     pdb.set_trace()
        if not mdic:
            dm = np.delete(dm, np.argmin(np.abs(dm)), axis=0)
            try:
                dummy = np.min(np.abs(dm))> 15
            except ValueError:
                continue
            if dummy:
                print('Date missing')
                continue
            ind = (np.abs(dm)).argmin()
            dt0 = dm[ind]
            ndate = date + dt.timedelta(minutes=int(dt0))
            m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
            mdic = m.get_data(llbox=[tdic['lon'].values.min(),  tdic['lon'].values.max(), tdic['lat'].values.min(),tdic['lat'].values.max()])

            if not mdic:
                print('Date missing')
                continue

        print('TRMM:', date, 'MSG:', ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute )

        lon1 = mdic['lon'].values
        lat1 = mdic['lat'].values
        mdic['t'].values[mdic['t'].values >= -40] = 0  # T threshold -10
        labels, numL = label(mdic['t'].values)

        u, inv = np.unique(labels, return_inverse=True)
        n = np.bincount(inv)

        goodinds = u[n > 39]  # defines minimum MCS size e.g. 350 km2 = 39 pix at 3x3km res
        print(goodinds)
        if not sum(goodinds) > 0:
            continue

        for gi in goodinds:
            if gi == 0:  # index 0 is always background, ignore!
                continue

            inds = np.where(labels == gi)

            # cut a box for every single blob from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
            latmax, latmin = mdic['lat'].values[inds].max(), mdic['lat'].values[inds].min()
            lonmax, lonmin = mdic['lon'].values[inds].max(), mdic['lon'].values[inds].min()
            mmeans = np.percentile(mdic['t'].values[inds], 90)
            td = t.get_ddata(date, cut=[latmin - 1, latmax + 1])

            # ensure minimum trmm rainfall in area
            # if len(np.where(td['p'].values > 0.1)[0]) < 1:  # at least 1 pixel with rainfall
            #     print('Kickout: TRMM min pixel < 1')
            #     continue

            dt0 = dm[ind]
            ndate = date + dt.timedelta(minutes=int(dt0))

            # if (ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute) == (2006, 6, 6, 5, 0):
            #     ipdb.set_trace()

            ml0 = m.get_data(llbox=[lonmin - 1,  lonmax + 1, latmin - 1, latmax + 1])
            if not ml0:
                continue

            #make salem grid
            grid = u_grid.make(ml0['lon'].values, ml0['lat'].values,5000)
            lon, lat = grid.ll_coordinates

            # interpolate TRM and MSG to salem grid
            inter, mpoints = u_grid.griddata_input(ml0['lon'].values, ml0['lat'].values,grid)
            inter, tpoints = u_grid.griddata_input(td['lon'].values, td['lat'].values, grid)

            # Interpolate TRMM using delaunay triangularization
            try:
                dummyt = griddata(tpoints, td['p'].values.flatten(), inter, method='linear')
            except ValueError:
                continue
            outt = dummyt.reshape((grid.ny, grid.nx))
            # if len(np.where(outt > 0.1)[0]) < 2:  # at least 2 pixel with rainfall
            #     print('Kickout: TRMM wavelet min pixel pcp < 2')
            #     continue

            if np.sum(np.isfinite(outt)) < 5:  # at least 2 valid pixel
                print('Kickout: TRMM wavelet min pixel  < 2')
                continue

            # Interpolate TRMM flags using nearest
            dummyf = griddata(tpoints, td['flags'].values.flatten(), inter, method='nearest')
            outf = dummyf.reshape((grid.ny, grid.nx))
            outf=outf.astype(np.float)
            isnot = np.isnan(outt)
            outf[isnot]=np.nan

            ##remove edges of interpolated TRMM
            for nb in range(5):
                boole = np.isnan(outt)
                outt[boole] = -1000
                grad = np.gradient(outt)
                outt[boole] = np.nan
                outt[abs(grad[1]) > 300] = np.nan
                outt[abs(grad[0]) > 300] = np.nan
                outf[abs(grad[1]) > 300] = np.nan
                outf[abs(grad[0]) > 300] = np.nan

            #get convective rainfall only
            outff = tm_utils.getTRMMconv(outf)
            outk = outt.copy()*0
            outk[np.where(outff)]=outt[np.where(outff)]


            # Interpolate MSG using delaunay triangularization
            dummy = griddata(mpoints, ml0['t'].values.flatten(), inter, method='linear')
            dummy = dummy.reshape((grid.ny, grid.nx))
            outl = np.full_like(dummy, np.nan)
            xl, yl = grid.transform(lon1[inds], lat1[inds], crs=salem.wgs84, nearest=True, maskout=True)
            outl[yl.compressed(), xl.compressed()] = dummy[yl.compressed(), xl.compressed()]

            # #### SHIFTING WITH RESPECT TO MIN T / MAX P - search for Pmax within 20km from Tmin, shift TRMM image
            #
            # tmin = np.argmin(outl)
            # pmax =
            #
            # dist =
            #

            tmask = np.isfinite(outt)
            mmask = np.isfinite(outl)
            mask2 = np.isfinite(outl[tmask])

            if (sum(mmask.flatten())*25 < 350) | (outt.max()>200):# or (sum(mmask.flatten())*25 > 1500000): #or (outt.max()<0.1)
                continue

            if sum(mask2.flatten()) < 5:  # sum(mmask.flatten())*0.3:
                print('Kickout: TRMM MSG overlap less than 3pix of cloud area')
                continue

            print('Hit:', gi)

            da = xr.Dataset({'p': (['x', 'y'], outt),
                             'pconv': (['x', 'y'], outk),
                             't_lag0': (['x', 'y'], dummy),
                             'tc_lag0': (['x', 'y'], outl),
                             },
                            coords={'lon': (['x', 'y'], lon),
                                    'lat': (['x', 'y'], lat),
                                    'time': date})
            da.attrs['lag0'] = dt0
            da.attrs['meanT'] = np.mean(outl[mmask])
            da.attrs['T90perc'] = mmeans
            da.attrs['meanT_cut'] = np.mean(outl[tmask][mask2])
            da.attrs['area'] = sum(mmask.flatten())
            da.attrs['area_cut'] = sum(mask2)
            da.close()
            savefile = '/users/global/cornkle/MCSfiles/WA350_4-8N_14W-10E_-40/' + date.strftime('%Y-%m-%d_%H:%M:%S') + '_' + str(gi) + '.nc'
            try:
                os.remove(savefile)
            except OSError:
                pass
            da.to_netcdf(path=savefile, mode='w')
            print('Saved ' + savefile)

            cnt = cnt + 1

    print('Saved ' + str(cnt) + ' MCSs as netcdf.')
