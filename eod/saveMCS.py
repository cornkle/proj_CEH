# -*- coding: utf-8 -*-


import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage.measurements import label
import datetime as dt
from eod import msg, trmm, tm_utils
import xarray as xr
import os
import ipdb
import matplotlib.pyplot as plt
from utils import u_grid

HOD = range(24)  # hours of day
YRANGE = range(2004, 2014)


def saveMCS_WA15():
    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_WA30'

    # make a salem grid
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')

    t = trmm.ReadWA(trmm_folder, yrange=YRANGE, area=[-15, 4, 20, 25])   # [-15, 15, 4, 21], [-10, 10, 10, 20]
    m = msg.ReadMsg(msg_folder)

    cnt = 0

    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered      
    for _y, _m, _d, _h, _mi in zip(t.dates.y, t.dates.m, t.dates.d, t.dates.h, t.dates.mi):

        tdic = t.get_ddata(_y, _m, _d, _h, _mi, cut=[3, 26])

        # define the "0 lag" frist
        arr = np.array([15, 30, 45, 60, 0])
        dm = arr - _mi
        ind = (np.abs(dm)).argmin()

        # set zero shift time for msg
        date = dt.datetime(_y, _m, _d, _h, _mi)

        dt0 = dm[ind]
        ndate = date + dt.timedelta(minutes=int(dt0))
        m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
        mdic = m.get_data(llbox=[tdic['lon'].values.min(),  tdic['lat'].values.min(), tdic['lon'].values.max(),tdic['lat'].values.max()])

        # check whether date is completely missing or just 30mins interval exists
        if not mdic:
            dm = np.delete(dm, np.argmin(np.abs(dm)), axis=0)
            if np.min(np.abs(dm))> 15:
                print('Date missing')
                continue
            ind = (np.abs(dm)).argmin()
            dt0 = dm[ind]
            ndate = date + dt.timedelta(minutes=int(dt0))
            m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
            mdic = m.get_data(llbox=[tdic['lon'].values.min(), tdic['lat'].values.min(), tdic['lon'].values.max(),
                                     tdic['lat'].values.max()])

            if not mdic:
                print('Date missing')
                continue

        # if (ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute) != (2006, 6, 6, 5, 0):
        #     continue

        print('TRMM:', date, 'MSG:', ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute )

        lon1 = mdic['lon'].values
        lat1 = mdic['lat'].values
        mdic['t'].values[mdic['t'].values >= -40] = 0  # T threshold -10
        labels, numL = label(mdic['t'].values)

        u, inv = np.unique(labels, return_inverse=True)
        n = np.bincount(inv)

        goodinds = u[n > 36]  # all blobs with more than 36 pixels = 18 km x*y = 324 km2 (meteosat ca. 3km)
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
            td = t.get_ddata(_y, _m, _d, _h, _mi, cut=[latmin - 1, latmax + 1])

            # ensure minimum trmm rainfall in area
            if len(np.where(td['p'].values > 0.1)[0]) < 1:  # at least 1 pixel with rainfall
                print('Kickout: TRMM min pixel < 1')
                continue

            dt0 = dm[ind]
            ndate = date + dt.timedelta(minutes=int(dt0))

            # if (ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute) == (2006, 6, 6, 5, 0):
            #     ipdb.set_trace()

            ml0 = m.get_data(llbox=[lonmin - 1, latmin - 1, lonmax + 1, latmax + 1])
            if not ml0:
                continue

            #make salem grid
            grid = u_grid.make(ml0['lon'].values, ml0['lat'].values,proj,5000)
            lon, lat = grid.ll_coordinates

            # interpolate TRM and MSG to salem grid
            inter, mpoints = u_grid.griddata_input(ml0['lon'].values, ml0['lat'].values,grid)
            inter, tpoints = u_grid.griddata_input(td['lon'].values, td['lat'].values, grid)

            # Interpolate TRMM using delaunay triangularization
            dummyt = griddata(tpoints, td['p'].values.flatten(), inter, method='linear')
            outt = dummyt.reshape((grid.ny, grid.nx))
            if len(np.where(outt > 0.1)[0]) < 2:  # at least 2 pixel with rainfall
                print('Kickout: TRMM wavelet min pixel pcp < 2')
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
            #
            # f = plt.figure()
            # ax = f.add_subplot(1, 3, 1)
            # plt.imshow(outt, cmap='jet')
            #
            # ax = f.add_subplot(1, 3, 2)
            # plt.imshow(outff, cmap='jet')
            #
            # ax = f.add_subplot(1, 3, 3)
            # plt.imshow(outk, cmap='jet')
            #
            # plt.show()
            #
            # return

            if cnt > 5:
               return

            # Interpolate MSG using delaunay triangularization
            dummy = griddata(mpoints, ml0['t'].values.flatten(), inter, method='linear')
            dummy = dummy.reshape((grid.ny, grid.nx))
            outl = np.full_like(dummy, np.nan)
            xl, yl = grid.transform(lon1[inds], lat1[inds], crs=salem.wgs84, nearest=True, maskout=True)
            outl[yl.compressed(), xl.compressed()] = dummy[yl.compressed(), xl.compressed()]

            tmask = np.isfinite(outt)
            mmask = np.isfinite(outl)
            mask2 = np.isfinite(outl[tmask])

            if (sum(mmask.flatten())*25 < 15000) or (outt.max()<0.1 or (outt.max()>200) or (sum(mmask.flatten())*25 > 1500000)):
                continue

            if sum(mask2.flatten()) < 3:  # sum(mmask.flatten())*0.3:
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
            #da.attrs['lag1'] = dt1
            da.attrs['meanT'] = np.mean(outl[mmask])
            da.attrs['T90perc'] = mmeans
            da.attrs['meanT_cut'] = np.mean(outl[tmask][mask2])
            da.attrs['area'] = sum(mmask.flatten())
            da.attrs['area_cut'] = sum(mask2)
            da.close()
            savefile = '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size/' + date.strftime('%Y-%m-%d_%H:%M:%S') + '_' + str(gi) + '.nc'
            try:
                os.remove(savefile)
            except OSError:
                pass
            da.to_netcdf(path=savefile, mode='w')
            print('Saved ' + savefile)

            cnt = cnt + 1

    print('Saved ' + str(cnt) + ' MCSs as netcdf.')
