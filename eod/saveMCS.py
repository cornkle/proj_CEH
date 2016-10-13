# -*- coding: utf-8 -*-


import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage.measurements import label
import datetime as dt
from eod import msg, trmm
import xarray as xr
import os

HOD = range(24)  # hours of day
YRANGE = range(2004, 2015)


def saveMCS_WA15():
    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_SA15'

    # make a salem grid
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')

    t = trmm.ReadWA(trmm_folder, yrange=YRANGE, area=[-15, 15, 10, 20])
    m = msg.ReadMsg(msg_folder)

    cnt = 0

    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered      
    for _y, _m, _d, _h, _mi in zip(t.dates.y, t.dates.m, t.dates.d, t.dates.h, t.dates.mi):

        tdic = t.get_ddata(_y, _m, _d, _h, _mi, cut=[9, 21])

        # define the "0 lag" frist
        arr = np.array([15, 30, 45, 60, 0])
        dm = arr - _mi
        ind = (np.abs(dm)).argmin()

        # set zero shift time for msg
        date = dt.datetime(_y, _m, _d, _h, _mi)
        print(date)

        dt0 = dm[ind]
        ndate = date + dt.timedelta(minutes=int(dt0))
        m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
        mdic = m.get_data(llbox=[tdic['lon'].min(), tdic['lon'].max(), tdic['lat'].min(), tdic['lat'].max()])
        if not mdic:
            print('Date missing')
            continue
        lon1 = mdic['lon']
        lat1 = mdic['lat']
        mdic['t'][mdic['t'] > -40] = 0
        labels, numL = label(mdic['t'])

        u, inv = np.unique(labels, return_inverse=True)
        n = np.bincount(inv)

        goodinds = u[n > 2500]  # all blobs with more than 2500 pixels - size threshold
        print(goodinds)
        if not sum(goodinds) > 0:
            continue

        for gi in goodinds:
            if gi == 0:  # index 0 is always background, ignore!
                continue

            inds = np.where(labels == gi)

            # cut a box for every single blob from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
            latmax, latmin = mdic['lat'][inds].max(), mdic['lat'][inds].min()
            lonmax, lonmin = mdic['lon'][inds].max(), mdic['lon'][inds].min()
            mmeans = np.percentile(mdic['t'][inds], 90)
            td = t.get_ddata(_y, _m, _d, _h, _mi, cut=[latmin - 1, latmax + 1])

            # ensure minimum trmm rainfall in area
            if len(np.where(td['p'] > 0)[0]) < 100:  # at least 100 pixel with rainfall
                print('Kickout: TRMM min pixel = 100')
                continue

            dt0 = dm[ind]
            ndate = date + dt.timedelta(minutes=int(dt0))
            # print('Date1', ndate)
            ml0 = m.get_data(llbox=[lonmin - 1, lonmax + 1, latmin - 1, latmax + 1])
            if not ml0:
                continue

            dt1 = dm[ind] - 15
            ndate = date + dt.timedelta(minutes=int(dt1))
            #   print('Date2', ndate)
            m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
            ml1 = m.get_data(llbox=[lonmin - 1, lonmax + 1, latmin - 1, latmax + 1])
            if not ml1:
                continue

            dt2 = dm[ind] - 30
            ndate = date + dt.timedelta(minutes=int(dt2))
            #   print('Date2', ndate)
            m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
            ml2 = m.get_data(llbox=[lonmin - 1, lonmax + 1, latmin - 1, latmax + 1])
            if not ml1:
                continue

            dt3 = dm[ind] - 45
            ndate = date + dt.timedelta(minutes=int(dt3))
            #   print('Date2', ndate)
            m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
            ml3 = m.get_data(llbox=[lonmin - 0.3, lonmax + 0.3, latmin - 0.25, latmax + 0.25])
            if not ml1:
                continue

            dtx = dm[ind] + 45
            ndate = date + dt.timedelta(minutes=int(dtx))
            #   print('Date2', ndate)
            m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
            mlx = m.get_data(llbox=[lonmin - 0.3, lonmax + 0.3, latmin - 0.25, latmax + 0.25])
            if not ml1:
                continue

                # create grid surrounding the blob
            # Transform lon, lats to the mercator projection
            x, y = pyproj.transform(salem.wgs84, proj, ml0['lon'], ml0['lat'])
            # take the min and max
            xmax, xmin = np.max(x), np.min(x)
            ymax, ymin = np.max(y), np.min(y)
            # Count the number of pixels
            dx = 5000
            nx, r = divmod(xmax - xmin, dx)
            ny, r = divmod(ymax - ymin, dx)
            # Here one could add + 1 to be sure that the last pixel is always included
            grid = salem.Grid(nxny=(nx, ny), dxdy=(dx, dx), ll_corner=(xmin, ymin), proj=proj)

            # interpolate TRM and MSG to salem grid
            xi, yi = grid.ij_coordinates
            lon, lat = grid.ll_coordinates

            # Transform lons, lats to grid
            xm, ym = grid.transform(ml0['lon'].flatten(), ml0['lat'].flatten(), crs=salem.wgs84)
            xt, yt = grid.transform(td['lon'].flatten(), td['lat'].flatten(), crs=salem.wgs84)

            # Convert for griddata input
            mpoints = np.array((ym, xm)).T
            tpoints = np.array((yt, xt)).T
            inter = np.array((np.ravel(yi), np.ravel(xi))).T

            # Interpolate using delaunay triangularization
            dummy = griddata(mpoints, ml0['t'].flatten(), inter, method='linear')
            dummy = dummy.reshape((grid.ny, grid.nx))
            outl = np.full_like(dummy, np.nan)
            xl, yl = grid.transform(lon1[inds], lat1[inds], crs=salem.wgs84, nearest=True, maskout=True)

            print('Lag0')
            outl[yl.compressed(), xl.compressed()] = dummy[yl.compressed(), xl.compressed()]

            # Interpolate using delaunay triangularization
            dummyt = griddata(tpoints, td['p'].flatten(), inter, method='linear')
            outt = dummyt.reshape((grid.ny, grid.nx))
            if len(np.where(outt > 0)[0]) < 100:  # at least 100 pixel with rainfall
                print('Kickout: TRMM wavelet min pixel pcp = 100')
                continue

                ##remove edges of interpolated TRMM
            for nb in range(5):
                boole = np.isnan(outt)
                outt[boole] = -1000
                grad = np.gradient(outt)
                outt[boole] = np.nan
                outt[abs(grad[1]) > 300] = np.nan
                outt[abs(grad[0]) > 300] = np.nan

            tmask = np.isfinite(outt)
            mmask = np.isfinite(outl)
            mask2 = np.isfinite(outl[tmask])

            if sum(mask2.flatten()) < 200:  # sum(mmask.flatten())*0.3:
                print('Kickout: TRMM MSG overlap less than 0.3 of cloud area')
                continue

            print('Hit:', gi)

            # lag -1

            # Interpolate using delaunay triangularization
            dummy1 = griddata(mpoints, ml1['t'].flatten(), inter, method='linear')
            dummy1 = dummy1.reshape((grid.ny, grid.nx))
            outl1 = np.full_like(dummy1, np.nan)
            print('Lag1')
            outl1[yl.compressed(), xl.compressed()] = dummy1[yl.compressed(), xl.compressed()]

            # lag -2

            # Interpolate using delaunay triangularization
            dummy2 = griddata(mpoints, ml2['t'].flatten(), inter, method='linear')
            dummy2 = dummy2.reshape((grid.ny, grid.nx))
            outl2 = np.full_like(dummy2, np.nan)
            print('Lag2')
            outl2[yl.compressed(), xl.compressed()] = dummy2[yl.compressed(), xl.compressed()]

            # lag -3

            # Interpolate using delaunay triangularization
            dummy3 = griddata(mpoints, ml3['t'].flatten(), inter, method='linear')
            dummy3 = dummy3.reshape((grid.ny, grid.nx))
            outl3 = np.full_like(dummy3, np.nan)
            print('Lag3')
            outl3[yl.compressed(), xl.compressed()] = dummy3[yl.compressed(), xl.compressed()]

            # lag x

            # Interpolate using delaunay triangularization
            dummyx = griddata(mpoints, mlx['t'].flatten(), inter, method='linear')
            dummyx = dummyx.reshape((grid.ny, grid.nx))
            outlx = np.full_like(dummyx, np.nan)
            print('Lagx')
            outlx[yl.compressed(), xl.compressed()] = dummyx[yl.compressed(), xl.compressed()]

            da = xr.Dataset({'p': (['x', 'y'], outt),
                             't_lag0': (['x', 'y'], dummy),
                             'tc_lag0': (['x', 'y'], outl),
                             't_lag1': (['x', 'y'], dummy1),
                             'tc_lag1': (['x', 'y'], outl1),
                             't_lag2': (['x', 'y'], dummy2),
                             'tc_lag2': (['x', 'y'], outl2),
                             't_lag3': (['x', 'y'], dummy3),
                             'tc_lag3': (['x', 'y'], outl3),
                             't_lagx': (['x', 'y'], dummyx),
                             'tc_lagx': (['x', 'y'], outlx),
                             'tmask': (['x', 'y'], mmask.astype(int)),
                             'pmask': (['x', 'y'], tmask.astype(int))},
                            coords={'lon': (['x', 'y'], lon),
                                    'lat': (['x', 'y'], lat),
                                    'time': date})
            da.attrs['lag0'] = dt0
            da.attrs['lag1'] = dt1
            da.attrs['meanT'] = np.mean(outl[mmask])
            da.attrs['T90perc'] = mmeans
            da.attrs['meanT_cut'] = np.mean(outl[tmask][mask2])
            da.attrs['area'] = sum(mmask.flatten())
            da.attrs['area_cut'] = sum(mask2)
            da.close()
            savefile = '/users/global/cornkle/MCSfiles/SA15_big/' + date.strftime('%Y-%m-%d_%H:%M:%S') + '_' + str(gi) + '.nc'
            try:
                os.remove(savefile)
            except OSError:
                pass
            da.to_netcdf(path=savefile, mode='w')
            print('Saved ' + savefile)

            cnt = cnt + 1

    print('Saved ' + str(cnt) + ' MCSs as netcdf.')


def saveMCS_WA30():
    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_WA30'

    # make a salem grid
    proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')

    t = re.trmm(trmm_folder, yrange=YRANGE, area=[-15, 15, 10, 20])
    m = re.msg(msg_folder)

    cnt = 0

    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered      
    for _y, _m, _d, _h, _mi in zip(t.dates.y, t.dates.m, t.dates.d, t.dates.h, t.dates.mi):

        tdic = t.getDData(_y, _m, _d, _h, _mi, cut=[9, 21])
        if not tdic:
            print('TRMM problem')
            continue

        # define the "0 lag" frist
        arr = np.array([30, 60, 0])
        dm = arr - _mi
        ind = (np.abs(dm)).argmin()

        # set zero shift time for msg
        date = dt.datetime(_y, _m, _d, _h, _mi)
        print(date)

        dt0 = dm[ind]
        ndate = date + dt.timedelta(minutes=int(dt0))
        mdic = m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute,
                         llbox=[tdic['lon'].min(), tdic['lon'].max(), tdic['lat'].min(), tdic['lat'].max()])
        if not mdic:
            print('Date missing')
            continue
        lon1 = mdic['lon']
        lat1 = mdic['lat']
        mdic['t'][mdic['t'] > -40] = 0
        labels, numL = label(mdic['t'])

        u, inv = np.unique(labels, return_inverse=True)
        n = np.bincount(inv)

        goodinds = u[n > 2500]  # all blobs with more than 2500 pixels - size threshold
        print(goodinds)
        if not sum(goodinds) > 0:
            continue

        for gi in goodinds:
            if gi == 0:  # index 0 is always background, ignore!
                continue

            inds = np.where(labels == gi)

            # cut a box for every single blob from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
            latmax, latmin = mdic['lat'][inds].max(), mdic['lat'][inds].min()
            lonmax, lonmin = mdic['lon'][inds].max(), mdic['lon'][inds].min()
            mmeans = np.percentile(mdic['t'][inds], 90)
            td = t.getDData(_y, _m, _d, _h, _mi, cut=[latmin - 0.2, latmax + 0.2])
            if not td:
                print('TRMM problem')
                continue

            # ensure minimum trmm rainfall in area
            if len(np.where(td['p'] > 0)[0]) < 100:  # at least 100 pixel with rainfall
                print('Kickout: TRMM min pixel = 100')
                continue

            dt0 = dm[ind]
            ndate = date + dt.timedelta(minutes=int(dt0))
            # print('Date1', ndate)
            ml0 = m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute,
                            llbox=[lonmin - 0.3, lonmax + 0.3, latmin - 0.25, latmax + 0.25])
            if not ml0:
                continue

            dt2 = dm[ind] - 30
            ndate = date + dt.timedelta(minutes=int(dt2))
            #   print('Date2', ndate)
            ml2 = m.getData(y=ndate.year, m=ndate.month, d=ndate.day, h=ndate.hour, mi=ndate.minute,
                            llbox=[lonmin - 0.3, lonmax + 0.3, latmin - 0.25, latmax + 0.25])
            if not ml2:
                continue

            # create grid surrounding the blob
            # Transform lon, lats to the mercator projection
            x, y = pyproj.transform(salem.wgs84, proj, ml0['lon'], ml0['lat'])
            # take the min and max
            xmax, xmin = np.max(x), np.min(x)
            ymax, ymin = np.max(y), np.min(y)
            # Count the number of pixels
            dx = 5000
            nx, r = divmod(xmax - xmin, dx)
            ny, r = divmod(ymax - ymin, dx)
            # Here one could add + 1 to be sure that the last pixel is always included
            grid = salem.Grid(nxny=(nx, ny), dxdy=(dx, dx), ll_corner=(xmin, ymin), proj=proj)

            # interpolate TRM and MSG to salem grid
            xi, yi = grid.ij_coordinates
            lon, lat = grid.ll_coordinates

            # Transform lons, lats to grid
            xm, ym = grid.transform(ml0['lon'].flatten(), ml0['lat'].flatten(), crs=salem.wgs84)
            xt, yt = grid.transform(td['lon'].flatten(), td['lat'].flatten(), crs=salem.wgs84)

            # Convert for griddata input
            mpoints = np.array((ym, xm)).T
            tpoints = np.array((yt, xt)).T
            inter = np.array((np.ravel(yi), np.ravel(xi))).T

            # Interpolate using delaunay triangularization
            dummy = griddata(mpoints, ml0['t'].flatten(), inter, method='linear')
            dummy = dummy.reshape((grid.ny, grid.nx))
            outl = np.full_like(dummy, np.nan)
            xl, yl = grid.transform(lon1[inds], lat1[inds], crs=salem.wgs84, nearest=True, maskout=True)

            print('Lag0')
            outl[yl.compressed(), xl.compressed()] = dummy[yl.compressed(), xl.compressed()]

            # Interpolate using delaunay triangularization
            dummyt = griddata(tpoints, td['p'].flatten(), inter, method='linear')
            outt = dummyt.reshape((grid.ny, grid.nx))
            if len(np.where(outt > 0)[0]) < 100:  # at least 100 pixel with rainfall
                print('Kickout: TRMM wavelet min pixel pcp = 100')
                continue

                ##remove edges of interpolated TRMM
            for nb in range(5):
                boole = np.isnan(outt)
                outt[boole] = -1000
                grad = np.gradient(outt)
                outt[boole] = np.nan
                outt[abs(grad[1]) > 300] = np.nan
                outt[abs(grad[0]) > 300] = np.nan

            tmask = np.isfinite(outt)
            mmask = np.isfinite(outl)
            mask2 = np.isfinite(outl[tmask])

            if sum(mask2.flatten()) < 200:  # sum(mmask.flatten())*0.3:
                print('Kickout: TRMM MSG overlap less than 0.3 of cloud area')
                continue

            print('Hit:', gi)

            # lag -2

            # Interpolate using delaunay triangularization
            dummy2 = griddata(mpoints, ml2['t'].flatten(), inter, method='linear')
            dummy2 = dummy2.reshape((grid.ny, grid.nx))
            outl2 = np.full_like(dummy2, np.nan)
            print('Lag2')
            outl2[yl.compressed(), xl.compressed()] = dummy2[yl.compressed(), xl.compressed()]

            da = xr.Dataset({'p': (['x', 'y'], outt),
                             't_lag0': (['x', 'y'], dummy),
                             'tc_lag0': (['x', 'y'], outl),
                             't_lag2': (['x', 'y'], dummy2),
                             'tc_lag2': (['x', 'y'], outl2),
                             'tmask': (['x', 'y'], mmask.astype(int)),
                             'pmask': (['x', 'y'], tmask.astype(int))},
                            coords={'lon': (['x', 'y'], lon),
                                    'lat': (['x', 'y'], lat),
                                    'time': date})
            da.attrs['lag0'] = dt0
            da.attrs['lag1'] = dt2
            da.attrs['meanT'] = np.mean(outl[mmask])
            da.attrs['T90perc'] = mmeans
            da.attrs['meanT_cut'] = np.mean(outl[tmask][mask2])
            da.attrs['area'] = sum(mmask.flatten())
            da.attrs['area_cut'] = sum(mask2)
            da.close()
            savefile = '/users/global/cornkle/MCSfiles/WA30/' + date.strftime('%Y-%m-%d_%H:%M:%S') + '_' + str(
                gi) + '.nc'
            try:
                os.remove(savefile)
            except OSError:
                pass
            da.to_netcdf(path=savefile, mode='w')
            print('Saved ' + savefile)

            cnt = cnt + 1

    print('Saved ' + str(cnt) + ' MCSs as netcdf.')
