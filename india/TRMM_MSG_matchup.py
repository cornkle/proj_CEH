import salem
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage.measurements import label
import datetime as dt
from eod import msg, trmm, tm_utils
import xarray as xr
import os
import matplotlib.pyplot as plt
from utils import u_grid
import pdb

HOD = range(24)  # hours of day
YRANGE = range(2004, 2015)


def saveMCS():
    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_WA30'

    t = trmm.ReadWA(trmm_folder, yrange=YRANGE, area=[-15, 4, 20, 25])  # (ll_lon, ll_lat, ur_lon, ur_lat) define initial TRMM box and scan for swaths in that box
    m = msg.ReadMsg(msg_folder)

    cnt = 0

    # minute array to find closest MSG minute
    arr = np.array([15, 30, 45, 60, 0])

    # loop through TRMM dates - only dates that have a certain number of pixels in llbox are considered
    for _y, _m, _d, _h, _mi in zip(t.dates.y, t.dates.m, t.dates.d, t.dates.h, t.dates.mi):

        tdic = t.get_ddata(_y, _m, _d, _h, _mi, cut=[3,26]) # cut TRMM data at lower/upper lat
        #get value of closest minute
        dm = arr - _mi
        dm = dm[dm<0]
        try:
            ind = (np.abs(dm)).argmin()
        except ValueError:
            continue

        # set smallest lag time for msg
        date = dt.datetime(_y, _m, _d, _h, _mi)

        dt0 = dm[ind]
        ndate = date + dt.timedelta(minutes=int(dt0))
        m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
        mdic = m.get_data(llbox=[tdic['lon'].values.min(),  tdic['lat'].values.min(), tdic['lon'].values.max(),tdic['lat'].values.max()])

        # check whether date is completely missing or just 30mins interval exists
        if not mdic:
            dm = np.delete(dm, np.argmin(np.abs(dm)), axis=0)
            # try second closest minute
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
            mdic = m.get_data(llbox=[tdic['lon'].values.min(), tdic['lat'].values.min(), tdic['lon'].values.max(),
                                     tdic['lat'].values.max()])
            if not mdic:
                print('Date missing')
                continue

        print('TRMM:', date, 'MSG:', ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute )

        lon1 = mdic['lon'].values # MSG coords
        lat1 = mdic['lat'].values
        mdic['t'].values[mdic['t'].values >= -10] = 0  # T threshold -10 for clouds
        ### filter minimum cloud size
        labels, numL = label(mdic['t'].values)
        u, inv = np.unique(labels, return_inverse=True)
        n = np.bincount(inv)
        goodinds = u[n > 39]  # defines minimum MCS size e.g. 9x39 ~ 350km2
        print(goodinds) # indices of clouds of "good size"

        if not sum(goodinds) > 0:
            continue

        for gi in goodinds:
            if gi == 0:  # index 0 is always background, ignore!
                continue

            inds = np.where(labels == gi) # position of cloud

            # cut a box for every single blob (cloud) from msg - get min max lat lon of the blob, cut upper lower from TRMM to match blob
            latmax, latmin = lat1[inds].max(), lat1[inds].min()
            lonmax, lonmin = lon1.values[inds].max(), lon1[inds].min()
            mmeans = np.percentile(mdic['t'].values[inds], 90)
            td = t.get_ddata(_y, _m, _d, _h, _mi, cut=[latmin - 1, latmax + 1]) # for each cloud, cut TRMM swath

            dt0 = dm[ind]

            ml0 = m.get_data(llbox=[lonmin - 1, latmin - 1, lonmax + 1, latmax + 1]) # cut cloud box in MSG
            if not ml0:
                continue

            #make salem grid
            grid = u_grid.make(ml0['lon'].values, ml0['lat'].values,5000)  # 5km regular grid from lat/lon coords
            lon, lat = grid.ll_coordinates # 5km grid lat/lon coordinates

            # interpolate TRMM and MSG to 5km common grid
            inter, mpoints = u_grid.griddata_input(ml0['lon'].values, ml0['lat'].values,grid)
            inter, tpoints = u_grid.griddata_input(td['lon'].values, td['lat'].values, grid)

            # Interpolate TRMM using delaunay triangularization
            try:
                dummyt = griddata(tpoints, td['p'].values.flatten(), inter, method='linear')
            except ValueError:
                continue
            outt = dummyt.reshape((grid.ny, grid.nx))

            if np.sum(np.isfinite(outt)) < 5:  # at least 5 valid pixel
                print('Kickout: TRMM min pixel  < 5')
                continue

            # Interpolate TRMM flags USING NEAREST
            dummyf = griddata(tpoints, td['flags'].values.flatten(), inter, method='nearest')
            outf = dummyf.reshape((grid.ny, grid.nx))
            outf=outf.astype(np.float)
            isnot = np.isnan(outt)
            outf[isnot]=np.nan

            ##remove artefact edges of interpolated TRMM
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
            outff = tm_utils.getTRMMconv(outf) ## from TRMM flags, get positions of convective rain
            outk = np.zeros_like(outt)
            outk[np.where(outff)]=outt[np.where(outff)]

            # Interpolate MSG using delaunay triangularization
            dummy = griddata(mpoints, ml0['t'].values.flatten(), inter, method='linear')
            dummy = dummy.reshape((grid.ny, grid.nx))
            outl = np.full_like(dummy, np.nan)
            xl, yl = grid.transform(lon1[inds], lat1[inds], crs=salem.wgs84, nearest=True, maskout=True)
            outl[yl.compressed(), xl.compressed()] = dummy[yl.compressed(), xl.compressed()]

            # TODO #### SHIFTING WITH RESPECT TO MIN T / MAX P - search for Pmax within 20km from Tmin, shift TRMM image
            #
            # tmin = np.argmin(outl)
            # pmax =
            #
            # dist =
            #

            tmask = np.isfinite(outt)
            mmask = np.isfinite(outl)
            mask2 = np.isfinite(outl[tmask])

            #last check for min area, crazy rainfall or crazy cloud size
            if (sum(mmask.flatten())*25 < 350) or (outt.max()>200) or (sum(mmask.flatten())*25 > 1500000):
                continue

            if sum(mask2.flatten()) < 5:  # Check minimum overlap between TRMM swath and MSG cloud
                print('Kickout: TRMM MSG overlap less than 3pix of cloud area')
                continue

            print('Hit:', gi)

            da = xr.Dataset({'p': (['x', 'y'], outt),  # rainfall field
                             'pconv': (['x', 'y'], outk), # convective rainfall
                             't_lag0': (['x', 'y'], dummy), # full T image in cutout region
                             'tc_lag0': (['x', 'y'], outl), # cloud area only
                             },
                            coords={'lon': (['x', 'y'], lon),
                                    'lat': (['x', 'y'], lat),
                                    'time': date})
            da.attrs['lag0'] = dt0  # lag in minutes between TRMM / MSG
            da.attrs['meanT'] = np.mean(outl[mmask])  # cloud mean T
            da.attrs['T90perc'] = mmeans # cloud 90perc T
            da.attrs['meanT_cut'] = np.mean(outl[tmask][mask2]) # cloud mean T in TRMM region
            da.attrs['area'] = sum(mmask.flatten()) # total cloud area
            da.attrs['area_cut'] = sum(mask2)  # cloud area overlapping with TRMM
            da.close()
            savefile = '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_zR/' + date.strftime('%Y-%m-%d_%H:%M:%S') + '_' + str(gi) + '.nc'
            try:
                os.remove(savefile)
            except OSError:
                pass
            da.to_netcdf(path=savefile, mode='w')
            print('Saved ' + savefile)

            cnt = cnt + 1

    print('Saved ' + str(cnt) + ' TRMM/MSG merged MCSs as netcdf.')
