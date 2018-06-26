import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
import datetime as dt
from eod import trmm, msg
import xarray as xr
import pandas as pd
import os
from eod.old import msg_old
import glob
import pdb
from utils import constants, u_arrays as ua
HOD = range(24)  # hours of day
YRANGE = range(2004, 2010)#2015)


def run(tthreshold):
    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_WA30'

    t = trmm.ReadWA(trmm_folder, yrange=YRANGE, area=[-15, 20, 5, 14])   # [-15, 15, 4, 21], [-10, 10, 10, 20]
    m = msg.ReadMsg(msg_folder)

    msg_zipped = zip(m.lon.flatten(), m.lat.flatten()) # (x,y)
    msg_tuples = np.array(list(msg_zipped))

    # define the "0 lag" frist
    arr = np.array([15, 30, 45, 60, 0])

    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered
    cnt = 0
    cold_cnt = 0
    for _y, _m, _d, _h, _mi in zip(t.dates.y, t.dates.m, t.dates.d, t.dates.h, t.dates.mi):

        tdic = t.get_ddata(_y, _m, _d, _h, _mi, cut=[5, 14])
        if tdic['p'].max() < 30:
            continue
        #get closest minute
        dm = arr - _mi
        dm = dm[dm<0]
        try:
            ind = (np.abs(dm)).argmin()
        except ValueError:
            continue

        # set zero shift time for msg
        date = dt.datetime(_y, _m, _d, _h, _mi)

        dt0 = dm[ind]
        ndate = date + dt.timedelta(minutes=int(dt0))
        m.set_date(ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute)
        if not m.dpath:
            continue
        mdic = m.get_data(llbox=[tdic['lon'].values.min(),  tdic['lon'].values.max(), tdic['lat'].values.min(),tdic['lat'].values.max()])

        # check whether date is completely missing or just 30mins interval exists
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
            if not m.dpath:
                print('Date missing')
                continue
            mdic = m.get_data(llbox=[tdic['lon'].values.min(),  tdic['lon'].values.max(), tdic['lat'].values.min(),
                                     tdic['lat'].values.max()])


        print('TRMM:', date, 'MSG:', ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute )

        cnt += np.nansum(tdic['p'].values >= 30)

        gt_lon = tdic['lon'].values[tdic['p'].values >= 30]
        gt_lat = tdic['lat'].values[tdic['p'].values >= 30]


        for tpoint in zip(gt_lon, gt_lat):

            isnearest = ua.closest_point(tpoint, msg_tuples)

            pos = np.where((mdic.lat==msg_tuples[isnearest][1]) & (mdic.lon==msg_tuples[isnearest][0]))
            try:
                ttest = mdic['t'].isel(y=slice(pos[0][0]-2, pos[0][0]+3 ), x=slice(pos[1][0]-2, pos[1][0]+3 ))
            except IndexError:
                continue
            if np.sum(ttest <= tthreshold) > 0:
                cold_cnt += 1

    print('Number of extreme pixels TRMM: ', cnt)
    print('Number of cold MSG <='+str(tthreshold)+'C overlapping extreme pixels TRMM: ', cold_cnt)
    print('Fraction below thresh:', cold_cnt/cnt)



def checkBig40():

    files = glob.glob('/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size_zR/*.nc')
    cnt = 0
    for f in files:
        ds = xr.open_dataset(f)
        p = ds['p']
        cnt += np.nansum(p.values >= 30)
        print(cnt)
    print('Number of extreme pixels clouds -40', cnt)