import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
import datetime as dt
from eod import trmm, msg
import xarray as xr
import pandas as pd
import os
import glob
import pdb
from utils import constants
HOD = range(24)  # hours of day
YRANGE = range(2004, 2015)


def run():
    trmm_folder = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    msg_folder = '/users/global/cornkle/data/OBS/meteosat_WA30'

    t = trmm.ReadWA(trmm_folder, yrange=YRANGE, area=[-15, 20, 4, 25])   # [-15, 15, 4, 21], [-10, 10, 10, 20]
    m = msg.ReadMsg(msg_folder)


    # define the "0 lag" frist
    arr = np.array([15, 30, 45, 60, 0])

    # cycle through TRMM dates - only dates tat have a certain number of pixels in llbox are considered
    cnt = 0
    for _y, _m, _d, _h, _mi in zip(t.dates.y, t.dates.m, t.dates.d, t.dates.h, t.dates.mi):

        tdic = t.get_ddata(_y, _m, _d, _h, _mi, cut=[4, 25])
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
            mdic = m.get_data(llbox=[tdic['lon'].values.min(),  tdic['lon'].values.max(), tdic['lat'].values.min(),
                                     tdic['lat'].values.max()])

            if not mdic:
                print('Date missing')
                continue

        print('TRMM:', date, 'MSG:', ndate.year, ndate.month, ndate.day, ndate.hour, ndate.minute )

        cnt += np.nansum(tdic['p'].values >= 30)

    print('Number of extreme pixels TRMM: ', cnt)


def checkBig40():

    files = glob.glob('/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size_zR/*.nc')
    cnt = 0
    for f in files:
        ds = xr.open_dataset(f)
        p = ds['p']
        cnt += np.nansum(p.values >= 30)
        print(cnt)
    print('Number of extreme pixels clouds -40', cnt)