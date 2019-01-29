import numpy as np
import os
from utils import u_arrays as uarr
import pandas as pd
import xarray as xr
from utils import u_time as ut
from utils import u_lists as ul
import datetime
# import glob
import itertools
import ipdb
import glob




"""
MFG folder:
on wllf012, 1983-2005, 1094x463,
"""
y1 = 1999#2006
y2 = 2006 #
class ReadMfg(object):
    def __init__(self, msg_folder, y1=y1, y2=y2, months=None):

        yrange = range(y1, y2+1)  # 1998, 2014
        if months is None:
            mrange = range(1,13)
        else:
            if len(months) > 1:
                mrange = months
            else:
                mrange = range(months[0],months[0]+1)

        try:
            lpath = uarr.locate('lon.npz', msg_folder, exclude = None)
        except:
            print('Not a directory or no msg lat/lon found')
            return

        mpath = os.path.join(msg_folder, 'mfg_raw_binary')

        try:
            os.path.isdir(mpath)
        except:
            print('No mfg_raw_binary')
            quit()

        rfiles = []

        for yr, mo in itertools.product(yrange, mrange):  # rain_f4 files only available for 6 to 10

            filepath = os.path.join(mpath, str(yr), str(mo).zfill(2))
            try:
                files = glob.glob(mpath + os.sep + str(yr) + str(mo).zfill(2) + '*')

            except OSError:
                continue
            #print(rfiles)
            rfiles.extend(files)

        rfiles.sort(key=ul.natural_keys)

        msg_latlon = np.load(lpath[0])
        mlon = msg_latlon['lon']
        mlat = msg_latlon['lat']

        self.lat = mlat
        self.lon = mlon
        self.nx = mlon.shape[1]
        self.ny = mlon.shape[0]

        years = []
        outfiles = []
        for r in rfiles:

            years.append(os.path.basename(r)[0:4])
            outfiles.append(r + os.sep + 'tir.gra')

        self.years = years
        self.root = msg_folder
        self.fpath = outfiles

    def set_date(self, yr, mon, day, hr, mins):

        self.dpath = os.path.join(self.root, 'mfg_raw_binary', str(yr) + str(mon).zfill(2) +
         str(day).zfill(2) , 'tir.gra')

        if not os.path.isfile(self.dpath):  # if the date is not found, it's silently omitted. Not perfect but allows loop without massive print
            self.dpath = False              # self.dpath = False can be caught as non existant file
            return

        self.date = [pd.datetime(yr, mon, day, hr, mins)]


    def get_data(self, llbox=None, netcdf_path=None):

        if not self.dpath:
            print('I found no mfg files in my dpath')
            return False

        rrShape = (self.ny, self.nx)  # msg shape
        rrMDI = np.uint8(255)
        rr = np.fromfile(self.dpath, dtype=rrMDI.dtype)
        rr.shape = rrShape

        rr = rr.astype(np.int32) - 173

        if llbox:
            i, j = np.where(
                (self.lon > llbox[0]) & (self.lon < llbox[1]) & (self.lat > llbox[2]) & (self.lat < llbox[3]))
            blat = self.lat[i.min():i.max() + 1, j.min():j.max() + 1]
            blon = self.lon[i.min():i.max() + 1, j.min():j.max() + 1]
            rr = rr[i.min():i.max() + 1, j.min():j.max() + 1]
        else:
            blat = self.lat
            blon = self.lon
            rr = rr

        date = self.date  # or np.atleast_1d(dt.datetime())

        da = xr.DataArray(rr[None, ...], coords={'time': (('time'), date),
                                                   'lat': (('y', 'x'), blat),
                                                   'lon': (('y', 'x'), blon)},
                          dims=['time', 'y', 'x']).isel(time=0)

        ds = xr.Dataset({'t': da})

        if netcdf_path:
            savefile = netcdf_path

            try:
                os.remove(savefile)
            except OSError:
                pass
            try:
                ds.to_netcdf(path=savefile, mode='w')
            except OSError:
                print('File cannot be saved. Maybe directory wrong. Check your permissions')
                raise

            print('Saved ' + savefile)

        return ds


    def read_data(self, file, llbox=None, netcdf_path=None):

        if not self.fpath:
            print('I found no msg files in my fpath')
            return False

        rrShape = (48, self.ny, self.nx)  # msg shape, 48 time steps per day
        rrMDI = np.uint8(255)

        rr = np.fromfile(file, dtype=rrMDI.dtype)
        rr.shape = rrShape
        rr = rr.astype(np.int32) #- 173

        if llbox:
            i, j = np.where(
                (self.lon > llbox[0]) & (self.lon < llbox[1]) & (self.lat > llbox[2]) & (self.lat < llbox[3]))
            blat = self.lat[i.min():i.max() + 1, j.min():j.max() + 1]
            blon = self.lon[i.min():i.max() + 1, j.min():j.max() + 1]
            rr = rr[i.min():i.max() + 1, j.min():j.max() + 1]
        else:
            blat = self.lat
            blon = self.lon
        rr = rr

        str = file.split(os.sep)[-2]
        curr_date = []
        for hour in range(24):
            for mins in np.array(['00','30']).astype(int):
                curr_date.append(
                pd.datetime(np.int(str[0:4]), np.int(str[4:6]), np.int(str[6:8]), hour, mins))

        date = curr_date  # or np.atleast_1d(dt.datetime())

        da = xr.DataArray(rr, coords={'time': (('time'), date),
                                                   'lat': (('y', 'x'), blat),
                                                   'lon': (('y', 'x'), blon)},
                          dims=['time', 'y', 'x'])#.isel(time=0)

        # if llbox:
        #     da = da.where((da.lon > llbox[0]) & (da.lon < llbox[1]) & (da.lat > llbox[2]) & (da.lat < llbox[3]), drop=True)

        ds = xr.Dataset({'t': da})

        if netcdf_path:
            savefile = netcdf_path

            try:
                os.remove(savefile)
            except OSError:
                pass
            try:
                ds.to_netcdf(path=savefile, mode='w')
            except OSError:
                print('File cannot be saved. Maybe directory wrong. Check your permissions')
                raise

            print('Saved ' + savefile)

        return ds
