import numpy as np
import os
from utils import u_arrays as uarr
import pandas as pd
import xarray as xr
from utils import u_time as ut
from utils import u_lists as ul
# import datetime as dt
# import glob
import itertools







"""
MSG folder:
meteosat_WA30: 2004 - 2014, JJAS, 4-21N, 18W - 32E, 580 x 1640 pixel, ~ 3-4km res, ev. 30 mins -> downloading 15mins!
meteosat_SA15: 2006 - 2010, JJAS, 10-20N, 10W - 10E, 350 x 728 pixel, ~ 3-4km, ev. 15 mins
"""

class ReadMsg(object):
    def __init__(self, msg_folder):

        try:
            lpath = uarr.locate('lon.npz', msg_folder)
        except:
            print('Not a directory or no msg lat/lon found')
            quit()

        mpath = os.path.join(msg_folder, 'msg_raw_binary')

        try:
            os.path.isdir(mpath)
        except:
            print('No msg_raw_binary')
            quit()

        msg_latlon = np.load(lpath[0])
        mlon = msg_latlon['lon']
        mlat = msg_latlon['lat']

        self.lat = mlat
        self.lon = mlon
        self.nx = mlon.shape[1]
        self.ny = mlon.shape[0]

        self.years = os.listdir(mpath)
        self.root = msg_folder

    def set_date(self, yr, mon, day, hr, mins):

        self.dpath = os.path.join(self.root, 'msg_raw_binary', str(yr), str(mon).zfill(2),
                                  str(yr) + str(mon).zfill(2) + str(day).zfill(2) + str(hr).zfill(2) + str(
                                      mins).zfill(2) + '.gra')

        if not os.path.isfile(self.dpath):  # if the date is not found, it's silently omitted. Not perfect but allows loop without massive print
            self.dpath = False              # self.dpath = False can be caught as non existant file
            return

        root = os.path.join(self.root, 'cell_blob_files')
        file = os.path.join(root, str(yr), str(mon).zfill(2),
                            str(yr) + str(mon).zfill(2) + str(day).zfill(2) + str(hr).zfill(2) + str(
                                mins).zfill(2) + '.gra')
        if os.path.isfile(file):
            self.bpath = file
        else:
            print('No blob file dir found!')
            self.bpath = False

        root = os.path.join(self.root, 'bigcell_area_table')
        file = os.path.join(root, 'rewrite',
                                'cell_40c_' + str(hr).zfill(2) + str(mins).zfill(2) + '_JJAS.txt')
        if os.path.isfile(file):
            self.tpath = file
        else:
            print('No table file dir found!')
            self.tpath = False

        self.date = [pd.datetime(yr, mon, day, hr, mins)]

    def get_data(self, llbox=None, netcdf_path=None):

        if not self.dpath:
            print('No data for date or date not set. Please set set_date first')
            return False

        rrShape = (self.ny, self.nx)  # msg shape
        rrMDI = np.uint8(255)
        rr = np.fromfile(self.dpath, dtype=rrMDI.dtype)
        rr.shape = rrShape
        rr = rr.astype(np.int32) - 173

        if llbox:
            i, j = np.where(
                (self.lon > llbox[0]) & (self.lon < llbox[2]) & (self.lat > llbox[1]) & (self.lat < llbox[3]))
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


    def get_blob(self, llbox=None):

        if not self.bpath:
            print('No blob file dir found!')
            return False

        rrShape = (self.ny, self.nx)  # msg shape
        rrMDI = np.uint16()
        rr = np.fromfile(self.bpath, dtype=rrMDI.dtype)
        rr.shape = rrShape
        if llbox:
            i, j = np.where(
                (self.lon > llbox[0]) & (self.lon < llbox[2]) & (self.lat > llbox[1]) & (self.lat < llbox[3]))
            blat = self.lat[i.min():i.max() + 1, j.min():j.max() + 1]
            blon = self.lon[i.min():i.max() + 1, j.min():j.max() + 1]
            rr = rr[i.min():i.max() + 1, j.min():j.max() + 1]
        else:
            blat = self.lat
            blon = self.lon
            rr = rr

        date = self.date

        da = xr.DataArray(rr[None, ...], coords={'time': (('time'), date),
                                                 'lat': (('y', 'x'), blat),
                                                 'lon': (('y', 'x'), blon)},
                          dims=['time', 'y', 'x']).isel(time=0)

        return da


    def get_table(self):

        if not self.tpath:
            print('No table file dir found!')
            return False

        tab = pd.read_csv(self.tpath)

        return tab
