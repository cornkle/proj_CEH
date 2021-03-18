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
import datetime




"""
MSG folder:
meteosat_WA30: 2004 - 2014, JJAS, 4-21N, 18W - 32E, 580 x 1640 pixel, ~ 3-4km res, ev. 30 mins -> downloading 15mins!
meteosat_SA15: 2006 - 2010, May-Oct, 10-20N, 10W - 10E, 350 x 728 pixel, ~ 3-4km, ev. 15 mins

meteosat_tropWA: 2004 - 2015, whole year, 4-10N, 14W - 25E, 350 x 728 pixel, ~ 3-4km, ev. 15 mins, daytime only
"""
y1 = 2004 #2006
y2 = 2020 #

class ReadMsg(object):
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
            lpath = uarr.locate('MSGbigDomain_1804_580_lat_lon.npz', msg_folder, exclude = None)
        except:
            print('Not a directory or no msg lat/lon found')
            return

        try:
            spath = uarr.locate('MSGstrip_164_580_lat_lon.npz', msg_folder, exclude = None)
        except:
            print('Not a directory or no msg lat/lon found')
            return


        mpath = os.path.join(msg_folder, 'msg_raw_binary')

        try:
            os.path.isdir(mpath)
        except:
            print('No msg_raw_binary')
            quit()

        rfiles = []
        for yr, mo in itertools.product(yrange, mrange):  # rain_f4 files only available for 6 to 10

            filepath = os.path.join(mpath, str(yr), str(mo).zfill(2))
            try:
                files = uarr.locate('.gra', filepath, exclude = ['_182', '_167'])
            except OSError:
                continue

            rfiles.extend(files)

        rfiles.sort(key=ul.natural_keys)

        msg_latlon = np.load(lpath[0])
        mlon = msg_latlon['lon']
        mlat = msg_latlon['lat']

        strip_shape = np.load(spath[0])['lon'].shape

        self.lat = mlat
        self.lon = mlon
        self.nx = mlon.shape[1]-strip_shape[1] # mlon.shape[1]
        self.nx_big = mlon.shape[1]
        self.ny = mlon.shape[0] #mlon.shape[0]
        self.s_shape = strip_shape
        self.l_shape = (mlon.shape[0]-strip_shape[0], mlon.shape[1])
        #ipdb.set_trace()
        self.years = os.listdir(mpath)
        self.root = msg_folder
        self.fpath = rfiles

    def set_date(self, yr, mon, day, hr, mins):

        self.dpath = os.path.join(self.root, 'msg_raw_binary', str(yr), str(mon).zfill(2),
                                  str(yr) + str(mon).zfill(2) + str(day).zfill(2) + str(hr).zfill(2) + str(
                                      mins).zfill(2) + '.gra')

        if not os.path.isfile(self.dpath):  # if the date is not found, it's silently omitted. Not perfect but allows loop without massive print
            self.dpath = False              # self.dpath = False can be caught as non existant file
            return


        root = os.path.join(self.root, 'msg_raw_strip')
        file = os.path.join(root, str(yr), str(mon).zfill(2),
                            str(yr) + str(mon).zfill(2) + str(day).zfill(2) + str(hr).zfill(2) + str(
                                mins).zfill(2) + '.gra')
        if os.path.isfile(file):
            self.bpath = file
        else:
            print('No strip file dir found!')
            self.bpath = None
            #return # only run through if strip exists

        self.date = [datetime.datetime(yr, mon, day, hr, mins)]


    def get_data(self, llbox=None, netcdf_path=None):

        if not self.dpath:
            print('I found no msg files in my dpath')
            return False

        rrShape = (self.ny, self.nx)  # msg shape
        rrMDI = np.uint8(255)
        rr = np.fromfile(self.dpath, dtype=rrMDI.dtype)
        rr.shape = rrShape

        rr = rr.astype(np.int32) - 173

        if (self.nx < self.lon.shape[0]) & (self.bpath is None):
            print('Stitch file missing for big domain. Returning.', self.date)
            return


        ssShape = self.s_shape  # msg shape
        ssMDI = np.uint8(255)
        ss = np.fromfile(self.dpath, dtype=ssMDI.dtype)
        ss.shape = ssShape
        try:
            ss = ss.astype(np.int32) - 173
        except:
            print('Something wrong with strip file')
            return

        stitched = np.flip(np.concatenate((ss,rr), axis=1), axis=0) #np.flip


        if llbox:
            i, j = np.where(
                (self.lon > llbox[0]) & (self.lon < llbox[1]) & (self.lat > llbox[2]) & (self.lat < llbox[3]))
            blat = self.lat[i.min():i.max() + 1, j.min():j.max() + 1]
            blon = self.lon[i.min():i.max() + 1, j.min():j.max() + 1]
            stitched = stitched[i.min():i.max() + 1, j.min():j.max() + 1]
        else:
            blat = self.lat
            blon = self.lon

        date = self.date  # or np.atleast_1d(dt.datetime())

        da = xr.DataArray(stitched[None, ...], coords={'time': (('time'), date),
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

        str = file.split(os.sep)[-1]
        curr_date = [
            datetime.datetime(np.int(str[0:4]), np.int(str[4:6]), np.int(str[6:8]), np.int(str[8:10]), np.int(str[10:12]))]
        date = curr_date  # or np.atleast_1d(dt.datetime())

        self.set_date(date[0].year, date[0].month, date[0].day, date[0].hour, date[0].minute)

        rrShape = (self.ny, self.nx)  # msg shape
        rrMDI = np.uint8(255)

        rr = np.fromfile(file, dtype=rrMDI.dtype)

        try:
            rr.shape = rrShape
        except ValueError:
            print('Got big domain, no slice needed')
            rrShape = (self.ny, self.nx_big)
            rr.shape = rrShape


        #rr.shape = self.l_shape
        rr = rr.astype(np.int32) - 173  #self.s_shape

        if (self.nx < self.lon.shape[0]) & (self.bpath is None):
            print('Stitch file missing for big domain or something else wrong with dimensions. Returning.', self.date)
            return

        if (self.bpath!=None) & (rrShape[1] != self.lon.shape[0]):
            ssShape = self.s_shape  # msg shape
            ssMDI = np.uint8(255)
            try:
                ss = np.fromfile(self.bpath, dtype=ssMDI.dtype)
            except:
                print('Something wrong with strip file')
                return

            ss.shape = ssShape

            ss = ss.astype(np.int32) - 173

            stitched = np.flip(np.concatenate((ss,rr), axis=1), axis=0)

        else:
            stitched = rr



        if llbox:
            i, j = np.where(
                (self.lon > llbox[0]) & (self.lon < llbox[1]) & (self.lat > llbox[2]) & (self.lat < llbox[3]))
            blat = self.lat[i.min():i.max() + 1, j.min():j.max() + 1]
            blon = self.lon[i.min():i.max() + 1, j.min():j.max() + 1]
            stitched = stitched[i.min():i.max() + 1, j.min():j.max() + 1]
        else:
            blat = self.lat
            blon = self.lon



        da = xr.DataArray(stitched[None, ...], coords={'time': (('time'), date),
                                                   'lat': (('y', 'x'),  blat),
                                                   'lon': (('y', 'x'), blon)},
                          dims=['time', 'y', 'x']).isel(time=0)

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
