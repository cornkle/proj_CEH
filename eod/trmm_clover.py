# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 13:50:12 2016

@author: cornkle
"""

import numpy as np
import os
from utils import u_arrays as uarr
from utils import u_time as ut
from utils import u_lists as ul
import itertools
import xarray as xr
import pandas as pd
import pandas as pd
import datetime as dt
import pdb

HOD = list(range(24))
YRANGE = range(2004, 2015)
MRANGE = range(3, 12)  # Jun - Sep
MTRESH = 0

"""
# =======================================================================================
# Reads TRMM 2A25 files to get a list with valid file names.
# Files with _rain_f4 in name mark files where rainfall is available over West Africa.
# However, we use the files with less floating points (without the _rain_f4 )
#
# Output: Filelist of TRMM files over West Africa
#
# Keywords:
# area: box [llon,llat,ulon,ulat], list
        filter for TRMM to have at least 5000 pixels in box
  hod:  list of hours, hours of day to be searched for, default range(24)
  yrange: list of years, default 2004-2012
  mrange: list of months, default 6-9
# =======================================================================================

# =======================================================================================
# Reads the data of TRMM 2A25 swaths binary files with lon, lat, rain, flag in it
# Every swath is 49 pixel wide. 2 bytes integer
# Lon = lon/100, Lat = lat/100, pcp = rain/10 (mm/h), flag
#
# =======================================================================================
"""

class ReadWA(object):
    def __init__(self, trmm_folder, yrange=YRANGE, mrange=MRANGE, hod=HOD, area=None):


        min_rain_swath = 200
        min_rain_box = 200
        min_tpixel = 2500
        rain_thresh = 0.1

        if not os.path.isdir(trmm_folder):
            print('Not a directory')
            return

        fdic = {'fpath': [], 'tmins': [], 'date': []}
        rfiles = []

        for yr, mo in itertools.product(yrange, mrange):  # rain_f4 files only available for 6 to 10

            tpath = os.path.join(trmm_folder, str(yr), str(mo).zfill(2))
            try:
                files = uarr.locate('.7.gra', tpath)
            except OSError:
                continue

            rfiles.extend(files)

        rfiles.sort(key=ul.natural_keys)

        if not rfiles:
            print('No trmm files found')
            return

            #  self.fpath=fdic['fpath']
            #  return
        for eachfile in rfiles:
            rain_str = eachfile
            time_str = eachfile.replace('.7.', '.7_time.')
            try:
                rr = np.fromfile(time_str, dtype=np.float32)  # seconds of day
            except FileNotFoundError:
                print(time_str+' missing, continue')
                continue

            secmean = rr.mean()
            try:
                t = ut.sec_to_time(secmean)
            except ValueError:
                print('ValueError sec to time')
                continue
            if not t.hour in hod:
                continue

            rr = np.fromfile(rain_str, dtype=np.int16)
            x = 49  # trmm swath is always 49 wide
            nb = rr.size
            single = int(nb / 4)  # variables lon lat rainrate flag

            lons = rr[0:single]
            lats = rr[single:2 * single]
            rainrs = rr[2 * single:3 * single]
            y = int(lons.size / x)
            lons = np.resize(lons, (y, x))
            lats = np.resize(lats, (y, x))
            rainrs = np.resize(rainrs, (y, x))
            lont = lons / 100.
            latt = lats / 100.
            rain = rainrs / 10.

            if np.sum(rain>rain_thresh) < min_rain_swath:  # minimum TRMM rainfall > 0.1 in swath
                continue
            if area:
                box = np.where((lont > area[0]) & (lont < area[1]) & (latt > area[2]) & (latt < area[3]))

                if not box[0].any():
                    continue
                    #       print(len(box[0]))
                if len(box[0]) < min_tpixel:  # minimum pixel overlap with TRMM and box (50000km2)
                    continue
                if np.sum(rain[box]>rain_thresh) < min_rain_box:  # minimum rainfall in defined box
                    continue


            fdic['fpath'].append(rain_str)
            # fdic['date'].add(int(rain_str[-20:-16]), int(rain_str[-16:-14]), int(rain_str[-14:-12]), t.hour, t.minute,
            #                  0)

            fdic['date'].append(pd.datetime(int(rain_str[-20:-16]), int(rain_str[-16:-14]), int(rain_str[-14:-12]), t.hour, t.minute,
                             0))

        self.fpaths = fdic['fpath']
        self.dates = pd.Series(fdic['date'])
        self.__area = area



    def get_ddata(self, date, cut=None, netcdf_path=None):
        """
        Gets a file with a certain date out of the initialised TRMM file list
        Keywords:
        cut: [lower,upper], list, gets the data and cuts swath at the given latitude upper and lower bondaries
        """

        pos = np.where(self.dates.dt.strftime('%Y-%m-%d_%H:%M') == date.strftime('%Y-%m-%d_%H:%M'))

        # print('Ind:', ind)

        if len(pos[0]) == 0:
            print('No data for date')
            return False

        tfile = self.fpaths[pos[0][0]]
        da = self.get_data(tfile, cut=cut, netcdf_path=netcdf_path)

        return da



    def get_data(self, path, cut=None, netcdf_path=None):
        """
           Gets TRMM data given the path to the file.
           Automatically crops TRMM to the initialised box (see ReadWA) or, if not given, to 3 - 22N, even if no "cut" is given!!
           Keywords:
           cut: [lower,upper], list, gets the data and cuts swath at the given latitude upper and lower bondaries. Can be used
           cut out smaller North-South ranges than initialised box (or 3-22N).
           """
        tfile = path
        if not os.path.isfile(tfile):
            print('File does not exist. Error')
            quit()

        rr = np.fromfile(tfile, dtype=np.int16)
        x = 49  # trmm swath is always 49 wide
        nb = rr.size
        single = int(nb / 4)  # variables lon lat rainrate flag

        lont = rr[0:single] / 100.
        latt = rr[single:2 * single] / 100.
        rain = rr[2 * single:3 * single] / 10.
        flags = rr[3 * single:4 * single]

        latpath = tfile.replace('.gra', '_lat_f4.gra')
        if os.path.isfile(latpath):
            latt = np.fromfile(latpath, dtype=np.float32)
            lonpath = tfile.replace('.gra', '_lon_f4.gra')
            lont = np.fromfile(lonpath, dtype=np.float32)


        y = int(lont.size / x)
        lont = np.resize(lont, (y, x))
        latt = np.resize(latt, (y, x))
        rain = np.resize(rain, (y, x))
        flags = np.resize(flags, (y, x))

        laty = latt[:, 0]

        if cut:
            lower = cut[0]
            upper = cut[1]
            cutt = np.where((laty <= upper) & (laty >= lower))
            cutt = np.array(cutt)

        try:
            rain = rain[cutt, :]
            lont = lont[cutt, :]
            latt = latt[cutt, :]
            flags = flags[cutt, :]
            yy = int(rain.size / x)

            rain = rain.reshape(yy, x)
            lont = lont.reshape(yy, x)
            latt = latt.reshape(yy, x)
            flags = flags.reshape(yy, x)
        except UnboundLocalError:
            pass  # don't cut if nothing to cut. Stay with former ranges

        dind = self.fpaths.index(tfile)
        date = self.dates[dind].to_pydatetime()

        da = xr.DataArray(rain[None, ...], coords={'time': (('time'), np.atleast_1d(date)),
                                                     'lat': (('y', 'x'), latt),
                                                     'lon': (('y', 'x'), lont)},
                            dims=['time', 'y', 'x']).isel(time=0)

        da2 = xr.DataArray(flags[None, ...], coords={'time': (('time'), np.atleast_1d(date)),
                                                   'lat': (('y', 'x'), latt),
                                                   'lon': (('y', 'x'), lont)},
                          dims=['time', 'y', 'x']).isel(time=0)

        ds = xr.Dataset({'p': da, 'flags': da2})


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
