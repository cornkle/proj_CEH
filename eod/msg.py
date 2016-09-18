import numpy as np
import os
from utils import u_arrays as uarr
import pandas as pd
from utils import u_time as ut
from utils import u_lists as ul
# import datetime as dt
# import glob
import itertools


class msg(object):
    def __init__(self, msg_folder):

        if not os.path.isdir(msg_folder):
            print('Not a directory')
            quit()

        lpath = uarr.locate('lon.npz', msg_folder)
        mpath = os.path.join(msg_folder, 'msg_raw_binary')

        if not os.path.isdir(mpath):
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
        if not os.path.isfile(self.dpath):
            self.dpath = False

        if os.path.isdir(os.path.join(self.root, 'cell_blob_files')):
            self.bpath = os.path.join(self.root, 'cell_blob_files', str(yr), str(mon).zfill(2),
                                      str(yr) + str(mon).zfill(2) + str(day).zfill(2) + str(hr).zfill(2) + str(
                                          mins).zfill(2) + '.gra')
        else:
            print('No blob file dir found!')
            self.bpath = False
        if os.path.isdir(os.path.join(self.root, 'bigcell_area_table')):
            self.tpath = os.path.join(self.root, 'bigcell_area_table', 'rewrite',
                                      'cell_40c_' + str(hr).zfill(2) + str(mins).zfill(2) + '_JJAS.txt')
        else:
            print('No table file dir found!')
            self.tpath = False
        self.date = {'year': yr, 'month': mon, 'day': day, 'hour': hr,
                     'minute': mins}  # could be turned into python date
