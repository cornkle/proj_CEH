# -*- coding: utf-8 -*-


import numpy as np
from saveCore_standalone_v2 import util
import xarray as xr
from saveCore_standalone_v2 import powerBlob_utils, util, run_powerBlobs
import multiprocessing
from utils import u_grid, u_interpolate as u_int
import glob
import ipdb
import os
import pandas as pd


def run(dataset='METEOSAT5K_vera', detection='sum'):

    met_folder = '/media/ck/Elements/Africa/WestAfrica/NFLICS/MCS_TIR/real_time_data/2020/09/06/'

    files = glob.glob(met_folder + 'IR_108_BT_*.nc')

    data_resolution, dist, start, nb = util.read_dic(util.DATASETS[dataset])

    dummy = xr.open_dataset(files[0], decode_times=False)

    #data_resolution = 5 # in km
    # make salem grid
    grid = u_grid.make(np.arange(-19,0), np.arange(4,20), data_resolution*1000)
    dlon = dummy['lon_2d'].squeeze().values.T
    dlat = dummy['lat_2d'].squeeze().values.T
    inds, weights, shape = u_int.interpolation_weights_grid(dlon, dlat, grid)

    for f in files:

        print('Doing ', f)

        outfile = f.replace('108', 'wavelet')
        outfile = outfile.replace('real_time_data', 'real_time_wavelet_old')

        if os.path.isfile(outfile):
            print('File exists, continue')
            continue

        try:
            ds = xr.open_dataset(f, decode_times=False)
        except:
            print('File error for file', f)
            continue

        ff = os.path.basename(f)
        year = ff[10:14]
        month = ff[14:16]
        day = ff[16:18]
        hour = ff[19:21]
        minute = ff[21:23]


        date = [pd.datetime(int(year), int(month), int(day), int(hour), int(minute))]

        data = ds['IR108_BT'].squeeze().values.T
        try:
            outt = u_int.interpolate_data(data, inds, weights, shape)
        except IndexError:
            print('Interpolation problem, continue')
        lon, lat = grid.ll_coordinates

        run_powerBlobs.wavelet_analysis(outt, lon, lat, date, outfile, data_resolution=data_resolution, detection=detection, dataset=dataset)

