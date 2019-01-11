# -*- coding: utf-8 -*-


import numpy as np
import datetime as dt
import xarray as xr
import os
import ipdb
import matplotlib.pyplot as plt
from utils import u_grid, u_arrays as ua
import ipdb
import pandas as pd
from utils import constants as cnst, u_met
import multiprocessing
import glob
import pickle as pkl


def dictionary():

    dic = {}
    vars = ['date', 'month', 'year', 'area', '70area',
            'minlon', 'minlat', 'maxlon', 'maxlat', 'clon', 'clat',
            'tmin', 'tmean', 't10', 't90']

    for v in vars:
        dic[v] = []
    return dic

def perSys():

    pool = multiprocessing.Pool(processes=6)
    tthresh = '-50'

    files = glob.glob(cnst.CP4_PATH + 'CLOVER/CP4hist/lw_out_PBLtop/*.nc')

    print('Nb files', len(files))
    #mdic = dictionary() #defaultdict(list)
    res = pool.map(file_loop, files)
    pool.close()

    # for f in files:
    #     res = file_loop(f)
    #     ipdb.set_trace()

    merged = {}
    for r in res:
        if r is not None:
            merged.update(r)

    test = pd.DataFrame.from_dict(merged, orient='index')

    pkl.dump(test, open(cnst.CLOVER_SAVES + 'StormLoc_CP4hist_-50_5000km_WA.p',
                           'wb'))


def file_loop(f):

    print('Doing ' + f)

    ydic = {}

    df = xr.open_dataset(f)
    df = df.sel(latitude=slice(-4,13), longitude=slice(-18,20))
    if (int(df['time.month'][0]) < 9) & (int(df['time.month'][0]) > 6):
         print('return')
         return

    df = df['lw_out_PBLtop'][df['time.hour'] == 18]

    for d in df:

        print('Read data')

        d.values = u_met.olr_to_bt(d.values)
        #ipdb.set_trace()
        labels, goodinds = ua.blob_define(d.values, -50, minmax_area=[258, 40000], max_area=None) # 4.4*4.4km: 258 pixel = 5000km2, 40000 pixel = 350000km2
        for g in goodinds:

            if g==0:
                continue

            pos = np.where(labels==g)

            dic = dictionary()

            ts = pd.to_datetime(d['time'].values)
            date = ts.strftime('%Y-%m-%d_%H:%M:%S')

            dic['date'] = ts


            dic['month'] = int(d['time.month'])
            dic['year'] = int(d['time.year'])

            storm = d.values[pos]

            dic['area'] = storm.size
            dic['70area'] = np.sum(storm<=-70)
            dic['minlon'] = np.min(d.longitude.values[pos[1]])
            dic['minlat'] = np.min(d.latitude.values[pos[0]])
            dic['maxlon'] = np.max(d.longitude.values[pos[1]])
            dic['maxlat'] = np.max(d.latitude.values[pos[0]])
            dic['clon'] = dic['minlon'] + (dic['maxlon'] - dic['minlon'])/2
            dic['clat'] = dic['minlat'] + (dic['maxlat'] - dic['minlat'])/2
            dic['tmin'] = np.min(storm)
            dic['tmean'] = np.mean(storm)
            dic['t10'] = np.percentile(storm, 10)
            dic['t90'] = np.percentile(storm, 10)

            ydic[date + '_' + str(g)] = dic

    return ydic
