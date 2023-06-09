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
from utils import constants as cnst
import multiprocessing
import glob
import pickle as pkl


def dictionary():

    dic = {}
    vars = ['date', 'month', 'hour', 'year', 'area', '70area',
            'minlon', 'minlat', 'maxlon', 'maxlat', 'clon', 'clat', 'tminlon', 'tminlat'
            'tmin', 'tmean', 't10', 't90', 'img_stormID']


    for v in vars:
        dic[v] = []
    return dic


def mcs_define(array, thresh, min_area=None, max_area=None, minmax_area=None):
    """

    :param array: 2d input array
    :param thresh: cloud threshold
    :param min_area: minimum area of the cloud
    :param max_area: maximum area of the cloud
    :param minmax_area: tuple indicating only clouds bigger than tuple[0] and smaller than tuple[1]
    :return: 2d array with labelled blobs
    """
    array[array >= thresh] = 0  # T threshold maskout
    array[np.isnan(array)] = 0  # set ocean nans to 0

    labels, numL = label(array)

    u, inv = np.unique(labels, return_inverse=True)
    n = np.bincount(inv)

    goodinds = u[u!=0]

    if min_area != None:
        goodinds = u[(n>=min_area) & (u!=0)]
        badinds = u[n<min_area]

        # for b in badinds:
        #     pos = np.where(labels==b)
        #     labels[pos]=0

    if max_area != None:
        goodinds = u[(n<=max_area)  & (u!=0)]
        badinds = u[n>max_area]

    if minmax_area != None:
        goodinds = u[(n <= minmax_area[1]) & (u != 0) & (n>=minmax_area[0])]
        badinds = u[(n > minmax_area[1]) | (n < minmax_area[0])]

    if (min_area is not None) | (max_area is not None) | (minmax_area is not None):
        for b in badinds:
            pos = np.where(labels==b)
            labels[pos]=0

    return labels, goodinds


def process_image(ctt, data_res, t_tresh=-50, min_mcs_size=5000):

    min_pix_nb = min_mcs_size / data_res**2
    max_pix_nb = 500000 / data_res**2  # this is to capture satellite artefacts that come in large contiguous stripes.
    labels, goodinds = ua.mcs_define(da.values, t_thresh, minmax_area=[min_pix_nb, max_pix_nb]) # 7.7x7.7km = 64km2 per pix in gridsat? 83 pix is 5000km2

    for g in goodinds:

        if g==0:
            continue

        pos = np.where(labels==g)
        dic = dictionary()

        ts = pd.to_datetime(ctt['time'].values)
        date = ts.strftime('%Y-%m-%d_%H:%M:%S')

        dic['date'] = ts
        dic['month'] = int(d['time.month'])
        dic['hour'] = int(d['time.hour'])
        dic['year'] = int(d['time.year'])

        storm = d.values[pos]
        tmin_pos = np.argmin(storm)

        dic['area'] = storm.size
        dic['70area'] = np.sum(storm<=-70)
        dic['minlon'] = np.nanmin(d.lon.values[pos[1]])
        dic['minlat'] = np.nanmin(d.lat.values[pos[0]])
        dic['maxlon'] = np.nanmax(d.lon.values[pos[1]])
        dic['maxlat'] = np.nanmax(d.lat.values[pos[0]])
        dic['clon'] = dic['minlon'] + (dic['maxlon'] - dic['minlon'])/2
        dic['clat'] = dic['minlat'] + (dic['maxlat'] - dic['minlat'])/2
        dic['tmin'] = np.nanmin(storm)
        dic['tminlat'] = d.lat[tmin_pos[0]]
        dic['tminlon'] = d.lon[tmin_pos[1]]
        dic['tmean'] = np.nanmean(storm)
        dic['tp5'] = np.percentile(storm, 5)
        dic['tp95'] = np.percentile(storm, 95)
        dic['img_stormID'] = date + '_' + str(g)

        ydic[date + '_' + str(g)] = dic

    return pd.DataFrame.from_dict(dic)
