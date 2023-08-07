import numpy as np
import xarray as xr
from utils import u_arrays as ua, u_darrays as uda
import matplotlib.pyplot as plt
import multiprocessing
import pickle as pkl
from collections import defaultdict
from utils import constants as cnst, u_met


import pdb
import glob
import ipdb
import pandas as pd

def dictionary():

    dic = {}
    vars = ['hour', 'month', 'year', 'area',
            'lon', 'lat', 'clon', 'clat',
            'tmin', 'tmean', 'tcw', 'shear', #'ef',
            'pmax', 'pmean',
             't', 'p' ]

    for v in vars:
        dic[v] = []
    return dic

def perSys():

    pool = multiprocessing.Pool(processes=3)
    tthresh = '-50'
    files = glob.glob('/home/ck/DIR/cornkle/data/DYAMOND/UM-5km/MCS/WA/UM-5km/*.nc')
    #ipdb.set_trace()


    pool = multiprocessing.Pool(processes=3)
    mdic = dictionary() #defaultdict(list)
    res = pool.map(file_loop, files)
    pool.close()

    # res = []
    # for f in files[500:1000]:
    #     out = file_loop(f)
    #     res.append(out)
    #
    #
    # ipdb.set_trace()

    print('Back from multiproc')
    keys = mdic.keys()
    for v in res:
        for k in keys:
            try:
                mdic[k].append(v[k])
            except TypeError:
                continue


    pkl.dump(mdic, open('/home/ck/DIR/cornkle/data/DYAMOND/UM-5km/MCS/WA/UM-5k_5000km2_-40_tir_newShear.p',
                           'wb'))


def file_loop(f):
    print('Doing file: ' + f)

    out = dictionary()

    ds = xr.open_dataset(f)

    #print(ds['time.hour'])
    out['hour'] = ds['time.hour'].item()

    #edate = pd.Timestamp(ds.time.values)


    tmask = (ds['prcp'].values>0.1)&(ds['tir'].values<-40)


    if np.sum(tmask) ==0:
        return
    #ipdb.set_trace()
    tt = ds['tir'].values#.copy()
    tt[~tmask] = np.nan
    pp = ds['prcp'].values#.copy()
    pp[~tmask] = np.nan

    maxpos = np.unravel_index(np.nanargmax(pp),pp.shape)
    minpos = np.unravel_index(np.nanargmin(tt), tt.shape)

    out['tmean'] = ds.attrs['meanT']
    out['tmin'] = np.nanmean(ua.cut_kernel(tt,minpos[1], minpos[0],1)) #ds.attrs['minT']
    out['pmean'] = ds.attrs['meanP']
    out['pmax'] = np.nanmean(ua.cut_kernel(pp,maxpos[1], maxpos[0],1))#ds.attrs['maxP']
    out['area'] = ds.attrs['area']
    out['t'] = ds['tir'].values[tmask]
    out['p'] = ds['prcp'].values[tmask]

    out['tcw'] = np.nanmean(ua.cut_kernel(ds['tcw'].values,minpos[1], minpos[0],1))
   # ipdb.set_trace()
    try:
        out['shear'] = np.nanmean(ua.cut_kernel(ds['shear'].values, minpos[1], minpos[0], 1))
    except:
        return
    #out['ef'] = np.nanmean(ua.cut_kernel(ds['ef'].values, minpos[1], minpos[0], 1))
    out['lon'] = ds['longitude'].values
    out['lat'] = ds['latitude'].values

    out['month'] = ds['time.month'].item()
    out['year'] = ds['time.year'].item()
    out['date'] = ds['time'].values

    if out['pmax'] < 1:
        return

    out['clat'] = np.min(out['lat'])+((np.max(out['lat'])-np.min(out['lat']))*0.5)
    out['clon'] = np.min(out['lon']) + ((np.max(out['lon']) - np.min(out['lon'])) * 0.5)


    return out
