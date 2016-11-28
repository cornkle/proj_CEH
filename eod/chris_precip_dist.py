# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
from wavelet import util
from utils import u_arrays as ua
from scipy import ndimage
import matplotlib.pyplot as plt
from eod import tm_utils
import multiprocessing
import ipdb
import pandas as pd
import pickle as pkl


def run():
    pool = multiprocessing.Pool(processes=7)
    files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E/')   # /WA30/
    out = '/users/global/cornkle/C_paper/chris2016/'
    #files = files[0:1000]
    print('Nb files', len(files))

    res = pool.map(file_loop, files)
    pool.close()
    res = [x for x in res if x is not None]

    nb_sys = len(res)

    print('Number systems: ', nb_sys)

    res = [item for sublist in res for item in sublist] # flatten list of lists



    dic = {'year' : [], 'month' : [], 'hour' : [], 'precip':[], 'sum30' : [],  'sum20' :[],  'sum':[], 'valid':[],
           'nz':[], 'clon' : [], 'clat' : [], 'cent': [] }  #  big , fin,shape, sum, sumvalid, tmin

    for v in res:

        dic['year'].append(v[0])
        dic['month'].append(v[1])
        dic['hour'].append(v[2])
        dic['precip'].append(v[3])
        dic['sum30'].append(v[4])
        dic['sum20'].append(v[5])
        dic['sum'].append(v[6])
        dic['valid'].append(v[7])
        dic['nz'].append(v[8])
        dic['clat'].append(v[9])
        dic['clon'].append(v[10])
        dic['cent'].append(v[11])


    #df = pd.DataFrame(dic)
    #df.to_pickle(out+'3dmax_gt15000_fakeprecip_-70.pkl')

    pkl.dump(dic, open(out+'chris_mcs_-40_gt1000.p',
                           'wb'))




def file_loop(f):
    ret = []

    print('Doing file: ' + f)

    dic = xr.open_dataset(f)

    outt = dic['tc_lag0'].values
    outp = dic['p'].values  #*0+1   ## ATTENTION CHANGED RAINFALL

    tmean = np.nanmean(outt)

    outp[np.isnan(outt)]=np.nan

    area = np.sum(outt < -40)

    pp = np.max(outp[(np.isfinite(outp)) & (np.isfinite(outt))])
    clat = np.min(dic.lat.values) + ((np.max(dic.lat.values) - np.min(dic.lat.values)) * 0.5)
    clon = np.min(dic.lon.values) + ((np.max(dic.lon.values) - np.min(dic.lon.values)) * 0.5)

    #area gt 3000km2 cause that's 30km radius if circular
    if (area * 25 < 1000) or (area * 25 > 500000) or (pp < 1) or (pp > 200): #or (np.sum(np.isfinite(outp)) < (np.sum(np.isfinite(outt))*0.1)):
        print(area*25)
        print('throw out')
        return


    sum30 = np.sum(outp>=30)
    sum20 = np.sum(outp>=20)
    try:
        cent = np.percentile(outp[(outp>0.1)],90)
    except IndexError:
        cent=np.nan
    sum = np.nansum(outp)
    valid = np.sum([np.isfinite(outp)])
    nz = np.sum([(outp>0.1)] )
    rval = outp[(outp>0.1)]

   # ipdb.set_trace()


    year= dic['time.year'].values.tolist()
    month=dic['time.month'].values.tolist()
    hour=dic['time.hour'].values.tolist()

    #### HOW TO GIVE BACK THE MAX SCALE PER SYSTEM??

    ret.append((year, month, hour, rval, sum30, sum20, sum, valid, nz, clat, clon, cent))

    dic.close()

    return ret


if __name__ == "__main__":
    run()
