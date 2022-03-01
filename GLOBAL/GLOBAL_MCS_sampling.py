# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import multiprocessing
import ipdb
import glob
import pandas as pd
import pickle as pkl
import salem
import os
### SIMILAR VERSION IN NOTEBOOK

mregions = {'WAf' : [[-18,25,4,25], 'spac', 0], # last is hourly offset to UCT # 12
 'SAf' : [[20,35, -35,-15], 'spac', 2], # 10
 'india' : [[70,90, 5,30], 'asia', 5], # 7
 'china' : [[105,115,25,40], 'asia', 8 ], # 4
 'australia' : [[120,140,-23, -11], 'asia', 9], # 3
 'sub_SA' : [[-68,-47, -40, -20.5], 'spac', -4] , # 16
 'trop_SA' : [[-75, -50, -20, -5], 'spac', -5], # 17
 'GPlains' : [[-100,-90,32,47], 'nam', -6] # # 18

}

def run_loop(region):

    dtag = mregions[region][1]
    box = mregions[region][0]
    files = glob.glob('/media/ck/Elements/global/MCS_Feng/tracks/' + dtag + '/*_201*.nc')

    out = '/home/ck/DIR/cornkle/data/GLOBAL_MCS/save_files/'

    dumpkeys = ['datetimestring', 'movement_r', 'movement_theta', 'movement_r_meters_per_second',
                'movement_time_lag', 'movement_storm_x', 'movement_storm_y', 'pf_nuniqpix', 'location_idx',
                'pixel_duration',
                'pixel_pcp', 'pf_skewness']

    for ff in files:

        fname = os.path.basename(ff)

        outname = fname[0:-3].replace('robust', region + '_noinit_')

        outfilename = out + outname + '.p'

        if os.path.isfile(outfilename):
            print('File exists, continue')
            #del ds
            continue

        ds = xr.open_dataset(ff).load()
        # ds = ds.sel(tracks=slice(0,5))

        # ipdb.set_trace()
        pfdic = {}

        for dv in ds.sel(times=0, tracks=0).data_vars:
            if dv in dumpkeys:
                continue
            pfdic[dv] = []
        pfdic['year'] = []
        pfdic['month'] = []
        pfdic['hour'] = []
        pfdic['minute'] = []
        pfdic['tracktime'] = []

        for ids, ai in enumerate(ds.tracks):
            track = ds.sel(tracks=ai)

            init_lat = track.sel(times=0)['meanlat']
            init_lon = track.sel(times=0)['meanlon']

            for tids in track.times:
                tt = track.sel(times=tids)

                # ipdb.set_trace()

                if np.isnan(tt['meanlat']):
                    continue
                #print('Doing', tt['base_time'].values)
                if (tt['meanlat'] < box[2]) | (tt['meanlat'] > box[3]) | (tt['meanlon'] < box[0]) | (
                        tt['meanlon'] > box[1]):
                    continue

                if (np.abs(tt['meanlat'].values-init_lat.values) < 1.5) & (np.abs(tt['meanlon'].values-init_lon.values) < 1.5):
                    print('Init too close, time number ', tids.values)

                print('Location ', tt['meanlat'].values, tt['meanlon'].values)

                print('Writing ', tt['base_time'].values)

                pfdic['tracktime'].append(tids.values)

                for dv in tt.data_vars:
                    if dv in dumpkeys:
                        continue

                    if dv == 'base_time':
                        dtime = pd.Timestamp(tt[dv].values)

                        pfdic['year'].append(dtime.year)
                        pfdic['month'].append(dtime.month)
                        pfdic['hour'].append(dtime.hour)
                        pfdic['minute'].append(dtime.minute)

                        # ipdb.set_trace()

                    pfdic[dv].append(tt[dv].values)

        pkl.dump(pfdic, open('/home/ck/DIR/cornkle/data/GLOBAL_MCS/save_files/' + outname + '.p', "wb"))
        del ds