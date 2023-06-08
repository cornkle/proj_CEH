# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import multiprocessing
import ipdb
import glob
import pandas as pd
import os
from utils import constants as cnst
from GLOBAL import glob_util
### SIMILAR VERSION IN NOTEBOOK

mregions = glob_util.MREGIONS

def multi():
    pool = multiprocessing.Pool(processes=5)
    yy = range(2000,2021)
    res = pool.map(run, yy)
    pool.close()


def run(year):
    for rr in ['china', 'WAf', 'india', 'sub_SA', 'australia', 'SAf', 'GPlains']:
        extract_box(rr, year)


def extract_box(region, year):
    dtag = region
    #files = glob.glob(cnst.lmcs_drive + 'MCS_Feng/tracks/custom/' + dtag + '/*_'+str(year)+'*.nc')
    files = glob.glob(cnst.lmcs_drive + 'MCS_Feng/global_v2/tracks/custom/' + dtag + '/*_'+str(year)+'*.nc')
    out = cnst.lmcs_drive + 'save_files_v2/'

    dumpkeys = ['datetimestring', 'movement_r', 'movement_theta', 'movement_r_meters_per_second',
                'movement_time_lag', 'movement_storm_x', 'movement_storm_y', 'pf_nuniqpix', 'location_idx',
                'pixel_duration', 'julian_day',
                'pixel_pcp', 'pf_skewness', 'mergecloudnumber', 'splitcloudnumber']

    for ff in files:

        ds = xr.open_dataset(ff)
        fname = os.path.basename(ff)

        outname = fname[0:-3].replace('robust', region + '_initTime_')

        outfilename = out + outname + '.csv'

        if os.path.isfile(outfilename):
            print('File exists, continue')
            continue

        pfdic = {}

        for dv in ds.isel(times=0, tracks=0).data_vars:
            if dv in dumpkeys:
                continue

            if "nmaxpf" in ds[dv].dims:

                for pftag in ['1','2','3']:
                    pfdic[str(dv)+str(pftag)] = []
            elif dv=="direction":
                for pftag in ['1','-1','-2','0']:
                    pfdic[str(dv)+str(pftag)] = []
            else:
                pfdic[dv] = []

        pfdic['utc_time'] = []
        pfdic['utc_date'] = []
        pfdic['utc_year'] = []
        pfdic['utc_month'] = []
        pfdic['utc_hour'] = []
        pfdic['utc_day'] = []
        pfdic['utc_minute'] = []
        pfdic['tracktime'] = []
        pfdic['trackid'] = []
        pfdic['londiff_loc-init'] = []
        pfdic['latdiff_loc-init'] = []
        pfdic['init_lon'] = []
        pfdic['init_lat'] = []
        pfdic['utc_init_time'] = []
        pfdic['utc_init_hour'] = []
        pfdic['lt_init_time'] = []
        pfdic['lt_init_hour'] = []
        pfdic['lt_time'] = []
        pfdic['lt_date'] = []
        pfdic['lt_hour'] = []

        for ids, ai in enumerate(ds.tracks):
            track = ds.sel(tracks=ai)

            init_lon = track.sel(times=0)['meanlon'].values
            init_lat = track.sel(times=0)['meanlat'].values

            init_date = pd.to_datetime(track.sel(times=0)['base_time'])
            init_date_lt = glob_util.UTC_to_LT_date(init_date, region)

            for tids in track.times:
                tt = track.sel(times=tids)
                tm1 = tids-1
                tm2 = tids-2
               
                tm0 = tids+1

                if np.isnan(tt['meanlat']):
                    continue
                print('Doing', tt['base_time'].values)

                print('Location ', tt['meanlat'].values, tt['meanlon'].values)

                print('Writing ', tt['base_time'].values)

                pfdic['tracktime'].append(int(tids.values))
                pfdic['trackid'].append(int(track.tracks.values))
                pfdic['londiff_loc-init'].append(float(tt['meanlon'].values - init_lon))
                pfdic['latdiff_loc-init'].append(float(tt['meanlat'].values - init_lat))
                pfdic['init_lon'].append(float(init_lon))
                pfdic['init_lat'].append(float(init_lat))

                pfdic['utc_init_hour'].append(int(init_date.hour))
                pfdic['utc_init_time'].append(init_date.values)
                pfdic['lt_init_hour'].append(int(init_date_lt.hour))
                pfdic['lt_init_time'].append(init_date_lt.values)


                for dv in tt.data_vars:
                    if dv in dumpkeys:
                        continue

                    if dv == 'base_time':
                        dtime = pd.to_datetime(tt[dv].values)
                        lt_dtime = glob_util.UTC_to_LT_date(dtime, region)

                        pfdic['utc_year'].append(dtime.year)
                        pfdic['utc_month'].append(dtime.month)
                        pfdic['utc_hour'].append(dtime.hour)
                        pfdic['utc_minute'].append(dtime.minute)
                        pfdic['utc_day'].append(dtime.day)
                        pfdic['utc_time'].append(tt[dv].values)
                        pfdic['utc_date'].append(tt[dv].replace(hour=0, minute=0).values)

                        pfdic['lt_time'].append(lt_dtime.values)
                        pfdic['lt_date'].append(lt_dtime.replace(hour=0, minute=0).values)
                        pfdic['lt_hour'].append(lt_dtime.hour)

                    elif dv == 'direction':
                        pfdic['direction0'].append(tt[dv].values)
                        if tm1>=0:
                            pfdic['direction-1'].append(track.sel(times=tm1)['direction'].values)
                        else:
                            pfdic['direction-1'].append(np.nan)
                        if tm2>=0:
                            pfdic['direction-2'].append(track.sel(times=tm2)['direction'].values)
                        else:
                            pfdic['direction-2'].append(np.nan)
                        try:
                            pfdic['direction1'].append(track.sel(times=tm0)['direction'].values)
                        except:
                            pfdic['direction1'].append(np.nan)

                    elif (tt[dv].size==3) & ("pf" in dv):
                        for pfids, pftag in enumerate(['1', '2', '3']):
                            pfdic[str(dv) + str(pftag)].append(tt[dv].values[pfids])
                    else:
                        pfdic[dv].append(tt[dv].values)


        # for kk in pfdic.keys():
        #     try:
        #         print(kk, len(pfdic[kk]))
        #     except:
        #         print(kk, pfdic[kk])
        df = pd.DataFrame.from_dict(pfdic)
        df.to_csv(outfilename, index=False)
        print('Saved ', outfilename)
