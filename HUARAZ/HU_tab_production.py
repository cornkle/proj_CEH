# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import ipdb
import pandas as pd
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_arrays as ua, constants as cnst, u_darrays
from scipy.interpolate import griddata
import multiprocessing
import os
import glob

import pickle as pkl
from wavelet import util as wutil


bigbox = [-79,-65,-17,-3]

fname = '/home/ck/DIR/cornkle/data/HUARAZ/shapes/riosan_sel_one.shp'
sdf = salem.read_shapefile(fname)
sdf = salem.transform_geopandas(sdf, to_crs=salem.wgs84)

file = cnst.GRIDSAT_PERU + 'daily_5000km2_LT_afternoon/gridsat_WA_-40Min_5000km2_13-19UTCperDay_'
file1 = cnst.GRIDSAT_PERU + 'daily_5000km2_UTC_afternoon/gridsat_WA_-40Min_5000km2_UTCDay_'
file2 = cnst.GRIDSAT_PERU + 'daily_ALLkm2_UTC_DAY/gridsat_WA_-40Min_ALLkm2_UTCDay_'

gbox = bigbox
isbuffer = gbox

#topo = xr.open_dataarray(cnst.TOPO_1MIN)
#topo = topo.sel(lon=slice(isbuffer[0], isbuffer[1]), lat=slice(isbuffer[2], isbuffer[3]))

chirps = xr.open_dataset(cnst.elements_drive + 'SouthAmerica/CHIRPS/chirps-v2.0.daily.peru.nc').chunk({'time':365})
chirps = chirps['precip'].sel(longitude=slice(isbuffer[0], isbuffer[1]), latitude=slice(isbuffer[2], isbuffer[3]))


def tab():
    
    wet_tab = pd.read_csv('/home/ck/DIR/cornkle/data/HUARAZ/dry_wet_schlump/wet_gl_sp.csv', names=['date', 'tag'])
    #dry_tab = pd.read_csv('/home/ck/DIR/cornkle/data/HUARAZ/dry_wet_schlump/dry_gl_sp.csv', names=['date', 'tag'])
    wet_tab = wet_tab[wet_tab['date'] > '1984-12-31']
    wet_tab = wet_tab
    wet_tab = setup_dic(wet_tab)
    wet_tab['dateM'] = wet_tab.date.map(lambda x: x[0:-3])

    chunk, chunk_ind, chunk_count = np.unique(wet_tab.dateM, return_index=True, return_counts=True)

    #ipdb.set_trace()

    chunks = [wet_tab.loc[wet_tab.index[ci:ci + cc]] for ci, cc in zip(chunk_ind, chunk_count)]  # daily chunks

    # res = []
    # for m in chunks[100:101]:
    #     out = file_loop(m)
    #     res.append(out)
    #
    # ipdb.set_trace()
    # return
    pool = multiprocessing.Pool(processes=5)

    res = pool.map(file_loop, chunks)
    pool.close()

    print('Returned from parallel')
    res = [x for x in res if x is not None]

    df_concat = pd.concat(res)

    #dic = df_concat.to_dict()

    df_concat.to_csv('/home/ck/DIR/cornkle/data/HUARAZ/dry_wet_schlump/wet_only_ERA5_GRIDSAT_UTCDay_allDates_q200.csv',
                   na_rep=-999, index_label='id')

    print('Save file written!')
    print('Dumped file')

def setup_dic(tab):

    #Huaraz
    tab['MCSmean'] = np.nan
    tab['MCSmin'] = np.nan
    tab['MCScover'] = np.nan
    
    tab['MCSmeanUTC'] = np.nan
    tab['MCSminUTC'] = np.nan
    tab['MCScoverUTC'] = np.nan
    
    tab['CloudcoverUTC'] = np.nan
    tab['CloudmeanUTC'] = np.nan
    tab['CloudminUTC'] = np.nan
    
    tab['CHIRPSmean'] = np.nan
    tab['CHIRPSmax'] = np.nan
    
    ##############
    
    tab['MCSmeanA'] = np.nan
    tab['MCSminA'] = np.nan
    tab['MCScoverA'] = np.nan
    
    tab['MCSmeanUTCA'] = np.nan
    tab['MCSminUTCA'] = np.nan
    tab['MCScoverUTCA'] = np.nan
    
    tab['CloudcoverUTCA'] = np.nan
    tab['CloudmeanUTCA'] = np.nan
    tab['CloudminUTCA'] = np.nan
    
    
    tab['CHIRPSmeanA'] = np.nan
    tab['CHIRPSmaxA'] = np.nan
    
    ##############
    
    tab['MCSmeanV'] = np.nan
    tab['MCSminV'] = np.nan
    tab['MCScoverV'] = np.nan
    
    tab['MCSmeanUTCV'] = np.nan
    tab['MCSminUTCV'] = np.nan
    tab['MCScoverUTCV'] = np.nan
    
    tab['CloudcoverUTCV'] = np.nan
    tab['CloudmeanUTCV'] = np.nan
    tab['CloudminUTCV'] = np.nan
    
    
    tab['CHIRPSmeanV'] = np.nan
    tab['CHIRPSmaxV'] = np.nan
    
    # Huaraz
    tab['q550'] = np.nan
    tab['q650'] = np.nan
    tab['q200'] = np.nan
    tab['t550'] = np.nan
    tab['t650'] = np.nan
    tab['u550'] = np.nan
    tab['r650'] = np.nan
    tab['tcwv'] = np.nan
    tab['shear'] = np.nan # 200-550
    tab['t2'] = np.nan 
    tab['cape'] = np.nan 
    tab['cin'] = np.nan
    tab['u200'] = np.nan 
    tab['v200'] = np.nan 
    tab['u550'] = np.nan
    tab['v550'] = np.nan 
    
    # Amazon [-73.6, -72.6, -10.2, -8.6]
    tab['tcwvA'] = np.nan
    tab['q550A'] = np.nan
    tab['q850A'] = np.nan
    tab['q200A'] = np.nan
    tab['t550A'] = np.nan
    tab['t850A'] = np.nan
    tab['shearA'] = np.nan  # 550-850
    tab['u200A'] = np.nan 
    tab['v200A'] = np.nan 
    tab['u550A'] = np.nan
    tab['v550A'] = np.nan 
    tab['u850A'] = np.nan
    tab['v850A'] = np.nan 
    tab['t2A'] = np.nan
    tab['capeA'] = np.nan
    tab['cinA'] = np.nan
    tab['r850A'] = np.nan
    
    # Valley [-75.6, -74.6, -10.2, -8.6]
    tab['tcwvV'] = np.nan
    tab['q550V'] = np.nan
    tab['q850V'] = np.nan
    tab['q200V'] = np.nan
    tab['t550V'] = np.nan
    tab['t850V'] = np.nan
    tab['shearV'] = np.nan # 550=850
    tab['u200V'] = np.nan 
    tab['v200V'] = np.nan 
    tab['u550V'] = np.nan
    tab['v550V'] = np.nan 
    tab['u850V'] = np.nan
    tab['v850V'] = np.nan 
    tab['t2V'] = np.nan 
    tab['capeV'] = np.nan
    tab['cinV'] = np.nan
    tab['r850V'] = np.nan
    
    tab['month'] = np.nan
    tab['year'] = np.nan
    tab['hour'] = np.nan

    return tab



def file_loop(df):

    monthdate = pd.to_datetime(df['date'].iloc[0])
    ff = file + str(monthdate.year) + '-' + str(monthdate.month).zfill(2) + '.nc'
    try:
        dat = xr.open_dataset(ff)
    except:
        return

    ff1 = file1 + str(monthdate.year) + '-' + str(monthdate.month).zfill(2) + '.nc'
    dat1 = xr.open_dataset(ff1)

    ff2 = file2 + str(monthdate.year) + '-' + str(monthdate.month).zfill(2) + '.nc'
    dat2 = xr.open_dataset(ff2)


    for dids, dit in df.iterrows():

        date = dit['date']
        dt = pd.to_datetime(date)

        print('Doing', date)

        gridsat = dat['tir'].sel(time=date, lon=slice(gbox[0], gbox[1]), lat=slice(gbox[2], gbox[3])).squeeze() / 100

        if gridsat.values.size == 0:
            continue

        gridsat1 = dat1['tir'].sel(time=date, lon=slice(gbox[0], gbox[1]), lat=slice(gbox[2], gbox[3])).squeeze() / 100
        gridsat2 = dat2['tir'].sel(time=date, lon=slice(gbox[0], gbox[1]), lat=slice(gbox[2], gbox[3])).squeeze() / 100

        try:
            era5 = xr.open_dataset(
                cnst.ERA5_HOURLY_PL_HU + '/ERA5_' + str(dt.year) + '_' + str(dt.month).zfill(2) + '_' + str(
                    dt.day).zfill(2) + '_pl.nc')
        except:
            continue
        # era5 = uda.flip_lat(era5)
        erah = 15
        eraH = era5.sel(time=era5['time.hour'] == erah, latitude=-9.51, longitude=-77.55,
                        method='nearest').squeeze()  # 10LT  -77.55,-9.51
        eraA = era5.sel(time=era5['time.hour'] == erah, latitude=slice(-8.5, -9), longitude=slice(-73.6, -73.1)).mean(
            ['latitude', 'longitude']).squeeze()
        eraV = era5.sel(time=era5['time.hour'] == erah, latitude=slice(-8.5, -9), longitude=slice(-75.6, -75.1)).mean(
            ['latitude', 'longitude']).squeeze()

        try:
            era5s = xr.open_dataset(
                cnst.ERA5_HOURLY_SRFC_HU + '/ERA5_' + str(dt.year) + '_' + str(dt.month).zfill(2) + '_' + str(
                    dt.day).zfill(2) + '_srfc.nc')
        except:
            continue
        # era5s = uda.flip_lat(era5s)
        eraHs = era5s.sel(time=era5s['time.hour'] == erah, latitude=-9.51, longitude=-77.55,
                          method='nearest').squeeze()  # 10LT  -77.55,-9.51
        eraAs = era5s.sel(time=era5s['time.hour'] == erah, latitude=slice(-8.5, -9),
                          longitude=slice(-73.6, -73.1)).mean(['latitude', 'longitude']).squeeze()
        eraVs = era5s.sel(time=era5s['time.hour'] == erah, latitude=slice(-8.5, -9),
                          longitude=slice(-75.6 - 75.1)).mean(['latitude', 'longitude']).squeeze()

        chirp = chirps.sel(time=date).squeeze()
        # This masks out the data which is not in the region
        gvalley = gridsat.salem.roi(shape=sdf)
        gvalley1 = gridsat1.salem.roi(shape=sdf)
        gvalley2 = gridsat2.salem.roi(shape=sdf)

        cvalley = chirp.salem.roi(shape=sdf)
        gA = gridsat.sel(lat=slice(-9, -8.5), lon=slice(-73.6, -73.1))
        gV = gridsat.sel(lat=slice(-9, -8.5), lon=slice(-75.6, -75.1))

        gA1 = gridsat1.sel(lat=slice(-9, -8.5), lon=slice(-73.6, -73.1))
        gV1 = gridsat1.sel(lat=slice(-9, -8.5), lon=slice(-75.6, -75.1))

        gA2 = gridsat2.sel(lat=slice(-9, -8.5), lon=slice(-73.6, -73.1))
        gV2 = gridsat2.sel(lat=slice(-9, -8.5), lon=slice(-75.6, -75.1))

        cA = chirp.sel(latitude=slice(-9, -8.5), longitude=slice(-73.6, -73.1))
        cV = chirp.sel(latitude=slice(-9, -8.5), longitude=slice(-75.6, -75.1))

        # ipdb.set_trace()

        # with warnings.catch_warnings():
        # warnings.simplefilter("ignore", category=RuntimeWarning)

        ####MCS LT
        MCScover = np.sum(gvalley.values < -10) / np.sum(np.isfinite(gvalley.values))
        MCSmean = np.nanmean(gvalley.values[gvalley.values < 0])
        MCSmin = np.nanmin(gvalley.values)

        MCScoverA = np.sum(gA.values < -10) / np.sum(np.isfinite(gA.values))
        MCSmeanA = np.nanmean(gA.values[gA.values < 0])
        MCSminA = np.nanmin(gA.values)

        MCScoverV = np.sum(gV.values < -10) / np.sum(np.isfinite(gV.values))
        MCSmeanV = np.nanmean(gV.values[gV.values < 0])
        MCSminV = np.nanmin(gV.values)

        CHIRPSmean = np.round(np.nanmean(cvalley.values[cvalley.values > 0.1]), decimals=2)
        CHIRPSmax = np.round(np.nanmax(cvalley.values), decimals=2)

        CHIRPSmeanA = np.round(np.nanmean(cA.values[cA.values > 0.1]), decimals=2)
        CHIRPSmaxA = np.round(np.nanmax(cA.values), decimals=2)

        CHIRPSmeanV = np.round(np.nanmean(cV.values[cV.values > 0.1]), decimals=2)
        CHIRPSmaxV = np.round(np.nanmax(cV.values), decimals=2)

        ####MCS UTC
        MCScover1 = np.sum(gvalley1.values < -10) / np.sum(np.isfinite(gvalley1.values))
        MCSmean1 = np.nanmean(gvalley1.values[gvalley1.values < 0])
        MCSmin1 = np.nanmin(gvalley1.values)

        MCScoverA1 = np.sum(gA1.values < -10) / np.sum(np.isfinite(gA1.values))
        MCSmeanA1 = np.nanmean(gA1.values[gA1.values < 0])
        MCSminA1 = np.nanmin(gA1.values)

        MCScoverV1 = np.sum(gV1.values < -10) / np.sum(np.isfinite(gV1.values))
        MCSmeanV1 = np.nanmean(gV1.values[gV1.values < 0])
        MCSminV1 = np.nanmin(gV1.values)

        ####Cloud UTC
        MCScover2 = np.sum(gvalley2.values < -10) / np.sum(np.isfinite(gvalley2.values))
        MCSmean2 = np.nanmean(gvalley2.values[gvalley2.values < 0])
        MCSmin2 = np.nanmin(gvalley2.values)

        MCScoverA2 = np.sum(gA2.values < -10) / np.sum(np.isfinite(gA2.values))
        MCSmeanA2 = np.nanmean(gA2.values[gA2.values < 0])
        MCSminA2 = np.nanmin(gA2.values)

        MCScoverV2 = np.sum(gV2.values < -10) / np.sum(np.isfinite(gV2.values))
        MCSmeanV2 = np.nanmean(gV2.values[gV2.values < 0])
        MCSminV2 = np.nanmin(gV2.values)

        #     plt.figure()
        #     plt.imshow(cvalley)

        # Huaraz
        # ipdb.set_trace()
        df.loc[dids, 'q550'] = eraH['q'].sel(level=550).values
        df.loc[dids, 'q650'] = eraH['q'].sel(level=650).values
        df.loc[dids, 'q200'] = eraH['q'].sel(level=200).values
        df.loc[dids, 't550'] = eraH['t'].sel(level=550).values
        df.loc[dids, 't650'] = eraH['t'].sel(level=650).values
        df.loc[dids, 'u550'] = eraH['u'].sel(level=550).values
        df.loc[dids, 'r650'] = eraH['r'].sel(level=650).values
        df.loc[dids, 'tcwv'] = eraHs['tcwv'].values
        df.loc[dids, 'shear'] = eraH['u'].sel(level=200).values - eraH['u'].sel(level=550).values  # 200-550
        df.loc[dids, 't2'] = eraHs['t2m'].values
        df.loc[dids, 'cape'] = eraHs['cape'].values
        df.loc[dids, 'cin'] = eraHs['cin'].values
        df.loc[dids, 'u200'] = eraH['u'].sel(level=200).values
        df.loc[dids, 'v200'] = eraH['v'].sel(level=200).values
        df.loc[dids, 'u550'] = eraH['u'].sel(level=550).values
        df.loc[dids, 'v550'] = eraH['v'].sel(level=550).values

        # Amazon [-73.6, -72.6, -10.2, -8.6]
        df.loc[dids, 'tcwvA'] = eraAs['tcwv'].values
        df.loc[dids, 'q550A'] = eraA['q'].sel(level=550).values
        df.loc[dids, 'q850A'] = eraA['q'].sel(level=850).values
        df.loc[dids, 'q200A'] = eraA['q'].sel(level=200).values
        df.loc[dids, 't550A'] = eraA['t'].sel(level=550).values
        df.loc[dids, 't850A'] = eraA['t'].sel(level=850).values
        df.loc[dids, 'shearA'] = eraA['u'].sel(level=550).values - eraA['u'].sel(level=850).values  # 550-850
        df.loc[dids, 'u200A'] = eraA['u'].sel(level=200).values
        df.loc[dids, 'v200A'] = eraA['v'].sel(level=200).values
        df.loc[dids, 'u550A'] = eraA['u'].sel(level=550).values
        df.loc[dids, 'v550A'] = eraA['v'].sel(level=550).values
        df.loc[dids, 'u850A'] = eraA['u'].sel(level=850).values
        df.loc[dids, 'v850A'] = eraA['v'].sel(level=850).values
        df.loc[dids, 't2A'] = eraAs['t2m'].values
        df.loc[dids, 'capeA'] = eraAs['cape'].values
        df.loc[dids, 'cinA'] = eraAs['cin'].values
        df.loc[dids, 'r850A'] = eraA['r'].sel(level=850).values

        # Valley [-75.6, -74.6, -10.2, -8.6]
        df.loc[dids, 'tcwvV'] = eraVs['tcwv'].values
        df.loc[dids, 'q550V'] = eraV['q'].sel(level=550).values
        df.loc[dids, 'q850V'] = eraV['q'].sel(level=850).values
        df.loc[dids, 'q200V'] = eraV['q'].sel(level=200).values
        df.loc[dids, 't550V'] = eraV['t'].sel(level=550).values
        df.loc[dids, 't850V'] = eraV['t'].sel(level=850).values
        df.loc[dids, 'shearV'] = eraV['u'].sel(level=550).values - eraV['u'].sel(level=850).values  # 550=850
        df.loc[dids, 'u200V'] = eraV['u'].sel(level=200).values
        df.loc[dids, 'v200V'] = eraV['v'].sel(level=200).values
        df.loc[dids, 'u850V'] = eraV['u'].sel(level=850).values
        df.loc[dids, 'v850V'] = eraV['v'].sel(level=850).values
        df.loc[dids, 'u550V'] = eraV['u'].sel(level=550).values
        df.loc[dids, 'v550V'] = eraV['v'].sel(level=550).values
        df.loc[dids, 't2V'] = eraVs['t2m'].values
        df.loc[dids, 'capeV'] = eraVs['cape'].values
        df.loc[dids, 'cinV'] = eraVs['cin'].values
        df.loc[dids, 'r850V'] = eraV['r'].sel(level=850).values

        df.loc[dids, 'MCSmean'] = np.round(MCSmean, decimals=2)
        df.loc[dids, 'MCSmin'] = np.round(MCSmin, decimals=2)
        df.loc[dids, 'MCScover'] = np.round(MCScover, decimals=2)
        df.loc[dids, 'MCSmeanUTC'] = np.round(MCSmean1, decimals=2)
        df.loc[dids, 'MCSminUTC'] = np.round(MCSmin1, decimals=2)
        df.loc[dids, 'MCScoverUTC'] = np.round(MCScover1, decimals=2)
        df.loc[dids, 'CloudmeanUTC'] = np.round(MCSmean2, decimals=2)
        df.loc[dids, 'CloudminUTC'] = np.round(MCSmin2, decimals=2)
        df.loc[dids, 'CloudcoverUTC'] = np.round(MCScover2, decimals=2)
        df.loc[dids, 'CHIRPSmean'] = np.round(CHIRPSmean, decimals=2)
        df.loc[dids, 'CHIRPSmax'] = np.round(CHIRPSmax, decimals=2)

        df.loc[dids, 'MCSmeanA'] = np.round(MCSmeanA, decimals=2)
        df.loc[dids, 'MCSminA'] = np.round(MCSminA, decimals=2)
        df.loc[dids, 'MCScoverA'] = np.round(MCScoverA, decimals=2)
        df.loc[dids, 'MCSmeanUTCA'] = np.round(MCSmeanA1, decimals=2)
        df.loc[dids, 'MCSminUTCA'] = np.round(MCSminA1, decimals=2)
        df.loc[dids, 'MCScoverUTCA'] = np.round(MCScoverA1, decimals=2)
        df.loc[dids, 'CloudmeanUTCA'] = np.round(MCSmeanA2, decimals=2)
        df.loc[dids, 'CloudminUTCA'] = np.round(MCSminA2, decimals=2)
        df.loc[dids, 'CloudcoverUTCA'] = np.round(MCScoverA2, decimals=2)
        df.loc[dids, 'CHIRPSmeanA'] = np.round(CHIRPSmeanA, decimals=2)
        df.loc[dids, 'CHIRPSmaxA'] = np.round(CHIRPSmaxA, decimals=2)

        df.loc[dids, 'MCSmeanV'] = np.round(MCSmeanV, decimals=2)
        df.loc[dids, 'MCSminV'] = np.round(MCSminV, decimals=2)
        df.loc[dids, 'MCScoverV'] = np.round(MCScoverV, decimals=2)
        df.loc[dids, 'MCSmeanUTCV'] = np.round(MCSmeanV1, decimals=2)
        df.loc[dids, 'MCSminUTCV'] = np.round(MCSminV1, decimals=2)
        df.loc[dids, 'MCScoverUTCV'] = np.round(MCScoverV1, decimals=2)
        df.loc[dids, 'CloudmeanUTCV'] = np.round(MCSmeanV2, decimals=2)
        df.loc[dids, 'CloudminUTCV'] = np.round(MCSminV2, decimals=2)
        df.loc[dids, 'CloudcoverUTCV'] = np.round(MCScoverV2, decimals=2)
        df.loc[dids, 'CHIRPSmeanV'] = np.round(CHIRPSmeanV, decimals=2)
        df.loc[dids, 'CHIRPSmaxV'] = np.round(CHIRPSmaxV, decimals=2)

        df.loc[dids, 'month'] = dt.month
        df.loc[dids, 'year'] = dt.year
        df.loc[dids, 'hour'] = dt.hour

        # ipdb.set_trace()

        del era5
        del era5s
        del eraH
        del eraV
        del eraA
        del eraHs
        del eraVs
        del eraAs
        del gridsat
        del gridsat1
        del gridsat2
        del chirp

    return df
