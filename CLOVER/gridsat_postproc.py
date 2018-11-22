import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy
from utils import u_plot as up
import cartopy.crs as ccrs
import os
import matplotlib as mpl
import pickle as pkl
from utils import constants as cnst
from scipy.ndimage.measurements import label


def month():
    y1 = 1982
    y2 =2017#2017
    years = list(range(1983,1985)) #+ list(range(2004,2014))

    msg_folder = cnst.GRIDSAT
    fname='aggs/gridsat_WA_-70_monthly.nc'

    if not os.path.isfile(msg_folder + fname):
        da = None
        for y in years:
            y = str(y)
            da1 = xr.open_dataset(cnst.GRIDSAT+'gridsat_WA_' + y + '.nc')
            print('Doing '+y)
            da1['tir'] = da1['tir'].where(da1['tir'] <= -70)

            da1 = da1.resample(time='m').mean('time')
            try:
                da = xr.concat([da, da1], 'time')
            except TypeError:
                da = da1.copy()


        enc = {'tir': {'complevel': 5,  'zlib': True}}
        da.to_netcdf(msg_folder + fname, encoding=enc)


    else:
        ds = xr.open_dataset(msg_folder + fname)
        da = ds['tir']
    da.values[da.values==0]=np.nan
    da.sel(lat=slice(11, 20))
    mean = da['tir'].mean(dim=['lat', 'lon'])

    mean.plot()

    plt.savefig('/users/global/cornkle/figs/CLOVER/GRIDSAT_cold_clouds/tests/trend_mcs.png', dpi=300)


def month_mean():

    years = list(range(1983,2018))

    msg_folder = cnst.GRIDSAT
    fname = 'aggs/gridsat_WA_-50_monthly_mean_SA.nc'

    if not os.path.isfile(msg_folder + fname):
        da = None
        da_box = None
        for y in years:
            y = str(y)
            da1 = xr.open_dataset(cnst.GRIDSAT + 'gridsat_WA_-50_' + y + '.nc')
            print('Doing ' + y)
            da1['tir'] = da1['tir'].where((da1['tir'] <= -50) & (da1['tir'] >= -108) )
            #da1['tir'].values[da1['tir'].values < -70] = 1

            da_res = da1.resample(time='m').mean('time')
            WA_box = [-13,13,4.5,8]
            SA_box = [20,30,-30,-10]
            boxed = da1['tir'].sel(lat=slice(SA_box[2],SA_box[3]), lon=slice(SA_box[0],SA_box[1])).resample(time='m').mean()

            try:
                da = xr.concat([da, da_res], 'time')
            except TypeError:
                da = da_res.copy()

            try:
                da_box = xr.concat([da_box, boxed], 'time')
            except TypeError:
                da_box = boxed.copy()


        # pkl.dump(np.array(boxed),
        #          open('/users/global/cornkle/data/CLOVER/saves/box_13W-13E-4-8N_meanT-50_from5000km2.p',
        #               'wb'))

        enc = {'tir': {'complevel': 5, 'zlib': True}}
        pdb.set_trace()
        da_box.to_netcdf(msg_folder + 'box_13W-13E-4-8N_meanT-50_from5000km2_SA.nc')

        da.to_netcdf(msg_folder + fname, encoding=enc)



def month_mean_climatology():

    years = list(range(1983,2018))

    msg_folder = cnst.GRIDSAT
    fname = 'aggs/gridsat_WA_-50_monthly_mean.nc'

    if not os.path.isfile(msg_folder + fname):
        da = None
        da_box = None
        for y in years:
            y = str(y)
            da1 = xr.open_dataset(cnst.GRIDSAT + 'gridsat_WA_-50_' + y + '.nc')
            print('Doing ' + y)
            da1['tir'] = da1['tir'].where((da1['tir'] <= -50) & (da1['tir'] >= -108) )
            #da1['tir'].values[da1['tir'].values < -70] = 1




            da_res = da1.resample(time='m').mean('time')

            boxed = da1['tir'].sel(lat=slice(4.5,8), lon=slice(-13,13)).resample(time='m').mean()

            try:
                da = xr.concat([da, da_res], 'time')
            except TypeError:
                da = da_res.copy()

            try:
                da_box = xr.concat([da_box, boxed], 'time')
            except TypeError:
                da_box = boxed.copy()


        # pkl.dump(np.array(boxed),
        #          open('/users/global/cornkle/data/CLOVER/saves/box_13W-13E-4-8N_meanT-50_from5000km2.p',
        #               'wb'))
        pdb.set_trace()
        enc = {'tir': {'complevel': 5, 'zlib': True}}

        da_box.to_netcdf(msg_folder + 'box_13W-13E-4-8N_meanT-50_from5000km2.nc')

        da.to_netcdf(msg_folder + fname, encoding=enc)


def month_count():

    years = list(range(1983, 2018))

    msg_folder = cnst.GRIDSAT
    fname = 'aggs/gridsat_WA_-70_monthly_count.nc'

    if not os.path.isfile(msg_folder + fname):
        da = None
        for y in years:
            y = str(y)
            da1 = xr.open_dataset(cnst.GRIDSAT + 'gridsat_WA_' + y + '.nc')
            print('Doing ' + y)
            da1['tir'] = da1['tir'].where((da1['tir'] <= -70) & (da1['tir'] >= -108))
            da1['tir'].values[da1['tir'].values < -70] = 1

            da1 = da1.resample(time='m').sum('time')
            try:
                da = xr.concat([da, da1], 'time')
            except TypeError:
                da = da1.copy()

        enc = {'tir': {'complevel': 5, 'zlib': True}}
        da.to_netcdf(msg_folder + fname, encoding=enc)


    # else:
    #     ds = xr.open_dataset(msg_folder + fname)
    #     da = ds['tir']
    # pdb.set_trace()
    # da.values[da.values == 0] = np.nan
    # da.sel(lat=slice(11, 20))
    # mean = da['t'].mean(dim=['lat', 'lon'])
    #
    # mean.plot()
    #
    # plt.savefig('/users/global/cornkle/figs/CLOVER/GRIDSAT_cold_clouds/tests/trend_mcs-50.png', dpi=300)
