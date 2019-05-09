import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy
from utils import u_plot as up, u_darrays
import cartopy.crs as ccrs
import os
import matplotlib as mpl
import pickle as pkl
from utils import constants as cnst, u_darrays, u_grid, u_plot
from scipy.ndimage.measurements import label
import ipdb


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

            da_res = da1.resample(time='m').count('time')
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


def test_rainfall():

    ts = xr.open_mfdataset(cnst.CHIRPS + "*.nc")

    res = ts.resample(time='m').max('time') # ts.where(ts['rainfall']>=20).notnull().resample(time='m').sum('time') #
    resmean = ts.resample(time='m').mean('time') #ts.where(ts['rainfall']>=0.1).notnull().resample(time='m').sum('time') #.count(ts['rainfall']>1,dim='time')

    trend = res['rainfall'][(res['time.month'] == 10) & (res['time.year'] >= 1990)].sel(longitude=slice(-12,12), latitude=slice(4.3,8.5))
    trendmean = resmean['rainfall'][resmean['time.month'] == 10 & (res['time.year'] >= 1990)].sel(longitude=slice(-12,12), latitude=slice(4.3,8.5))

    mapp = plt.imshow(trend[10,:,:], origin='lower', vmin=0, vmax=20)

    mtrend = trend.max(['latitude', 'longitude'])
    mmtrend = trendmean.max(['latitude', 'longitude'])


    plt.plot(trend['time.year'], mtrend)


def trend_rainfall():
    ts = xr.open_dataset(cnst.CHIRPS + "/global/chirps-v2.0.monthly.nc")
    WA = [-15,15,4,10]
    SA = [5,55,-40,0]
    mean_years = ts.sel(longitude=slice(SA[0], SA[1]), latitude=slice(SA[2], SA[3]))

    #mean_years = u_darrays.flip_lat(mean_years)
    mean_years = (mean_years['precip'])[(mean_years['time.month']>=11) | (mean_years['time.month']<=1)]


    mean_years = mean_years.groupby('time.year').mean('time')


    mask = np.where(np.isnan(mean_years.values))
    mean_y = u_grid.refactor_da(mean_years, 0.25)
    clim = mean_y.mean('time')
    mean_y.values[mask] = 0

    # stack lat and lon into a single dimension called allpoints
    datastacked = mean_y.stack(allpoints=['latitude', 'longitude'])

    # apply the function over allpoints to calculate the trend at each point
    print('Entering trend calc')

    alpha = 0.05
    # NaNs means there is not enough data, slope = 0 means there is no significant trend.

    dtrend = datastacked.groupby('allpoints').apply(u_darrays.linear_trend_mk, alpha=alpha, eps=0.0001,nb_missing=10)
    dtrend = dtrend.unstack('allpoints')
    sig=True
    if sig:
            (dtrend['slope'].values)[dtrend['ind'].values==0] = 0

    contour = {'data': clim, 'x': mean_y['latitude'], 'y': mean_y['longitude'], 'levels': np.arange(0, 100, 10),
               'cmap': 'PuOr'}

    u_plot.draw_map(dtrend['slope'], mean_y['longitude'], mean_y['latitude'], cmap='RdBu_r',contour=contour)