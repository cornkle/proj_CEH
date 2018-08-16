import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import salem
from utils import u_plot as up, u_mann_kendall as mk


def t_trend():
    #file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = '/localscratch/wllf030/cornkle/ERA-I/monthly/old/ERA-Int-Monthly-2mTemp.nc'

    fpath = '/users/global/cornkle/figs/CLOVER/months/'


    dam = xr.open_dataarray(file)
    months = np.arange(1,13)

    for m in months:

        da = dam[(dam['time.month']==m)]
        da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))
        da = da.groupby('time.year').mean(axis=0)

        lons = da.longitude
        lats = np.flip(da.latitude.values, axis=0)

        # define a function to compute a linear trend of a timeseries
        def linear_trend(x):

            #pf = np.polyfit(np.arange(len(x)), x, 1)
            pf, slope, int, p, ind = mk.test(np.arange(len(x)),x.squeeze().values, eps=0.001, alpha=0.05, Ha='upordown')

            # we need to return a dataarray or else xarray's groupby won't be happy

            if ind == 1:
                issig = slope
            else:
                issig = np.nan

            return xr.DataArray(issig, )

        # stack lat and lon into a single dimension called allpoints
        stacked = da.stack(allpoints=['latitude','longitude'])
        # apply the function over allpoints to calculate the trend at each point
        trend = stacked.groupby('allpoints').apply(linear_trend)
        # unstack back to lat lon coordinates
        trend_unstacked = trend.unstack('allpoints')

        trend_unstacked = trend_unstacked*10. # warming over decade
        da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

        fp = fpath + 'ttrend_'+str(m).zfill(2)+'.png'

        up.quick_map_salem(da2, vmin=-0.4, vmax=0.4, cmap='RdBu_r', save=fp)  #

        plt.close('all')


def t_trend_slice():
    #file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = '/localscratch/wllf030/cornkle/ERA-I/monthly/old/ERA-Int-Monthly-2mTemp.nc'

    fpath = '/users/global/cornkle/figs/CLOVER/months/'


    dam = xr.open_dataarray(file)
    lower = 9
    higher = 11

    da = dam[(dam['time.month']>=lower) & (dam['time.month']<=higher)]
    da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))
    da = da.groupby('time.year').mean(axis=0)

    lons = da.longitude
    lats = np.flip(da.latitude.values, axis=0)

    # define a function to compute a linear trend of a timeseries
    def linear_trend(x):

        #pf = np.polyfit(np.arange(len(x)), x, 1)
        pf, slope, int, p, ind = mk.test(np.arange(len(x)),x.squeeze().values, eps=0.001, alpha=0.01, Ha='upordown')

        # we need to return a dataarray or else xarray's groupby won't be happy

        if ind == 1:
            issig = slope
        else:
            issig = np.nan

        return xr.DataArray(issig, )

    # stack lat and lon into a single dimension called allpoints
    stacked = da.stack(allpoints=['latitude','longitude'])
    # apply the function over allpoints to calculate the trend at each point
    trend = stacked.groupby('allpoints').apply(linear_trend)
    # unstack back to lat lon coordinates
    trend_unstacked = trend.unstack('allpoints')

    trend_unstacked = trend_unstacked*10. # warming over decade
    da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

    fp = fpath + 'ttrend_'+str(lower).zfill(2)+'-'+str(higher).zfill(2)+'.png'

    up.quick_map_salem(da2, vmin=-0.4, vmax=0.4, cmap='RdBu_r', save=fp)  #

    plt.close('all')




def t_mean():
    # file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = '/users/global/cornkle/data/ERA-I monthly/ERA-Int-Monthly-2mTemp.nc'

    fpath = '/users/global/cornkle/figs/gap_filling_Tgrad/months/'

    dam = xr.open_dataarray(file)
    months = np.arange(1, 13)

    for m in months:
        da = dam[(dam['time.month'] == m)]
        da = da.sel(longitude=slice(-18, 51), latitude=slice(36, -37))
        da = da.mean(axis=0)-273.15

        fp = fpath + 'tmean_' + str(m).zfill(2) + '.png'

        up.quick_map_salem(da, levels=np.arange(20,41,2), cmap='jet', save=fp)



def gridsat_70trend():

    fpath = '/users/global/cornkle/mymachine/GRIDSAT/MCS18/'
    out = '/users/global/cornkle/figs/CLOVER/GRIDSAT_cold_clouds/tests'


    dam = xr.open_mfdataset(fpath+'*.nc')
    dam = dam['tir']
    months = np.arange(1,13)

    for m in months:

        da = dam[(dam['time.month']==m)]
        #da = da.sel(longitude=slice(-18,51), latitude=slice(-37,36))
        da.values[da.values<-40] = 1
        da = da.groupby('time.year').mean(axis=0)

        lons = da.longitude
        lats = np.flip(da.latitude.values, axis=0)

        # define a function to compute a linear trend of a timeseries
        def linear_trend(x):

            #pf = np.polyfit(np.arange(len(x)), x, 1)
            pf, slope, int, p, ind = mk.test(np.arange(len(x)),x.squeeze().values, eps=0.001, alpha=0.05, Ha='upordown')

            # we need to return a dataarray or else xarray's groupby won't be happy

            if ind == 1:
                issig = slope
            else:
                issig = np.nan

            return xr.DataArray(issig, )

        # stack lat and lon into a single dimension called allpoints
        stacked = da.stack(allpoints=['lat','lon'])
        # apply the function over allpoints to calculate the trend at each point
        trend = stacked.groupby('allpoints').apply(linear_trend)
        # unstack back to lat lon coordinates
        trend_unstacked = trend.unstack('allpoints')

        trend_unstacked = trend_unstacked*10. # warming over decade
        da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['lat', 'lon'])

        fp = fpath + 'ttrend_'+str(m).zfill(2)+'.png'

        up.quick_map_salem(da2, cmap='RdBu_r', save=fp)  #

        plt.close('all')