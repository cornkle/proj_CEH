import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import salem
from utils import u_plot as up, u_mann_kendall as mk
from utils import u_darrays


def t_trend():
    #file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'
    file = '/localscratch/wllf030/cornkle/ERA-I/monthly/monthly_1979-2017_srfc.nc'

    fpath = '/users/global/cornkle/figs/CLOVER/months/'


    dam = xr.open_dataset(file)
    dam = dam['t2m']
    months = np.arange(1,13)

    for m in months:

        da = dam[(dam['time.month']==m)]
        da = da.sel(longitude=slice(-18,51), latitude=slice(36, -37))
        da = da.groupby('time.year').mean(axis=0)

        lons = da.longitude
        lats = np.flip(da.latitude.values, axis=0)

        # stack lat and lon into a single dimension called allpoints
        stacked = da.stack(allpoints=['latitude','longitude'])
        # apply the function over allpoints to calculate the trend at each point
        trend = stacked.groupby('allpoints').apply(u_darrays.linear_trend)
        # unstack back to lat lon coordinates
        trend_unstacked = trend.unstack('allpoints')

        trend_unstacked = trend_unstacked*10. # warming over decade
        da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

        fp = fpath + 'ttrend_WA'+str(m).zfill(2)+'.png'

        up.quick_map_salem(da2, levels=np.arange(-0.5,0.5,0.1), cmap='RdBu_r', save=fp)  #

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




def t_trend_polyfit():
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



