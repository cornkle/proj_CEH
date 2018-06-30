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
from scipy.ndimage.measurements import label


def monthly(years):

    for y in years:
        y = str(y)
        da = xr.open_dataset('/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gridsat_WA_' + y + '.nc')
        da = da['t']
        da = da.where(da <= -80)

        month = da.groupby('time.month').mean(dim='time')

        simple = month.plot(x='lon', y='lat', col='month', col_wrap=3, cmap='viridis_r', transform=ccrs.PlateCarree(),
                            subplot_kws={'projection': ccrs.PlateCarree()}, vmin=-85)

        for ax in simple.axes.flat:
            ax.coastlines()
            ax.gridlines()
            ax.set_extent([-17.5, 30, -6, 20])
            ax.set_aspect('equal', 'box-forced')

        plt.suptitle(y)
        plt.savefig('/users/global/cornkle/figs/cold_clouds/monthly_Tmean<-80/monthly_'+y+'.png')
        plt.close('all')



def climatology_month():
    years = np.arange(2005,2016)#2017)

    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'
    fname='gridsat_WA_cold_climatology_mean.nc'

    if not os.path.isfile(msg_folder + fname):
        da = xr.open_dataset('/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gridsat_WA_' + str(2004) + '.nc')
        da = da['t']
        da = da.where(da <= -60)

        month = da.groupby('time.month').mean(dim='time')
        for y in years:
            y = str(y)
            da = xr.open_dataset('/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gridsat_WA_' + y + '.nc')
            da = da['t']
            da = da.where(da <= -60)

            month = month + da.groupby('time.month').mean(dim='time')

        month = month / (len(years)+1)

        enc = {'t': {'complevel': 5,  'zlib': True}}
        month.to_netcdf(msg_folder + fname, encoding=enc)


    else:
        ds = xr.open_dataset(msg_folder + fname)
        month = ds['t']
    month.values[month.values==0]=np.nan

    simple = month.plot(x='lon', y='lat', col='month', col_wrap=4, cmap='inferno', transform=ccrs.PlateCarree(),
                        subplot_kws={'projection': ccrs.PlateCarree()}, vmax=-65, vmin=-75, levels = 6)

    for ax in simple.axes.flat:
        ax.coastlines()
        ax.gridlines()
        ax.set_extent([-17.5, 30, -6, 20])
        ax.set_aspect('equal', 'box-forced')


    plt.savefig('/users/global/cornkle/VERA/plots/leeds_june_2017/mean_t.png', dpi=300)


def month():
    y1 = 1982
    y2 =2017#2017
    years = np.arange(y1+1,y2)#2017)

    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'
    fname='gridsat_WA_-70_monthly.nc'

    if not os.path.isfile(msg_folder + fname):
        da = xr.open_dataset('/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gridsat_WA_' + str(y1) + '.nc')
        da['t'] = da['t'].where(da['t'] <= -70)


        da = da.resample('m', dim='time', how='mean')

        for y in years:
            y = str(y)
            da1 = xr.open_dataset('/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gridsat_WA_' + y + '.nc')
            print('Doing '+y)
            da1['t'] = da1['t'].where(da1['t'] <= -70)

            da1 = da1.resample('m', dim='time', how='mean')

            da = xr.concat([da, da1], 'time')
            da1.close()

        enc = {'t': {'complevel': 5,  'zlib': True}}
        da.to_netcdf(msg_folder + fname, encoding=enc)


    else:
        ds = xr.open_dataset(msg_folder + fname)
        da = ds['t']
    da.values[da.values==0]=np.nan
    da.sel(lat=slice(11, 20))
    mean = da['t'].mean(dim=['lat', 'lon'])

    mean.plot()

    plt.savefig('/users/global/cornkle/VERA/plots/leeds_june_2017/trend_mcs.png', dpi=300)


def month_count():
    y1 = 1982
    y2 = 2017  # 2017
    years = np.arange(y1 + 1, y2)  # 2017)

    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'
    fname = 'gridsat_WA_-40_monthly_count.nc'

    if not os.path.isfile(msg_folder + fname):
        da = xr.open_dataset('/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gridsat_WA_' + str(y1) + '.nc')

        da['t'] = da['t'].where(da['t'] <= -40)
        da['t'].values[da['t'].values <= -40] = 1

        da = da.resample('m', dim='time', how='sum')

        for y in years:
            y = str(y)
            da1 = xr.open_dataset('/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gridsat_WA_' + y + '.nc')
            print('Doing ' + y)
            da1['t'] = da1['t'].where(da1['t'] <= -40)
            da1['t'].values[da1['t'].values <= -40] = 1

            da1 = da1.resample('m', dim='time', how='sum')

            da = xr.concat([da, da1], 'time')
            da1.close()

        enc = {'t': {'complevel': 5, 'zlib': True}}
        da.to_netcdf(msg_folder + fname, encoding=enc)


    else:
        ds = xr.open_dataset(msg_folder + fname)
        da = ds['t']
    da.values[da.values == 0] = np.nan
    da.sel(lat=slice(11, 20))
    mean = da['t'].mean(dim=['lat', 'lon'])

    mean.plot()

    plt.savefig('/users/global/cornkle/VERA/plots/leeds_june_2017/trend_mcs.png', dpi=300)


def hourly_count():
    y1 = 1982
    y2 = 1984  # 2017
    years = np.arange(y1 + 1, y2)  # 2017)

    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/yearly_files/'
    fname = 'gridsat_WA_-70_hourly_count.nc'

    if not os.path.isfile(msg_folder + fname):
        da = xr.open_dataset(msg_folder+'gridsat_WA_' + str(y1) + '.nc')

        da['t'] = da['t'].where(da['t'] <= -70)
        da['t'].values[da['t'].values <= -70] = 1
        da = da['t']

        da = da[(da['time.month']>=6) & (da['time.month']<=9)]

        da= da.groupby('time.hour').sum(dim='time.hour')

        for y in years:
            y = str(y)
            da1 = xr.open_dataset('/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/gridsat_WA_hours_' + y + '.nc')
            print('Doing ' + y)
            da1['t'] = da1['t'].where(da['t'] <= -70)
            da1['t'].values[da1['t'].values <= -70] = 1

            da1 = da1[(da1['time.month'] >= 6) & (da1['time.month'] <= 9)]

            da1 = da1.groupby('time.hour').sum(dim='time.hour')

            pdb.set_trace()

            da = xr.concat([da, da1], 'time')
            da1.close()

        enc = {'t': {'complevel': 5, 'zlib': True}}
        da.to_netcdf(msg_folder + fname, encoding=enc)


    else:
        ds = xr.open_dataset(msg_folder + fname)
        da = ds['t']
    da.values[da.values == 0] = np.nan
    da.sel(lat=slice(11, 20))
    mean = da['t'].mean(dim=['lat', 'lon'])

    mean.plot()

    plt.savefig('/users/global/cornkle/VERA/plots/leeds_june_2017/trend_mcs.png', dpi=300)


def timeline_trend_count():
    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'
    fname = 'gridsat_WA_-70_monthly_count.nc'

    da = xr.open_dataarray(msg_folder + fname)
    da = da.sel(lat=slice(8,12), lon=slice(-17, 20))
    #da=da.sel(lat=slice(5,10))
    #da[da==0]=np.nan
    mean = da.mean(dim=['lat', 'lon'])
    #mean = mean[(mean['time.month']==8)]
    f= plt.figure(figsize=(10,6))
    for i in range(6,9):
        bla = mean[(mean['time.month'] == i)]
        bla.plot(label=str(i), marker='o')
    plt.title('Average number of pixels <= -70C, 11-18N')
    plt.legend()
   # plt.ylim(-76,-72)

def timeline_trend_mean():
    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'
    fname = 'gridsat_WA_-70_monthly.nc'

    da = xr.open_dataarray(msg_folder + fname)
    da=da.sel(lat=slice(5,7), lon=slice(-17,20))
    da[da==0]=np.nan
    mean = da.mean(dim=['lat', 'lon'])
    #mean = mean[(mean['time.month']==8)]
    f= plt.figure(figsize=(10,6))
    for i in range(4,6):
        bla = mean[(mean['time.month'] == i)]
        bla.plot(label=str(i), marker='o')
    plt.title('Monthly mean temperature of pixels <= -40C, 11-18N')
    plt.legend()
    plt.ylim(-78,-71)


def trend_map():

    def linear_trend(x):
        pf = np.polyfit(np.arange(len(x)), x, 1)
        # we need to return a dataarray or else xarray's groupby won't be happy
        return xr.DataArray(pf[0])

    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'
    fpath = '/users/global/cornkle/figs/gap_filling_Tgrad/gridsat/'

    fname = 'gridsat_WA_-70_monthly_count.nc'

    dam = xr.open_dataarray(msg_folder + fname)
    dam = dam.sel(lon=slice(-18,51), lat=slice(-37, 36))
    lons = dam.lon
    lats = dam.lat
    months = np.arange(1, 13)
    #da[da == 0] = np.nan
    for m in months:
        da = dam[(dam['time.month']==m) & (dam['time.year']>1982)]
        da = da.groupby('time.day').sum(axis=0)
        da = da.groupby('time.year').mean(axis=0)  # average daily frequency per year
        # stack lat and lon into a single dimension called allpoints
        stacked = da.stack(allpoints=['lat', 'lon'])
        # apply the function over allpoints to calculate the trend at each point

        trend = stacked.groupby('allpoints').apply(linear_trend)
        # unstack back to lat lon coordinates
        trend_unstacked = trend.unstack('allpoints')

        trend_unstacked = trend_unstacked * 10.
        da2 = xr.DataArray(trend_unstacked, coords=[lats, lons], dims=['latitude', 'longitude'])

        fp = fpath + 'ttrend_' + str(m).zfill(2) + '.png'

        up.quick_map_salem(da2, vmin=-0.4, vmax=0.4, cmap='RdBu_r', save=fp)



def t_ratio():

    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/'
    fname = 'gridsat_WA_-70_monthly_count.nc'
    da70 = xr.open_dataarray(msg_folder + fname)
    fname = 'gridsat_WA_-40_monthly_count.nc'
    da40 = xr.open_dataarray(msg_folder + fname)

    da70.values[da70.values == 0] = np.nan
    da40.values[da40.values == 0] = np.nan
    da70.sel(lat=slice(11, 20))
    da40.sel(lat=slice(11, 20))
    #ratio = da70/da40

    msg40 = da40[(da40['time.year'] >= 2007) & (da40['time.year'] <= 2015)]
    msg70 = da70[(da70['time.year'] >= 2007) & (da70['time.year'] <= 2015)]

    mfg40 = da40[(da40['time.year'] >= 1990) & (da40['time.year'] <= 2005)]
    mfg70 = da70[(da70['time.year'] >= 1990) & (da70['time.year'] <= 2005)]

    # f = plt.figure(figsize=(10, 6))
    # ax = f.add_subplot(111)

    msg40 = msg40.groupby('time.season').sum(dim='time')
    msg70 = msg70.groupby('time.season').sum(dim='time')

    mfg40 = mfg40.groupby('time.season').sum(dim='time')
    mfg70 = mfg70.groupby('time.season').sum(dim='time')

    msg_ratio = msg70/msg40*100
    mfg_ratio = mfg70 / mfg40*100

    simple = msg_ratio.plot(x='lon', y='lat', col='season', col_wrap=2, cmap='viridis', transform=ccrs.PlateCarree(),
                        subplot_kws={'projection': ccrs.PlateCarree()}, vmax=20)

    for ax in simple.axes.flat:
        ax.coastlines()
        ax.gridlines()
        ax.set_extent([-17.5, 30, -6, 20])
        ax.set_aspect('equal', 'box-forced')
        ax.set_title('MSG')

    simple = mfg_ratio.plot(x='lon', y='lat', col='season', col_wrap=2, cmap='viridis', transform=ccrs.PlateCarree(),
                        subplot_kws={'projection': ccrs.PlateCarree()}, vmax=20)

    for ax in simple.axes.flat:
        ax.coastlines()
        ax.gridlines()
        ax.set_extent([-17.5, 30, -6, 20])
        ax.set_aspect('equal', 'box-forced')
        ax.set_title('MFG')

    ratio = (msg_ratio-mfg_ratio)


    simple = ratio.plot(x='lon', y='lat', col='season', col_wrap=2, cmap='RdBu', transform=ccrs.PlateCarree(),
                        subplot_kws={'projection': ccrs.PlateCarree()}, vmax=10, vmin=-10, levels=[-15,-10,-5,5,10,15])

    for ax in simple.axes.flat:
        ax.coastlines()
        ax.gridlines()
        ax.set_extent([-17.5, 30, -6, 20])
        ax.set_aspect('equal', 'box-forced')



def size_trend():

    msg_folder = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/yearly_files/'
    data = xr.open_mfdataset(msg_folder + 'gridsat*.nc')

    cut = data.sel(lat=slice(10,17), lon=slice(-17,-10))
    cut = cut.isel(time= ((cut['time.year']>1984) & (cut['time.month']==8)))
    cut=cut['t']

    dic= {}
    for p in np.arange(1985,2017,1):
        dic[p] = []

    def mcs_find(image, thresh=None):
        if not thresh:
            print('Give threshold')
            return

        image[image > thresh] = 0
        image[image <= thresh] = 1
        image[np.isnan(image)] = 0

        if np.sum(image<10):
            return []

        labels, numL = label(image)

        ret = []

        for l in np.unique(labels):
            if l == 0:
                continue

            blob = np.sum(labels == l)

            pdb.set_trace()

            if np.sum(len(blob[0])) < 100:  # at least 1000m2
                continue

            ret.append(blob*49)

        return ret

    for i in np.arange(cut.shape[0]):

        ret = mcs_find(cut[i,:,:].values, thresh=-40)
        if ret == []:
            continue
        pdb.set_trace()
        dic[cut['time.year']].append(ret)

    pdb.set_trace()

    for d in dic:
        d = [item for sublist in d for item in sublist]




