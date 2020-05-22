import pandas as pd
import numpy as np
import xarray as xr
import ipdb
import matplotlib.pyplot as plt
import cartopy
from utils import u_plot as up, u_darrays as uda
import cartopy.crs as ccrs
import os
import matplotlib as mpl
from utils import constants as cnst
from scipy.ndimage.measurements import label
import salem
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature


def monthly(years):

    for y in years:
        y = str(y)
        da = xr.open_dataset(cnst+'gridsat_WA_' + y + '.nc')
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
        plt.savefig('/users/global/cornkle/figs/CLOVER/GRIDSAT_cold_clouds/tests/monthly_'+y+'.png')
        plt.close('all')



def climatology_month(h):
    years = np.arange(1986,2018)#2017)

    htag = str(h).zfill(2)

    msg_folder = cnst.GRIDSAT_PERU
    fname='aggs/gridsat_WA_-40C_climatology_mean_'+htag+'UTC.nc'
    cnt = 1
    if not os.path.isfile(msg_folder + fname):
        da = xr.open_dataset(msg_folder+'gridsat_WA_-40_5000km2_13-19UTC' + str(1985) + '.nc')
        da = da.sel(time=(da['time.hour'] == 21))
        da = da['tir']/100
        da = da.where(da <= -40)
        #ipdb.set_trace()
        month = da.groupby('time.month').mean(dim='time')
        for y in years:
            y = str(y)
            try:
                da = xr.open_dataset(msg_folder+'gridsat_WA_-40_5000km2_13-19UTC' + y + '.nc')
            except:
                print('Error continue!', y)
                continue
            da = da.sel(time=(da['time.hour'] == 21))
            da = da['tir']/100

            da = da.where(da <= -40)

            month = month + da.groupby('time.month').mean(dim='time', skipna=True)
            cnt +=1
            del da
        #ipdb.set_trace()
        month = month / cnt

        enc = {'tir': {'complevel': 5,  'zlib': True}}
        month.to_netcdf(msg_folder + fname, encoding=enc)


    else:
        da = xr.open_dataset(msg_folder + fname)
        #ipdb.set_trace()
        month = da['tir']
        del da
    #month.values[month.values==0]=np.nan

    simple = month.plot(x='lon', y='lat', col='month', col_wrap=4, cmap='jet', transform=ccrs.PlateCarree(),
                        subplot_kws={'projection': ccrs.PlateCarree()}, vmax=-45, vmin=-60, levels = 6, extend='both')

    for ax in simple.axes.flat:
        ax.coastlines()
        ax.gridlines()
        #ax.set_extent([-17.5, 30, -6, 20])
        ax.set_aspect('equal', 'box-forced')


    plt.savefig(cnst.network_data + '/figs/HUARAZ/mean_t_'+htag+'UTC.png', dpi=300)


def climatology_month_mf(h):

    htag = str(h).zfill(2)+'UTC_1985-1993'

    msg_folder = cnst.GRIDSAT_PERU
    fname='aggs/gridsat_WA_-40C_climatology_mean_mf_'+htag+'_coarse.nc'
    cnt = 1
    if not os.path.isfile(msg_folder + fname):
        da = xr.open_mfdataset(msg_folder+'gridsat_WA_-40_5000km2_13-19UTC*.nc')
        da = da.sel(time=((da['time.hour'] == 21)&(da['time.year'] >= 1985)&(da['time.year'] <= 1993)))
        #da = da['tir']/100
        da = da.where(da <= -4000)
        da['ymonth'] = (
        'time', [str(y) + '-' + str(m) for (y, m) in zip(da['time.year'].values, da['time.month'].values)])
        #ipdb.set_trace()
        #trendmonth = da.resample(time='1M')

        mean = da['tir'].groupby('time.month').mean('time')

        grid = mean.salem.grid.regrid(factor=0.25)
        coarse = grid.lookup_transform(mean)
        grid = grid.to_dataset()
        mean = xr.DataArray(coarse, coords=[mean['month'], grid['y'], grid['x']], dims=['month', 'lat', 'lon'])

        month = mean/100
        month.name='tir'

       # ipdb.set_trace()

        enc = {'tir': {'complevel': 5,  'zlib': True}}
        month.to_netcdf(msg_folder + fname, encoding=enc)


    else:
        da = xr.open_dataset(msg_folder + fname)
        #ipdb.set_trace()
        month = da['tir']
        del da
    #month.values[month.values==0]=np.nan

    simple = month.plot(x='lon', y='lat', col='month', col_wrap=4, cmap='jet', transform=ccrs.PlateCarree(),
                        subplot_kws={'projection': ccrs.PlateCarree()}, vmax=-40, vmin=-65, levels = 10, extend='both')

    for ax in simple.axes.flat:
        ax.coastlines()
        ax.gridlines()
        #ax.set_extent([-17.5, 30, -6, 20])
        ax.set_aspect('equal', 'box-forced')


    plt.savefig(cnst.network_data + '/figs/HUARAZ/mean_t_coarse_'+htag+'.png', dpi=300)


def plot(h):
    msg_folder = cnst.GRIDSAT_PERU
    htag1 = str(h).zfill(2) + 'UTC_2010-2017'
    p1 = 'aggs/gridsat_WA_-40C_climatology_mean_mf_'+htag1+'.nc'
    dat1 = xr.open_dataset(msg_folder + p1)
    #ipdb.set_trace()
    month1 = dat1['tir']

    htag2 = str(h).zfill(2) + 'UTC_1985-1993'
    p2 = 'aggs/gridsat_WA_-40C_climatology_mean_mf_'+htag2+'.nc'
    dat2 = xr.open_dataset(msg_folder + p2)
    month2 = dat2['tir']


    for d, tag in [(month1, htag1), (month2,htag2)]:
        d = d.sel(lon=slice(-79, -74), lat=slice(-12, -7))
        simple = d.plot(x='lon', y='lat', col='month', col_wrap=4, cmap='jet', transform=ccrs.PlateCarree(),
                            subplot_kws={'projection': ccrs.PlateCarree()}, vmax=-40, vmin=-65, levels=10,
                            extend='both')

        for ax in simple.axes.flat:
            ax.coastlines()
            ax.gridlines()
            ax.plot(-77.55, -9.51, 'ko')
            # ax.set_extent([-17.5, 30, -6, 20])
            ax.set_aspect('equal', 'box-forced')

        plt.savefig(cnst.network_data + '/figs/HUARAZ/mean_t_' + tag + '_small.png', dpi=300)

    fname = '/home/ck/DIR/cornkle/data/HUARAZ/shapes/riosan_sel_one.shp'

    shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                   ccrs.PlateCarree(), facecolor='none')


    diff = month1-month2
    diff = diff.sel(lon=slice(-79,-74), lat=slice(-12,-7))
    simple = diff.plot(x='lon', y='lat', col='month', col_wrap=4, cmap='RdBu_r', transform=ccrs.PlateCarree(),
                    subplot_kws={'projection': ccrs.PlateCarree()}, vmax=5, vmin=-5, levels=10,
                    extend='both')

    for ax in simple.axes.flat:
        ax.coastlines()
        ax.gridlines()
        # ax.set_extent([-17.5, 30, -6, 20])
        ax.set_aspect('equal', 'box-forced')
        ax.add_feature(shape_feature)
        ax.add_geometries(Reader(fname).geometries(),
                          ccrs.PlateCarree(),
                          facecolor='white', hatch='xxxx')

    plt.savefig(cnst.network_data + '/figs/HUARAZ/mean_t_diff_fistlast8years_16UTC_small.png', dpi=300)
    plt.close('all')



def climatology_month_mf_count(h):

    htag = str(h).zfill(2)+'UTC_1985-1993'

    msg_folder = cnst.GRIDSAT_PERU
    fname='aggs/gridsat_WA_-40C_climatology_mean_mf_'+htag+'_count.nc'
    cnt = 1
    if not os.path.isfile(msg_folder + fname):
        da = xr.open_mfdataset(msg_folder+'gridsat_WA_-40_5000km2_13-19UTC*.nc')
        da = da.sel(time=((da['time.hour'] == 21)&(da['time.year'] >= 1985)&(da['time.year'] <= 1993)))
        #da = da['tir']/100
        da = da.where(da <= -4000)
        da['ymonth'] = (
        'time', [str(y) + '-' + str(m) for (y, m) in zip(da['time.year'].values, da['time.month'].values)])
        #ipdb.set_trace()
        #trendmonth = da.resample(time='1M')

        mean = da['tir'].groupby('time.month').count('time')

        # grid = mean.salem.grid.regrid(factor=0.25)
        # coarse = grid.lookup_transform(mean)
        # grid = grid.to_dataset()
        # mean = xr.DataArray(coarse, coords=[mean['month'], grid['y'], grid['x']], dims=['month', 'lat', 'lon'])

        month = mean/100
        month.name='tir'

       # ipdb.set_trace()

        enc = {'tir': {'complevel': 5,  'zlib': True}}
        month.to_netcdf(msg_folder + fname, encoding=enc)


    else:
        da = xr.open_dataset(msg_folder + fname)
        #ipdb.set_trace()
        month = da['tir']
        del da
    #month.values[month.values==0]=np.nan

    simple = month.plot(x='lon', y='lat', col='month', col_wrap=4, cmap='jet', transform=ccrs.PlateCarree(),
                        subplot_kws={'projection': ccrs.PlateCarree()}, vmax=1, vmin=0.1, levels = 10, extend='both')

    for ax in simple.axes.flat:
        ax.coastlines()
        ax.gridlines()
        #ax.set_extent([-17.5, 30, -6, 20])
        ax.set_aspect('equal', 'box-forced')


    plt.savefig(cnst.network_data + '/figs/HUARAZ/mean_t_count_'+htag+'.png', dpi=300)



def climatology_month_mf_count(h):

    htag = str(h).zfill(2)+'UTC_1985-1993'

    msg_folder = cnst.GRIDSAT_PERU
    fname='aggs/gridsat_WA_-40C_climatology_mean_mf_'+htag+'_count.nc'
    cnt = 1
    if not os.path.isfile(msg_folder + fname):
        da = xr.open_mfdataset(msg_folder+'gridsat_WA_-40_5000km2_13-19UTC*.nc')
        da = da.sel(time=((da['time.hour'] == 21)&(da['time.year'] >= 1985)&(da['time.year'] <= 1993)))
        #da = da['tir']/100
        da = da.where(da <= -4000)
        da['ymonth'] = (
        'time', [str(y) + '-' + str(m) for (y, m) in zip(da['time.year'].values, da['time.month'].values)])
        #ipdb.set_trace()
        #trendmonth = da.resample(time='1M')

        mean = da['tir'].groupby('time.month').count('time')

        # grid = mean.salem.grid.regrid(factor=0.25)
        # coarse = grid.lookup_transform(mean)
        # grid = grid.to_dataset()
        # mean = xr.DataArray(coarse, coords=[mean['month'], grid['y'], grid['x']], dims=['month', 'lat', 'lon'])

        month = mean/100
        month.name='tir'

       # ipdb.set_trace()

        enc = {'tir': {'complevel': 5,  'zlib': True}}
        month.to_netcdf(msg_folder + fname, encoding=enc)


    else:
        da = xr.open_dataset(msg_folder + fname)
        #ipdb.set_trace()
        month = da['tir']
        del da
    #month.values[month.values==0]=np.nan

    simple = month.plot(x='lon', y='lat', col='month', col_wrap=4, cmap='jet', transform=ccrs.PlateCarree(),
                        subplot_kws={'projection': ccrs.PlateCarree()}, vmax=1, vmin=0.1, levels = 10, extend='both')

    for ax in simple.axes.flat:
        ax.coastlines()
        ax.gridlines()
        #ax.set_extent([-17.5, 30, -6, 20])
        ax.set_aspect('equal', 'box-forced')


    plt.savefig(cnst.network_data + '/figs/HUARAZ/mean_t_count_'+htag+'.png', dpi=300)
