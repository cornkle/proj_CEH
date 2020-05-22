# -*- coding: utf-8 -*-


import ipdb
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import cartopy
import cartopy.crs as ccrs
from utils import constants as cnst
from eod import msg
from utils import u_darrays as uda
from scipy import stats
import warnings
import multiprocessing


def corr(a, b, bsingle=None, c_box=None):
    ds = xr.Dataset()
    #ipdb.set_trace()
    ds['pval'] = a.copy(deep=True).sum('time') * np.nan
    ds['r'] = a.copy(deep=True).sum('time') * np.nan
    ds['slope'] = a.copy(deep=True).sum('time') * np.nan
    ds['intercept'] = a.copy(deep=True).sum('time') * np.nan

    corr_box = c_box
    perPixel = False
    if bsingle:
        bb = b
    elif c_box:
        bb = b.sel(latitude=slice(corr_box[2], corr_box[3]), longitude=slice(corr_box[0], corr_box[1])).mean(dim=['latitude', 'longitude'])
    else:
        perPixel=True


    for lat in a.latitude.values:
        for lon in a.longitude.values:
            aa = a.sel(latitude=lat, longitude=lon)
            if bsingle:
                r, p = stats.pearsonr(aa.values, bb)

                #pf = np.polyfit(aa.values, bb, 1)
                pf, intercept, r, p, std_err = stats.linregress(aa.values, bb)
            elif c_box:
                # r, p = stats.pearsonr(aa.values, bb.values)
                # pf = np.polyfit(aa.values, bb.values, 1)
                pf, intercept, r, p, std_err = stats.linregress(aa.values, bb.values)
            elif perPixel:
                bb = b.sel(latitude=lat, longitude=lon)
                pf, intercept, r, p, std_err = stats.linregress(aa.values, bb.values)


            slope = pf#[0]

#                 if (np.nansum(aa.values == 0) >= 10):
#                     p = np.nan
#                     r = np.nan

            ds['r'].loc[{'latitude': lat, 'longitude': lon}] = r
            ds['pval'].loc[{'latitude': lat, 'longitude': lon}] = p
            ds['slope'].loc[{'latitude': lat, 'longitude': lon}] = slope
            ds['intercept'].loc[{'latitude': lat, 'longitude': lon}] = intercept

    return ds

def readERA():

    u200orig = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/u_15UTC_1981-2019_peru_big.nc')
    u200orig = uda.flip_lat(u200orig)
    datetimes = pd.to_datetime(u200orig.time.values)
    newtimes = []
    for t in datetimes:
        newtimes.append(t.replace(hour=0))
    u200orig['time'] = ('time', newtimes)
    return u200orig


def saveCHIRPS():

    chirpsbox = [-81, -68, -17, 0]

    u200orig = readERA()
    u200 = u200orig.sel(longitude=slice(chirpsbox[0], chirpsbox[1]), latitude=slice(chirpsbox[3], chirpsbox[2]))

    chirps_all = xr.open_dataset('/home/ck/DIR/mymachine/CHIRPS/peru/CHIRPS_peru_onERA5.nc')
    pos = np.intersect1d(chirps_all.time, u200.time)

    chirps = chirps_all.sel(time=pos)
    u200 = u200.sel(time=pos)

    cdoy = chirps['precip'].rolling(time=3, min_periods=1, center=True).mean(dim='time')
    udoy = u200['u'].rolling(time=3, min_periods=1, center=True).mean(dim='time')

    dslist = []

    for doy in np.arange(1,366): #366
        dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))

    pool = multiprocessing.Pool(processes=6)

    res = pool.map(run_doy, dslist)
    pool.close()

    chirps_ds = xr.concat(res, dim='dayofyear')
    chirps_ds.attrs['years'] = np.unique(chirps['time.year'])
    chirps_ds.to_netcdf('/home/ck/DIR/mymachine/CHIRPS/peru/CHIRPS_u200_correlation_peru.nc')


def saveGRIDSAT():

    u200orig = readERA()
    u200orig = u200orig.sel(time=(u200orig['time.year'] <= 1998))

    gridsat = xr.open_mfdataset('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/daily_ALLkm2_UTC_DAY_onBIGERA/*.nc',
                                combine='nested', concat_dim='time').load()
    posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

    u200orig = u200orig.sel(time=posgrid)
    gridsat = gridsat.sel(time=posgrid)



    udoy = u200orig['u'].rolling(time=3, min_periods=1, center=True).mean(dim='time')
    cdoy = gridsat['tir'].rolling(time=3, min_periods=1, center=True).mean(dim='time')

    dslist = []
    print('Did rolling aggregation')

    for doy in np.arange(1,366): #366
        dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))

    pool = multiprocessing.Pool(processes=6)

    res = pool.map(run_doy, dslist)
    pool.close()

    # #for d in np.arange(1, 6):  # 366
    #     ccdoy = cdoy['tir'].sel(time=(cdoy['time.dayofyear'] == d))
    #     uudoy = udoy['u'].sel(time=(udoy['time.dayofyear'] == d))
    #
    #     # ipdb.set_trace()
    #
    #     diff1 = xr.DataArray(ccdoy.values[1::, :, :] - ccdoy.values[0:-1, :, :],
    #                          coords=[ccdoy.time[1::], ccdoy.latitude, ccdoy.longitude],
    #                          dims=['time', 'latitude', 'longitude'])
    #     diff2 = xr.DataArray(uudoy.values[1::, :, :] - uudoy.values[0:-1, :, :],
    #                          coords=[uudoy.time[1::], uudoy.latitude, uudoy.longitude],
    #                          dims=['time', 'latitude', 'longitude'])

        # outarr = corr(diff1, diff2)
        # dslist.append(outarr)
    chirps_ds = xr.concat(res, dim='dayofyear')
    chirps_ds.attrs['years'] = np.unique(gridsat['time.year'])
    chirps_ds.to_netcdf('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/correlations/GRIDSAT_u_correlation_SouthAmerica_1985-2000.nc')


def saveGPM():

    u200orig = readERA()
    #u200orig = u200orig.sel(time=(u200orig['time.year'] > 1998))

    gridsat = xr.open_mfdataset('/media/ck/Elements/SouthAmerica/GPM/daily_onERA/*.nc',
                                combine='nested', concat_dim='time').load()
    posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

    u200orig = u200orig.sel(time=posgrid)
    gridsat = gridsat.sel(time=posgrid)

    grid = gridsat.isel(time=0)
    grid = grid.where(np.isfinite(grid['precip']), drop=True)


    u200orig = u200orig.sel(latitude=slice(grid.latitude.min(), grid.latitude.max()), longitude=slice(grid.longitude.min(), grid.longitude.max()))
    gridsat = gridsat.sel(latitude=slice(grid.latitude.min(), grid.latitude.max()), longitude=slice(grid.longitude.min(), grid.longitude.max()))

    ipdb.set_trace()

    udoy = u200orig['u'].rolling(time=3, min_periods=1, center=True).mean(dim='time')
    cdoy = gridsat['precip'].rolling(time=3, min_periods=1, center=True).mean(dim='time')

    dslist = []
    print('Did rolling aggregation')

    for doy in np.arange(1,366): #366
        dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))

    pool = multiprocessing.Pool(processes=6)

    res = pool.map(run_doy, dslist)
    pool.close()

    # #for d in np.arange(1, 6):  # 366
    #     ccdoy = cdoy['tir'].sel(time=(cdoy['time.dayofyear'] == d))
    #     uudoy = udoy['u'].sel(time=(udoy['time.dayofyear'] == d))
    #
    #     # ipdb.set_trace()
    #
    #     diff1 = xr.DataArray(ccdoy.values[1::, :, :] - ccdoy.values[0:-1, :, :],
    #                          coords=[ccdoy.time[1::], ccdoy.latitude, ccdoy.longitude],
    #                          dims=['time', 'latitude', 'longitude'])
    #     diff2 = xr.DataArray(uudoy.values[1::, :, :] - uudoy.values[0:-1, :, :],
    #                          coords=[uudoy.time[1::], uudoy.latitude, uudoy.longitude],
    #                          dims=['time', 'latitude', 'longitude'])

        # outarr = corr(diff1, diff2)
        # dslist.append(outarr)
    chirps_ds = xr.concat(res, dim='dayofyear')
    chirps_ds.attrs['years'] = np.unique(gridsat['time.year'])
    chirps_ds.to_netcdf('/home/ck/DIR/mymachine/GPM/GPM_u_correlation_SouthAmerica_2000-2018.nc')


def saveGPM_GRIDSAT():

    u200orig =  xr.open_mfdataset('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/daily_ALLkm2_UTC_DAY_onBIGERA/*.nc',
                                combine='nested', concat_dim='time').load()
    #u200orig = u200orig.sel(time=(u200orig['time.year'] > 1998))

    gridsat = xr.open_mfdataset('/media/ck/Elements/SouthAmerica/GPM/daily_onERA/*.nc',
                                combine='nested', concat_dim='time').load()
    posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

    u200orig = u200orig.sel(time=posgrid)
    gridsat = gridsat.sel(time=posgrid)

    grid = gridsat.isel(time=0)
    grid = grid.where(np.isfinite(grid['precip']), drop=True)


    u200orig = u200orig.sel(latitude=slice(grid.latitude.min(), grid.latitude.max()), longitude=slice(grid.longitude.min(), grid.longitude.max()))
    gridsat = gridsat.sel(latitude=slice(grid.latitude.min(), grid.latitude.max()), longitude=slice(grid.longitude.min(), grid.longitude.max()))

    #ipdb.set_trace()

    udoy = u200orig['tir'].rolling(time=3, min_periods=1, center=True).mean(dim='time')
    cdoy = gridsat['precip'].rolling(time=3, min_periods=1, center=True).mean(dim='time')

    dslist = []
    print('Did rolling aggregation')

    for doy in np.arange(1,366): #366
        dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))

    pool = multiprocessing.Pool(processes=6)

    res = pool.map(run_doy, dslist)
    pool.close()

    # #for d in np.arange(1, 6):  # 366
    #     ccdoy = cdoy['tir'].sel(time=(cdoy['time.dayofyear'] == d))
    #     uudoy = udoy['u'].sel(time=(udoy['time.dayofyear'] == d))
    #
    #     # ipdb.set_trace()
    #
    #     diff1 = xr.DataArray(ccdoy.values[1::, :, :] - ccdoy.values[0:-1, :, :],
    #                          coords=[ccdoy.time[1::], ccdoy.latitude, ccdoy.longitude],
    #                          dims=['time', 'latitude', 'longitude'])
    #     diff2 = xr.DataArray(uudoy.values[1::, :, :] - uudoy.values[0:-1, :, :],
    #                          coords=[uudoy.time[1::], uudoy.latitude, uudoy.longitude],
    #                          dims=['time', 'latitude', 'longitude'])

        # outarr = corr(diff1, diff2)
        # dslist.append(outarr)
    chirps_ds = xr.concat(res, dim='dayofyear')
    chirps_ds.attrs['years'] = np.unique(gridsat['time.year'])
    chirps_ds.to_netcdf('/home/ck/DIR/mymachine/GPM/GPM_GRIDSAT_correlation_SouthAmerica_2000-2018.nc')




def run_doy(x):

    d = x[0]
    print('Doing dayofyear', d)

    udoy = x[1]
    cdoy = x[2]

    ccdoy = cdoy.sel(time=(cdoy['time.dayofyear'] == d))
    uudoy = udoy.sel(time=(udoy['time.dayofyear'] == d))

    # ipdb.set_trace()

    diff1 = xr.DataArray(ccdoy.values[1::, :, :] - ccdoy.values[0:-1, :, :],
                         coords=[ccdoy.time[1::], ccdoy.latitude, ccdoy.longitude],
                         dims=['time', 'latitude', 'longitude'])
    diff2 = xr.DataArray(uudoy.values[1::, :, :] - uudoy.values[0:-1, :, :],
                         coords=[uudoy.time[1::], uudoy.latitude, uudoy.longitude],
                         dims=['time', 'latitude', 'longitude'])

    outarr = corr(diff1, diff2)

    del diff1
    del diff2
    del ccdoy
    del uudoy

    return outarr



def saveERA5():
    #chirpsbox = [-81, -68, -18.5, 0]  # peru daily

    era5 = xr.open_mfdataset(cnst.ERA5_HOURLY_PL_HU + '/ERA5_*_pl.nc', concat_dim='time', combine='nested')
    u200 = era5['u'].sel(level=200, time=(era5[
                                              'time.hour'] == 15))  # .load()   #longitude=slice(bigbox[0], bigbox[1]), latitude=slice(bigbox[3], bigbox[2]),
    u200 = uda.flip_lat(u200)

    u200 = u200.sel(time=((u200['time.year'] > 1980) & (u200['time.year'] < 2019)))
    u200.name = 'u'
    comp = dict(zlib=True, complevel=5)
    encoding = {'u': comp}
    u200.to_netcdf('/media/ck/Elements/SouthAmerica/ERA5/hourly/pressure_levels/u_15UTC_1981-2018_peru.nc', mode='w',
                  encoding=encoding, format='NETCDF4')



def saveGRIDSAT_fullYear():

    u200orig = readERA()
    u200orig = u200orig.sel(time=(u200orig['time.year'] > 1998))

    gridsat = xr.open_mfdataset('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/daily_ALLkm2_UTC_DAY_onBIGERA/*.nc',
                                combine='nested', concat_dim='time').load()
    posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

    u200orig = u200orig.sel(time=posgrid)
    gridsat = gridsat.sel(time=posgrid)



    udoy = u200orig['u'].rolling(time=3, min_periods=1, center=True).mean(dim='time')
    cdoy = gridsat['tir'].rolling(time=3, min_periods=1, center=True).mean(dim='time')

    dslist = []
    print('Did rolling aggregation')

    # for doy in np.arange(1,366): #366
    #     dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))

    pool = multiprocessing.Pool(processes=6)

    res = pool.map(run_year, (udoy,cdoy))
    pool.close()

    # #for d in np.arange(1, 6):  # 366
    #     ccdoy = cdoy['tir'].sel(time=(cdoy['time.dayofyear'] == d))
    #     uudoy = udoy['u'].sel(time=(udoy['time.dayofyear'] == d))
    #
    #     # ipdb.set_trace()
    #
    #     diff1 = xr.DataArray(ccdoy.values[1::, :, :] - ccdoy.values[0:-1, :, :],
    #                          coords=[ccdoy.time[1::], ccdoy.latitude, ccdoy.longitude],
    #                          dims=['time', 'latitude', 'longitude'])
    #     diff2 = xr.DataArray(uudoy.values[1::, :, :] - uudoy.values[0:-1, :, :],
    #                          coords=[uudoy.time[1::], uudoy.latitude, uudoy.longitude],
    #                          dims=['time', 'latitude', 'longitude'])

        # outarr = corr(diff1, diff2)
        # dslist.append(outarr)
    chirps_ds = xr.concat(res, dim='dayofyear')
    chirps_ds.attrs['years'] = np.unique(gridsat['time.year'])
    chirps_ds.to_netcdf('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/correlations/GRIDSAT_u_correlation_SouthAmerica_2000-2018_fullyear.nc')

def run_year(x):

    #d = x[0]
    #print('Doing dayofyear', d)

    uudoy = x[0]
    ccdoy = x[1]

    # ccdoy = cdoy.sel(time=(cdoy['time.dayofyear'] == d))
    # uudoy = udoy.sel(time=(udoy['time.dayofyear'] == d))

    # ipdb.set_trace()

    diff1 = xr.DataArray(ccdoy.values[1::, :, :] - ccdoy.values[0:-1, :, :],
                         coords=[ccdoy.time[1::], ccdoy.latitude, ccdoy.longitude],
                         dims=['time', 'latitude', 'longitude'])
    diff2 = xr.DataArray(uudoy.values[1::, :, :] - uudoy.values[0:-1, :, :],
                         coords=[uudoy.time[1::], uudoy.latitude, uudoy.longitude],
                         dims=['time', 'latitude', 'longitude'])

    outarr = corr(diff1, diff2)

    del diff1
    del diff2
    del ccdoy
    del uudoy

    return outarr