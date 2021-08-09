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
import pandas as pd
import scipy


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

def lag_linregress_3D(x, y, lagx=0, lagy=0, stride=3):
    """
    Input: Two xr.Datarrays of any dimensions with the first dim being time.
    Thus the input data could be a 1D time series, or for example, have three
    dimensions (time,lat,lon).
    Datasets can be provided in any order, but note that the regression slope
    and intercept will be calculated for y with respect to x.
    Output: Covariance, correlation, regression slope and intercept, p-value,
    and standard error on regression between the two datasets along their
    aligned time dimension.
    Lag values can be assigned to either of the data, with lagx shifting x, and
    lagy shifting y, with the specified lag amount.
    """
    #1. Ensure that the data are properly alinged to each other.
    x,y = xr.align(x,y)

    #2. Add lag information if any, and shift the data accordingly
    if lagx!=0:

        # If x lags y by 1, x must be shifted 1 step backwards.
        # But as the 'zero-th' value is nonexistant, xr assigns it as invalid
        # (nan). Hence it needs to be dropped
        x   = x.shift(time = -lagx).dropna(dim='time')

        # Next important step is to re-align the two datasets so that y adjusts
        # to the changed coordinates of x
        x,y = xr.align(x,y)

    if lagy!=0:
        y   = y.shift(time = -lagy).dropna(dim='time')
        x,y = xr.align(x,y)

    #3. Compute data length, mean and standard deviation along time axis:
    n = y.notnull().sum(dim='time')
    xmean = x.mean(axis=0)
    ymean = y.mean(axis=0)
    xstd  = x.std(axis=0)
    ystd  = y.std(axis=0)

    #4. Compute covariance along time axis
    cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)

    #5. Compute correlation along time axis
    cor   = cov/(xstd*ystd)

    #6. Compute regression slope and intercept:
    slope     = cov/(xstd**2)
    intercept = ymean - xmean*slope

    #7. Compute P-value and standard error
    #Compute t-statistics

    #tstats = cor*np.sqrt(n-2)/np.sqrt(1-cor**2)  # Tstat with T normal distribution
    TINY = 1.0e-20
    tstats = cor * np.sqrt(n-2 / ((1.0 - cor + TINY)*(1.0 + cor + TINY)))

    stderr = slope/tstats

    from scipy.stats import t
    pval   = t.sf(np.abs(tstats), n-2)*2
    #pval = t.sf(np.abs(tstats), n - 2)*2
    pval = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)

    #ipdb.set_trace()
    # nullmask = np.sum((y.values == 0), axis=0)
    # pos = np.where(nullmask>=0.75*stride)

    return cov,cor,slope,intercept,pval,stderr



def readERA(var):

    u200orig = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/'+var+'_15UTC_1981-2019_peru_big.nc')
    #u200orig = xr.open_dataset('/media/ck/Elements/SouthAmerica/ERA5/hourly/v850_15UTC_1981-2019_peru_big.nc')
    u200orig = uda.flip_lat(u200orig)
    datetimes = pd.to_datetime(u200orig.time.values)
    newtimes = []
    for t in datetimes:
        newtimes.append(t.replace(hour=0))
    u200orig['time'] = ('time', newtimes)
    return u200orig


def saveCHIRPS_BIG(y1,y2):

    u200orig = readERA('u200')
    u200 = u200orig#.sel(longitude=slice(chirpsbox[0], chirpsbox[1]), latitude=slice(chirpsbox[2], chirpsbox[3]))
    u200 = u200.sel(time=((u200orig['time.year'] >= y1) & (u200orig['time.year'] < y2)))

    chirps_all = xr.open_mfdataset('/media/ck/Elements/SouthAmerica/CHIRPS/SA_daily_onERA/*.nc')
    pos = np.intersect1d(chirps_all.time, u200.time)

    chirps = chirps_all.sel(time=pos).load()
    u200 = u200.sel(time=pos)

    cdoy = chirps['precip'].where((chirps['precip']>1) | np.isnan(chirps['precip']), other=0)
    #ipdb.set_trace()
    cdoy = cdoy.rolling(time=5, min_periods=3, center=True).mean(dim='time')
    udoy = u200['u'].rolling(time=5, min_periods=3, center=True).mean(dim='time')

    dslist = []

    for doy in np.arange(1,366): #366
        dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))

    pool = multiprocessing.Pool(processes=4)

    # res = pool.map(run_doy_correlation, dslist)
    # pool.close()

    res = []
    for d in dslist:
        out = run_doy_correlation(d)
        res.append(out)

    chirps_ds = xr.concat(res, dim='dayofyear')
    chirps_ds.attrs['years'] = np.unique(chirps['time.year'])

    chirps_ds.to_netcdf('/home/ck/DIR/mymachine/CHIRPS/peru/CHIRPS_u200_correlation_5dRolling_1mm_peruBIG_'+str(y1)+'-'+str(y2-1)+'_diffs.nc')



def runCHIRPS():
    ylist = [(1985,2004), (2000,2019)]#[(1985,2003), (2002,2019), (1985,2019)] #[(1985,2005), (1998,2019), (1985,2019), (1985,1995), (1985,2002), (2002,2019), (1995,2005), (2005,2019)]
    for y in ylist:
        saveCHIRPS_BIG(y[0], y[1])

def saveCHIRPS(y1,y2):

    chirpsbox = [-81, -68, -17, 0]

    u200orig = readERA('u200')
    u200 = u200orig.sel(longitude=slice(chirpsbox[0], chirpsbox[1]), latitude=slice(chirpsbox[2], chirpsbox[3]))
    u200 = u200.sel(time=((u200orig['time.year'] >= y1) & (u200orig['time.year'] < y2)))

    chirps_all = xr.open_dataset('/home/ck/DIR/mymachine/CHIRPS/peru/CHIRPS_peru_onERA5.nc')
    pos = np.intersect1d(chirps_all.time, u200.time)

    chirps = chirps_all.sel(time=pos)
    u200 = u200.sel(time=pos)

    cdoy = chirps['precip'].rolling(time=5, min_periods=3, center=True).mean(dim='time')
    udoy = u200['u'].rolling(time=5, min_periods=3, center=True).mean(dim='time')

    dslist = []

    for doy in np.arange(1,366): #366
        dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))

    pool = multiprocessing.Pool(processes=4)

    # res = pool.map(run_doy_correlation, dslist)
    # pool.close()

    res = []
    for d in dslist:
        out = run_doy_correlation(d)
        res.append(out)

    chirps_ds = xr.concat(res, dim='dayofyear')
    chirps_ds.attrs['years'] = np.unique(chirps['time.year'])
  #  ipdb.set_trace()
    chirps_ds.to_netcdf('/home/ck/DIR/mymachine/CHIRPS/peru/CHIRPS_u200_correlation_5dRolling_peru_'+str(y1)+'-'+str(y2-1)+'_diffs.nc')


def runGridsat():
    ylist = [(1985,2004), (2000,2019)]#[(1985,2003), (2002,2019), (1985,2019)]
    for y in ylist:
        saveGRIDSAT(y[0], y[1])

def saveGRIDSAT(y1,y2):

    u200orig = readERA('u200')
    u200orig = u200orig.sel(time=((u200orig['time.year'] >= y1) & (u200orig['time.year'] < y2)))

    gridsat = xr.open_mfdataset('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/daily_-15ALLkm2_UTC_DAY_onBIGERA/*.nc',
                                combine='nested', concat_dim='time')
    posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

    u200orig = u200orig.sel(time=posgrid)
    gridsat = gridsat.sel(time=posgrid).load()

    udoy = u200orig['u'].rolling(time=3, min_periods=1, center=True).mean(dim='time') # time=3, min_periods=1
    cdoy = gridsat['tir'].where(gridsat['tir']<=-3000, other=0).rolling(time=5, min_periods=3, center=True).mean(dim='time')

    dslist = []
    print('Did rolling aggregation')

    for doy in np.arange(1,366): #366
        dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))
    #ipdb.set_trace()
    # pool = multiprocessing.Pool(processes=4)
    #
    # res = pool.map(run_doy_correlation, dslist)
    # pool.close()
    res = []
    for d in dslist:
        out = run_doy_correlation(d)
        res.append(out)

    chirps_ds = xr.concat(res, dim='dayofyear')
    chirps_ds.attrs['years'] = np.unique(gridsat['time.year'])
    chirps_ds.to_netcdf('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/correlations/GRIDSAT-30_u_correlation_3dRolling_SouthAmerica_'+str(y1)+'-'+str(y2-1)+'_diffs.nc')


def saveGRIDSAT_DJF(y1,y2):

    u200orig = readERA('u200')
    u200orig = u200orig.sel(time=((u200orig['time.year'] >= y1) & (u200orig['time.year'] < y2)))

    gridsat = xr.open_mfdataset('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/daily_-40ALLkm2_UTC_DAY_onBIGERA/*.nc',
                                combine='nested', concat_dim='time')
    posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

    u200orig = u200orig.sel(time=posgrid)
    gridsat = gridsat.sel(time=posgrid).load()

    udoy = uda.season_mean(u200orig['u']) # time=3
    cdoy = uda.season_mean(gridsat['tir'])

    print('Did rolling aggregation')

    res = corr(udoy, cdoy)
    #ipdb.set_trace()
    #chirps_ds = xr.concat(res, dim='dayofyear')
    #chirps_ds.attrs['years'] = np.unique(gridsat['time.year'])
    res.to_netcdf('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/correlations/GRIDSAT-40_u_correlation_SouthAmerica_'+str(y1)+'-'+str(y2-1)+'_diffs_DJF.nc')


def saveGRIDSAT_NDVI_DJF(y1,y2):


    u200orig =  xr.open_mfdataset('/media/ck/Elements/SouthAmerica/NDVI/onERA/*.nc')

    u200orig = u200orig['precip'].load()
    u200orig = u200orig.sel(time=((u200orig['time.year'] >= y1) & (u200orig['time.year'] < y2)))

    gridsat = xr.open_mfdataset('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/daily_-40ALLkm2_UTC_DAY_onBIGERA/*.nc',
                                combine='nested', concat_dim='time')
    posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

    u200orig = u200orig.sel(time=posgrid)
    gridsat = gridsat.sel(time=posgrid).load()

    udoy = uda.season_mean(u200orig) # time=3
    cdoy = uda.season_mean(gridsat['tir'])

    print('Did rolling aggregation')

    res = corr(udoy, cdoy)
    #ipdb.set_trace()
    #chirps_ds = xr.concat(res, dim='dayofyear')
    #chirps_ds.attrs['years'] = np.unique(gridsat['time.year'])
    res.to_netcdf('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/correlations/GRIDSAT-40_NDVI_correlation_SouthAmerica_'+str(y1)+'-'+str(y2-1)+'_diffs_DJF.nc')



def saveCHIRPS_DJF(y1,y2):

    chirpsbox = [-81, -68, -17, 0]

    u200orig = readERA('u200')
    u200 = u200orig.sel(longitude=slice(chirpsbox[0], chirpsbox[1]), latitude=slice(chirpsbox[2], chirpsbox[3]))
    u200 = u200.sel(time=((u200orig['time.year'] >= y1) & (u200orig['time.year'] < y2)))

    chirps_all = xr.open_dataset('/home/ck/DIR/mymachine/CHIRPS/peru/CHIRPS_peru_onERA5.nc')
    pos = np.intersect1d(chirps_all.time, u200.time)

    chirps = chirps_all.sel(time=pos)
    u200 = u200.sel(time=pos)


    udoy = uda.season_mean(u200['u']) # time=3
    cdoy = uda.season_mean(chirps['precip'])

    print('Did rolling aggregation')

    res = corr(udoy, cdoy)

  #  ipdb.set_trace()
    res.to_netcdf('/home/ck/DIR/mymachine/CHIRPS/peru/CHIRPS_u200_correlation_peru_'+str(y1)+'-'+str(y2-1)+'_diffs_DJF.nc')


def saveCHIRPS_NDVI_DJF(y1,y2):

    chirpsbox = [-81, -68, -17, 0]

    u200orig = xr.open_mfdataset('/media/ck/Elements/SouthAmerica/NDVI/onCHIRPSbox/*')
    #u200orig = uda.flip_lat(u200orig)
    u200 = u200orig['precip'].sel(longitude=slice(chirpsbox[0], chirpsbox[1]), latitude=slice(chirpsbox[2], chirpsbox[3])).load()

    u200 = u200.sel(time=((u200orig['time.year'] >= y1) & (u200orig['time.year'] < y2)))

    chirps_all = xr.open_dataset('/home/ck/DIR/mymachine/CHIRPS/peru/CHIRPS_peru_onERA5.nc')
    pos = np.intersect1d(chirps_all.time, u200.time)

    chirps = chirps_all.sel(time=pos)
    u200 = u200.sel(time=pos)


    udoy = uda.season_mean(u200) # time=3
    cdoy = uda.season_mean(chirps['precip'])

    print('Did rolling aggregation')

    res = corr(udoy, cdoy)

    #ipdb.set_trace()
    res.to_netcdf('/home/ck/DIR/mymachine/CHIRPS/peru/CHIRPS_NDVI_correlation_peru_'+str(y1)+'-'+str(y2-1)+'_diffs_DJF.nc')


def saveGRIDSAT_tropicalJet(y1,y2):

    u200orig = readERA('v850')
    isjet = [-75.5, -74.5, -8.5, -6.5]
    u200orig = u200orig.sel(time=((u200orig['time.year'] >= y1) & (u200orig['time.year'] < y2)))

    gridsat = xr.open_mfdataset('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/daily_-40ALLkm2_UTC_DAY_onBIGERA/*.nc',
                                combine='nested', concat_dim='time')
    posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

    u200orig = u200orig.sel(time=posgrid)
    gridsat = gridsat.sel(time=posgrid).load()

    udoy = u200orig['v'].rolling(time=3, min_periods=1, center=True).mean(dim='time') # time=3
    cdoy = gridsat['tir'].rolling(time=3, min_periods=1, center=True).mean(dim='time')

    dslist = []
    print('Did rolling aggregation')

    for doy in np.arange(1,366): #366
        dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))
    #ipdb.set_trace()
    # pool = multiprocessing.Pool(processes=4)
    #
    # res = pool.map(run_doy_correlation, dslist)
    # pool.close()
    res = []
    for d in dslist:
        out = run_doy(d, c_box=isjet)
        res.append(out)

    chirps_ds = xr.concat(res, dim='dayofyear')
    chirps_ds.attrs['years'] = np.unique(gridsat['time.year'])
    chirps_ds.to_netcdf('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/correlations/GRIDSAT_v850_correlation_SouthAmerica_'+str(y1)+'-'+str(y2-1)+'_diffs.nc')


def saveGPM():

    u200orig = readERA()
    u200orig = u200orig.sel(time=((u200orig['time.year'] >= 2000) & (u200orig['time.year'] <= 2018)))

    gridsat = xr.open_mfdataset('/media/ck/Elements/SouthAmerica/GPM/daily_onERA/*.nc',
                                combine='nested', concat_dim='time')
    posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

    u200orig = u200orig.sel(time=posgrid)
    gridsat = gridsat.sel(time=posgrid).load()

    grid = gridsat.isel(time=0)
    grid = grid.where(np.isfinite(grid['precip']), drop=True)


    u200orig = u200orig.sel(latitude=slice(grid.latitude.min(), grid.latitude.max()), longitude=slice(grid.longitude.min(), grid.longitude.max()))
    gridsat = gridsat.sel(latitude=slice(grid.latitude.min(), grid.latitude.max()), longitude=slice(grid.longitude.min(), grid.longitude.max()))

    #ipdb.set_trace()

    udoy = u200orig['u'].rolling(time=3, min_periods=1, center=True).mean(dim='time')
    cdoy = gridsat['precip'].rolling(time=3, min_periods=1, center=True).mean(dim='time')

    dslist = []
    print('Did rolling aggregation')

    for doy in np.arange(1,366): #366
        dslist.append((doy, udoy.sel(time=(udoy['time.dayofyear'] == doy)), cdoy.sel(time=(cdoy['time.dayofyear'] == doy))))

    # pool = multiprocessing.Pool(processes=4)
    #
    # res = pool.map(run_doy_correlation, dslist)
    # pool.close()
    #
    res = []
    for d in dslist:
        out = run_doy_correlation(d)
        res.append(out)


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

    res = pool.map(run_doy_correlation_optim, dslist)
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




def run_doy(x, c_box=None):

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

    outarr = corr(diff1, diff2, c_box=c_box)

    del diff1
    del diff2
    del ccdoy
    del uudoy

    return outarr



def saveERA5():
    #chirpsbox = [-81, -68, -18.5, 0]  # peru daily

    era5 = xr.open_mfdataset(cnst.ERA5_HOURLY_PL_HU + '/ERA5_*_pl.nc', concat_dim='time', combine='nested')
    u200 = era5['w'].sel(level=500, time=(era5[
                                              'time.hour'] == 15))  # .load()   #longitude=slice(bigbox[0], bigbox[1]), latitude=slice(bigbox[3], bigbox[2]),
    u200 = uda.flip_lat(u200)

    u200 = u200.sel(time=((u200['time.year'] > 1980) & (u200['time.year'] < 2019)))
    u200.name = 'w'
    comp = dict(zlib=True, complevel=5)
    encoding = {'w': comp}
    u200.to_netcdf('/media/ck/Elements/SouthAmerica/ERA5/hourly/pressure_levels/w500_15UTC_1981-2018_peru.nc', mode='w',
                  encoding=encoding, format='NETCDF4')



def saveGRIDSATPerYear_optim():

    u200 = readERA()

    for year in range(1985,2019):
        u200orig = u200.sel(time=(u200['time.year'] == year))

        gridsat = xr.open_mfdataset('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/daily_-40ALLkm2_UTC_DAY_onBIGERA/*.nc',
                                    combine='nested', concat_dim='time')

        posgrid = np.intersect1d(u200orig.time.values, gridsat.time.values)

        u200orig = u200orig.sel(time=posgrid)
        gridsat = gridsat.sel(time=posgrid).load()

        udoy = u200orig['u']
        cdoy = gridsat['tir']

        dslist = []
        print('Did rolling aggregation')

        # rw = rolling_window(np.arange(1,366), 15)
        res = []

        #for doy in np.arange(11,16): #366
        stride=7

        for doy in np.arange(1, 366):  # 366

            d1 = doy - stride
            d2 = doy + stride

            try:
                updoy = udoy.sel(time=((udoy['time.dayofyear'] >= d1) & (udoy['time.dayofyear'] <= d2)))
            except:
                continue
            try:
                cpdoy = cdoy.sel(time=((udoy['time.dayofyear'] >= d1) & (udoy['time.dayofyear'] <= d2)))
            except:
                continue

            #ipdb.set_trace()
            if (len(cpdoy.time) < stride) | (len(updoy.time) < stride) | (np.sum(udoy['time.dayofyear'] == doy)==0) | (np.sum(cdoy['time.dayofyear'] == doy)==0) :
                continue

            dslist.append(
                (doy, updoy, cpdoy))

        #ipdb.set_trace()
        pool = multiprocessing.Pool(processes=4)

        res = pool.map(run_doy_correlation_optim, dslist)
        pool.close()

        # for d in dslist:
        #     run_doy_correlation_optim(d)
        #
        # ipdb.set_trace()

        concat_ds = xr.concat(res, dim='time')
        concat_ds.to_netcdf('/home/ck/DIR/mymachine/GRIDSAT/MCS18_peru/correlations/rollingCorr/GRIDSAT_u200_rollingCorr_SouthAmerica_'+str(year)+'.nc')





def run_doy_correlation_optim(x):

    doy = x[0]

    udoy = x[1]
    cdoy = x[2]

    print('Doing ', doy)

    ds = xr.Dataset()
    try:
        ds['pval'] = udoy.sel(time=((udoy['time.dayofyear'] == doy))) * np.nan
    except:
        ipdb.set_trace()
    ds['r'] = udoy.sel(time=((udoy['time.dayofyear'] == doy))) * np.nan
    ds['slope'] = udoy.sel(time=((udoy['time.dayofyear'] == doy))) * np.nan

    if (np.isnan(udoy.all())) | (np.isnan(cdoy.all())):
        return ds

    uudoy = udoy#.sel(time=((udoy['time.dayofyear'] >= d1) & (udoy['time.dayofyear'] <= d2)))
    ccdoy = cdoy

    diff1 = xr.DataArray(ccdoy.values[1::, :, :] - ccdoy.values[0:-1, :, :],
                         coords=[ccdoy.time[1::], ccdoy.latitude, ccdoy.longitude],
                         dims=['time', 'latitude', 'longitude'])
    diff2 = xr.DataArray(uudoy.values[1::, :, :] - uudoy.values[0:-1, :, :],
                         coords=[uudoy.time[1::], uudoy.latitude, uudoy.longitude],
                         dims=['time', 'latitude', 'longitude'])

    cov, cor, slope, intercept, pval, stderr = lag_linregress_3D(diff2, diff1, lagx=0, lagy=0)

    try:
        ds['r'].values = cor.values[np.newaxis, :]
    except:
        ipdb.set_trace()
    try:
        ds['pval'].values = pval.values[np.newaxis, :]
    except:
        ipdb.set_trace()
    try:
        ds['slope'].values = slope.values[np.newaxis, :]
    except:
        ipdb.set_trace()
    #ipdb.set_trace()

    return ds


def run_doy_correlation(x):

    doy = x[0]

    udoy = x[1]
    cdoy = x[2]

    print('Doing ', doy)

    ds = xr.Dataset()

    ds['pval'] = udoy.copy(deep=True).sum('time') * np.nan
    ds['r'] = udoy.copy(deep=True).sum('time') * np.nan
    ds['slope'] = udoy.copy(deep=True).sum('time') * np.nan
    ds['intercept'] = udoy.copy(deep=True).sum('time') * np.nan

    if (np.isnan(udoy.all())) | (np.isnan(cdoy.all())):
        return ds

    uudoy = udoy#.sel(time=((udoy['time.dayofyear'] >= d1) & (udoy['time.dayofyear'] <= d2)))
    ccdoy = cdoy

    diff1 = xr.DataArray(ccdoy.values[1::, :, :] - ccdoy.values[0:-1, :, :],
                         coords=[ccdoy.time[1::], ccdoy.latitude, ccdoy.longitude],
                         dims=['time', 'latitude', 'longitude'])
    diff2 = xr.DataArray(uudoy.values[1::, :, :] - uudoy.values[0:-1, :, :],
                         coords=[uudoy.time[1::], uudoy.latitude, uudoy.longitude],
                         dims=['time', 'latitude', 'longitude'])

    cov, cor, slope, intercept, pval, stderr = lag_linregress_3D(diff2, diff1, lagx=0, lagy=0)

    try:
        ds['r'].values = cor.values
    except:
        ipdb.set_trace()
    try:
        ds['pval'].values = pval.values
    except:
        ipdb.set_trace()
    try:
        ds['slope'].values = slope.values
    except:
        ipdb.set_trace()
    try:
        ds['intercept'].values = intercept.values
    except:
        ipdb.set_trace()
    #ipdb.set_trace()

    return ds