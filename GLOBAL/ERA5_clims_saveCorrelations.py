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
from utils import u_darrays
from scipy import stats
import warnings
import multiprocessing
import pandas as pd
import scipy
import glob
import os

def calc_grad(da, tcoord='time', dim='y', dist=1.5):
    out = da.copy(deep=True) * np.nan
    for latid, lat in enumerate(da.latitude.values):
        for lonid, lon in enumerate(da.longitude.values):
            if dim == 'y':
                aa = da.sel(latitude=slice(lat - dist, lat + dist), longitude=lon).polyfit(dim='latitude', deg=1)
                out.values[:, latid, lonid] = aa['polyfit_coefficients'][0].values

            if dim == 'x':
                aa = da.sel(longitude=slice(lon - dist, lon + dist), latitude=lat).polyfit(dim='longitude', deg=1)
                out.values[:, latid, lonid] = aa['polyfit_coefficients'][0].values

            # out.loc[{'year' : tt, 'latitude': lat, 'longitude': lon}] = aa['polyfit_coefficients'][0]
    return out



def corr(a, b, bsingle=None, c_box=None, tcoord='time'):
    ds = xr.Dataset()
    #ipdb.set_trace()
    ds['pval'] = a.copy(deep=True).sum(tcoord) * np.nan
    ds['r'] = a.copy(deep=True).sum(tcoord) * np.nan
    ds['slope'] = a.copy(deep=True).sum(tcoord) * np.nan
    ds['intercept'] = a.copy(deep=True).sum(tcoord) * np.nan

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


def write_grads(dist=3, meridional=True):

    mainpath = cnst.lmcs_drive + '/ERA5_global_0.7/monthly/'

    era5f_pl = xr.open_mfdataset(glob.glob(mainpath + 'pressure_levels/*.nc'))
    era5f_pl = u_darrays.flip_lat(era5f_pl)

    era5f_srfc = xr.open_mfdataset(glob.glob(mainpath + 'surface/*.nc'))
    era5f_srfc = u_darrays.flip_lat(era5f_srfc)

    for mm in range(1,13):
        print('Doing month', mm)

        if meridional:
            tag = 'meridional'
            direct = 'y'
        else:
            tag = 'zonal'
            direct = 'x'

        if len(glob.glob(mainpath + 'gradients/ef_polyGrad_plusMinus'+str(dist)+'deg_'+tag+'_'+str(mm).zfill(2)+'.nc')) > 0:
            print('Month exists, continue')
            continue

        sh = era5f_srfc['sshf'].sel(time=era5f_pl['time.month']==mm).groupby('time.year').mean(['time']).squeeze().load()
        sh = sh / -86400
        t2 = era5f_srfc['t2m'].sel(time=era5f_pl['time.month']==mm).groupby('time.year').mean(['time']).squeeze().load()
        lh = era5f_srfc['slhf'].sel(time=era5f_pl['time.month']==mm).groupby('time.year').mean(['time']).squeeze().load()
        lh = lh / -86400

        era_t = era5f_pl['t'].sel(time=era5f_pl['time.month']==mm, level=925).groupby('time.year').mean(['time']).squeeze().load()

        ef = lh / (sh+lh)

        vnames = ['sh', 'ef']#, 'lh']#, 't2', 't925']

        for ids, das in enumerate([sh, ef]):#, lh, t2, era_t]):
            print('Doing var', vnames[ids])

            shpoly = calc_grad(das, tcoord='year', dim=direct, dist=dist)
            shpoly.to_netcdf(mainpath + 'gradients/'+vnames[ids]+'_polyGrad_plusMinus'+str(dist)+'deg_'+tag+'_'+str(mm).zfill(2)+'.nc')





def write_corr(dist=3, meridional=True):

    mainpath = cnst.lmcs_drive + '/ERA5_global_0.7/monthly/'
    #era5f_pl = xr.open_mfdataset('/media/ck/LStorage/global_water/ERA5_global_0.7/monthly/pressure_levels/*.nc')
    era5f_pl = xr.open_mfdataset(glob.glob(mainpath + 'pressure_levels/*.nc'))
    era5f_pl = u_darrays.flip_lat(era5f_pl)

    #era5f_pl = xr.open_mfdataset('/media/ck/LStorage/global_water/ERA5_global_0.7/monthly/pressure_levels/*.nc')
    era5f_srfc = xr.open_mfdataset(glob.glob(mainpath + 'surface/*.nc'))
    era5f_srfc = u_darrays.flip_lat(era5f_srfc)

    # era5u_mean = era5f_pl['u'].load()
    # era5v_mean = era5f_pl['v'].load()

    if meridional:
        tag = 'meridional'
    else:
        tag = 'zonal'

    vnames = ['sh', 'ef'] #, 'lh', 't2', 't925']

    for mm in range(1,13):
        try:
            era_uh = era5f_pl['u'].sel(time=era5f_pl['time.month']==mm, level=650).groupby('time.year').mean(['time']).squeeze().load()
        except:
            continue
        era_ul = era5f_pl['u'].sel(time=era5f_pl['time.month']==mm, level=925).groupby('time.year').mean(['time']).squeeze().load()
        era_vh = era5f_pl['v'].sel(time=era5f_pl['time.month']==mm, level=650).groupby('time.year').mean(['time']).squeeze().load()
        era_vl = era5f_pl['v'].sel(time=era5f_pl['time.month']==mm, level=925).groupby('time.year').mean(['time']).squeeze().load()

        era_ushear = era_uh - era_ul
        era_vshear = era_vh - era_vl
        era_shear = np.sqrt(era_ushear ** 2 + era_vshear ** 2)

        wnames = ['ushear', 'vshear', 'shear']#, 'ul', 'uh', 'vl', 'vh']

        for ids, das in enumerate(vnames):

            for idx, was in enumerate([era_ushear, era_shear, era_vshear]):  #, era_ul, era_uh, era_vl, era_vh
                outf = mainpath + 'correlations/'+vnames[ids]+'_versus_'+wnames[idx]+'_plusMinus'+str(dist)+'deg_'+tag+'_'+str(mm).zfill(2)+'.nc'
                if os.path.isfile(outf):
                    print('File exists')
                    continue

                print('Doing', outf)
                grad_da = xr.open_dataarray(mainpath + 'gradients/'+vnames[ids]+'_polyGrad_plusMinus'+str(dist)+'deg_'+tag+'_'+str(mm).zfill(2)+'.nc')
                corr_ds = corr(was, grad_da, tcoord='year')

                corr_ds.to_netcdf(outf)

                del grad_da
                del was
                del corr_ds
