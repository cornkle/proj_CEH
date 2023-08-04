import numpy as np
import xarray as xr
import pandas as pd
import multiprocessing
import glob
import pickle as pkl
from scipy.ndimage.measurements import label
from GLOBAL import glob_util
from utils import u_darrays, constants as cnst
import datetime


def dictionary():

    dic = {}
    vars = ['date', 'month', 'hour', 'minute', 'year', 'day', 'area', '70area', 'tmin',
            'minlon', 'minlat', 'maxlon', 'maxlat', 'clon', 'clat', 'tminlon', 'tminlat',
            'tmin', 'tmean', 'tp1', 'tp99', 'stormID', 'cloudMask', 'tir']


    for v in vars:
        dic[v] = []
    return dic


def mcs_define(array, thresh, min_area=None, max_area=None, minmax_area=None):
    """

    :param array: 2d input array
    :param thresh: cloud threshold
    :param min_area: minimum area of the cloud
    :param max_area: maximum area of the cloud
    :param minmax_area: tuple indicating only clouds bigger than tuple[0] and smaller than tuple[1]
    :return: 2d array with labelled blobs
    """
    array[array >= thresh] = 0  # T threshold maskout
    array[np.isnan(array)] = 0  # set ocean nans to 0

    labels, numL = label(array)

    u, inv = np.unique(labels, return_inverse=True)
    n = np.bincount(inv)

    goodinds = u[u!=0]

    if min_area != None:
        goodinds = u[(n>=min_area) & (u!=0)]
        badinds = u[n<min_area]

        # for b in badinds:
        #     pos = np.where(labels==b)
        #     labels[pos]=0

    if max_area != None:
        goodinds = u[(n<=max_area)  & (u!=0)]
        badinds = u[n>max_area]

    if minmax_area != None:
        goodinds = u[(n <= minmax_area[1]) & (u != 0) & (n>=minmax_area[0])]
        badinds = u[(n > minmax_area[1]) | (n < minmax_area[0])]

    if (min_area is not None) | (max_area is not None) | (minmax_area is not None):
        for b in badinds:
            pos = np.where(labels==b)
            labels[pos]=0

    return labels, goodinds


def process_tir_image(ds, data_res, t_thresh=-50, min_mcs_size=5000):
    """
    This function cuts out MCSs. By default, an MCS is defined as contiguous brightness temperature area at <=-50 degC over >= 5000km2.
    :param ctt: brightness temperature image (in degC)
    :param data_res: spatial resolution of input image (approximately, in km - this defines how many pixel are needed to define an MCS)
    :param t_thresh: temperature threshold (in degC) for considered contiguous cloud / MCS check.
    :param min_mcs_size: minimum size of contiguous cloud to be considered an MCS (in km2)
    :return: dictionary with dates and MCS characteristics from cloud top information only.
    """
    ctt = (ds['tb']).squeeze()-273.15
    min_pix_nb = min_mcs_size / data_res**2

    max_pix_nb = 300000 / data_res**2  # this is to capture satellite artefacts that come in large contiguous stripes.
    labels, goodinds = mcs_define(ctt.values, t_thresh, minmax_area=[min_pix_nb, max_pix_nb]) # 7.7x7.7km = 64km2 per pix in gridsat? 83 pix is 5000km2
    dic = dictionary()
    #plt.figure()
    #plt.pcolormesh(labels)
    #plt.colorbar()
    #plt.show()
    for g in goodinds:

        if g==0:
            continue

        pos = np.where(labels==g)
        npos = np.where(labels!=g)
        datestr = str(int(ctt['time.year'].values))+'-'+str(int(ctt['time.month'].values)).zfill(2)+'-'+str(int(ctt['time.day'].values)).zfill(2)+'_'+\
                      str(int(ctt['time.hour'].values)).zfill(2)+':'+str(int(ctt['time.minute'].values)).zfill(2)
        
        dic['date'].append(datestr)
        dic['month'].append(int(ctt['time.month']))
        dic['hour'].append(int(ctt['time.hour']))
        dic['year'].append(int(ctt['time.year']))
        dic['day'].append(int(ctt['time.day']))
        dic['minute'].append(int(ctt['time.minute']))

        storm = ctt.copy()
        storm.values[npos] = np.nan
        tmin_pos = np.nanargmin(storm.values)
        tpos_2d = np.unravel_index(tmin_pos, storm.shape)
        
        latmin = np.nanmin(ctt.lat.values[pos[0]])
        latmax = np.nanmax(ctt.lat.values[pos[0]])
        lonmin = np.nanmin(ctt.lon.values[pos[1]])
        lonmax = np.nanmax(ctt.lon.values[pos[1]])
        dic['area'].append(np.sum(np.isfinite(storm.values))*data_res**2)
        dic['70area'].append(np.sum(storm.values<=-70)*data_res**2)
        dic['minlon'].append(lonmin)
        dic['minlat'].append(latmin)
        dic['maxlon'].append(lonmax)
        dic['maxlat'].append(latmax)
        dic['clon'].append(lonmin + (lonmax - lonmin)/2)
        dic['clat'].append(latmin + (latmax - latmin)/2)
        dic['tmin'].append(np.nanmin(storm))
        dic['tminlat'].append(float(ctt.lat[tpos_2d[0]].values))
        dic['tminlon'].append(float(ctt.lon[tpos_2d[1]].values))
        dic['tmean'].append(float(np.nanmean(storm)))
        dic['tp1'].append(float(np.nanpercentile(storm, 1)))
        dic['tp99'].append(float(np.nanpercentile(storm, 99)))
        dic['stormID'].append(datestr + '_' + str(g))
        dic['cloudMask'].append(labels==g)
        dic['tir'].append(storm.values)

    # for k in dic.keys():
    #     print(k, len(dic[k]))
    return dic


def add_environment_toTable(tab, in_ds, envvar_take=[],tabvar_skip=[], rainvar_name=None, env_tformat="%Y-%m-%d %H:%M:%S", env_hour=12):
    """
    NEEDS TO BE IN SAME FILE (as in Zhe Fengs TIR/PRECIP files)
    This function saves rainfall and MCS environment variables. Rainfall is saved as mean across contiguous MCS area and max. rainfall at ~15km resolution (0.15deg)
    centred on the maximum rainfall pixel. All other variables are sampled as ~80km (0.75deg) averages around the pixel of minimum MCS cloud top temperature.

    :param file: string to netcdf file, path to file with environmental variables to be added to MCS table
    :param tab: dictionary, the associated MCS table, output from "process tir image"
    :param envvar_take: list, variable names to extract from environment netcdf file
    :param tabvar_skip: list, old table variable names to skip (remove) from new merged MCS/environment table.
    :param rainvar_name: string, rainfall variable name if rainfall is to be extracted.
    :param env_tformat: string, time format of environment netcdf file
    :param env_hour: int, time of day at which MCS environment should be extracted
    :return: copied MCS list - with optional variables removed - including saved MCS environment for provided variables.
    """

    tab_tformat = "%Y-%m-%d_%H:%M"
    dic = {}
    for k in tab.keys():
        if k in tabvar_skip:   # option to add variables to be excluded from environment table
            continue
        dic[k] = tab[k]

    ds = in_ds
    envdates = pd.to_datetime(ds.time.dt.floor('T'), format=env_tformat)
    ###### sample variables
    for tlat, tlon, date, mask, tir in zip(dic['tminlat'], dic['tminlon'], dic['date'], dic['cloudMask'], dic['tir']):


        #save cloud-wide rainfall stats, and rainfall maximum at ~0.15deg
        if rainvar_name is not None:
            tabdate = pd.to_datetime(date, format=tab_tformat)  # rainfall sampling same time as TIR

            pos = envdates == tabdate
            rain = ds[rainvar_name].isel(time=pos).where(mask).squeeze() # to mm/h
            pmax_pos = np.nanargmax(rain.values)
            ppos_2d = np.unravel_index(pmax_pos, rain.shape)
            pmax_lon = rain.lon[ppos_2d[1]]
            pmax_lat = rain.lat[ppos_2d[0]]
            pmax = rain.sel(lon=slice(pmax_lon - 0.075, pmax_lon + 0.075), lat=slice(pmax_lat - 0.075, pmax_lat + 0.075)).mean().values

            if (rainvar_name + '_mean') not in dic.keys():
                for tag in ['_mean', '_max', '_p95', '_p99']:
                    dic[rainvar_name + tag] = []
            dic[rainvar_name + '_mean'].append(float(rain.mean().values)) # full cloud mean
            dic[rainvar_name + '_max'].append(float(pmax)) # ~0.15deg rain max
            dic[rainvar_name + '_p95'].append(float(rain.quantile(0.95).values))
            dic[rainvar_name + '_p99'].append(float(rain.quantile(0.99).values))

        ## save mean environments at ~0.7deg centred on location of minimum storm temperature
        if len(envvar_take) > 0:
            tabdate = pd.to_datetime(date, format=tab_tformat).replace(hour=env_hour, minute=0)  # hour of environment sampling
            pos = envdates == tabdate
            single = ds.isel(time=pos).sel(longitude=slice(tlon - 0.375, tlon + 0.375), latitude=slice(tlat - 0.375, tlat + 0.375)).mean()
            for vt in envvar_take:
                if vt in dic.keys():
                    dic[vt].append(single[vt].values)
                else:
                    dic[vt] = [single[vt].values]
    return dic


def ERA_dictionary(tab):
    dic = {}
    vars = [
        'direction', 'tminlon', 'tminlat', 'tmin_calc',
        'tmin', 'tmean_core', 'tmean_ccs', 'tcwv', 'tgrad2m', 'tgrad925', 'smgrad', 'efgrad', 'shgrad', 'lhgrad',
        'pmax', 'pmean', 'ptot',
        'q925', 'q650', 'q850', 'era_precip', 'sm', 'ef',
        'u925', 'u650',
        'v925', 'v650',
        'w925', 'w650',
        'rh925', 'rh650',
        't925', 't650',
        'div925', 'div850', 'div650',
        'pv925', 'pv650',
        'ushear925_650', 'ushear850_650', 'ushear925_850', 'ushear100m_650',
        'vshear925_650', 'vshear850_650', 'vshear925_850', 'vshear100m_650',
        'shear925_650', 'shear850_650', 'shear925_850', 'shear100m_650',
        'cape', 't2m']

    for v in vars:
        dic[v] = []

    dummy = tab
    for v in dummy.keys():
        dic[v] = []

def run_ERA5_regional(tab):
    #run per daily chunks
    inpath = cnst.lmcs_drive + '/ERA5/hourly/'
    outtab = ERA_dictionary(tab)

    mcs_local_time = tab['lt_date']
    era_hour = 10
    etime_local = mcs_local_time.replace(hour=era_hour, minute=0).floor('S')
    edate = glob_util.LT_to_UTC_date(etime_local)

    try:
        era_pl = xr.open_dataset(
            inpath + 'pressure_levels/' + mreg + '/ERA5_' + str(edate.year) + '_' + str(edate.month).zfill(
                2) + '_' + str(edate.day).zfill(2) + '_' + mreg + '_pl.nc')
    except:
        print('ERA5 missing')
        return
    try:
        era_srfc = xr.open_dataset(
            inpath + 'surface/' + mreg + '/ERA5_' + str(edate.year) + '_' + str(edate.month).zfill(
                2) + '_' + str(edate.day).zfill(2) + '_' + mreg + '_srfc.nc')
    except:
        print('ERA5 srfc missing')
        return

    try:
        era_wi100 = xr.open_dataset(
            inpath + 'surface/'+mreg+'/100mWind_ERA5_' + str(edate.year) + '_' + str(edate.month).zfill(
                2) + '_'+ str(edate.day).zfill(2) + '_'+mreg+'_srfc.nc')
    except:
        print('ERA5 missing')
        return


    era_pl = u_darrays.flip_lat(era_pl)
    era_srfc = u_darrays.flip_lat(era_srfc)
    era_wi100 = u_darrays.flip_lat(era_wi100)

    era_pl_day = era_pl.sel(time=edate)
    era_srfc_day = era_srfc.sel(time=edate)
    era_wi100_day = era_wi100.sel(time=edate)


for regs in REGIONS:
    REGION = regs
    MONTHS = (MREGIONS[REGION])[4]
    MHOUR_SLICE = (15,21)
    composite(MHOUR_SLICE)
