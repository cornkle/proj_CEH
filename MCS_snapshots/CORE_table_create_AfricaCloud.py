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
import ipdb


def dictionary():
    dic = {}
    vars = ['date', 'month', 'hour', 'minute', 'year', 'day', 'area', '70area', '60area', '50area',
            'minlon', 'minlat', 'maxlon', 'maxlat', 'clon', 'clat', 'tminlon', 'tminlat',
            'tmin', 'tmean', 'tp1', 'tp99', 'coreID', 'coreMask', 'tir',
            'powmean', 'powmax', 'MCSID'
            ]

    for v in vars:
        dic[v] = []
    return dic


def core_define(array, thresh, min_area=None, max_area=None, minmax_area=None):
    """

    :param array: 2d input array
    :param thresh: cloud threshold
    :param min_area: minimum area of the cloud
    :param max_area: maximum area of the cloud
    :param minmax_area: tuple indicating only clouds bigger than tuple[0] and smaller than tuple[1]
    :return: 2d array with labelled blobs
    """
    array[array <= thresh] = 0  # T threshold maskout
    array[np.isnan(array)] = 0  # set ocean nans to 0

    labels, numL = label(array)

    u, inv = np.unique(labels, return_inverse=True)
    n = np.bincount(inv)

    goodinds = u[u != 0]

    if min_area != None:
        goodinds = u[(n >= min_area) & (u != 0)]
        badinds = u[n < min_area]

        # for b in badinds:
        #     pos = np.where(labels==b)
        #     labels[pos]=0

    if max_area != None:
        goodinds = u[(n <= max_area) & (u != 0)]
        badinds = u[n > max_area]

    if minmax_area != None:
        goodinds = u[(n <= minmax_area[1]) & (u != 0) & (n >= minmax_area[0])]
        badinds = u[(n > minmax_area[1]) | (n < minmax_area[0])]

    if (min_area is not None) | (max_area is not None) | (minmax_area is not None):
        for b in badinds:
            pos = np.where(labels == b)
            labels[pos] = 0

    return labels, goodinds


def process_core_image(ds, data_res, t_thresh=0.1, min_mcs_size=9):
    """
    This function cuts out MCSs. By default, an MCS is defined as contiguous brightness temperature area at <=-50 degC over >= 5000km2.
    :param ctt: brightness temperature image (in degC)
    :param data_res: spatial resolution of input image (approximately, in km - this defines how many pixel are needed to define an MCS)
    :param t_thresh: temperature threshold (in degC) for considered contiguous cloud / MCS check.
    :param min_mcs_size: minimum size of contiguous cloud to be considered an MCS (in km2)
    :return: dictionary with dates and MCS characteristics from cloud top information only.
    """
    ctt = (ds['tir']).squeeze()
    ctt.values[ctt>-40] = 0
    ctt.values[np.isnan(ctt)]=0
    mcs_label, mL = label(ctt.values)

    cores = (ds['cores']).squeeze()
    cores.values[cores<0.1] = 0
    cores.values[np.isnan(cores)]=0
    labels, goodinds = label(cores.values)
    

    min_pix_nb = min_mcs_size / data_res ** 2

    max_pix_nb = 300000 / data_res ** 2  # this is to capture satellite artefacts that come in large contiguous stripes.
    labels, goodinds = core_define(cores.values, t_thresh, minmax_area=[min_pix_nb, max_pix_nb])  # 7.7x7.7km = 64km2 per pix in gridsat? 83 pix is 5000km2
    
    dic = dictionary()
    # plt.figure()
    # plt.pcolormesh(labels)
    # plt.colorbar()
    # plt.show()
    ipdb.set_trace()
    for glon, glat in zip(ds.maxlon, ds.maxlat):
        
        point = (glon, glat)
        points = np.array(list(zip(ds['dlon'].values, ds['dlat'].values)))
        dist_2 = np.sum((points - point) * (points - point), axis=1)
        pmax_pos = np.argmin(dist_2)

        ipdb.set_trace()
        g = labels[pmax_pos[0],pmax_pos[1]]
        pos = np.where(labels == g)
        npos = np.where(labels != g)
        datestr = str(int(ctt['time.year'].values)) + '-' + str(int(ctt['time.month'].values)).zfill(2) + '-' + str(int(ctt['time.day'].values)).zfill(2) + '_' + \
                  str(int(ctt['time.hour'].values)).zfill(2) + ':' + str(int(ctt['time.minute'].values)).zfill(2)

        dic['date'].append(datestr)
        dic['month'].append(int(ctt['time.month']))
        dic['hour'].append(int(ctt['time.hour']))
        dic['year'].append(int(ctt['time.year']))
        dic['day'].append(int(ctt['time.day']))
        dic['minute'].append(int(ctt['time.minute']))

        storm = ctt.copy().astype('float')
        storm.values[npos] = np.nan
        tmin_pos = np.nanargmin(storm.values)
        tpos_2d = np.unravel_index(tmin_pos, storm.shape)

        power = cores.copy().astype('float')
        power.values[npos] = np.nan
        #ipdb.set_trace()
        mcs_id = mcs_label[tpos_2d[0], tpos_2d[1]]

        latmin = np.nanmin(ds['dlat'].values[pos[0]])
        latmax = np.nanmax(ds['dlat'].values[pos[0]])
        lonmin = np.nanmin(ds['dlon'].values[pos[1]])
        lonmax = np.nanmax(ds['dlon'].values[pos[1]])
        dic['area'].append(np.sum(np.isfinite(storm.values)) * data_res ** 2)
        dic['70area'].append(np.sum(storm.values <= -70) * data_res ** 2)
        dic['60area'].append(np.sum(storm.values <= -60) * data_res ** 2)
        dic['50area'].append(np.sum(storm.values <= -50) * data_res ** 2)
        dic['minlon'].append(lonmin)
        dic['minlat'].append(latmin)
        dic['maxlon'].append(lonmax)
        dic['maxlat'].append(latmax)
        dic['clon'].append(lonmin + (lonmax - lonmin) / 2)
        dic['clat'].append(latmin + (latmax - latmin) / 2)
        dic['tmin'].append(np.nanmin(storm))
        dic['tminlat'].append(float(ctt.lat[tpos_2d[0]].values))
        dic['tminlon'].append(float(ctt.lon[tpos_2d[1]].values))
        dic['tmean'].append(float(np.nanmean(storm)))
        dic['powmean'].append(float(np.nanmean(power)))
        dic['powmax'].append(float(np.nanmax(power)))
        dic['tp1'].append(float(np.nanpercentile(storm, 1)))
        dic['tp99'].append(float(np.nanpercentile(storm, 99)))
        dic['coreID'].append(datestr + '_' + str(g))
        dic['MCSID'].append(datestr + '_' + str(mcs_id))
        dic['coreMask'].append(labels == g)
        dic['tir'].append(storm.values)

    # for k in dic.keys():
    #     print(k, len(dic[k]))
    return dic


def add_environment_toTable(tab, in_ds, envvar_take=[], tabvar_skip=[], rainvar_name=None, env_tformat="%Y-%m-%d %H:%M:%S", env_hour=12):
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
        if k in tabvar_skip:  # option to add variables to be excluded from environment table
            continue
        dic[k] = tab[k]

    ds = in_ds
    envdates = pd.to_datetime(ds.time.dt.floor('T'), format=env_tformat)
    ###### sample variables
    for tlat, tlon, date, mask, tir in zip(dic['tminlat'], dic['tminlon'], dic['date'], dic['cloudMask'], dic['tir']):

        # save cloud-wide rainfall stats, and rainfall maximum at ~0.15deg
        if rainvar_name is not None:
            tabdate = pd.to_datetime(date, format=tab_tformat)  # rainfall sampling same time as TIR

            pos = envdates == tabdate
            rain = ds[rainvar_name].isel(time=pos).where(mask).squeeze()  # to mm/h
            pmax_pos = np.nanargmax(rain.values)
            ppos_2d = np.unravel_index(pmax_pos, rain.shape)
            pmax_lon = rain.lon[ppos_2d[1]]
            pmax_lat = rain.lat[ppos_2d[0]]
            pmax = rain.sel(lon=slice(pmax_lon - 0.075, pmax_lon + 0.075), lat=slice(pmax_lat - 0.075, pmax_lat + 0.075))
            #  ipdb.set_trace()
            pmax = pmax.mean().values
            if (rainvar_name + '_mean') not in dic.keys():
                for tag in ['_mean', '_max', '_p95', '_p99']:
                    dic[rainvar_name + tag] = []
            dic[rainvar_name + '_mean'].append(float(rain.mean().values))  # full cloud mean
            dic[rainvar_name + '_max'].append(float(pmax))  # ~0.15deg rain max
            dic[rainvar_name + '_p95'].append(float(rain.quantile(0.95).values))
            dic[rainvar_name + '_p99'].append(float(rain.quantile(0.99).values))

        ## save mean environments at ~0.7deg centred on location of minimum storm temperature
        if len(envvar_take) > 0:
            tabdate = pd.to_datetime(date, format=tab_tformat).replace(hour=env_hour, minute=0)  # hour of environment sampling
            pos = envdates == tabdate
            single = ds.isel(time=pos).sel(longitude=slice(tlon - 0.375, tlon + 0.375), latitude=slice(tlat - 0.375, tlat + 0.375))
            #  ipdb.set_trace()
            single = single.mean()
            for vt in envvar_take:
                if vt in dic.keys():
                    dic[vt].append(single[vt].values)
                else:
                    dic[vt] = [single[vt].values]
    return dic
