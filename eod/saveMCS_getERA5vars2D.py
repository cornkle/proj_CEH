# -*- coding: utf-8 -*-


import numpy as np
import datetime as dt
import xarray as xr
import os
import ipdb
import matplotlib.pyplot as plt
from utils import u_grid, u_arrays as ua, u_darrays as uda, u_arrays as ua
import ipdb
import pandas as pd
from utils import constants as cnst
import multiprocessing
import glob
import pickle as pkl


def dictionary_storm():

    dic = {}
    vars = [
            'u925_s','u650_s','q925_s','q700_s','t2_s', 'divMoist_s','slp_s', 'd925_s','tcwv_s','CAPE_s',
            'tcwv','CAPE','u925','u650','q925','q700', 't2','divMoist', 'slp', 'd925',
            'v925', 'v925_s', 'u10', 'v10', 'u10_s', 'v10_s']

    for v in vars:
        dic[v] = []
    return dic



def perSys(clim=False):

    pool = multiprocessing.Pool(processes=5)

    pdf = pkl.load(open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA.p', 'rb'))
    pdf_all = pdf.where((pdf.clat >= 5) & (pdf.clat <= 8) & (pdf.clon >= -12) & (pdf.clon <= 10) &
                        (pdf.year >= 2000) & (pdf.year <= 2014))
    pdf_all = pdf_all.dropna()

    #ipdb.set_trace()

    dates = np.unique(np.array(pdf_all.date))
    print(dates)

    inputs = []
    for d in dates:
        pos = np.where(pdf_all.date==d)

        inputs.append((pd.Timestamp(d),pdf_all,pos[0],clim)) #,era_pl,era_srfc,

    #
    res = pool.map(get_ERA5, inputs)
    pool.close()
    # res = []
    # for ins in inputs[0:20]:
    #     res.append(get_ERA5(ins))

    #merged = ua.merge_dicts(res, merge_lists=True)

    merged = xr.concat(res, 'id')

    #test = pd.DataFrame.from_dict(merged, orient='index')
    # if clim:
    #     pkl.dump(merged, open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA_ERA5_allmonth_5-8N_2000-2014_18UTC_front_CLIM.p',
    #                        'wb'))
    # else:
    #     pkl.dump(merged, open(cnst.CLOVER_SAVES + 'StormLoc_-50_5000km_WA_ERA5_allmonth_5-8N_2000-2014_18UTC_front.p',
    #                        'wb'))

    merged.to_netcdf(cnst.CLOVER_SAVES + '2d_ERA5_comp1000x1000km_-50_5000km_WA_5-8N_12W-10E.nc')

def get_ERA5(inputs):

    date = inputs[0]
    indic = inputs[1]
    # era_pl = inputs[2]
    # era_srfc = inputs[3]
    ids = inputs[2]
    clim = inputs[3]

    dic = dictionary_storm()

    if clim:
        time = str(date.month).zfill(2) + '-12'
        stormtime = str(date.month).zfill(2) + '-18'

        pl_str = 'ERA5/monthly/synop_selfmade/CLIM_2000-2014/pressure_levels/ERA5_2000-2014_CLIM_'
        srfc_str = 'ERA5/monthly/synop_selfmade/CLIM_2000-2014/surface/ERA5_2000-2014_CLIM_'

        try:
            print('Open '+ cnst.local_data + pl_str + time + '_pl.nc')
            era_day_pl = xr.open_dataset(
                cnst.local_data + pl_str + time + '_pl.nc')
        except (TypeError, IndexError, KeyError):
            print('Era missing:', date)
            #             for k in dic.keys():
            #                 dic[k].append(np.nan)
            return

        era_day_sf = xr.open_dataset(
            cnst.local_data + srfc_str + time + '_srfc.nc')
        try:
            print('Open ' + cnst.local_data + srfc_str + stormtime + '_srfc.nc')
            era_day_sft = xr.open_dataset(
                cnst.local_data + srfc_str + stormtime + '_srfc.nc')
        except IndexError:
            return
        era_day_plt = xr.open_dataset(
            cnst.local_data + pl_str + stormtime + '_pl.nc')

    else:

        era_pl = xr.open_dataset(cnst.local_data + 'ERA5/hourly/pressure_levels/ERA5_' +str(date.year) + '_' + str(date.month).zfill(2) + '_pl.nc')
        era_srfc = xr.open_dataset(
            cnst.local_data + 'ERA5/hourly/surface/ERA5_' + str(date.year) + '_' + str(date.month).zfill(2) + '_srfc.nc')

        era_pl = era_pl.rename({'longitude' : 'lon', 'latitude' : 'lat'})
        era_srfc = era_srfc.rename({'longitude' : 'lon', 'latitude' : 'lat'})

        time = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + 'T12'
        stormtime = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + 'T18'

        try:
            era_day_pl = era_pl.sel(time=time).isel(time=0)
        except (TypeError, IndexError, KeyError):
            print('Era missing:', date)
            #             for k in dic.keys():
            #                 dic[k].append(np.nan)
            return

        era_day_sf = era_srfc.sel(time=time).isel(time=0)
        try:
            era_day_sft = era_srfc.sel(time=stormtime).isel(time=0)
        except IndexError:
            return
        era_day_plt = era_pl.sel(time=stormtime).isel(time=0)

    era_day_sf = uda.flip_lat(era_day_sf)
    era_day_pl = uda.flip_lat(era_day_pl)

    for id in ids:

        print('Doing', date)


        # elat = indic.clat[id]
        # elon = indic.clon[id]

        elat = indic.clat[id]
        elon = indic.minlon[id]

        # dic['dates'].append(date)
        # dic['lat'].append(elat)
        # dic['lon'].append(elon)
        # ipdb.set_trace()
        point = era_day_pl.sel(lat=elat, lon=elon, method='nearest')

        posx = int(np.where(era_day_sf.lon == point.lon)[0])
        posy = int(np.where(era_day_sf.lat == point.lat)[0])

        posxx = int(np.where(era_day_pl.lon == point.lon)[0])
        posyy = int(np.where(era_day_pl.lat == point.lat)[0])

        dist = 18  # ca 200km i.e. 45 *4.4km

        ds_pl = era_day_pl.apply(uda.cut_box, xpos=posxx, ypos=posyy, dist=dist)
        ds_srfc = era_day_sf.apply(uda.cut_box, xpos=posx, ypos=posy, dist=dist)

        ds_plt = era_day_plt.apply(uda.cut_box, xpos=posxx, ypos=posyy, dist=dist)
        ds_srfct = era_day_sft.apply(uda.cut_box, xpos=posx, ypos=posy, dist=dist)

        try:
            #ipdb.set_trace()
            dic['u925'].append(xr.DataArray(ds_pl['u'].sel(level=925).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id], 'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        except ValueError:
            print('Value problem')
            continue
        # ipdb.set_trace()
        dic['u650'].append(xr.DataArray(ds_pl['u'].sel(level=650).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id], 'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['v925'].append(xr.DataArray(ds_pl['v'].sel(level=925).values,
                                        coords={'tmean': indic.tmean[id], 'tmin': indic.tmin[id],
                                                'area': indic.area[id], 'time': date, 'lat': ds_pl['y'].values,
                                                'lon': ds_pl['x'].values}, dims=['lat', 'lon']))


        dic['q925'].append(xr.DataArray(ds_pl['q'].sel(level=925).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['q700'].append(xr.DataArray(ds_pl['q'].sel(level=700).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['d925'].append(xr.DataArray(ds_pl['d'].sel(level=925).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['CAPE'].append(xr.DataArray(ds_srfc['cape'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['tcwv'].append(xr.DataArray(ds_srfc['tcwv'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['t2'].append(xr.DataArray(ds_srfc['t2m'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['divMoist'].append(xr.DataArray(ds_srfc['p84.162'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['slp'].append(xr.DataArray(ds_srfc['msl'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['u10'].append(xr.DataArray(ds_srfc['u10'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['v10'].append(xr.DataArray(ds_srfc['v10'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))


        dic['u925_s'].append(
            xr.DataArray(ds_plt['u'].sel(level=925).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['v925_s'].append(
            xr.DataArray(ds_plt['v'].sel(level=925).values,
                         coords={'tmean': indic.tmean[id], 'tmin': indic.tmin[id], 'area': indic.area[id], 'time': date,
                                 'lat': ds_pl['y'].values, 'lon': ds_pl['x'].values}, dims=['lat', 'lon']))

        dic['u650_s'].append(
            xr.DataArray(ds_plt['u'].sel(level=650).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['q925_s'].append(
            xr.DataArray(ds_plt['q'].sel(level=925).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['q700_s'].append(
            xr.DataArray(ds_plt['q'].sel(level=700).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['d925_s'].append(
            xr.DataArray(ds_plt['d'].sel(level=925).values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['CAPE_s'].append(xr.DataArray(ds_srfct['cape'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['tcwv_s'].append(xr.DataArray(ds_srfct['tcwv'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['t2_s'].append(xr.DataArray(ds_srfct['t2m'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['divMoist_s'].append(xr.DataArray(ds_srfct['p84.162'].values, coords={'tmean':indic.tmean[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['slp_s'].append(xr.DataArray(ds_srfct['msl'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['u10_s'].append(xr.DataArray(ds_srfct['u10'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))
        dic['v10_s'].append(xr.DataArray(ds_srfct['v10'].values, coords={'tmean':indic.tmean[id],'tmin':indic.tmin[id],'area':indic.area[id],'time': date, 'lat': ds_pl['y'].values, 'lon':ds_pl['x'].values}, dims=['lat', 'lon']))

        print('DID', date)

    ds = xr.Dataset()

    try:
        for k in dic.keys():
            print('doing', k)

            ds[k] = xr.concat(dic[k], 'id')

    except ValueError:
        return


    return ds

