# -*- coding: utf-8 -*-
"""
Created on Thu May 19 15:50:34 2016

Routines rewriting stuff!

@author: cornkle
"""

import numpy as np
import os
from utils import u_arrays as uarr
import pandas as pd
from utils import u_time as ut, u_interpolate as uint, constants as cnst, u_darrays
import datetime as dt
from utils import u_grid
import xarray as xr
import ipdb
from scipy.ndimage.measurements import label

#========================================================================================
# Rewrites 580x1640 msg lat lon to something nice (lat lon from blobs)
#========================================================================================
def rewriteMsgLonLat_WA():
    llFile = cnst.network_data + 'data/OBS/MSG_WA30/MSG_1640_580_lat_lon.gra'

    llShape = (580,1640)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)

    ipdb.set_trace()
    lon = ll[0:580*1640]
    lat = ll[580*1640:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = cnst.network_data + 'data/OBS/MSG_WA30/MSG_1640_580_lat_lon'
    np.savez(llsavefile,lon=lon,lat=lat)


#========================================================================================
# Rewrites 580x1640 msg lat lon to something nice (lat lon from blobs)
#========================================================================================
def rewriteMFGLonLat_WA():
    llFile = cnst.network_data + 'data/OBS/MFG_JJAS/msat_grid_463x1094.gra'

    llShape = (463,1094)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)
    lon = ll[0:463*1094]
    lat = ll[463*1094:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = cnst.network_data + 'data/OBS/MFG_JJAS/MFG_1094_463_lat_lon'
    np.savez(llsavefile,lon=lon,lat=lat)

#========================================================================================
# Rewrites 504x1368 msg lat lon to something nice (lat lon from blobs)
#========================================================================================
def rewriteMSG_MAMON_LonLat_WA():
    llFile = cnst.network_data + 'data/OBS/MSG_MAMON/lat_lon_1368_504.gra'

    llShape = (504,1368)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)

    lat = ll[0:504*1368]
    lon = ll[504*1368:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = cnst.network_data + 'data/OBS/MSG_MAMON/MSG_1368_504_lat_lon'
    np.savez(llsavefile,lon=lon,lat=lat)


#========================================================================================
# Rewrites 580x1640 msg lat lon to something nice (lat lon from blobs)
#========================================================================================
def rewriteMFG_MAMON_LonLat_WA():
    llFile = cnst.network_data + 'data/OBS/MFG_MAMON/msat_grid_336x914.gra'

    llShape = (336,914)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)

    lon = ll[0:336*914]
    lat = ll[336*914:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = cnst.network_data + 'data/OBS/MFG_MAMON/MFG_336_914_lat_lon'
    np.savez(llsavefile,lon=lon,lat=lat)


#========================================================================================
# Reads Chris' METEOSAT cell table files, gives them a header and just reads columns
# of interest for comparison with TRMM. Creates an easy to read object.

# Table files: every 15/30 minutes, identified single systems with a random number
# Columns: 5-area (km2), 6/7-lat/lon, 8-min col number (start at 1),
# 12 - mean T of cell, 13 -Temp treshold to identify cell
#========================================================================================
def parseCellTab(tab):

    dt.datetime.strptime('2004_0601_1315', '%Y_%m%d_%H%M')
  #  parser = lambda date: pd.datetime.strptime(date, '%Y_%m%d_%h%M')

    df = pd.read_csv(tab, sep='\s+', header=None, converters={'Mday': lambda x: str(x)}, names=["Year", "Mday", "Slot", "Pixel", "Area", "Lat", "Lon", "Mincol", "a", "b", "c", "Temp", "Tresh"])

    sec = df["Slot"]*30.*60.
    t = ut.sec_to_time(sec[0])
    df["Hour"] = df["Slot"]*0+t.hour
    df["Minute"] = df["Slot"]*0+t.minute
    df["Hour"] = df.Hour.map("{:02}".format)
    df["Minute"] = df.Minute.map("{:02}".format)

    small=df.loc[:, ["Pixel", "Area", "Lat", "Lon", "Mincol", "Temp", "Tresh"]]

    small["Date"] = df.Year.astype(str).str.cat(df.Mday.astype(str), sep='_')
    small["Date"] = small.Date.astype(str).str.cat(df.Hour.astype(str), sep='_')
    small["Date"] = small.Date.astype(str).str.cat(df.Minute.astype(str), sep='')
    small["Date"] = pd.to_datetime(small["Date"], format='%Y_%m%d_%H%M')

    return small


#=========================================================================================
# Rewrites the METEOSAT cell table files from parseCellTables() object to nice *.txt files.
#=========================================================================================
def rewriteCellTab():

    path = "/users/global/cornkle/data/OBS/meteosat/bigcell_area_table/"
    out = path + 'rewrite/'
    print(out)
    os.system('rm '+out+'*.txt')
    ok=uarr.locate("*.txt", path)

    for a in ok:
        print('Doing '+a)
        tab = parseCellTab(a)
        minute=tab["Date"][0].minute
        hour=tab["Date"][0].hour
        tab.to_csv(out+'cell_40c_'+str(hour).zfill(2)+str(minute).zfill(2)+'_JJAS.txt')

#========================================================================================
# Rewrites 350x728 msg lat lon to something nice (lat lon from blobs)
#========================================================================================
def rewriteMsgLonLat_Sahel():
    llFile = '/users/global/cornkle/data/OBS/meteosat_SA15/Sahel_728_350_lat_lon.gra'

    llShape = (350,728)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)
    lon = ll[0:350*728]
    lat = ll[350*728:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = '/users/global/cornkle/data/OBS/meteosat_SA15/MSG_728_350_lat_lon'
    np.savez(llsavefile,lon=lon,lat=lat)


#========================================================================================
# Rewrites 201x326 msg lat lon to something nice (lat lon from blobs)
#  file: lat lon grads file
#  ny : pixel in y direction
#  nx : pixel in x direction
#========================================================================================
def rewriteMsgLonLat(file, nx, ny):
    llFile = file

    llShape = (ny,nx)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)
    lon = ll[0:ny*nx]
    lat = ll[ny*nx:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = file.replace('.gra', '')
    np.savez(llsavefile,lon=lon,lat=lat)

#========================================================================================
# Rewrites modis lat lon to something nice (lat lon from blobs)
#  file: lat lon grads file
#  ny : pixel in y direction
#  nx : pixel in x direction
#========================================================================================
def rewriteMODISLstLonLat(file, nx, ny):
    llFile = file

    llShape = (ny,nx)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)
    lat = ll[0:ny*nx]
    lon = ll[ny*nx:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = file.replace('.gra', '')
    np.savez(llsavefile,lon=lon,lat=lat)

#========================================================================================
# Rewrites modis lat lon to something nice (lat lon from blobs)
#  file: lat lon grads file
#  ny : pixel in y direction
#  nx : pixel in x direction
#========================================================================================

def rewriteLSTA_toNetcdf(file, interp, write=None):

    out = file.replace('lsta_raw_binary_new', 'lsta_netcdf_new')
    if '.gz' in out:
        out = out.replace('.gra.gz', '.nc')
    else:
        out = out.replace('.gra', '.nc')

    print('Doing '+file)

    # if os.path.isfile(out):
    #     print('File exists, continue: ', out)
    #     return

    ll = np.load(cnst.network_data + 'data/OBS/MSG_LSTA/lsta_728_348_lat_lon.npz')

    blat = ll['lat']
    blon = ll['lon']

    ffile = os.path.split(file)[-1]

    llist = ffile.split('.')
    yr = int(llist[0][-8:-4])
    mon = int(llist[0][-4:-2])
    day = int(llist[0][-2::])

    date = [pd.datetime(yr, mon, day, 0, 0)]

    rrShape = (blat.shape[0],blat.shape[1])

    rr = np.fromfile(file, dtype=np.float32())#dtype='>f')
    addVar = False

    latmin = np.min(blat)
    latmax = np.max(blat)
    lonmin = -9.98 #np.min(blon)
    lonmax = 9.98 #np.max(blon)
    dist = np.round(np.float(np.mean(blon[0,:][1::] - blon[0,:][0:-1])), decimals=4)

    lat_regular = np.arange(latmin + 10*dist, latmax - 10*dist , dist)
    lon_regular = np.arange(lonmin , lonmax  , dist)

    ds = xr.Dataset(coords={'time':  date,'lat':  interp['y'],'lon': interp['x']})

    if rrShape[0]*rrShape[1]*2 == rr.size:
        rr2 = rr[int(rr.size/2)::].copy()
        rr2.shape = rrShape
        #rr2_regridded = uint.interpolate_data(rr2, interp['inds'], interp['weights'], interp['shape'])
        rr2_regridded = uint.griddata_int(rr2, blon, blat, lon_regular, lat_regular, isll=False, method='nearest')


        ds['NbSlot'] = (('time', 'lat', 'lon'), rr2_regridded[None, ...])

        rr3 = rr[0:int(rr.size/2)].copy()
        rr3.shape = rrShape
        rr3[rr3 == -999] = np.nan

        #rr3_regridded = uint.interpolate_data(rr3, interp['inds'], interp['weights'], interp['shape'])
        rr3_regridded = uint.griddata_int(rr3, blon, blat, lon_regular, lat_regular, isll=False, method='nearest')

        ds['LSTA'] = (('time', 'lat', 'lon'), rr3_regridded[None, ...])

    else:
        rr.shape = rrShape
        rr[rr == -999] = np.nan
        #rr_regridded = uint.interpolate_data(rr, interp['inds'], interp['weights'], interp['shape'])
        rr_regridded = uint.griddata_int(rr, blon, blat, lon_regular, lat_regular, isll=False, method='nearest')

        ds['LSTA'] = (('time', 'lat', 'lon'), rr_regridded[None, ...])


    # rr.shape = rrShape
    #
    #
    #
    # da = xr.DataArray(rr[None, ...], ,
    #                   dims=['time', 'lat', 'lon'])#.isel(time=0)
    #
    # ds = xr.Dataset({'LSTA': da})
    # if addVar:
    #     ds['NbSlot'] = (('time', 'lat', 'lon'), rr2[None, ...])

    if write is not None:
        try:
            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in ds.data_vars}
            ds.to_netcdf(path=out, mode='w', encoding=encoding, format='NETCDF4')

        except OSError:
            print('Did not find ' + out)
            print('Out directory not found')
        print('Wrote ' + out)
    return ds, out

#========================================================================================
# Rewrites modis lat lon to something nice (lat lon from blobs)
#  file: lat lon grads file
#  ny : pixel in y direction
#  nx : pixel in x direction
#========================================================================================
def rewriteLSTAClim_toNetcdf(file):

    out = file.replace('modis_raw_binary', 'modis_netcdf')
    if '.gz' in out:
        out = out.replace('.gra.gz', '.nc')
    else:
        out = out.replace('.gra', '.nc')

    print('Doing '+file)

    if os.path.isfile(out):
        print('File exists, continue: ', out)
        return

    ll = np.load('/users/global/cornkle/data/OBS/MSG_LSTA/lsta_728_348_lat_lon.npz')

    blat = ll['lat']
    blon = ll['lon']
    ffile = os.path.split(file)[-1]

    llist = ffile.split('.')
    yr = int(llist[0][-4::])


    rr = np.fromfile(file, dtype='>f')
    single_image = blat.size
    nb_image = int(rr.size / single_image)

    rr = np.reshape(rr, (nb_image,blat.shape[0],blat.shape[1]))
    date = pd.date_range(str(yr)+'-01-01 06:30:00', freq='15min', periods=nb_image)


    da = xr.DataArray(rr, coords={'time':  date,
                                             'lat':  blat[:,0],
                                             'lon': blon[0,:]},
                      dims=['time', 'lat', 'lon'])#.isel(time=0)

    ds = xr.Dataset({'LST': da})

    try:
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(path=out, mode='w', encoding=encoding, format='NETCDF4')

    except OSError:
        print('Did not find ' + out)
        print('Out directory not found')
    print('Wrote ' + out)
    return


def rewrite_AMSRE(file):

    out = file.replace('raw', 'nc')
    out = out.replace('.gra', '.nc')

    if os.path.isfile(out):
        return

    day=True
    if day:
        hour = 13
    else:
        hour = 1

    print('Doing ' + file)
    #full_globe, 0.25deg
    blat = np.arange(-89.875,90, 0.25)
    blon = np.arange(-179.875,180, 0.25)

    ffile = os.path.split(file)[-1]

    llist = ffile.split('.')
    yr = int(llist[0][-8:-4])
    mon = int(llist[0][-4:-2])
    day = int(llist[0][-2::])

    date = [pd.datetime(yr, mon, day, hour, 0)]

    rrShape = (blat.shape[0], blon.shape[0])

    rr = np.array(np.fromfile(file, dtype=np.int8()), dtype=float)
    rr.shape = rrShape
    rr = np.flip(rr,axis=0)

    da = xr.DataArray(rr[None, ...], coords={'time': date,
                                             'lat': blat,
                                             'lon': blon},
                      dims=['time', 'lat', 'lon'])  # .isel(time=0)

    da.values[da.values==-1] = np.nan
    if np.sum(da.values) == np.nan:
        return

    ds = xr.Dataset({'SM': da})

    ds = ds.sel(lon=slice(-18,30), lat=slice(0,27))  # ds.sel(lon=slice(-17,20), lat=slice(0,20)) #lon=slice(-20,55), lat=slice(-40,40) AFRICA

    try:
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(path=out, mode='w', encoding=encoding, format='NETCDF4')

    except OSError:
        print('Did not find ' + out)
        print('Out directory not found')
    print('Wrote ' + out)
    return ds

def rewrite_AMSR2(file):

    out = file.replace('raw', 'nc')
    #out = out.replace('.nc4', '.nc')
    out = out.replace('LPRM-AMSR2', 'AMSR2')

    if '_A_' in file:
        day=True
    else:
        day = False

    if day:
        hour = 13
    else:
        hour = 1

    cut = os.path.basename(out)
    path = os.path.dirname(out)
    path = path.replace(path[-6::], '')
    pieces = cut.split('_')
    time = (pieces[5])[0:8]
    out = path + os.sep + pieces[0] + '_' +  pieces[1] +'_LPRMv05_' +  pieces[2] + '_'  + time + '.nc'

    #ipdb.set_trace()
    if os.path.isfile(out):
        return

    yr = int(time[0:4])
    mon = int(time[4:6])
    day = int(time[6:8])
    date = [pd.datetime(yr, mon, day, hour, 0)]

    ds = xr.open_dataset(file)
    ds = u_darrays.flip_lat(ds)

    da = xr.DataArray((ds['soil_moisture_c1']).values.T[None, ...], coords={'time': date,
                                             'lat': ds.Latitude.values,
                                             'lon': ds.Longitude.values},
                      dims=['time', 'lat', 'lon'])  # .isel(time=0)

    da.values[da.values<-1] = np.nan
    if np.sum(da.values) == np.nan:
        return

    ds = xr.Dataset({'SM': da})

    ds = ds.sel(lon=slice(-18,30), lat=slice(0,27)) #lon=slice(-20,55), lat=slice(-40,40) AFRICA

    try:
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(path=out, mode='w', encoding=encoding, format='NETCDF4')

    except OSError:
        print('Did not find ' + out)
        print('Out directory not found')
    print('Wrote ' + out)
    return ds


def rewrite_CMORPH(file):

    out = file.replace('raw', 'nc')
    out = out.replace('.gra', '.nc')
    day=True

    hours = [0,3,6,9,12,15,18,21]

    print('Doing ' + file)
    #full_globe, 0.25deg
    blat = np.arange(-59.875,60, 0.25)
    blon = np.arange(0.125, 360, 0.25)

    ffile = os.path.split(file)[-1]

    rrShape = (8, blat.shape[0], blon.shape[0])

    rr = np.array(np.fromfile(file, dtype=np.int16()))

    rr.shape = rrShape
    rr = np.array(rr/100., dtype=np.float)
    #pdb.set_trace()
    #rr = np.flip(rr,axis=0)


    llist = ffile.split('.')

    yr = int(llist[2][-8:-4])
    mon = int(llist[2][-4:-2])
    day = int(llist[2][-2::])
    date = []
    for h in hours:
        datel = pd.datetime(yr, mon, day, h, 0)
        date.append(datel)

    date = pd.to_datetime(date)
    dates = pd.date_range(date[0], date[-1], freq='3H')

    blon[blon>180] = blon[blon>180]-360

    da = xr.DataArray(rr, coords={'time': dates,
                                             'lat': blat,
                                             'lon': blon},
                      dims=['time', 'lat', 'lon'])  # .isel(time=0)

    #
    da.values[da.values==-1] = np.nan
    # if np.sum(da.values) == np.nan:
    #     return

    da = da.roll(lon=np.sum(blon<0))
    ds = xr.Dataset({'pr': da})

    ds = ds.sel(lon=slice(-17.5,20), lat=slice(2,20)) #lon=slice(-20,55), lat=slice(-40,40) AFRICA

    try:
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(path=out, mode='w', encoding=encoding, format='NETCDF4')

    except OSError:
        print('Did not find ' + out)
        print('Out directory not found')
    print('Wrote ' + out)
    return ds
