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
from utils import u_time as ut
import datetime as dt

#========================================================================================
# Rewrites 580x1640 msg lat lon to something nice (lat lon from blobs)
#========================================================================================
def rewriteMsgLonLat_WA():
    llFile = '/users/global/cornkle/data/OBS/meteosat_WA30/MSG_1640_580_lat_lon.gra'

    llShape = (580,1640)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)
    lon = ll[0:580*1640]
    lat = ll[580*1640:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = '/users/global/cornkle/data/OBS/meteosat_WA30/MSG_1640_580_lat_lon'
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
#========================================================================================
def rewriteMsgLonLat_IODC_TIR():
    llFile = '/users/global/cornkle/Emma/TIR/lat_lon_nx326_ny201.gra'

    llShape = (201,326)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)
    lon = ll[0:201*326]
    lat = ll[201*326:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = '/users/global/cornkle/Emma/TIR/lat_lon_nx326_ny201'
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