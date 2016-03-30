# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 16:56:03 2016

@author: cornkle
"""

import numpy as np
import glob
import itertools
from utils import u_time as ut
from utils import u_arrays as uarr 
import pandas as pd
from pandas import *
import datetime as dt
from os import remove

#==============================================================================
# Reads TRMM 2A25 files to get a list with valid file names (files of interest 
# where time is close to full or half hour to be potentially correspondant to METEOSAT).
#
# Output: Filelist of TRMM files over West Africa at full/half hour; 
# corrsponding time rounded to half or full hour
#==============================================================================
def read_TRMMswath():
    
    path = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
    pattern = path + '{0:d}/{1:02d}/2A25.{0:d}{1:02d}{2:02d}.*.7_rain_f4.gra' 
    
    yrange=[2013]
    
    fdic = {'fpath' : [], 'date' : ut.date_list()}
    
    for yr,mo,dy in itertools.product(yrange,range(6,10),range(1,32)):
        
        a=''
        date = np.array([yr,mo,dy])
        
        
        extpath = pattern.format(date[0],date[1],date[2]) 
        a=glob.glob(extpath)
        
        if a: 
            for eachfile in a:
                
                rain_str = eachfile.replace('_rain_f4', '')
                time_str = eachfile.replace('_rain_f4', '_time')
                
                rr = np.fromfile(time_str,dtype=np.float32) # seconds of day
                
                secmean = rr.mean()
                t = ut.sec_to_time(secmean)
                
                # test whether close to 30mins or full hour
 
                minute = 0   # guessing that t.minute is shortly after full
                
                if t.minute > 15 and t.minute < 45: 
                    minute = 30 
                 
                fdic['fpath'].append(rain_str)
                fdic['date'].add(yr,mo,dy,t.hour,minute,0)
    
    return fdic
    

#==============================================================================
# Reads the METEOSAT cell table files, gives them a header and just reads columns
# of interest for comparison with TRMM. Creates an easy to read object. 
#==============================================================================
def parseCellTables(tab):
     
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
    

#==============================================================================
# Rewrites the METEOSAT cell table files from the nice object. 
#==============================================================================
def rewriteBigcellTab():

    path = "/users/global/cornkle/data/OBS/meteosat/bigcell_area_table/"    
    out = path + 'rewrite/'
    print out
    os.system('rm '+out+'*.txt')
    ok=uarr.locate("*.txt", path)
    
    for a in ok:
        print 'Doing '+a
        tab = parseCellTables(a)
        minute=tab["Date"][0].minute
        hour=tab["Date"][0].hour
        tab.to_csv(out+'cell_40c_'+str(hour).zfill(2)+str(minute).zfill(2)+'_JJAS.txt')
  
#==============================================================================
# Reads the METEOSAT blob files       
#==============================================================================
def readBlobFile(bfile):
    blob_path = "/users/global/cornkle/data/OBS/meteosat/cell_blob_files/"
    bfile = blob_path+str(yr)+'/'+str(mon).zfill(2)+'/'+str(yr)+str(mon).zfill(2)+str(day).zfill(2)+str(hr).zfill(2)+str(mins).zfill(2)
    rr = np.fromfile(bfile,dtype=rrMDI.dtype) 
    rr.shape = rrShape      
        
#==============================================================================
# Compares TRMM to METEOSAT
#==============================================================================
def compare_TRMMmsg():
   msg  = "/users/global/cornkle/data/OBS/meteosat/bigcell_area_table/rewrite/"
   
   trmm = read_TRMMswath()
   rrShape = (580,1640)
   rrMDI = np.uint16()
   
   for hr, mins, yr, mon, day in zip(trmm['date'].hours, trmm['date'].minutes,trmm['date'].years, trmm['date'].months, trmm['date'].days):
        # pick hour/ minute of msg where TRMM swath exists
       df = pd.read_csv(msg+'cell_40c_'+str(hr).zfill(2)+str(mins).zfill(2)+'_JJAS.txt')
       df.set_index('Date', inplace=True, drop=True)
       sel = df.loc[str(yr)+'-'+str(mon).zfill(2)+'-'+str(day).zfill(2)+' '+str(hr).zfill(2)+':'+str(mins).zfill(2)+':'+str(00)]
       big = sel.loc[sel['Area'] >= 25000]
       
       blobs = readBlobFile
       
       get blob hhere or extra??
       
       
       
       
       
    
        