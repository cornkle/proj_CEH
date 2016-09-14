# -*- coding: utf-8 -*-
"""
Created on Thu May 19 15:11:49 2016

@author: cornkle
"""

import numpy as np
import os
from utils import u_arrays as uarr
import pandas as pd
from utils import u_time as ut
from utils import u_lists as ul
#import datetime as dt
#import glob
import itertools

HOD=range(24)
YRANGE=range(2004,2013)
MTRESH=0

#=======================================================================================
# Reads TRMM 2A25 files to get a list with valid file names (files of interest 
# where time is close to full or half hour to be potentially correspondant to METEOSAT).

# Files with _rain_f4 in name mark files where rainfall is available over West Africa. 
# However, we use the files with less floating points (without the _rain_f4 )
#
# Output: Filelist of TRMM files over West Africa at full/half hour; 
# corrsponding time rounded to half or full hour
#=======================================================================================
        
#=======================================================================================
# Reads the data of TRMM 2A25 swaths binary files with lon, lat, rain, flag in it
# Every swath is 49 pixel wide. 2 bytes integer
# Lon = lon/100, Lat = lat/100, pcp = rain/10 (mm/h), flag
#     
#=======================================================================================  
    
class trmm(object):
     def __init__(self, trmm_folder, hod=HOD, yrange=YRANGE, area=None):
        
       if not os.path.isdir(trmm_folder):
            print('Not a directory')    
            quit()
       
       fdic = {'fpath' : [], 'tmins' : [], 'date' : ut.date_list()}
       rfiles = []                         
       
       for yr,mo in itertools.product(list(yrange),list(range(6,10))):  # rain_f4 files only available for 6 to 10
           
           tpath= os.path.join(trmm_folder, str(yr), str(mo).zfill(2))
           files=uarr.locate('_rain_f4.gra', tpath)         
           rfiles.extend(files)
       rfiles.sort(key=ul.natural_keys)  
       
       if not rfiles:
           print('Not trmm files found')
          
     #  self.fpath=fdic['fpath']     
     #  return
       for eachfile in rfiles:           
            rain_str = eachfile.replace('_rain_f4', '')
            time_str = eachfile.replace('_rain_f4', '_time')   
            rr = np.fromfile(time_str,dtype=np.float32) # seconds of day
             
            secmean = rr.mean()
            t = ut.sec_to_time(secmean)
                
            if not t.hour in hod: 
                continue    
            
            rr = np.fromfile(rain_str,dtype=np.int16) 
            x = 49   # trmm swath is always 49 wide
            nb = rr.size
            single = int(nb/4) # variables lon lat rainrate flag

            lons = rr[0:single]
            lats = rr[single:2*single]             
            y = int(lons.size/x)
            lons = np.resize(lons, (y,x))
            lats = np.resize(lats, (y,x))       
            lont=lons/100.
            latt=lats/100.
            if area:          
                 box=np.where((lont > area[0]) & (lont < area[1]) & (latt > area[2]) & (latt < area[3]))
                 if not box[0].any():
                    continue
         #       print(len(box[0]))
                 if not len(box[0]) > 5000:   # minimum pixel overlap with TRMM and box
                    continue
                 
            
            fdic['fpath'].append(rain_str)        
            fdic['date'].add(int(rain_str[-20:-16]),int(rain_str[-16:-14]),int(rain_str[-14:-12]),t.hour,t.minute,0)                            
            
       self.fpaths =fdic['fpath']  
       self.dates = fdic['date']           
      
                                    
        
     def getDData(self, yr, mo, dy, hr, mi, cut=None):                 
        
         ind=self.dates.getInd(yr,mo,dy,hr,mi)
         #print('Ind:', ind)
         if not ind:
            print('No data for date')
            return False
         
 
         tfile = self.fpaths[ind[0]]  
         dic = self.getData(tfile, cut=cut)    
    
         return dic 
         
         
     def getData(self, path, cut=None):
           
         tfile = path
         if not os.path.isfile(tfile):
             print('Ind found but file does not exist. Error')
             quit()
   
         rr = np.fromfile(tfile,dtype=np.int16) 
         x = 49   # trmm swath is always 49 wide
         nb = rr.size
         single = int(nb/4) # variables lon lat rainrate flag

         lons = rr[0:single]
         lats = rr[single:2*single]
         rainrs = rr[2*single:3*single]
         flags = rr[3*single:4*single]
    
         y = int(lons.size/x)
         lons = np.resize(lons, (y,x))
         lats = np.resize(lats, (y,x))
         rainrs = np.resize(rainrs, (y,x))
         flags = np.resize(flags, (y,x))
         dlont=lons/100.
         dlatt=lats/100.
         rain=rainrs/10.

         laty=dlatt[:,0]
         
         lower=3
         upper=22
         
         if cut:
             lower = cut[0] 
             upper = cut[1] 
    
         cutt=np.where((laty<=upper) & (laty>=lower))   # cut TRMM down to 3 - 22 degN for smaller array
         cutt = np.array(cutt)
         if cutt.size:
             rain=rain[cutt, :]
             lont=dlont[cutt, :]
             latt=dlatt[cutt, :]
             flags=flags[cutt, :]
             yy=int(rain.size/49)
    
             rain=rain.reshape(yy,49)
             lont=lont.reshape(yy,49)
             latt=latt.reshape(yy,49)
             flags=flags.reshape(yy,49)   
         else:                           
             print('Cut problem')
             return 
         
         self.lon = dlont
         self.lat = dlatt
         self.nx = 49
         self.ny = y  
         
         dic = {'p' : rain, 'lon' : lont, 'lat' : latt, 'flags' : flags}
           
    
         return dic
                
    
    
class msg(object):
    
    def __init__(self, msg_folder):
        
       if not os.path.isdir(msg_folder):
            print('Not a directory')    
            quit()
       lpath=uarr.locate('lon.npz', msg_folder)
       
       mpath=os.path.join(msg_folder, 'msg_raw_binary')        
       if not os.path.isdir(mpath):
           print('No msg_raw_binary')    
           quit()
       
       msg_latlon=np.load(lpath[0])
       mlon = msg_latlon['lon']
       mlat = msg_latlon['lat']      

       self.lat = mlat
       self.lon = mlon
       self.nx = mlon.shape[1]
       self.ny = mlon.shape[0]       
                               
       self.years = os.listdir(mpath)
       self.root =  msg_folder
       

    def setDate(self, yr,mon,day,hr,mins):
        self.dpath = os.path.join(self.root, 'msg_raw_binary', str(yr) , str(mon).zfill(2), str(yr)+str(mon).zfill(2)+str(day).zfill(2)+str(hr).zfill(2)+str(mins).zfill(2)+'.gra' )      
        if not os.path.isfile(self.dpath):
            self.dpath = False
        
        if  os.path.isdir(os.path.join(self.root, 'cell_blob_files')):       
            self.bpath = os.path.join(self.root, 'cell_blob_files', str(yr) , str(mon).zfill(2), str(yr)+str(mon).zfill(2)+str(day).zfill(2)+str(hr).zfill(2)+str(mins).zfill(2)+'.gra' )                  
        else:
            print('No blob file dir found!')
            self.bpath = False
        if  os.path.isdir(os.path.join(self.root, 'bigcell_area_table')):     
            self.tpath = os.path.join(self.root, 'bigcell_area_table', 'rewrite',  'cell_40c_'+str(hr).zfill(2)+str(mins).zfill(2)+'_JJAS.txt')  
        else:
            print('No table file dir found!')  
            self.tpath = False              
        self.date = {'year' : yr, 'month': mon, 'day' : day, 'hour' : hr, 'minute' : mins}      # could be turned into python date
        
        
    def getData(self, y=None,m=None,d=None,h=None,mi=None, llbox=None):
        
        if y:
            self.dpath=os.path.join(self.root, 'msg_raw_binary', str(y) , str(m).zfill(2), str(y)+str(m).zfill(2)+str(d).zfill(2)+str(h).zfill(2)+str(mi).zfill(2)+'.gra' ) 
        
        if not os.path.isfile(self.dpath):
            print('No data for date')
            return False
                
                
        rrShape = (self.ny, self.nx) # msg shape
        rrMDI = np.uint8(255)
        rr = np.fromfile(self.dpath,dtype=rrMDI.dtype) 
        rr.shape = rrShape
        rr=rr.astype(np.int32) - 173
        
        if llbox:
           i,j=np.where((self.lon > llbox[0]) & (self.lon < llbox[1]) & (self.lat > llbox[2]) & (self.lat < llbox[3]))           
           blat=self.lat[i.min():i.max()+1, j.min():j.max()+1 ]
           blon=self.lon[i.min():i.max()+1, j.min():j.max()+1 ]    
           rr=rr[i.min():i.max()+1, j.min():j.max()+1 ] 
        else: 
           blat=self.lat
           blon=self.lon
           rr=rr
           
        dic = {'t' : rr, 'lon' : blon, 'lat' : blat}   
        
        return dic
        
        
    def getBlob(self, llbox=None):
        
       if not self.bpath:
           print('No blob file dir found!')
           return False
           
       rrShape = (self.ny, self.nx) # msg shape
       rrMDI = np.uint16() 
       rr = np.fromfile(self.bpath,dtype=rrMDI.dtype) 
       rr.shape = rrShape
       if llbox:
           i,j=np.where((self.lon > llbox[0]) & (self.lon < llbox[1]) & (self.lat > llbox[2]) & (self.lat < llbox[3]))           
           blat=self.lat[i.min():i.max()+1, j.min():j.max()+1 ]
           blon=self.lon[i.min():i.max()+1, j.min():j.max()+1 ]    
           rr=rr[i.min():i.max()+1, j.min():j.max()+1 ]           
           
       dic = {'t' : rr, 'lon' : blon, 'lat' : blat}   
        
       return dic
          
        
    def getTable(self):
        
         if not self.tpath:
           print('No table file dir found!')
           return False
         
         tab = pd.read_csv(self.tpath)
           
         return tab      
    
