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
import datetime as dt
import os
import pickle as pkl

HOD=range(24)   # hours of day
CFAC=-781648343
LFAC=-781648343
COFF=1856
LOFF=1856
YRANGE=range(2004,2014)
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
def extract_TRMMfile(tpath, hod=HOD, yrange=YRANGE, mtresh=MTRESH):
    
    path = tpath
    pattern = path + '{0:d}/{1:02d}/2A25.{0:d}{1:02d}{2:02d}.*.7_rain_f4.gra' 
    
    fdic = {'fpath' : [], 'tmins' : [], 'date' : ut.date_list()}
    
    for yr,mo,dy in itertools.product(yrange,range(6,10),range(1,31)):
     #   print yr
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
                
                if not t.hour in hod: 
                    continue
                
                # test whether close to 30mins or full hour
                if mtresh:
                    if (t.minute > 3) & (t.minute < 27):
                        continue
                    if (t.minute > 33) & (t.minute < 57):
                        continue
                
           #     print t.minute
                minute = 0   # guessing that t.minute is shortly after full
              
               # print(t.minute)
               # I could include a better minute filter here 
                if t.minute > 15 and t.minute < 45: 
                    minute = 30 
                  
                
           #     print 'MSG', minute
                fdic['fpath'].append(rain_str)
                fdic['date'].add(yr,mo,dy,t.hour,minute,0)
                fdic['tmins'].append(t.minute)
    
    return fdic
    

#========================================================================================
# Reads the METEOSAT cell table files, gives them a header and just reads columns
# of interest for comparison with TRMM. Creates an easy to read object. 

# Table files: every 30 minutes, identified single systems with a random number
# Columns: 5-area (km2), 6/7-lat/lon, 8-min col number (start at 1), 
# 12 - mean T of cell, 13 -Temp treshold to identify cell
#========================================================================================
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
    

#=========================================================================================
# Rewrites the METEOSAT cell table files from parseCellTables() object to nice *.txt files. 
#=========================================================================================
def rewriteBigcellTab():

    path = "/users/global/cornkle/data/OBS/meteosat/bigcell_area_table/"    
    out = path + 'rewrite/'
    print(out)
    os.system('rm '+out+'*.txt')
    ok=uarr.locate("*.txt", path)
    
    for a in ok:
        print('Doing '+a)
        tab = parseCellTables(a)
        minute=tab["Date"][0].minute
        hour=tab["Date"][0].hour
        tab.to_csv(out+'cell_40c_'+str(hour).zfill(2)+str(minute).zfill(2)+'_JJAS.txt')
        


#========================================================================================
# Rewrites 580x1640 msg lat lon to something nice (lat lon from blobs)
#========================================================================================
def rewriteMsgLonLat():
    llFile = '/users/global/cornkle/data/OBS/meteosat/MSG_1640_580_lat_lon.gra'

    llShape = (580,1640)
    llMDI = np.float32(13.5)
    ll = np.fromfile(llFile,dtype=llMDI.dtype)
    lon = ll[0:580*1640]
    lat = ll[580*1640:]
    lat.shape = llShape
    lon.shape = llShape

    llsavefile = '/users/global/cornkle/data/OBS/meteosat/MSG_1640_580_lat_lon'
    np.savez(llsavefile,lon=lon,lat=lat)        

  
#========================================================================================
# Reads the METEOSAT blob files       
#========================================================================================
def readBlobFile(bfile):
    if not os.path.isfile(bfile):
        return np.array([False])
    rrShape = (580,1640) # msg shape
    rrMDI = np.uint16() 
    rr = np.fromfile(bfile,dtype=rrMDI.dtype) 
    rr.shape = rrShape
    return rr
    
#========================================================================================
# Reads the METEOSAT raw files     
#========================================================================================
def readMSGraw(bfile):
    if not os.path.isfile(bfile):
        return np.array([False])
        
    rrShape = (580,1640) # msg shape
    rrMDI = np.uint8(255)
    rr = np.fromfile(bfile,dtype=rrMDI.dtype) 
    rr.shape = rrShape
    rr=rr.astype(np.int32) - 173
    msg_latlon=np.load('/users/global/cornkle/data/OBS/meteosat/MSG_1640_580_lat_lon.npz')
    mlon = msg_latlon['lon']
    mlat = msg_latlon['lat']
    
    msg_obj = {'t' : rr, 'lons' : mlon, 'lats' : mlat} # lats lons numpy arrays!      
    
    return msg_obj   
    
    
#=======================================================================================
# Reads the data of TRMM 2A25 swaths binary files with lon, lat, rain, flag in it
# Every swath is 49 pixel wide. 2 bytes integer
# Lon = lon/100, Lat = lat/100, pcp = rain/10 (mm/h), flag
#     
#=======================================================================================
def readTRMMswath(tpath):
    
    if not os.path.isfile(tpath):
        return np.array([False])
    rr = np.fromfile(tpath,dtype=np.int16) 
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
    lont=lons/100.
    latt=lats/100.
    rain=rainrs/10.
  #  lonmin, lonmax=np.amin(lon),np.amax(lon)
   # latmin, latmax=np.amin(lat),np.amax(lat)
 #   lonx=lon[0,:]
    laty=latt[:,0]
    
    cut=np.where((laty<=22) & (laty>=0))   # cut TRMM down to 0 - 22 degN for smaller array
    cut = np.array(cut)
    if cut.size:
        rain=rain[cut, :]
        lont=lont[cut, :]
        latt=latt[cut, :]
        flags=flags[cut, :]
        yy=int(rain.size/49)
    
        rain=rain.reshape(yy,49)
        lont=lont.reshape(yy,49)
        latt=latt.reshape(yy,49)
        flags=flags.reshape(yy,49)
    
    lattmin, lattmax = latt.min(), latt.max()
    lontmin, lontmax = lont.min(), lont.max()
    
    laty=latt[:,0]
    
    trmm_obj = {'pcp' : rain, 'lons' : lont, 'lats' : latt, 'flags' : flags, 'laty' : laty, 'latmin' : lattmin, 'latmax' : lattmax, 'lonmin' : lontmin, 'lonmax' : lontmax} # lats lons numpy arrays!
    
    return trmm_obj
    
    
#==============================================================================
# To MSG indices
#==============================================================================    
def ll_toMSG(lons,lats,cfac=CFAC,lfac=LFAC,coff=COFF,loff=LOFF): 
    
  
    
    lats = np.array(lats)
    lons = np.array(lons)
  
    if not lats.shape == lons.shape:
        print('Lats lons must have same dimensions')
        return
    if not lats.size == lons.size:
        print('Lats lons must have same size')
        return 
            
    if (np.min(lats) < -90.) | (np.max(lats) > 90.):
        print('Lats are out of range')
        return
    if (np.min(lons) < -180.) | (np.max(lons) > 180.):          
        print('Lons are out of range')
        return
        
    pi=3.14159265359 # Define as double precision.
    sat_height=42164.0 # Height of satellite.
    r_eq=6378.169    # Radius of Earth at equator.
    r_pol=6356.5838  # Radius of Earth at pole.
    sub_lon=0.0      # Longitude of sub-satellite point.
    
    # Convert lats and longs to radians.
    lats_r=lats*pi/180.
    lons_r=lons*pi/180.

    # Calculate geocentric latitude from the geographic one.
    c_lat=np.arctan(0.993243*np.sin(lats_r)/np.cos(lats_r))

    # Use c_lat to calculate the length from the Earth centre
    # to the surface of the Earth ellipsoid equations.
    re=r_pol/np.sqrt(1.-0.00675701*np.cos(c_lat)*np.cos(c_lat))

    # Calculate the forward projection.
    r1=sat_height-re*np.cos(c_lat)*np.cos(lons_r-sub_lon)
    r2=-re*np.cos(c_lat)*np.sin(lons_r-sub_lon)
    r3=re*np.sin(c_lat)
    rn=np.sqrt(r1*r1+r2*r2+r3*r3)

    # Create output arrays.
    cols=np.empty_like(lats)
    rows=np.empty_like(lats)


    # Check for visibility, whether the point is visible from the satellite.
    dotprod=np.array([r1*re*np.cos(c_lat)*np.cos(lons_r-sub_lon)-r2*r2-r3*r3*(r_eq/r_pol)**2.])

#    t=np.where(dotprod > 0., True, False)
#
#    if ~t.all():
#        cols[~t]=-999
#        rows[~t]=-999
#
#    if t.any():
#        cols[t]=math.atan(-r2[t]/r1[t])
#        rows[t]=math.asin(-r3[t]/rn[t])
#
#    cols[t]=np.round(coff+cols[t]*cfac/2.**16) # This seems to be incorrect in the example program.
#    rows[t]=np.round(loff+rows[t]*lfac/2.**16)

    cols=np.arctan(-r2/r1)
    rows=np.arcsin(-r3/rn)
    
    cols=np.round(coff+cols*cfac/2.**16) # This seems to be incorrect in the example program.
    rows=np.round(loff+rows*lfac/2.**16)
   
    dic = {'x' : cols, 'y' : rows}

    return dic      

        
   
   
def getTRMMconv(flags):
  
  bla=flags.astype(int)
  npfalse=[]

  for b, i in zip(np.nditer(bla), range(bla.size)): 
       bb='{0:016b}'.format(int(b))    
       npfalse.append(int(bb[-6]))

  return npfalse
  

#==============================================================================
# Compares TRMM to METEOSAT  | MSG indices
#==============================================================================
def compare_TRMMmsg_indices(hod=HOD, yrange=YRANGE):
   msg  = "/users/global/cornkle/data/OBS/meteosat/bigcell_area_table/rewrite/"
   blob_path = "/users/global/cornkle/data/OBS/meteosat/cell_blob_files/"
   msg_latlon=np.load('/users/global/cornkle/data/OBS/meteosat/MSG_1640_580_lat_lon.npz')
   ttpath = "/users/global/cornkle/data/OBS/TRMM/trmm_swaths_WA/"
   msg_rawp="/users/global/cornkle/data/OBS/meteosat/msg_raw_binary/"
   mlon = msg_latlon['lon']
   mlat = msg_latlon['lat']
   print('Start extract Trmm')
   trmm = extract_TRMMfile(ttpath, yrange=yrange, hod=hod)
   
   tpath = trmm['fpath']
   
   mdic = {'perc' : [], 'p' : [], 'pp' : [], 'pi' : [], 'pmax' : [],  't' : [], 'hod' : [],
           'lat' : [], 'lon' : [], 'yr' : [], 'mon' : [], 'tpixel' : []}
   mdic_f = {'perc' : [],  'pconv': [],'ppconv': [], 'piconv': [],  'pmaxconv' : [], 'tconv' : [], 
             'tnfconv' : [], 'hod' : [], 'lat' : [],'lon' : [], 'yr' : [], 'mon' : [], 'tpixel' : []}
             
   mll = ll_toMSG(mlon,mlat) 
   
   cnt=0      
   for hr, mins, yr, mon, day, tfile in zip(trmm['date'].hours, trmm['date'].minutes,trmm['date'].years, trmm['date'].months, trmm['date'].days, tpath):       
     
       print('Doing '+str(yr)+'-'+str(mon).zfill(2)+'-'+str(day).zfill(2)+' '+str(hr).zfill(2)+':'+str(mins).zfill(2))
    #   print 'TRMM:'+tfile 
       # pick hour/ minute of msg where TRMM swath exists, check blobfile!
       d=dict()
       df = pd.read_csv(msg+'cell_40c_'+str(hr).zfill(2)+str(mins).zfill(2)+'_JJAS.txt')
       dstring=str(yr)+'-'+str(mon).zfill(2)+'-'+str(day).zfill(2)+' '+str(hr).zfill(2)+':'+str(mins).zfill(2)+':'+str(00).zfill(2)
       
       if not dstring in df['Date'].as_matrix():
           continue
       
       sel=df.loc[df['Date'] == dstring]       
       big = sel.loc[sel['Area'] >= 25000] # only mcs over 25.000km2
       if big.empty:
           continue
       
       trs = readTRMMswath(tfile)
       if not trs:
           print('TRMM file not found')
           return
           
      # print ' went trough trmm'    
       
       
       fflags=trs['flags'][np.where(trs['lats']>4)]
       fpcp=trs['pcp'][np.where(trs['lats']>4)]
       fmask=getTRMMconv(fflags)  # filter for convective rain
       fmask=np.array(fmask)
       if fpcp.flat[np.where(fmask)].any(): perc= np.percentile(fpcp.flat[np.where(fmask)], 98) 
       mdic['perc'].append(perc)       

       mfile=msg_rawp+str(yr)+'/'+str(mon).zfill(2)+'/'+str(yr)+str(mon).zfill(2)+str(day).zfill(2)+str(hr).zfill(2)+str(mins).zfill(2)+'.gra'
       
       msg_raw=readMSGraw(mfile)          
       
       if not msg_raw:
           continue  
#       print 'MSG Went trough'
       msg_t=msg_raw['t']
       
       d['lat']=big['Lat'].values.tolist()
       d['lon']=big['Lon'].values.tolist()
       d['area']=big['Area'].values.tolist()
       d['temp']=big['Temp'].values.tolist()
       d['mincol']=big['Mincol'].values.tolist()
       
       bfile = blob_path+str(yr)+'/'+str(mon).zfill(2)+'/'+str(yr)+str(mon).zfill(2)+str(day).zfill(2)+str(hr).zfill(2)+str(mins).zfill(2)+'.gra'
       blobs = readBlobFile(bfile)               
       if not blobs.any():
           continue
       
       # loop for each blob center, to find the whole area
       for lon, lat, mt, mc in zip(d['lon'], d['lat'], d['temp'] , d['mincol']):
          
          pp = ll_toMSG(lon,lat)
          point=np.where((mll['x']==pp['x']) & (mll['y']==pp['y']))
          
          # blob number
          nb = blobs[point]
          # if we find a 0 instead of a blob, continue
          if not nb[0]:       
              continue          
     
          # find where else is blob number
          isblob = np.where(blobs == nb)          
          #  print min(isblob[1])
          #  print mc-1
          # check again whether found blob is approx 25.000km2
          if isblob[0].size < 2500:
              continue
          
          # lat lons of complete blob
          blats=mlat[isblob]
          blons=mlon[isblob]
          
          # msg indices of complete blob
          my=mll['y'][isblob]
          mx=mll['x'][isblob]          
          mpair = (mx+my)*(mx+my+1)/2+my
          
          blatmin, blatmax = blats.min(), blats.max()
          blonmin, blonmax = blons.min(), blons.max()          
                            
          # whole blob must be inside TRMM. Attention: This draws a rectangle. 
          # There is still a chance that blob is not in TRMM. Checked later!               
          if not (trs['lonmin'] < blonmin) & (trs['lonmax'] > blonmax):                
                 continue
          if not (trs['latmin'] < blatmin) & (trs['latmax'] > blatmax):                 
                 continue
                             
          ll_trmm = ll_toMSG(trs['lons'],trs['lats'])
          
          tx = ll_trmm['x']
          ty = ll_trmm['y']
          
          tpair = (tx+ty)*(tx+ty+1)/2+ty               
          inter=np.in1d(tpair, mpair)                   #returns false and true, whole grid
          inter_rev=np.in1d(mpair, tpair.flat[inter])   # Attention: this leaves out meteosat cells were no closest TRMM cell (since TRMM is coarser!)
          
          # No intersection between Trmm and MSG? Continue!
          # if not inter.any(): 
          #     continue  
          
          # have at least 500 pixels shared for MCS between TRMM and MSG
          if sum(inter) < 500:  
              continue
          
          bprcp=trs['pcp'].flat[inter]            
          flags=trs['flags'].flat[inter] 
          mtt=msg_t[isblob].flat[inter_rev]       
          
          #we need same number of TRMM and MSG per plot to do the masking
          if not bprcp.size==mtt.size: 
              print('Tprcp and MSGT not same, someting wrong!')
              continue
          
          # if TRMM zero for whole blob, continue
          #if not bprcp.any():
          #    continue
          
           # check for at least 500 TRMM pixels in MSG above 0 rain
          if np.count_nonzero(bprcp) < 20:
              continue
          
          mask=getTRMMconv(flags)  # filter for convective rain
          mask=np.array(mask)
                    
          
          # remove all these zero rainfall from blob
          bprcp=bprcp.flat[np.where(bprcp)]
          mtt=mtt.flat[np.where(bprcp)]   
          flags=flags.flat[np.where(bprcp)] 
          mask=getTRMMconv(flags)  #list of 0 and 1, flattened!
          
          if sum(mask) < 2:
              continue
          
          pm=np.nanmean(bprcp)  
          tm=np.nanmean(mtt)
          ppm=np.percentile(bprcp, 98)  
          pmaxm=bprcp.max()  
          pi=float(np.greater(bprcp, 30.).sum())/float(bprcp.size)
          
          mdic['p'].append(pm)
          mdic['pp'].append(ppm)
          mdic['pmax'].append(pmaxm)
          mdic['pi'].append(pi)
          mdic['t'].append(tm)
          mdic['hod'].append(hr)
          mdic['yr'].append(yr)
          mdic['mon'].append(mon)
          mdic['lat'].append(lat)
          mdic['lon'].append(lon)
          mdic['tpixel'].append(bprcp.size) 
                   
          # if less than 100 pixel convective in blob, continue
          if sum(mask) < 50:
              continue
          
                   
          # check for at least 500 TRMM pixels in MSG above 0 rain
          #if np.count_nonzero(bprcp) < 500:
          #    continue
                    
                                                                   
          pc=np.nanmean(bprcp.flat[np.where(mask)])
          tc=np.nanmean(mtt.flat[np.where(mask)])
          pic=float(np.greater(bprcp.flat[np.where(mask)], 30.).sum())/float(sum(mask))
          ppc=np.percentile(bprcp.flat[np.where(mask)], 98)
          pmaxc=bprcp.flat[np.where(mask)].max()
          
          
          
          # 30 is the treshold for intensive rainfall
                         
          

        #  print 'Nb', nb
          mdic_f['pconv'].append(pc)
          mdic_f['piconv'].append(pic)
          mdic_f['ppconv'].append(ppc)
          mdic_f['pmaxconv'].append(pmaxc)
          mdic_f['tconv'].append(tc)  
          mdic_f['tnfconv'].append(tm)                  
          mdic_f['hod'].append(hr)
          mdic_f['yr'].append(yr)
          mdic_f['mon'].append(mon)
          mdic_f['lat'].append(lat)
          mdic_f['lon'].append(lon)
          mdic_f['tpixel'].append(sum(mask)) 
          cnt=cnt+1
          
          myDicts=[mdic, mdic_f]
          
   pkl.dump(myDicts, open('/users/global/cornkle/data/OBS/MSG_TRMM_temp_pcp_300px'+str(yrange[0])+'-'+str(yrange[-1])+'.p', 'wb'))
   print('Saved '+'MSG_TRMM_temp_pcp_'+str(yrange[0])+'-'+str(yrange[-1])+'.p with '+str(cnt)+' MCSs')
          
  # savefile1 = '/users/global/cornkle/data/OBS/MSG_TRMM_temp_pcp_'+str(yrange[0])+'-'+str(yrange[-1])+'.npy'
  # savefile2 = '/users/global/cornkle/data/OBS/MSG_TRMM_temp_pcp_conv_'+str(yrange[0])+'-'+str(yrange[-1])+'.npy'
  # np.save(savefile1, mdic)  
  # np.save(savefile2, mdic_f)
  # print 'Saved '+'MSG_TRMM_temp_pcp_'+str(yrange[0])+'-'+str(yrange[-1])+'.npy'
   
       
NB=0
def quickreadTrmmMSG(tfile, nb=NB):
    
          
    msg_rawp="/users/global/cornkle/data/OBS/meteosat/msg_raw_binary/"

    yr=tfile['date'].years[nb]
    mon=tfile['date'].months[nb]
    hr=tfile['date'].hours[nb]
    day=tfile['date'].days[nb]
    mins=tfile['date'].minutes[nb]
    files = str(yr)+'/'+str(mon).zfill(2)+'/'+str(yr)+str(mon).zfill(2)+str(day).zfill(2)+str(hr).zfill(2)+str(mins).zfill(2)+'.gra' 
    bfile=msg_rawp+files
    msg=readMSGraw(bfile)
    trmm = readTRMMswath(tfile['fpath'][nb])
    
    return (msg, trmm)
          
    
        