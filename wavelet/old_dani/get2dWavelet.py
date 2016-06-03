# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 11:20:48 2016

@author: cornkle
"""

import numpy as np
from scipy import ndimage
import twod as w2d
from eod import msg
import itertools







YRANGE=range(2012,2013)
MRANGE=range(6,10)
DRANGE=range(32)
HRANGE=range(24)
MIRANGE=range(31)   # minute either 0 or 30

def get2dWavelet(yrange=YRANGE, drange=DRANGE, hrange=HRANGE, mrange=MRANGE, mirange=MIRANGE):
    
    msg_rawp="/users/global/cornkle/data/OBS/meteosat/msg_raw_binary/"
    savepath = "/users/global/cornkle/data/pythonWorkspace/proj_CEH/wavelet/saves/"
    msg_latlon=np.load('/users/global/cornkle/data/OBS/meteosat/MSG_1640_580_lat_lon.npz')
    mlon = msg_latlon['lon']
    mlat = msg_latlon['lat']
    lon_xx=2
    lon_x=-2
    lat_xx=14
    lat_x=12
    box=np.where((mlon < lon_xx) & (mlon > lon_x) & (mlat < lat_xx) & (mlat > lat_x))  
    milon,malon = min(box[1]), max(box[1])
    milat,malat = min(box[0]), max(box[0])            
                
    dt = 3.
    
   # dic= {'powers' : [], 'scales' : [], 'xlocs' : [], 'ylocs' : [], 'years' : [], 'months' : [], 
   #   'days' : [], 'hours' : [], 'tirmin' : [], 'tirmean' : [], 'tirloc' : [], }
    
    Powers=[]
    Scales=[]
    xLocs=[]
    yLocs=[]
    years=[]
    months=[]
    days=[]
    hours=[]
    tirmin=[]
    tirmean=[]
    tirloc=[]

    for hr, mins, yr, mon, day in itertools.product(hrange, mirange, yrange, mrange ,drange):    
 

        mfile=msg_rawp+str(yr)+'/'+str(mon).zfill(2)+'/'+str(yr)+str(mon).zfill(2)+str(day).zfill(2)+str(hr).zfill(2)+str(mins).zfill(2)+'.gra'
        print(mfile)
        tir=msg.readMSGraw(mfile)       
        if not tir.any():
           continue  
       
        print("Doing "+str(yr)+'-'+str(mon).zfill(2)+'-'+str(day).zfill(2)+' '+str(hr).zfill(2)+':'+str(mins).zfill(2))

        tir = tir.astype(np.float)
        tir = tir[milat:malat, milon:malon]
        tiror = np.copy(tir)
        tir[tir>20] = 20
        tir = tir - np.mean(tir)
        
        mother2d = w2d.Mexican_hat()
        
        wavel2d, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=1./12, s0=6./mother2d.flambda(), J=106)
        
        wavel2d[np.real(wavel2d>=0)] = 0.01
        power2d = (np.abs(wavel2d)) ** 2 # Normalized wavelet power spectrum
        period2d = 1. / freqs2d
        scales2d.shape = (len(scales2d),1,1)
        power2ds = power2d / (scales2d**2)
        
        #Find maxima of the 2D wavelet spectrum    
        igra = power2ds
        maxigra = (igra == ndimage.maximum_filter(igra,size=(5,5,5),
                   mode='constant',cval=np.amax(igra)+1))
        maxigra = maxigra.astype(int)
        #igrathresh = np.repeat(50*(scales2d/5.)**.5,igra.shape[1],axis=1) #NOTICE the scale-dependant threshold for power
        igrathresh = np.repeat(1*(scales2d/1.)**.5,igra.shape[1],axis=1) #NOTICE the scale-dependant threshold for power
        igrathresh = np.repeat(igrathresh,igra.shape[2],axis=2)
        igrabkg = (igra < igrathresh).astype(int) #Apply the threshold
        igrapeaks = maxigra - igrabkg  #Remove all smaller than the threshold
        zpks,ypks,xpks = np.argwhere(igrapeaks == 1).T
        
        Powers = np.concatenate((Powers,igra[zpks,ypks,xpks]),axis=0)
        Scales = np.concatenate((Scales,period2d[zpks]/2.),axis=0)
        xLocs = np.concatenate((xLocs,xpks),axis=0)
        yLocs = np.concatenate((yLocs,ypks),axis=0)
        years = np.concatenate((years,np.repeat(yr,len(xpks))),axis=0)
        months = np.concatenate((months,np.repeat(mon,len(xpks))),axis=0)
        days = np.concatenate((days,np.repeat(day,len(xpks))),axis=0)
        hours = np.concatenate((hours,np.repeat(hr,len(xpks))),axis=0)
        tirloc = np.concatenate((tirloc,tiror[ypks,xpks]),axis=0)
        iscale = (np.ceil(period2d[zpks]/4./dt)).astype(int)
        
        mintir = []
        meantir = []
       
        for imtir in range(len(xpks)):
            #Find all indices within the local circle of radius iscale...
            # ... Then average over those indices
            xloc1 = np.arange(xpks[imtir]-iscale[imtir],xpks[imtir]+iscale[imtir]+1)
            yloc1 = np.arange(ypks[imtir]-iscale[imtir],ypks[imtir]+iscale[imtir]+1)
            xloc,yloc = np.meshgrid(xloc1,yloc1)
            distloc = ( (xloc-xpks[imtir])**2 + (yloc-ypks[imtir])**2 ) ** .5
            indloc = (distloc <= iscale[imtir]).nonzero()
            yindtir = indloc[0] - iscale[imtir] + ypks[imtir]
            xindtir = indloc[1] - iscale[imtir] + xpks[imtir]
            #remove indices outside the whole study domain
            idel = np.nonzero( (yindtir<0) | (yindtir>tiror.shape[0]-1) | 
            (xindtir<0) | (xindtir>tiror.shape[1]-1) )
            xindtir = np.delete(xindtir,idel)
            yindtir = np.delete(yindtir,idel)
            #calculate mean and min temperature
            mintir = np.append(mintir,np.amin(tiror[yindtir,xindtir]))
            meantir = np.append(meantir,np.mean(tiror[yindtir,xindtir]))
            
        tirmin = np.concatenate((tirmin,mintir),axis=0)
        tirmean = np.concatenate((tirmean,meantir),axis=0)  
        
                
                            
    filepath = savepath+'wavelet_'+str(yrange[0])+'-'+str(yrange[-1])+'.npz'
    np.savez(filepath,tirmin=tirmin,tirmean=tirmean,tirloc=tirloc,Scales=Scales,xLocs=xLocs,yLocs=yLocs,Powers=Powers,years=years,months=months,days=days,hours=hours)
    print 'Saved '+'wavelet_'+str(yrange[0])+'-'+str(yrange[-1])+'.npz'