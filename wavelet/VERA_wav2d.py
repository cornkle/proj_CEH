#!/bin/env python

import numpy as np
import os.path
from scipy import ndimage
import twod as w2d
import itertools

def readCloud(date):
    filePattern = "/users/global/cmt/msg/tropWA/archive_tropWA/ch9/{0:d}/{1:02d}/{0:d}{1:02d}{2:02d}{3:02d}{4:02d}.gra"
    if (date[0] > 2012):
        filePattern = "/users/global/danbel/msg/tropWA/archive_tropWA/ch9/{0:d}/{1:02d}/{0:d}{1:02d}{2:02d}{3:02d}{4:02d}.gra"
    rrFile = filePattern.format(date[0],date[1],date[2],date[3],date[4])
    FileExists = os.path.isfile(rrFile)
    
    if FileExists:
        rrShape = (222,1384)
        rrMDI = np.uint8(255)
        rr = np.fromfile(rrFile,dtype=rrMDI.dtype)
        rr.shape = rrShape
        rr = rr.astype(np.int32) - 173
        return FileExists, rr
    else:
        return FileExists, np.zeros(1)


dt = 3. #approximate grid spacing in km
#Tai park area
#indxb = 25
#indxe= 106
#indyb = 170
#indye = 301
#NW of Tai park area
#indxb = 85
#indxe= 215
#indyb = 70
#indye = 200
#Ghana
#indxb = 25
#indxe= 165
#indyb = 335
#indye = 500
#Lake Volta
indxb = 65
indxe= 195
indyb = 452
indye = 538

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

yrange = range(2012,2013)
#yrange = range(2010,2016)
mn = 0
for yr,mo,dy,hr in itertools.product(yrange,range(1,13),range(1,32),range(10,19)):
    d = np.array([yr,mo,dy,hr,mn])
    print yr,mo,dy,hr
    filexist, tir = readCloud(d)
    if filexist:
        tir = tir.astype(np.float)
        tir = tir[indxb:indxe,indyb:indye]
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
        months = np.concatenate((months,np.repeat(mo,len(xpks))),axis=0)
        days = np.concatenate((days,np.repeat(dy,len(xpks))),axis=0)
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
                            
filePattern = "/users/global/cornkle/data/pythonWorkspace/wavelet/saves/VERAwav2DpeaksT20Volta{0:d}_{1:d}"
rrFile = filePattern.format(yrange[0],yrange[-1])
np.savez(rrFile,tirmin=tirmin,tirmean=tirmean,tirloc=tirloc,Scales=Scales,xLocs=xLocs,yLocs=yLocs,Powers=Powers,years=years,months=months,days=days,hours=hours)
