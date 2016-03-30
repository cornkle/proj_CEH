#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#from read_deforest_data import read_deforest_data

dt = 3.

#***Change this as needed (size thresholds)***
sizthresh = 6 #Min size of the main structures
upperthresh = 2000 #Max size of the main structure
powthresh0 = 0 #Initial value for the scale-dependant power threshold
#********************************************

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

# Read lat/lon
ll = np.load('/users/global/danbel/msg/tropWA/latlon.npz')
lon = ll['lon']
lat = ll['lat']
lat = lat[indxb:indxe,indyb:indye]
lon = lon[indxb:indxe,indyb:indye]

# Load deforestation data
#cover, loss = read_deforest_data()
#cover = cover[indxb:indxe,indyb:indye]

# Create map projection
#Tai park area
#LON_RANGE=[-9.0,-5.5]
#LAT_RANGE=[4.6,6.8]
#NW of Tai park area
#LON_RANGE=[-11.8,-8.3]
#LAT_RANGE=[6.3,9.8]
#General area
LON_RANGE=[np.amin(lon),np.amax(lon)]
LAT_RANGE=[np.amin(lat),np.amax(lat)]
M = Basemap(projection='geos',resolution='i', lon_0=0, \
            llcrnrlat =LAT_RANGE[0],urcrnrlat=LAT_RANGE[1], \
            llcrnrlon =LON_RANGE[0],urcrnrlon=LON_RANGE[1] )
xx,yy = M(lon,lat)


#Read data:
sizes1 = []
powers1 = []
xlocs1 = []
ylocs1 = []
years1 = []
months1 = []
days1 = []
hours1 = []
tirmin1 = []
tirmean1 = []
#yrange = ['2004_2005','2006_2007','2008_2009','2010_2011','2012_2013','2014_2015']
yrange = ['2012_2012']
for iyr in range(len(yrange)):
    wavfile = '/users/global/cornkle/data/pythonWorkspace/wavelet/saves/VERAwav2DpeaksT20Volta2012_2012.npz'
    wav = np.load(wavfile)
    sizes1.append(wav['Scales'])
    powers1.append(wav['Powers'])
    xlocs1.append(wav['xLocs'])
    ylocs1.append(wav['yLocs'])
    years1.append(wav['years'])
    months1.append(wav['months'])
    days1.append(wav['days'])
    hours1.append(wav['hours'])
    tirmin1.append(wav['tirmin'])
    tirmean1.append(wav['tirmean'])
sizes1 = np.hstack(sizes1)
powers1 = np.hstack(powers1)
xlocs1 = np.hstack(xlocs1)
ylocs1 = np.hstack(ylocs1)
years1 = np.hstack(years1)
months1 = np.hstack(months1)
days1 = np.hstack(days1)
hours1 = np.hstack(hours1)
tirmin1 = np.hstack(tirmin1)
tirmean1 = np.hstack(tirmean1)

figsy = min((indxe-indxb)/8, 14) 
figsx = figsy + abs((indye-indyb) - (indxe-indxb)) / 5
#Calculate:
for mo in range(1,13):
    print mo
    #Create figures
    F = plt.figure(figsize=(figsx,figsy))
    Fcount = plt.figure(figsize=(figsx,figsy))
    Fpow = plt.figure(figsize=(figsx,figsy))
    Fpowcum = plt.figure(figsize=(figsx,figsy))
    Ftmin = plt.figure(figsize=(figsx,figsy))
    Ftmean = plt.figure(figsize=(figsx,figsy))
#    Fs = plt.figure(figsize=(figsx,figsy))
#    Fscount = plt.figure(figsize=(figsx,figsy))
#    Fspow = plt.figure(figsize=(figsx,figsy))
#    Fspowcum = plt.figure(figsize=(figsx,figsy))
#    Fstmin = plt.figure(figsize=(figsx,figsy))
#    Fstmean = plt.figure(figsize=(figsx,figsy))
    ihr = 0
    for hr in range(11,19,2):
        ihr = ihr + 1
            #Initialise 2D climatology fields
        strsiz = np.zeros((indxe-indxb,indye-indyb))
        strcount = np.zeros((indxe-indxb,indye-indyb))
        strpow = np.zeros((indxe-indxb,indye-indyb))
        strtirmin = np.zeros((indxe-indxb,indye-indyb))
        strtirmean = np.zeros((indxe-indxb,indye-indyb))
#        sssiz = np.zeros((indxe-indxb,indye-indyb))
#        sscount = np.zeros((indxe-indxb,indye-indyb))
#        sspow = np.zeros((indxe-indxb,indye-indyb))
#        sstirmin = np.zeros((indxe-indxb,indye-indyb))
#        sstirmean = np.zeros((indxe-indxb,indye-indyb))
        #Extract only month/hour combinations of interest
        yrhr = (hours1==hr) & (months1==mo)
        siz4med = sizes1[yrhr]
        day4med = days1[yrhr]
        year4med = years1[yrhr]
        xloc4med = xlocs1[yrhr]
        yloc4med = ylocs1[yrhr]
        tirmin4med = tirmin1[yrhr]
        tirmean4med = tirmean1[yrhr]
        pow4med = powers1[yrhr]
        # Introduce power threshold
        powthresh = powthresh0 * (siz4med/2.)**.5
        siz4med = siz4med[pow4med>powthresh]
        day4med = day4med[pow4med>powthresh]
        year4med = year4med[pow4med>powthresh]
        xloc4med = xloc4med[pow4med>powthresh]
        yloc4med = yloc4med[pow4med>powthresh]
        tirmin4med = tirmin4med[pow4med>powthresh]
        tirmean4med = tirmean4med[pow4med>powthresh]
        pow4med = pow4med[pow4med>powthresh]
        # Search for max-scale structure & find substructures
        siz4str = []
        pow4str = []
        tirmin4str = []
        tirmean4str = []
        xloc4str = []
        yloc4str = []
        count4ss = []
        siz4ss = []
        pow4ss = []
        tirmin4ss = []
        tirmean4ss = []
        xloc4ss = []
        yloc4ss = []
        while True:
            if (siz4med.size==0):
                break
            imax = np.argmax(siz4med)
            if (siz4med[imax]<sizthresh):
                break
            # Main structure
            delind = imax
            iscale = (np.ceil(siz4med[imax]/2./dt)+1).astype(int)
            if (siz4med[imax]<upperthresh):
                siz4str = np.append(siz4str,siz4med[imax])
                pow4str = np.append(pow4str,pow4med[imax])
                tirmin4str = np.append(tirmin4str,tirmin4med[imax])
                tirmean4str = np.append(tirmean4str,tirmean4med[imax])
                xloc4str = np.append(xloc4str,xloc4med[imax])
                yloc4str = np.append(yloc4str,yloc4med[imax])
            # Substructures
            ssind = ( (year4med==year4med[imax]) & (day4med==day4med[imax]) &
                    (xloc4med>xloc4med[imax]-iscale) & (xloc4med<xloc4med[imax]+iscale) &
                    (siz4med<siz4med[imax]) & (yloc4med>yloc4med[imax]-iscale) & (yloc4med<yloc4med[imax]+iscale) )
            if np.count_nonzero(ssind):
#                if (siz4med[imax]<upperthresh):
                if 0:
                    siz4ss = np.append(siz4ss,siz4med[ssind])
                    pow4ss = np.append(pow4ss,pow4med[ssind])
                    tirmin4ss = np.append(tirmin4ss,tirmin4med[ssind])
                    tirmean4ss = np.append(tirmean4ss,tirmean4med[ssind])
                    xloc4ss = np.append(xloc4ss,xloc4med[ssind])
                    yloc4ss = np.append(yloc4ss,yloc4med[ssind])
                delind = np.append(np.nonzero(ssind),imax)
    
            # Remove the structure and substructures!
            siz4med = np.delete(siz4med,delind)
            day4med = np.delete(day4med,delind)
            year4med = np.delete(year4med,delind)
            pow4med = np.delete(pow4med,delind)
            tirmin4med = np.delete(tirmin4med,delind)
            tirmean4med = np.delete(tirmean4med,delind)
            xloc4med = np.delete(xloc4med,delind)
            yloc4med = np.delete(yloc4med,delind)
        
        # mo/hr statistics for each grid point
        # Structures
        for icel in range(len(xloc4str)):
            strsiz[yloc4str[icel],xloc4str[icel]] = strsiz[yloc4str[icel],xloc4str[icel]] + siz4str[icel]
            strpow[yloc4str[icel],xloc4str[icel]] = strpow[yloc4str[icel],xloc4str[icel]] + pow4str[icel]
            strcount[yloc4str[icel],xloc4str[icel]] = strcount[yloc4str[icel],xloc4str[icel]] + 1
            strtirmin[yloc4str[icel],xloc4str[icel]] = strtirmin[yloc4str[icel],xloc4str[icel]] + tirmin4str[icel]
            strtirmean[yloc4str[icel],xloc4str[icel]] = strtirmean[yloc4str[icel],xloc4str[icel]] + tirmean4str[icel]
        #Average for structures
        strsiz[strcount>0] = strsiz[strcount>0] / strcount[strcount>0]
        strpowcum = np.copy(strpow)
        strpow[strcount>0] = strpow[strcount>0] / strcount[strcount>0]
        strtirmin[strcount>0] = strtirmin[strcount>0] / strcount[strcount>0]
        strtirmean[strcount>0] = strtirmean[strcount>0] / strcount[strcount>0]
        # Substructures
#        for icel in range(len(xloc4ss)):
#            sssiz[yloc4ss[icel],xloc4ss[icel]] = sssiz[yloc4ss[icel],xloc4ss[icel]] + siz4ss[icel]
#            sspow[yloc4ss[icel],xloc4ss[icel]] = sspow[yloc4ss[icel],xloc4ss[icel]] + pow4ss[icel]
#            sscount[yloc4ss[icel],xloc4ss[icel]] = sscount[yloc4ss[icel],xloc4ss[icel]] + 1
#            sstirmin[yloc4ss[icel],xloc4ss[icel]] = sstirmin[yloc4ss[icel],xloc4ss[icel]] + tirmin4ss[icel]
#            sstirmean[yloc4ss[icel],xloc4ss[icel]] = sstirmean[yloc4ss[icel],xloc4ss[icel]] + tirmean4ss[icel]
#        #Average for substructures
#        sssiz[sscount>0] = sssiz[sscount>0] / sscount[sscount>0]
#        sspowcum = sspow[:]
#        sspow[sscount>0] = sspow[sscount>0] / sscount[sscount>0]
#        sstirmin[sscount>0] = sstirmin[sscount>0] / sscount[sscount>0]
#        sstirmean[sscount>0] = sstirmean[sscount>0] / sscount[sscount>0]
        
        #Plot:
        #Structure
        #Size
        plt.figure(F.number)
        ttl = "Size, Month: {0:02d}"
        title = ttl.format(mo)
        plt.suptitle(title,fontsize=14)
        A = plt.subplot(2,2,ihr)
        ttl = "{0:02d} UTC"
        title = ttl.format(hr)
        A.set_title(title,fontsize=14)
        IMAGE = M.contourf(xx,yy,strsiz,np.arange(6,100,4))
    #    IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='brown',linewidths=1)
        M.drawcountries(linewidth=0.4)
        M.drawcoastlines(linewidth=0.4)
        M.drawmapboundary()
        parallels = np.arange(-90.,90.,2)
        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
        meridians = np.arange(-180.,180.,2)
        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
        cbar = M.colorbar(IMAGE,location='right',pad="5%")
        #Count
        plt.figure(Fcount.number)
        ttl = "Number of systems, Month: {0:02d}"
        title = ttl.format(mo)
        plt.suptitle(title,fontsize=14)
        A = plt.subplot(2,2,ihr)
        ttl = "{0:02d} UTC"
        title = ttl.format(hr)
        A.set_title(title,fontsize=14)
        IMAGE = M.contourf(xx,yy,strcount,np.arange(0,np.amax(strcount)))
   #     IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='w',linewidths=1)
        M.drawcountries(linewidth=0.4)
        M.drawcoastlines(linewidth=0.4)
        M.drawmapboundary()
        parallels = np.arange(-90.,90.,2)
        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
        meridians = np.arange(-180.,180.,2)
        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
        cbar = M.colorbar(IMAGE,location='right',pad="5%")
        #Power
        plt.figure(Fpow.number)
        ttl = "Power, Month: {0:02d}"
        title = ttl.format(mo)
        plt.suptitle(title,fontsize=14)
        A = plt.subplot(2,2,ihr)
        ttl = "{0:02d} UTC"
        title = ttl.format(hr)
        A.set_title(title,fontsize=14)
#        IMAGE = M.contourf(xx,yy,strpow,np.arange(50,500,20))
        IMAGE = M.contourf(xx,yy,strpow,30)
     #   IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='brown',linewidths=1)
        M.drawcountries(linewidth=0.4)
        M.drawcoastlines(linewidth=0.4)
        M.drawmapboundary()
        parallels = np.arange(-90.,90.,2)
        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
        meridians = np.arange(-180.,180.,2)
        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
        cbar = M.colorbar(IMAGE,location='right',pad="5%")
        #Cumulative Power
        plt.figure(Fpowcum.number)
        ttl = "Cumulative Power, Month: {0:02d}"
        title = ttl.format(mo)
        plt.suptitle(title,fontsize=14)
        A = plt.subplot(2,2,ihr)
        ttl = "{0:02d} UTC"
        title = ttl.format(hr)
        A.set_title(title,fontsize=14)
#        IMAGE = M.contourf(xx,yy,strpowcum,np.arange(50,500,20))
        IMAGE = M.contourf(xx,yy,strpowcum,30)
   #     IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='brown',linewidths=1)
        M.drawcountries(linewidth=0.4)
        M.drawcoastlines(linewidth=0.4)
        M.drawmapboundary()
        parallels = np.arange(-90.,90.,2)
        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
        meridians = np.arange(-180.,180.,2)
        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
        cbar = M.colorbar(IMAGE,location='right',pad="5%")
        #Tmin
        plt.figure(Ftmin.number)
        ttl = "Tmin, Month: {0:02d}"
        title = ttl.format(mo)
        plt.suptitle(title,fontsize=14)
        A = plt.subplot(2,2,ihr)
        ttl = "{0:02d} UTC"
        title = ttl.format(hr)
        A.set_title(title,fontsize=14)
        IMAGE = M.contourf(xx,yy,strtirmin,30)
    #    IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='k',linewidths=1)
        M.drawcountries(linewidth=0.4)
        M.drawcoastlines(linewidth=0.4)
        M.drawmapboundary()
        parallels = np.arange(-90.,90.,2)
        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
        meridians = np.arange(-180.,180.,2)
        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
        cbar = M.colorbar(IMAGE,location='right',pad="5%")
        #Tmean
        plt.figure(Ftmean.number)
        ttl = "Tmean, Month: {0:02d}"
        title = ttl.format(mo)
        plt.suptitle(title,fontsize=14)
        A = plt.subplot(2,2,ihr)
        ttl = "{0:02d} UTC"
        title = ttl.format(hr)
        A.set_title(title,fontsize=14)
        IMAGE = M.contourf(xx,yy,strtirmean,30)
     #   IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='k',linewidths=1)
        M.drawcountries(linewidth=0.4)
        M.drawcoastlines(linewidth=0.4)
        M.drawmapboundary()
        parallels = np.arange(-90.,90.,2)
        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
        meridians = np.arange(-180.,180.,2)
        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
        cbar = M.colorbar(IMAGE,location='right',pad="5%")
        #Substructures
        #Size
#        plt.figure(Fs.number)
#        ttl = "Sub Size, Month: {0:02d}"
#        title = ttl.format(mo)
#        plt.suptitle(title,fontsize=14)
#        A = plt.subplot(2,2,ihr)
#        ttl = "{0:02d} UTC"
#        title = ttl.format(hr)
#        A.set_title(title,fontsize=14)
#        IMAGE = M.contourf(xx,yy,sssiz,np.arange(6,15,1))
#        IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='brown',linewidths=1)
#        M.drawcountries(linewidth=0.4)
#        M.drawcoastlines(linewidth=0.4)
#        M.drawmapboundary()
#        parallels = np.arange(-90.,90.,2)
#        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
#        meridians = np.arange(-180.,180.,2)
#        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
#        cbar = M.colorbar(IMAGE,location='right',pad="5%")
#        #Count
#        plt.figure(Fscount.number)
#        ttl = "Number of subsystems, Month: {0:02d}"
#        title = ttl.format(mo)
#        plt.suptitle(title,fontsize=14)
#        A = plt.subplot(2,2,ihr)
#        ttl = "{0:02d} UTC"
#        title = ttl.format(hr)
#        A.set_title(title,fontsize=14)
#        IMAGE = M.contourf(xx,yy,sscount,np.arange(20))
#        IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='w',linewidths=1)
#        M.drawcountries(linewidth=0.4)
#        M.drawcoastlines(linewidth=0.4)
#        M.drawmapboundary()
#        parallels = np.arange(-90.,90.,2)
#        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
#        meridians = np.arange(-180.,180.,2)
#        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
#        cbar = M.colorbar(IMAGE,location='right',pad="5%")
#        #Power
#        plt.figure(Fspow.number)
#        ttl = "Sub Power, Month: {0:02d}"
#        title = ttl.format(mo)
#        plt.suptitle(title,fontsize=14)
#        A = plt.subplot(2,2,ihr)
#        ttl = "{0:02d} UTC"
#        title = ttl.format(hr)
#        A.set_title(title,fontsize=14)
#        IMAGE = M.contourf(xx,yy,sspow,np.arange(50,300,10))
#        IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='brown',linewidths=1)
#        M.drawcountries(linewidth=0.4)
#        M.drawcoastlines(linewidth=0.4)
#        M.drawmapboundary()
#        parallels = np.arange(-90.,90.,2)
#        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
#        meridians = np.arange(-180.,180.,2)
#        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
#        cbar = M.colorbar(IMAGE,location='right',pad="5%")
#        #Cumulative Power
#        plt.figure(Fspowcum.number)
#        ttl = "Sub Cumulative Power, Month: {0:02d}"
#        title = ttl.format(mo)
#        plt.suptitle(title,fontsize=14)
#        A = plt.subplot(2,2,ihr)
#        ttl = "{0:02d} UTC"
#        title = ttl.format(hr)
#        A.set_title(title,fontsize=14)
#        IMAGE = M.contourf(xx,yy,sspowcum,np.arange(50,300,10))
#        IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='brown',linewidths=1)
#        M.drawcountries(linewidth=0.4)
#        M.drawcoastlines(linewidth=0.4)
#        M.drawmapboundary()
#        parallels = np.arange(-90.,90.,2)
#        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
#        meridians = np.arange(-180.,180.,2)
#        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
#        cbar = M.colorbar(IMAGE,location='right',pad="5%")
#        #Tmin
#        plt.figure(Fstmin.number)
#        ttl = "Sub Tmin, Month: {0:02d}"
#        title = ttl.format(mo)
#        plt.suptitle(title,fontsize=14)
#        A = plt.subplot(2,2,ihr)
#        ttl = "{0:02d} UTC"
#        title = ttl.format(hr)
#        A.set_title(title,fontsize=14)
#        IMAGE = M.contourf(xx,yy,sstirmin,30)
#        IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='k',linewidths=1)
#        M.drawcountries(linewidth=0.4)
#        M.drawcoastlines(linewidth=0.4)
#        M.drawmapboundary()
#        parallels = np.arange(-90.,90.,2)
#        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
#        meridians = np.arange(-180.,180.,2)
#        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
#        cbar = M.colorbar(IMAGE,location='right',pad="5%")
#        #Tmean
#        plt.figure(Fstmean.number)
#        ttl = "Sub Tmean, Month: {0:02d}"
#        title = ttl.format(mo)
#        plt.suptitle(title,fontsize=14)
#        A = plt.subplot(2,2,ihr)
#        ttl = "{0:02d} UTC"
#        title = ttl.format(hr)
#        A.set_title(title,fontsize=14)
#        IMAGE = M.contourf(xx,yy,sstirmean,30)
#        IMAGE1 = M.contour(xx,yy,cover,np.arange(60,61,25),colors='k',linewidths=1)
#        M.drawcountries(linewidth=0.4)
#        M.drawcoastlines(linewidth=0.4)
#        M.drawmapboundary()
#        parallels = np.arange(-90.,90.,2)
#        M.drawparallels(parallels, labels =[1,0,0,0], linewidth = 0.5, fontsize=13)
#        meridians = np.arange(-180.,180.,2)
#        M.drawmeridians(meridians, labels =[0,0,0,1], linewidth = 0.5, fontsize=13)
#        cbar = M.colorbar(IMAGE,location='right',pad="5%")
            
        
    savePattern = "/users/global/cornkle/data/pythonWorkspace/wavelet/figures/VERA2DT20Volta_Pks2DSizes_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
    F.savefig(savefile)
    plt.close(F)
    savePattern = "/users/global/cornkle/data/pythonWorkspace/wavelet/figures/VERA2DT20Volta_Pks2DCount_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
    Fcount.savefig(savefile)
    plt.close(Fcount)
    savePattern = "/users/global/cornkle/data/pythonWorkspace/wavelet/figures/VERA2DT20Volta_Pks2DPow_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
    Fpow.savefig(savefile)
    plt.close(Fpow)
    savePattern = "/users/global/cornkle/data/pythonWorkspace/wavelet/figures/VERA2DT20Volta_Pks2DPowCum_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
    Fpowcum.savefig(savefile)
    plt.close(Fpowcum)
    savePattern = "/users/global/cornkle/data/pythonWorkspace/wavelet/figures/VERA2DT20Volta_Pks2DTmin_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
    Ftmin.savefig(savefile)
    plt.close(Ftmin)
    savePattern = "/users/global/cornkle/data/pythonWorkspace/wavelet/figures/VERA2DT20Volta_Pks2DTmean_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
    Ftmean.savefig(savefile)
    plt.close(Ftmean)

#    savePattern = "../wavelet/figures/VERA2DTai_Pks2DSSizes_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
#    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
#    Fs.savefig(savefile)
#    plt.close(Fs)
#    savePattern = "../wavelet/figures/VERA2DTai_Pks2DSCount_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
#    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
#    Fscount.savefig(savefile)
#    plt.close(Fscount)
#    savePattern = "../wavelet/figures/VERA2DTai_Pks2DSPow_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
#    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
#    Fspow.savefig(savefile)
#    plt.close(Fspow)
#    savePattern = "../wavelet/figures/VERA2DTai_Pks2DSPowCum_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
#    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
#    Fspowcum.savefig(savefile)
#    plt.close(Fspowcum)    
#    savePattern = "../wavelet/figures/VERA2DTai_Pks2DSTmin_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
#    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
#    Fstmin.savefig(savefile)
#    plt.close(Fstmin)
#    savePattern = "../wavelet/figures/VERA2DTai_Pks2DSTmean_{0:d}_{1:d}km_{3:d}powthresh_{2:02d}.png"
#    savefile = savePattern.format(sizthresh,upperthresh,mo,powthresh0)
#    Fstmean.savefig(savefile)
#    plt.close(Fstmean)
    