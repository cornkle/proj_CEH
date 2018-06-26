#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import twod as w2d

#Read Conni's file
ccfile = "/users/global/cornkle/wtest/845_mt_wavelet_test.nc"
fh = Dataset(ccfile, mode='r')
pcp = fh.variables['trmm'][:]
tir = fh.variables['msg'][:]
lat = fh.variables['lat'][:]
lon = fh.variables['lon'][:]

dt = 5. #grid spacing in km

#2D continuous wavelet analysis:
#TIR
tiror = np.copy(tir)
tir[tir>0] = 0
tir = tir - np.mean(tir) 
mother2d = w2d.Mexican_hat()
powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=1./12, s0=30./mother2d.flambda(), J=45)
powerTIR[np.real(powerTIR>=0)] = 0.01
powerTIR = (np.abs(powerTIR)) ** 2 # Normalized wavelet power spectrum
period2d = 1. / freqs2d
scales2d.shape = (len(scales2d),1,1)
powerTIR = powerTIR / (scales2d**2)
#Precip
powerPCP, scales2d, freqs2d = w2d.cwt2d(pcp, dt, dt, dj=1./12, s0=30./mother2d.flambda(), J=45)
powerPCP[np.real(powerPCP<=0)] = 0.01
powerPCP = (np.abs(powerPCP)) ** 2 # Normalized wavelet power spectrum
scales2d.shape = (len(scales2d),1,1)
powerPCP = powerPCP / (scales2d**2)

#Plot:
klevels = 2. ** np.arange(-1,11)
# TIR
F = plt.figure(figsize=(24,14))
A = plt.subplot(3,4,1)
title = "IRT & Precip"
A.set_title(title,fontsize=12)
plt.contourf(lon,lat,tir,50)
#plt.colorbar()
plt.contour(lon,lat,pcp,np.arange(5,67,20),colors='w',linewidths=2)
#Plot wavelets
for isc in range(11):
    pltscale = isc*4
    B = plt.subplot(3,4,isc+2)
    ttl = "IRT CWT at {0:d} km"
    title = ttl.format(int(np.rint(period2d[pltscale]/2.)))
    B.set_title(title,fontsize=12)
    plt.contourf(lon,lat,np.log2(powerTIR[pltscale,:,:]), np.log2(klevels),   # here he uses log2
        cmap=plt.cm.jet)
    plt.contour(lon,lat,pcp,np.arange(5,67,20),colors='k',linewidths=2)
plt.tight_layout(pad=4.0)
plt.show()
# Precip
F = plt.figure(figsize=(24,14))
A = plt.subplot(3,4,1)
title = "IRT & Precip"
A.set_title(title,fontsize=12)
plt.contourf(lon,lat,tir,50)
#plt.colorbar()
plt.contour(lon,lat,pcp,np.arange(5,67,20),colors='w',linewidths=2)
#Plot wavelets
for isc in range(11):
    pltscale = isc*4
    B = plt.subplot(3,4,isc+2)
    ttl = "Precip CWT at {0:d} km"
    title = ttl.format(int(np.rint(period2d[pltscale]/2.)))
    B.set_title(title,fontsize=12)
    plt.contourf(lon,lat,np.log2(powerPCP[pltscale,:,:]), np.log2(klevels),
        cmap=plt.cm.jet)
    plt.contour(lon,lat,pcp,np.arange(5,67,20),colors='k',linewidths=2)
plt.tight_layout(pad=4.0)
plt.show()
