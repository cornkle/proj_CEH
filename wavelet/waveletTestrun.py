# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:29:48 2016

@author: cornkle
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from wavelet import twod as w2d
import matplotlib.pyplot as plt
from scipy import signal
from scipy import misc
from scipy.ndimage.measurements import label
from skimage import data
from skimage.feature import match_template


def readWaveletData():
    
    dic={'t' : [], 'p' : []}    
    
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
    print(period2d)
    scales2d.shape = (len(scales2d),1,1)
    powerTIR = powerTIR / (scales2d**2)
    #Precip
    powerPCP, scales2d, freqs2d = w2d.cwt2d(pcp, dt, dt, dj=1./12, s0=30./mother2d.flambda(), J=45)
    powerPCP[np.real(powerPCP<=0)] = 0.01
    powerPCP = (np.abs(powerPCP)) ** 2 # Normalized wavelet power spectrum
    scales2d.shape = (len(scales2d),1,1)
    powerPCP = powerPCP / (scales2d**2)
    
    print(powerTIR.shape)
    print(powerPCP.shape)
    
    dic['t']=powerTIR
    dic['p']=powerPCP
    
    return
    
    

def waveletTestrun():
    
    dic=readWaveletData()
    
    t=dic['t']
    p=dic['p']
    
    xx=[]
    yy=[]
    for k in range(10,26):
        t=t[k,:,:]
        p=p[k,:, :]
        corr2 = match_template(t, p[10:-10, 10:-10])
        y,x = np.unravel_index(np.argmax(corr2), corr2.shape)
        xx.append(x)
        yy.append(y)
    ttest=t[20,:,:]
    ptest=p[20,:, :]
    corr2 = match_template(ttest, ptest[10:-10, 10:-10])
    plt.pcolormesh(corr2)
    plt.colorbar()
    plt.plot(10,10, 'bo')
    for k in range(len(xx)):
        plt.plot(xx[k],yy[k], 'ro')
    plt.show()    
    