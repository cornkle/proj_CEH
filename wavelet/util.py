# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 14:08:18 2016

@author: cornkle
"""
import numpy as np
from wavelet import twod as w2d 
from scipy import ndimage
import pdb
def waveletTP(t, p, dt, max = False):
        
    dic= {}    
        
    #2D continuous wavelet analysis:
    #TIR   
    tir=t
    tir[tir>0] = 0
    tir = tir - np.mean(tir) 
    mother2d = w2d.Mexican_hat()
    
    powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=1./12, s0=30./mother2d.flambda(), J=45)  # s0=30./
    powerTIR[np.real(powerTIR>=0)] = 0.01
    powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR)) # Normalized wavelet power spectrum
    period2d = 1. / freqs2d    
    scales2d.shape = (len(scales2d),1,1)
    powerTIR = powerTIR / (scales2d*scales2d)
    
    #Precip
    powerPCP, scales2d, freqs2d = w2d.cwt2d(p, dt, dt, dj=1./12, s0=30./mother2d.flambda(), J=45)
    powerPCP[np.real(powerPCP<=0)] = 0.01
    powerPCP = (np.abs(powerPCP)) * (np.abs(powerPCP)) # Normalized wavelet power spectrum
    scales2d.shape = (len(scales2d),1,1)
    powerPCP = powerPCP / (scales2d*scales2d)
        
    dic['t']=powerTIR
    dic['p']=powerPCP
    dic['scales'] = (period2d/2.)

    if max:
        # Find maxima of the 2D wavelet spectrum
        # TIR
        maxpowerTIR = (powerTIR == ndimage.maximum_filter(powerTIR, size=(5, 5, 5), mode='constant', cval=np.amax(powerTIR) + 1))
        maxpowerTIR = maxpowerTIR.astype(int)
        zpks, ypks, xpks = np.where((maxpowerTIR == 1) & (powerTIR > 10))
        dic['z'] = zpks
        dic['y'] = ypks
        dic['x'] = xpks
    
    
    return dic


def waveletT(t, dt, max=False):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    tir = t
    tir[tir > 0] = 0
    tir = tir - np.mean(tir)
    mother2d = w2d.Mexican_hat()

    powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=1. / 12, s0=30. / mother2d.flambda(), J=45)  # s0=30./
    powerTIR[np.real(powerTIR >= 0)] = 0.01
    powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    dic['t'] = powerTIR
    dic['scales'] = (period2d / 2.)

    if max:
        # Find maxima of the 2D wavelet spectrum
        # TIR
        maxpowerTIR = (powerTIR == ndimage.maximum_filter(powerTIR, size=(5, 5, 5), mode='constant', cval=np.amax(powerTIR) + 1))
        maxpowerTIR = maxpowerTIR.astype(int)
        zpks, ypks, xpks = np.where((maxpowerTIR == 1) & (powerTIR > 10))
        dic['z'] = zpks
        dic['y'] = ypks
        dic['x'] = xpks

    return dic

def waveletT8(t, dt, max=False):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    tir = t
    tir[tir > 0] = 0
    tir = tir - np.mean(tir)
    mother2d = w2d.Mexican_hat()

    powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=1. / 12, s0=40. / mother2d.flambda(), J=45)  # s0=30./
    powerTIR[np.real(powerTIR >= 0)] = 0.01
    powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    dic['t'] = powerTIR
    dic['scales'] = (period2d / 2.)

    if max:
        # Find maxima of the 2D wavelet spectrum
        # TIR
        maxpowerTIR = (powerTIR == ndimage.maximum_filter(powerTIR, size=(5, 5, 5), mode='constant', cval=np.amax(powerTIR) + 1))
        maxpowerTIR = maxpowerTIR.astype(int)
        zpks, ypks, xpks = np.where((maxpowerTIR == 1) & (powerTIR > 10))
        dic['z'] = zpks
        dic['y'] = ypks
        dic['x'] = xpks

    return dic

def waveletSurface(t, dt):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    tir = t
    tir[tir < 0] = 100
    tir[np.isnan(tir)]= 100
    #tir = tir - np.mean(tir)
    mother2d = w2d.Mexican_hat()

    powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.75, s0=1500. / mother2d.flambda(), J=10)  # s0=30./
    #powerTIR[np.real(powerTIR >= 0)] = 0.01
    powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    dic['power'] = powerTIR
    dic['scales'] = (period2d / 2.)



    return dic

def waveletLSTA(t, dt, dry=None, wet=None):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    # dj: distance between scales
    # s0: start scale, approx 2*3*pixel scale (3 pix necessary for one wave)
    # j: number of scales
    if dry==wet:
        'Please set either dry or wet keyword. Either none or both given.'
    if wet:
        tir = t.copy()
        tir[np.isnan(tir)] = 0
        tir = tir*(-1)
    else:
        tir = t.copy()
        tir[np.isnan(tir)] = 0
    tir[tir < 0] = 0
    tir[np.isnan(tir)]= 0
    #tir = tir - np.mean(tir)
    mother2d = w2d.Mexican_hat()

    powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.3, s0=18/mother2d.flambda(), J=20)  # s0=30./
    #powerTIR[np.real(powerTIR >= 0)] = 0.01
    powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    dic['power'] = powerTIR
    dic['scales'] = (period2d / 2.)


    return dic




    return dic

def waveletSurfaceneg(t, dt):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    tir = t*(-1)
    tir[tir < 0] = -100
    tir[np.isnan(tir)]= -100
    #tir = tir - np.mean(tir)
    mother2d = w2d.Mexican_hat()

    powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.75, s0=1500. / mother2d.flambda(), J=10)  # s0=30./
    powerTIR[np.real(powerTIR <= 0)] = 0.01
    powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    dic['power'] = powerTIR
    dic['scales'] = (period2d / 2.)



    return dic

def waveletnegT(t, dt, max=False):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    tir = t*(-1)
    tir[tir < 0] = 0

    tir = tir - np.mean(tir)
    mother2d = w2d.Mexican_hat()

    powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=1. / 12, s0=30. / mother2d.flambda(), J=45)  # s0=30./
    powerTIR[np.real(powerTIR <= 0)] = 0.01
    powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    dic['t'] = powerTIR
    dic['scales'] = (period2d / 2.)

    if max:
        # Find maxima of the 2D wavelet spectrum
        # TIR
        maxpowerTIR = (powerTIR == ndimage.maximum_filter(powerTIR, size=(5, 5, 5), mode='constant', cval=np.amax(powerTIR) + 1))
        maxpowerTIR = maxpowerTIR.astype(int)
        zpks, ypks, xpks = np.where((maxpowerTIR == 1) & (powerTIR > 10))
        dic['z'] = zpks
        dic['y'] = ypks
        dic['x'] = xpks


    return dic