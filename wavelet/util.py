# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 14:08:18 2016

@author: cornkle
"""
import numpy as np
from wavelet import twod as w2d 
from scipy import ndimage
import pdb
import matplotlib.pyplot as plt
def waveletTP(t, p, dt):
        
    dic= {}    
        
    #2D continuous wavelet analysis:
    #TIR   
    tir=t.copy()
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
    
    return dic


def waveletT(t, dt):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    tir = t.copy()
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


    return dic

def waveletT8(t, dt):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    tir = t.copy()
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

    return dic

def waveletSurface(t, dt):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    tir = t.copy()
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


def waveletLSTA_dom(t, dt):
    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    # dj: distance between scales
    # s0: start scale, approx 2*3*pixel scale (3 pix necessary for one wave)
    # j: number of scales

    # 2D continuous wavelet analysis:
    # TIR
    tir = t.copy()
    tir = tir
    #
    #tir = tir - np.mean(tir)

    #tir = tir + np.min(tir)
    #tir[tir < 0] = 0

    mother2d = w2d.Mexican_hat()

    powerTIRR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.28, s0=18. / mother2d.flambda(), J=14)  # s0=30./
    #powerTIRR[np.real(powerTIRR <= 0)] = 0

    powerTIR = (np.abs(powerTIRR)) * (np.abs(powerTIRR))  # Normalized wavelet power spectrum
    #iszero = np.where(np.abs(tir) < 0.25)
    #powerTIR[:,iszero[0],iszero[1]] = 0
    #powerTIR[neg] = powerTIR[neg]*-1
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    maxperpix = np.argmax(powerTIR, axis=0)
    dom_scale = np.zeros_like(tir)
    scales = (period2d / 2.)
    for i in range(maxperpix.shape[0]):
        for j in range(maxperpix.shape[1]):
            max = maxperpix[i,j]
            pt = np.real(powerTIRR[:,i,j])
            ptt = powerTIR[:,i,j]
            scal = scales[max]
            ptest = pt[max]
            pttest = ptt[max]
            if ptest < 0:
                scal = scal*-1
            if pttest < 0.05:
                scal=np.nan

            dom_scale[i,j] = scal

    dic['power'] = powerTIR
    dic['scales'] = scales
    dic['dominant'] = dom_scale

    return dic


def waveletLSTA_domLocMax(t, dt):
    dic = {}

    tir = t.copy()
    tir = tir

    tir = tir - np.mean(tir)
    #tir[np.abs(tir)<0.25]=0
    mother2d = w2d.Mexican_hat()

    powerTIRR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.28, s0=18. / mother2d.flambda(), J=14)  # s0=30./
    # powerTIRR[np.real(powerTIRR <= 0)] = 0

    neg = np.where(np.real(powerTIRR) < 0)

    powerTIR = (np.abs(powerTIRR)) * (np.abs(powerTIRR))  # Normalized wavelet power spectrum
    # powerTIR[neg] = powerTIR[neg]*-1
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    dom_scale = np.zeros_like(tir)
    scales = (period2d / 2.)
    for i in range(powerTIR.shape[1]):
        for j in range(powerTIR.shape[2]):
            parr = powerTIR[:, i, j]

            maxoutt = (parr == ndimage.maximum_filter(parr, 15, mode='reflect'))

            try:
                x = np.where((maxoutt == 1))  # ((wl >= np.percentile(wl[wl >= 0.5], 90)) &
            except IndexError:
                continue
            ind = (x[0])[0]
            pt = np.real(powerTIRR[:, i, j])

            #print('Max scales', scales[x], x)

            scal = scales[ind]
            ptest = pt[ind]
            pttest = parr[ind]
            if ptest < 0:
                scal = scal * -1
            if pttest < 0.05:
                scal = np.nan

            dom_scale[i, j] = scal

    dic['power'] = powerTIR
    dic['scales'] = scales
    dic['dominant'] = dom_scale

    return dic

def waveletLSTA_domSmall(t, dt):
    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    # dj: distance between scales
    # s0: start scale, approx 2*3*pixel scale (3 pix necessary for one wave)
    # j: number of scales

    # 2D continuous wavelet analysis:
    # TIR
    tir = t.copy()
    tir = tir
    #tir[tir < 0] = 0
    tir = tir - np.mean(tir)
    mother2d = w2d.Mexican_hat()

    powerTIRR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.28, s0=18. / mother2d.flambda(), J=14)  # s0=30./
    #powerTIRR[np.real(powerTIRR <= 0)] = 0

    neg = np.where(powerTIRR<0)

    powerTIR = (np.abs(powerTIRR)) * (np.abs(powerTIRR))  # Normalized wavelet power spectrum
    #powerTIR[neg] = powerTIR[neg]*-1
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    maxperpix = np.argmax(powerTIR, axis=0)
    dom_scale = np.zeros_like(tir)*np.nan
    scales = (period2d / 2.)
    for i in range(len(scales)-1, -1, -1):

        scal = scales[i]
        parray = powerTIR[i, :, :]
        sign = np.real(powerTIRR[i, :, :]) < 0

        parray[parray <= np.percentile(parray[parray >= 0.1], 75)] = np.nan
        dom_scale[np.isfinite(parray)]=scal
        (dom_scale)[sign] =  (dom_scale)[sign]*-1

    dic['power'] = powerTIR
    dic['scales'] = scales
    dic['dominant'] = dom_scale

    return dic

def waveletLSTA_test(t, dt):

    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    tir = t.copy()

    tir = tir - np.nanmean(tir)
    tir[np.isnan(tir)] = 0
    #tir[tir < 0] = 0

    nanpos =  np.isnan(tir)
    mother2d = w2d.Mexican_hat()

    powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.28, s0=18. / mother2d.flambda(), J=14)  # s0=30./
    powerTIR[np.real(powerTIR <= 0)] = 0
    powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)
    powerTIR[:,nanpos[0], nanpos[1]]= np.nan
    dic['power'] = powerTIR
    dic['scales'] = (period2d / 2.)

    return dic


def waveletLSTA(t, dt, method=None):
    dic = {}

    # 2D continuous wavelet analysis:
    # TIR
    # dj: distance between scales
    # s0: start scale, approx 2*3*pixel scale (3 pix necessary for one wave)
    # j: number of scales

    # 2D continuous wavelet analysis:
    # TIR
    tir = t.copy()
    tir = tir - np.mean(tir)
    mother2d = w2d.Mexican_hat()
    #tir[tir<0] = 0
    powerTIRR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.28, s0=18. / mother2d.flambda(), J=14)

    isneg = np.where(powerTIRR<0)

    powerTIR = (np.abs(powerTIRR)) * (np.abs(powerTIRR))  # Normalized wavelet power spectrum

    powerTIR[isneg] = powerTIR[isneg]*-1
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    if method == 'dry':
        pos = np.where(tir<=0)
        powerTIR[:, pos[0], pos[1]] = np.nan

    if method == 'wet':
        pos = np.where(tir>=0)
        powerTIR[:, pos[0], pos[1]] = np.nan

    dic['power'] = powerTIR
    dic['scales'] = (period2d / 2.)

    return dic


def waveletSurfaceneg(t, dt):

    dic = {}
    tir = t.copy()
    # 2D continuous wavelet analysis:
    # TIR
    tir = tir*(-1)
    tir[tir < 0] = -100
    tir[np.isnan(tir)]= -100
    #tir = tir - np.mean(tir)
    mother2d = w2d.Mexican_hat()

    powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.75, s0=1500. / mother2d.flambda(), J=10)  # s0=30./
    powerTIR[np.real(powerTIR <= 0)] = 0.001
    powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
    period2d = 1. / freqs2d
    scales2d.shape = (len(scales2d), 1, 1)
    powerTIR = powerTIR / (scales2d * scales2d)

    dic['power'] = powerTIR
    dic['scales'] = (period2d / 2.)

    return dic