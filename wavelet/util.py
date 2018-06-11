# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 14:08:18 2016

@author: cornkle
"""
import numpy as np
#from wavelet import twod as w2d
from scipy import ndimage
from wavelet import wav
import pdb
import matplotlib.pyplot as plt

def read_dic(dic):
    dt = dic['dx']
    dist = dic['dist']
    start = dic['start']
    nb = dic['nb']

    return dt, dist, start, nb

def _create_dic(dx, dist, start, nb):

    dic = {}
    dic['dx'] = dx
    dic['dist'] = dist
    dic['start'] = start
    dic['nb'] = nb

    return dic

############ Frequently used datasets
DATASETS = {
    'METEOSAT5K': _create_dic(5, 1 / 12., 15, 45),
    'GRIDSAT': _create_dic(8, 1 / 12., 15, 45),
    'METSRFC': _create_dic(3, 0.45, 9, 10)
}


def waveletTP(t, p, dx=None, dist=None,start=None, nb=None, dataset=None):

    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx,dist,nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    #2D continuous wavelet analysis:
    #TIR   
    tir=t.copy()
    tir[tir>0] = 0
    tir = tir - np.mean(tir)

    obj = wav.wavelet(dx, dist, nb, start=start)

    #TIR
    coeffsTIR, powerTIR = obj.calc_coeffs(tir, ge_thresh=0, fill=0.01)
    #Precip
    coeffsPCP, powerPCP = obj.calc_coeffs(p, le_thresh=0, fill=0.01)
        
    dic['t']=powerTIR
    dic['p']=powerPCP
    dic['scales'] = obj.scales
    dic['res'] = obj.res
    return dic


def waveletT(t, dx=None, dist=None,start=None, nb=None, dataset=None):

    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx, dist, nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    tir = t.copy()
    tir[tir > 0] = 0
    tir = tir - np.mean(tir)

    obj = wav.wavelet(dx, dist, nb, start=start)
    # TIR
    coeffsTIR, powerTIR = obj.calc_coeffs(tir, ge_thresh=0, fill=0.01)

    dic['t'] = powerTIR
    dic['scales'] = obj.scales
    dic['res'] = obj.res

    return dic
#
# def waveletT8(t, dt):
#
#     dic = {}
#
#     # 2D continuous wavelet analysis:
#     # TIR
#     tir = t.copy()
#     tir[tir > 0] = 0
#     tir = tir - np.mean(tir)
#     mother2d = w2d.Mexican_hat()
#
#     powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=1. / 12, s0=40. / mother2d.flambda(), J=45)  # s0=30./
#     powerTIR[np.real(powerTIR >= 0)] = 0.01
#     powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
#     period2d = 1. / freqs2d
#     scales2d.shape = (len(scales2d), 1, 1)
#     powerTIR = powerTIR / (scales2d * scales2d)
#
#     dic['t'] = powerTIR
#     dic['scales'] = (period2d / 2.)
#
#     return dic
#
# def waveletSurface(t, dt):
#
#     dic = {}
#
#     # 2D continuous wavelet analysis:
#     # TIR
#     tir = t.copy()
#     tir[tir < 0] = 100
#     tir[np.isnan(tir)]= 100
#     #tir = tir - np.mean(tir)
#     mother2d = w2d.Mexican_hat()
#
#     powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.75, s0=1500. / mother2d.flambda(), J=10)  # s0=30./
#
#     #powerTIR[np.real(powerTIR >= 0)] = 0.01
#     powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
#     period2d = 1. / freqs2d
#     scales2d.shape = (len(scales2d), 1, 1)
#     powerTIR = powerTIR / (scales2d * scales2d)
#
#     dic['power'] = powerTIR
#     dic['scales'] = (period2d / 2.)
#
#     return dic
#
#
# def waveletLSTA_dom(t, dt):
#     dic = {}
#
#     # 2D continuous wavelet analysis:
#     # TIR
#     # dj: distance between scales
#     # s0: start scale, approx 2*3*pixel scale (3 pix necessary for one wave)
#     # j: number of scales
#
#     # 2D continuous wavelet analysis:
#     # TIR
#     tir = t.copy()
#     tir = tir
#
#     mother2d = w2d.Mexican_hat()
#
#     #powerTIRR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.28, s0=18. / mother2d.flambda(), J=14)  # s0=30./
#     powerTIRR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.45, s0=18. / mother2d.flambda(), J=10)
#     #powerTIRR[np.real(powerTIRR <= 0)] = 0
#
#     powerTIR = (np.abs(powerTIRR)) * (np.abs(powerTIRR))  # Normalized wavelet power spectrum
#     #iszero = np.where(np.abs(tir) < 0.25)
#     #powerTIR[:,iszero[0],iszero[1]] = 0
#     #powerTIR[neg] = powerTIR[neg]*-1
#     period2d = 1. / freqs2d
#     scales2d.shape = (len(scales2d), 1, 1)
#     powerTIR = powerTIR / (scales2d * scales2d)
#
#     maxperpix = np.argmax(powerTIR, axis=0)
#     perc = []
#     for i in range(powerTIR.shape[0]):
#         sp = np.percentile(powerTIR[i,:,:],20)
#         perc.append(sp)
#     perc = np.array(perc)
#
#     dom_scale = np.zeros_like(tir)
#     power_scale = np.zeros_like(tir)
#     scales = (period2d / 2.)
#     for i in range(maxperpix.shape[0]):
#         for j in range(maxperpix.shape[1]):
#             max = maxperpix[i,j]
#             pt = np.real(powerTIRR[:,i,j])
#             ptt = powerTIR[:,i,j]
#             scal = scales[max]
#             ptest = pt[max]
#             pttest = ptt[max]
#             if ptest < 0:
#                 scal = scal*-1
#                 pttest = pttest*-1
#             # if pttest < perc[max]:
#             #     scal=np.nan
#             power_scale[i,j] = pttest
#             dom_scale[i,j] = scal
#
#     dic['power'] = power_scale
#     dic['scales'] = scales
#     dic['dominant'] = dom_scale
#
#     return dic
#
#
#
# def waveletLSTA_domLocMax(t, dt):
#     dic = {}
#
#     tir = t.copy()
#     mother2d = w2d.Mexican_hat()
#
#     powerTIRR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.45, s0=18. / mother2d.flambda(), J=10)  # s0=30./
#     #powerTIRR[np.real(powerTIRR <= 0)] = 0
#
#     powerTIR = (np.abs(powerTIRR)) * (np.abs(powerTIRR))  # Normalized wavelet power spectrum
#     period2d = 1. / freqs2d
#     scales2d.shape = (len(scales2d), 1, 1)
#     powerTIR = powerTIR / (scales2d * scales2d)
#
#     dom_scale = np.zeros_like(tir)*np.nan
#     scales = (period2d / 2.)
#     for i in range(powerTIR.shape[1]):
#         for j in range(powerTIR.shape[2]):
#             parr = powerTIR[:, i, j]
#
#             maxoutt = (parr == ndimage.maximum_filter(parr, 15, mode='reflect'))
#
#             try:
#                 x = np.where((maxoutt == 1) )  #& (np.abs(parr) > scales ** .5)# ((wl >= np.percentile(wl[wl >= 0.5], 90)) &
#             except IndexError:
#                 continue
#             try:
#                 ind = (x[0])[0]
#             except IndexError:
#                 continue
#             pt = np.real(powerTIRR[:, i, j])
#
#            # print('Max scales', scales[x], x)
#
#             scal = scales[ind]
#             ptest = pt[ind]
#             tt = tir[i,j]
#             pttest = parr[ind]
#             if (ptest < 0) | (tt < -1.) :
#                 scal = scal * -1
#             if pttest < 0.02:
#                 scal = np.nan
#
#             dom_scale[i, j] = scal
#
#     dic['power'] = powerTIR
#     dic['scales'] = scales
#     dic['dominant'] = dom_scale
#
#     return dic
#
#
#
# def waveletLSTA_both(t, dt, dom=False):
#
#     def dom_get(wseries, tpoint, scales):
#
#         values = (np.abs(wseries)) * (np.abs(wseries))
#         values = values / (scales * scales)
#
#         maxoutt = (values == ndimage.maximum_filter(values, 15, mode='reflect'))
#
#         try:
#             x = np.where(
#                 (maxoutt == 1))  # & (np.abs(parr) > scales ** .5)# ((wl >= np.percentile(wl[wl >= 0.5], 90)) &
#         except IndexError:
#             return False
#         try:
#             ind = (x[0])[0]
#         except IndexError:
#             return False
#
#         scal = scales[ind]
#         wavelet_pos_neg = np.real(wseries)[ind]
#
#         if (wavelet_pos_neg < 0) : #| (tpoint < -1.)
#             scal = scal * -1
#         # if values[ind] < 0.02:
#         #     scal = np.nan
#         return scal
#
#     dic = {}
#
#     tir = t.copy()
#     mother2d = w2d.Mexican_hat()
#     #powerTIRR_dry, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.45, s0=18. / mother2d.flambda(), J=10)
#     powerTIRR_dry, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.28, s0=18. / mother2d.flambda(), J=14)
#     powerTIRR_dry[np.real(powerTIRR_dry <= 0)] = 0
#
#     tir = tir * -1
#     #powerTIRR_wet, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.45, s0=18. / mother2d.flambda(), J=10)
#     powerTIRR_wet, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.28, s0=18. / mother2d.flambda(), J=14)
#     powerTIRR_wet[np.real(powerTIRR_wet <= 0)] = 0
#
#     period2d = 1. / freqs2d
#     scales = (period2d / 2.)
#     scales2d.shape = (len(scales2d), 1, 1)
#
#     powerTIR_dry = (np.abs(powerTIRR_dry)) * (np.abs(powerTIRR_dry))
#     powerTIR_dry = powerTIR_dry / (scales2d * scales2d)
#
#     powerTIR_wet = (np.abs(powerTIRR_wet)) * (np.abs(powerTIRR_wet))
#     powerTIR_wet = powerTIR_wet / (scales2d * scales2d)
#
#     for id, s in enumerate(scales):
#         (powerTIR_wet[id, :, :])[powerTIR_wet[id, :, :] <= np.percentile(powerTIR_wet[id, :, :],50)] = 0
#         (powerTIR_dry[id, :, :])[powerTIR_dry[id, :, :] <= np.percentile(powerTIR_dry[id, :, :],50)] = 0
#
#     dic['power_dry'] = powerTIR_dry
#     dic['power_wet'] = powerTIR_wet
#     dic['scales'] = scales
#
#     if dom:
#         dom_scale = np.zeros_like(tir) * np.nan
#         for i in range(powerTIRR_dry.shape[1]):
#             for j in range(powerTIRR_dry.shape[2]):
#                 parr_dry = powerTIRR_dry[:, i, j]
#                 parr_wet = powerTIRR_wet[:, i, j]
#                 tpoint = t[i,j]
#
#                 scal_dry = dom_get(parr_dry, tpoint, scales)
#                 scal_wet = dom_get(parr_wet, tpoint, scales)
#
#                 scal_pos = np.argmin([np.abs(scal_dry), np.abs(scal_wet)])
#                 scal_arr = scal_dry # [scal_dry,scal_wet][scal_pos]
#
#                 dom_scale[i, j] = scal_arr
#         dic['dominant'] = dom_scale
#
#     return dic
#
#
# def waveletSurfaceneg(t, dt):
#
#     dic = {}
#     tir = t.copy()
#     # 2D continuous wavelet analysis:
#     # TIR
#     tir = tir*(-1)
#     tir[tir < 0] = -100
#     tir[np.isnan(tir)]= -100
#     #tir = tir - np.mean(tir)
#     mother2d = w2d.Mexican_hat()
#
#     powerTIR, scales2d, freqs2d = w2d.cwt2d(tir, dt, dt, dj=0.75, s0=1500. / mother2d.flambda(), J=10)  # s0=30./
#     powerTIR[np.real(powerTIR <= 0)] = 0.001
#     powerTIR = (np.abs(powerTIR)) * (np.abs(powerTIR))  # Normalized wavelet power spectrum
#     period2d = 1. / freqs2d
#     scales2d.shape = (len(scales2d), 1, 1)
#     powerTIR = powerTIR / (scales2d * scales2d)
#
#     dic['power'] = powerTIR
#     dic['scales'] = (period2d / 2.)
#
#     return dic