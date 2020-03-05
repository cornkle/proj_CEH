# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 14:08:18 2016

@author: cornkle
"""
import numpy as np
#from wavelet import twod as w2d
from scipy import ndimage
from wavelet import wav, wav1d
import pdb
import matplotlib.pyplot as plt
import ipdb

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
    'METEOSAT5K': _create_dic(5, 1 / 12., 15, 45),  # resolution, distance between scales,start scale, number of scales
    'METEOSAT5K_vera': _create_dic(5, 0.5, 25, 2),  #28     0.5
    'METEOSAT5K_veraTest': _create_dic(5, 0.5, 15, 2),  #28     0.5
    'METEOSAT8K': _create_dic(8, 1 / 12., 24, 40),
    'METEOSAT10K': _create_dic(10, 1 / 12., 30, 40),
    'GRIDSAT': _create_dic(8, 1 / 12., 24, 40),
    'METSRFC': _create_dic(3, 0.45, 9, 10), # nb =14 also tested
    'METSRFC_LS': _create_dic(3, 0.4, 9, 15), # larger scales
    'METSRFC24': _create_dic(3, 0.25, 50, 10), # nb =14 also tested
    'LSTATREND5K': _create_dic(5.55, 0.4, 16, 10)  # nb =14 also tested
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
    coeffsPCP, powerPCP = obj.calc_coeffs(p, le_thresh=0, fill=-0.01)

    dic['t']=powerTIR
    dic['p']=powerPCP
    dic['scales'] = obj.scales
    dic['res'] = obj.res
    return dic

def waveletT(t, dx=None, dist=None,start=None, nb=None, dataset=None):


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

    dic['t']=powerTIR
    dic['scales'] = obj.scales
    dic['res'] = obj.res
    dic['coeffs'] = coeffsTIR
    return dic



def waveletT_normalised(t, dx=None, dist=None,start=None, nb=None, dataset=None):

    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx,dist,nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    #2D continuous wavelet analysis:
    #TIR
    tir=t.copy()
    #tir[tir>0] = 0
    tir = tir - np.mean(tir)

    obj = wav.wavelet(dx, dist, nb, start=start)

    #TIR
    coeffsTIR, powerTIR = obj.calc_coeffs(tir, ge_thresh=0, fill=0.01)

    dic['t']=powerTIR
    dic['scales'] = obj.scales
    dic['res'] = obj.res
    return dic


def applyHat(t, dx=None, dist=None,start=None, nb=None, dataset=None):

    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx, dist, nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    tir = t.copy()
    #tir[tir > 0] = 0
    tir = tir #- np.mean(tir)

    obj = wav.wavelet(dx, dist, nb, start=start)
    # TIR
    coeffsTIR, powerTIR = obj.calc_coeffs(tir, ge_thresh=0, fill=0.01)

    dic['power'] = powerTIR
    dic['scales'] = obj.scales
    dic['res'] = obj.res
    dic['coeffs'] = coeffsTIR

    return dic

def applyHat_pure(t, dx=None, dist=None,start=None, nb=None, dataset=None):

    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx, dist, nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    tir = t.copy()
    #tir[tir > 0] = 0
    tir = tir #- np.mean(tir)

    obj = wav.wavelet(dx, dist, nb, start=start)
    # TIR
    coeffsTIR, powerTIR = obj.calc_coeffs(tir)

    dic['power'] = powerTIR
    dic['scales'] = obj.scales
    dic['res'] = obj.res
    dic['coeffs'] = coeffsTIR

    return dic
#
#
def LSTA_maxPowerScale(t, dx=None, dist=None, start=None, nb=None, dataset=None):
    """
    Calculates dominant scale per pixel in an image i.e. it scans each scale column in the 3d wavelet
     power spectrum and identifies maximum power.
    :param t: data , numpy array
    :param dx: resolution of the data in unit of output (e.g. 3 if it's 3km)
    :param dist: distance between scale disaggregation levels, will be logarithmic
    :param start: smallest scale for scale decomposition
    :param nb: number of scales to decompose into
    :param dataset: abbreviation for existing dataset dictionary for quick use
    :return:
    """
    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx, dist, nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    tir = t.copy()
    obj = wav.wavelet(dx, dist, nb, start=start)
    coeffsTIR, powerTIR = obj.calc_coeffs(tir)


    maxperpix = np.argmax(powerTIR, axis=0)
    perc = []
    for i in range(powerTIR.shape[0]):
        sp = np.percentile(powerTIR[i,:,:],20)
        perc.append(sp)
    perc = np.array(perc)  ## potential filtering of low powers

    dom_scale = np.zeros_like(tir)
    power_scale = np.zeros_like(tir)
    scales = obj.scales
    for i in range(maxperpix.shape[0]):
        for j in range(maxperpix.shape[1]):
            max = maxperpix[i,j]
            coeffs = np.real(coeffsTIR[:,i,j])
            power = powerTIR[:,i,j]
            scal = scales[max]
            coeffs_max = coeffs[max]
            power_max = power[max]
            if coeffs_max < 0:
                scal = scal*-1
                power_max = power_max*-1
            # if pttest < perc[max]:
            #     scal=np.nan
            power_scale[i,j] = power_max
            dom_scale[i,j] = scal

    dic['power'] = power_scale
    dic['scales'] = scales
    dic['dominant'] = dom_scale

    return dic
#
#
#
def LSTA_LocalMax(t, dx=None, dist=None, start=None, nb=None, dataset=None):
    """
    Calculates local power maxima in a 2d image within a scale column
    allowing superimposed maxima detection i.e. preferential pick
     of small scales
    :param t: data , numpy array
    :param dx: resolution of the data in unit of output (e.g. 3 if it's 3km)
    :param dist: distance between scale disaggregation levels, will be logarithmic
    :param start: smallest scale for scale decomposition
    :param nb: number of scales to decompose into
    :param dataset: abbreviation for existing dataset dictionary for quick use
    :return:
    """

    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx, dist, nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    tir = t.copy()
    obj = wav.wavelet(dx, dist, nb, start=start)
    coeffsTIR, powerTIR = obj.calc_coeffs(tir)

    dom_scale = np.zeros_like(tir)*np.nan
    scales = obj.scales
    for i in range(powerTIR.shape[1]):
        for j in range(powerTIR.shape[2]):
            power = powerTIR[:, i, j]

            maxoutt = (power == ndimage.maximum_filter(power, 15, mode='reflect')) # local max in 15x15 pixel distance

            try:
                x = np.where((maxoutt == 1) )  #& (np.abs(parr) > scales ** .5)# ((wl >= np.percentile(wl[wl >= 0.5], 90)) &
            except IndexError:
                continue
            try:
                ind = (x[0])[0]
            except IndexError:
                continue
            coeffs = np.real(coeffsTIR[:, i, j])

           # print('Max scales', scales[x], x)

            scal = scales[ind]
            coeff_locMax = coeffs[ind]
            tt = tir[i,j]
            power_locMax = power[ind]
            if (coeff_locMax < 0) | (tt < -1.) :  # coefficient negative or small negative T
                scal = scal * -1
            if power_locMax < 0.02:   # minimum power val, should find a more objective way, like standard deviation!
                scal = np.nan

            dom_scale[i, j] = scal

    dic['power'] = powerTIR
    dic['scales'] = scales
    dic['dominant'] = dom_scale

    return dic
#
#
#
def LSTA_bothSigns(t, dx=None, dist=None, start=None, nb=None, dataset=None, dom=False):
    """
    :param t: data , numpy array
    :param dx: resolution of the data in unit of output (e.g. 3 if it's 3km)
    :param dist: distance between scale disaggregation levels, will be logarithmic
    :param start: smallest scale for scale decomposition
    :param nb: number of scales to decompose into
    :param dataset: abbreviation for existing dataset dictionary for quick use
    :param dom: boolean, true if dominant scale should additionally be returned in the dictionary
    :return: wavelet dictionary
    """

    def dom_get(wseries, tpoint, scales):

        values = (np.abs(wseries)) * (np.abs(wseries))
        values = values / (scales * scales)

        maxoutt = (values == ndimage.maximum_filter(values, 15, mode='reflect'))

        try:
            x = np.where(
                (maxoutt == 1))  # & (np.abs(parr) > scales ** .5)# ((wl >= np.percentile(wl[wl >= 0.5], 90)) &
        except IndexError:
            return False
        try:
            ind = (x[0])[0]
        except IndexError:
            return False

        scal = scales[ind]
        wavelet_pos_neg = np.real(wseries)[ind]

        if (wavelet_pos_neg < 0) : #| (tpoint < -1.)
            scal = scal * -1
        # if values[ind] < 0.02:
        #     scal = np.nan
        return scal

    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx, dist, nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    tir = t.copy()
    obj = wav.wavelet(dx, dist, nb, start=start)

    coeffsTIR_dry, powerTIR_dry = obj.calc_coeffs(tir, le_thresh=0, fill=0.01)

    tir = tir * -1
    coeffsTIR_wet, powerTIR_wet = obj.calc_coeffs(tir, le_thresh=0, fill=0.01)


    for id, s in enumerate(obj.scales):
        (powerTIR_wet[id, :, :])[powerTIR_wet[id, :, :] <= np.percentile(powerTIR_wet[id, :, :],50)] = 0
        (powerTIR_dry[id, :, :])[powerTIR_dry[id, :, :] <= np.percentile(powerTIR_dry[id, :, :],50)] = 0

    dic['power_dry'] = powerTIR_dry
    dic['power_wet'] = powerTIR_wet
    dic['scales'] = obj.scales

    if dom:
        dom_scale = np.zeros_like(tir) * np.nan
        for i in range(coeffsTIR_dry.shape[1]):
            for j in range(coeffsTIR_dry.shape[2]):
                parr_dry = coeffsTIR_dry[:, i, j]
                parr_wet = coeffsTIR_wet[:, i, j]
                tpoint = t[i,j]

                scal_dry = dom_get(parr_dry, tpoint, obj.scales)
                scal_wet = dom_get(parr_wet, tpoint, obj.scales)

                scal_pos = np.argmin([np.abs(scal_dry), np.abs(scal_wet)])
                scal_arr = scal_dry # [scal_dry,scal_wet][scal_pos]

                dom_scale[i, j] = scal_arr
        dic['dominant'] = dom_scale

    return dic
#

def waveletT1D(t, dx=None, dist=None,start=None, nb=None, dataset=None, mask=None):


    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx,dist,nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    #1D continuous wavelet analysis:
    #TIR
    tir=t.copy()
    tir[tir<0] = 0
    tir = (tir - np.mean(tir)) / np.std(tir)

    obj = wav1d.wavelet(dx, dist, nb, start=start)

    #TIR
    coeffsTIRx, powerTIRx = obj.calc_coeffs(tir, le_thresh=0, fill=0.0001, direction='x')
    coeffsTIRy, powerTIRy = obj.calc_coeffs(tir, le_thresh=0, fill=0.0001, direction='y')


    tir_shuff = tir.copy()
    np.random.shuffle(tir_shuff)

    coeffsSIGx, powerSIGx = obj.calc_coeffs(tir_shuff, le_thresh=0, fill=0.0001, direction='x')
    coeffsSIGy, powerSIGy = obj.calc_coeffs(tir_shuff, le_thresh=0, fill=0.0001, direction='y')

    sigx = obj.significance(powerTIRx, powerSIGx, direction='x', mask=mask)
    sigy = obj.significance(powerTIRy, powerSIGy, direction='y', mask=mask)


    dic['powerx']=powerTIRx
    dic['powery'] = powerTIRy
    dic['scales'] = obj.scales
    dic['res'] = obj.res
    dic['coeffsx'] = coeffsTIRx
    dic['coeffsy'] = coeffsTIRy
    dic['sigx'] = sigx
    dic['sigy'] = sigy

    return dic