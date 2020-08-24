# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 14:08:18 2016

@author: cornkle
"""
import numpy as np
from saveCore_standalone_NFLICS import wav


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
    'METEOSAT8K': _create_dic(8, 1 / 12., 24, 40),
    'METEOSAT10K': _create_dic(10, 1 / 12., 30, 40),
    'GRIDSAT': _create_dic(8, 1 / 12., 24, 40),
}



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



def applyHat(t, dx=None, dist=None,start=None, nb=None, dataset=None):

    dic = {}

    if dataset in DATASETS:
        dx, dist, start, nb = read_dic(DATASETS[dataset])

    if not np.array([dx, dist, nb]).all():
        print('Information missing. Please provide either dataset or dx, dist and nb explicitly.')
        return

    tir = t.copy()

    obj = wav.wavelet(dx, dist, nb, start=start)
    print('Scales: ', obj.scales)
    # TIR
    coeffsTIR, powerTIR = obj.calc_coeffs(tir, ge_thresh=0, fill=0.01)

    dic['power'] = powerTIR
    dic['scales'] = obj.scales
    dic['res'] = obj.res
    dic['coeffs'] = coeffsTIR

    return dic