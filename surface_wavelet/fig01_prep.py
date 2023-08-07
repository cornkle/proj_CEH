# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import ipdb
import pandas as pd
from collections import OrderedDict
import salem
from utils import u_met, u_parallelise, u_arrays as ua, constants as cnst, u_darrays
from scipy.interpolate import griddata
import multiprocessing
import glob
from scipy import ndimage
from utils import u_statistics as ustats
import salem
from metpy import calc
from metpy.units import units


import pickle as pkl


matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)


h = 17
eh = -5

dic = {}
dic2 = {}
dic3={}
dic4 = {}


name2='ERA5_cores_2hOverlap_AMSRE_SMALL_MCSfilter_minusMean_smallDomain_'
name3 = 'ERA5_cores_2hOverlap_AMSRE_SMALL_AMSL_MCSfilter'
name="ERA5_cores_2hOverlap_AMSRE_SMALL_MCSfilter_minusMean_bigDomain_init400km_"
name4 ='ERA5_cores_NEWTRACKING_AMSRE_ALL222_'


def coll(dic2, h, eh, year, name):
    print(h)
    core = pkl.load(open(
        cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/core_txt/" + name + str(eh) + "UTCERA" + str(h).zfill(
            2) + '_' + str(year) + ".p", "rb"))

    # core = pkl.load(open(
    #     cnst.network_data + "figs/LSTA/corrected_LSTA/new/ERA5/" + name + str(eh) + "UTCERA" + str(h).zfill(
    #         2) + '_' + str(year) + "_small_cores.p", "rb"))

    for id, k in enumerate(core.keys()):
        try:
            dic2[k] = dic2[k] + core[k]
        except KeyError:
            dic2[k] = core[k]

for y in range(2006, 2010):
    coll(dic, h, eh, y, name)

for y in range(2006, 2010):
    coll(dic2, h, eh, y, name2)

for y in range(2006, 2010):
    coll(dic3, h, eh, y, name3)

for y in range(2006, 2010):
    coll(dic4, h, eh, y, name4)


outdic = dic
outdic['u925'] = dic2['u925'] / dic2['cnte']
outdic['v925'] = dic2['v925'] / dic2['cnte']

amsr = ndimage.gaussian_filter(dic2['lsta0'] / dic2['cnt0'], 2, mode='nearest')
outdic['lsta0'] = amsr
vals = np.percentile(amsr,[8,92])
mask = (amsr <= vals[0]) | (amsr >= vals[1])
outdic['lsta0_sig'] = mask

outdic['shear'] = dic3['shear'] / dic3['cnte']

div = ndimage.gaussian_filter(((dic['div'])/ dic['cnte'])*100, 3, mode='nearest')
vals = np.percentile(div,[9,91])
dmask =  (div <= vals[0])
outdic['div'] = div
outdic['div_sig'] = dmask

outdic['v925_orig'] = dic3['v925_orig'] / dic3['cnte']
outdic['itd'] = dic4['itd'] / dic4['cnte']
outdic['t'] = ndimage.gaussian_filter((dic['t']-dic['tclim']) / dic['cnte'], 6, mode='nearest')
qq = ndimage.gaussian_filter(((dic['q']-dic['qclim'])*1000/ dic['cnte']), 3, mode='nearest')
qqmask = np.percentile(qq,[9,95])
qmask =  qq>=vals[1]
outdic['q'] = ndimage.gaussian_filter(((dic['q']-dic['qclim'])*1000/ dic['cnte']), 3, mode='nearest')
outdic['q_sig'] = qmask



pkl.dump(outdic, open(cnst.network_data + "figs/LSTA/paper/saves/fig01.p", "wb"))
