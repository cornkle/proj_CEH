# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 10:15:40 2016

@author: cornkle
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
import pdb
import pandas as pd
from utils import u_met, u_parallelise, u_gis, u_arrays, constants, u_grid
import salem
from scipy.interpolate import griddata

import pickle as pkl


file = constants.MCS_15K
file2 = constants.MCS_CENTRE70

mcs = xr.open_dataarray(file)
centre = xr.open_dataarray(file2)

for id, time in enumerate(mcs.time.values):

    pos = np.where(centre.time == time)
    if len(pos[0]) > 0:
        image = centre[pos[0],:,:].squeeze()
        image2 = mcs[id, :, :].squeeze()
        pos2 = np.where(image == 2)
        if len(pos2[0]) > 0:

            print(image2.values[pos2])