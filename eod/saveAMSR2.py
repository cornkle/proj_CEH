# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import pdb
import os
import xarray as xr
from utils import constants as cnst
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import ipdb

def saveNetcdf(day=True):

    sm_folder = cnst.network_data + 'data/OBS/AMSRE/aqua/raw_day_AMSR2'
    pool = multiprocessing.Pool(processes=7)

    if day:
        daystring = '_A_'
    else:
        daystring = '_D_'

    files = glob.glob(sm_folder+'/LPRM-AMSR*'+daystring+'*.nc4')
    print('start loop')
    ipdb.set_trace()
    for f in files:
        ds = rewrite_data.rewrite_AMSR2(f)

    #res = pool.map(rewrite_data.rewrite_AMSR2, files)
