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

def saveNetcdf():

    path = cnst.lmcs_drive + 'ASCAT/raw/'
    outpath = cnst.lmcs_drive + 'ASCAT/raw/'
    for y in range(2007,2013):
        if not os.path.isdir(outpath + str(y)):
            os.mkdir(path+str(y))

        for ff in glob.glob(path+str(y)+'/*.gra'):
                rewrite_data.rewrite_ASCAT(ff)
