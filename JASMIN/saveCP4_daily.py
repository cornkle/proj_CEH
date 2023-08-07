# -*- coding: utf-8 -*-


import multiprocessing
import glob
from eod import rewrite_data
import ipdb
import os
import xarray as xr
import numpy as np
import pandas as pd
from utils import constants as cnst, u_met
from scipy.interpolate import griddata


def saveDaily():

    files = glob.glob('/prj/AMMA2050/CP4/historical/25km/precip/*.nc')

    for svar in files:

        out = svar.replace('A1hr_mean', 'Aday_mean')
        out = out.replace('precip', 'precip_day')

        ds = xr.open_dataset(svar)
        ds = ds.groupby('time.day').sum('time')

        ds.to_netcdf(out)