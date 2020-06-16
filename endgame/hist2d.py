import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
import matplotlib.pylab as pylab
import cartopy.feature as cfeature

import scipy.stats as stats
import xarray as xr
import ipdb
import glob
import itertools
import numpy.ma as ma
from utils import u_statistics as u_stat
from scipy.stats import gaussian_kde
import pickle as pkl
from utils import u_plot as uplot
import pandas as pd
from scipy.stats import gaussian_kde, linregress
import matplotlib.cm as cm
from utils import u_met
import seaborn
import metpy
from metpy import calc
from metpy.units import units

from numpy.polynomial import polynomial as P

##for regridding, install xesmf:
# conda install esmpy
# pip install xesmf
from utils import constants as cnst


def create_2dhist(xvar, yvar, xbins, ybins, vardic, varpick):
    """

    :param xvar: xvar of the 2dhist
    :param yvar: yvar of the 2d hist
    :param xbins: bins to use for the xvar
    :param ybins: bins to use for the yvar
    :param varlist: dictionary of variables to put into histogram
    :param varpick: list of variables in dic to calculate
    :return:
    """
    outdic = {}

    for varstr in varpick:
        outdic[varstr] = np.zeros((len(ybins), len(xbins)))
        outdic[varstr+'_val'] = np.zeros((len(ybins), len(xbins)))
        # outdic[varstr + '_xbin'] = np.zeros((len(ybins), len(xbins)))
        # outdic[varstr + '_ybin'] = np.zeros((len(ybins), len(xbins)))

        calcvar = vardic[varstr]

        for isq, qql in enumerate(ybins[0:-1]):

            for issh, shl in enumerate(xbins[0:-1]):

                poss_ds = (xvar >= shl) & (xvar < xbins[issh + 1]) & (yvar >= qql) & (yvar < ybins[isq + 1])
                #ipdb.set_trace()
                try:
                    ds_mmean = np.nansum(calcvar[poss_ds])  # np.percentile(ds.tmin[poss_ds], 50)
                    ds_val = np.sum(np.isfinite(calcvar[poss_ds]))
                    if ds_val < 5:
                        ds_mean = np.nan
                    else:
                        ds_mean = ds_mmean / ds_val
                except IndexError:
                    ds_mean = np.nan


                (outdic[varstr])[isq, issh] = ds_mean
                (outdic[varstr+'_val'])[isq, issh] = np.sum(poss_ds)

                # (outdic[varstr + '_xbin'])[isq, issh] = shl
                # (outdic[varstr + '_ybin'])[isq, issh] = qql

    outdic['xbins'] = xbins
    outdic['ybins'] = ybins

    return outdic



def create_2dhist_centile(xvar, yvar, xbins, ybins, vardic, varpick, percentile=90, valmin=10):
    """

    :param xvar: xvar of the 2dhist
    :param yvar: yvar of the 2d hist
    :param xbins: bins to use for the xvar
    :param ybins: bins to use for the yvar
    :param varlist: dictionary of variables to put into histogram
    :param varpick: list of variables in dic to calculate
    :return:
    """
    outdic = {}

    for varstr in varpick:
        outdic[varstr] = np.zeros((len(ybins), len(xbins)))
        outdic[varstr+'_val'] = np.zeros((len(ybins), len(xbins)))
        # outdic[varstr + '_xbin'] = np.zeros((len(ybins), len(xbins)))
        # outdic[varstr + '_ybin'] = np.zeros((len(ybins), len(xbins)))

        calcvar = vardic[varstr]

        for isq, qql in enumerate(ybins[0:-1]):

            for issh, shl in enumerate(xbins[0:-1]):

                poss_ds = (xvar >= shl) & (xvar < xbins[issh + 1]) & (yvar >= qql) & (yvar < ybins[isq + 1])

                try:
                    ds_mmean = np.percentile(calcvar[poss_ds], percentile) #np.nansum(calcvar[poss_ds])  # np.percentile(ds.tmin[poss_ds], 50)
                    ds_val = np.sum(np.isfinite(calcvar[poss_ds]))
                    if ds_val < valmin:
                        ds_mean = np.nan
                    else:
                        ds_mean = ds_mmean #/ ds_val
                except IndexError:
                    ds_mean = np.nan


                (outdic[varstr])[isq, issh] = ds_mean
                (outdic[varstr+'_val'])[isq, issh] = np.sum(poss_ds)

                # (outdic[varstr + '_xbin'])[isq, issh] = shl
                # (outdic[varstr + '_ybin'])[isq, issh] = qql

    outdic['xbins'] = xbins
    outdic['ybins'] = ybins

    return outdic