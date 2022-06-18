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

        varstrr = varstr
        if varstr == 'pall':
            varstrr = 'prcp'

        outdic[varstr] = np.zeros((len(ybins), len(xbins)))
        outdic[varstr+'_val'] = np.zeros((len(ybins), len(xbins)))
        # outdic[varstr + '_xbin'] = np.zeros((len(ybins), len(xbins)))
        # outdic[varstr + '_ybin'] = np.zeros((len(ybins), len(xbins)))

        calcvar = vardic[varstr]

        try:
            calcvar = np.concatenate(np.array(calcvar), axis=0)
        except:
            pass

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


                (outdic[varstrr])[isq, issh] = ds_mean
                (outdic[varstrr+'_val'])[isq, issh] = np.sum(poss_ds)

                # (outdic[varstr + '_xbin'])[isq, issh] = shl
                # (outdic[varstr + '_ybin'])[isq, issh] = qql

    outdic['xbins'] = xbins
    outdic['ybins'] = ybins

    return outdic



def create_2dhist_centile(xvar, yvar, xbins, ybins, vardic, varpick, percentile=95, valmin=10):
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

        varstrr = varstr
        if varstr == 'pall':
            varstrr = 'prcp'

        outdic[varstrr] = np.zeros((len(ybins), len(xbins)))
        outdic[varstrr+'_val'] = np.zeros((len(ybins), len(xbins)))
        # outdic[varstr + '_xbin'] = np.zeros((len(ybins), len(xbins)))
        # outdic[varstr + '_ybin'] = np.zeros((len(ybins), len(xbins)))

        calcvar = vardic[varstr]


        for isq, qql in enumerate(ybins[0:-1]):

            for issh, shl in enumerate(xbins[0:-1]):

                # if issh == 0:
                #     poss_ds =  (xvar < xbins[issh + 1]) & (yvar >= qql) & (yvar < ybins[isq + 1])
                #
                # elif issh == len(xbins[0:-1])-1:
                #     poss_ds =  (xvar >= shl) & (yvar >= qql) & (yvar < ybins[isq + 1])
                # elif isq == 0:
                #     poss_ds = (xvar >= shl) & (xvar < xbins[issh + 1]) & (yvar < ybins[isq + 1])
                # elif isq == len(ybins[0:-1])-1:
                #     poss_ds = (xvar >= shl) & (xvar < xbins[issh + 1]) & (yvar >= qql)
                # else:
                poss_ds = (xvar >= shl) & (xvar < xbins[issh + 1]) & (yvar >= qql) & (yvar < ybins[isq + 1])

                try:
                    isdata = (calcvar[poss_ds])
                    if varstr == 'pall':
                        try:
                            isdata = np.concatenate(np.array(isdata), axis=0)
                            #isdata = isdata[isdata>1]
                        except:
                            pass

                    if len(isdata) == 0:
                        ds_mean = np.nan
                    else:

                        try:
                            isdata = isdata[np.isfinite(isdata)]
                        except:
                            ipdb.set_trace()

                        ds_mmean = np.percentile(isdata, percentile) #np.nansum(calcvar[poss_ds])  # np.percentile(ds.tmin[poss_ds], 50)

                        ds_val = np.sum(np.isfinite(isdata))
                        if ds_val < valmin:
                            ds_mean = np.nan
                        else:
                            ds_mean = ds_mmean #/ ds_val
                        #print(np.sum(poss_ds), ds_val)
                except IndexError:
                    ds_mean = np.nan


                (outdic[varstrr])[isq, issh] = ds_mean
                (outdic[varstrr+'_val'])[isq, issh] = np.sum(poss_ds)

                # (outdic[varstr + '_xbin'])[isq, issh] = shl
                # (outdic[varstr + '_ybin'])[isq, issh] = qql

    outdic['xbins'] = xbins
    outdic['ybins'] = ybins

    return outdic


def create_2dhist_maxYear(xvar, yvar, xbins, ybins, vardic, varpick, valmin=2):
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

        varstrr = varstr

        outdic[varstrr] = np.zeros((len(ybins)-1, len(xbins)-1))*np.nan
        outdic[varstrr+'_val'] = np.zeros((len(ybins)-1, len(xbins)-1))*np.nan
        outdic[varstrr + '_std'] = np.zeros((len(ybins) - 1, len(xbins) - 1))*np.nan

        calcvar = vardic#[varstr]


        for isq, qql in enumerate(ybins[0:-1]):

            for issh, shl in enumerate(xbins[0:-1]):


                poss_ds = (xvar >= shl) & (xvar < xbins[issh + 1]) & (yvar >= qql) & (yvar < ybins[isq + 1])

                try:
                    isdata = (calcvar[poss_ds].groupby('year').max())
                    valdata = (calcvar[poss_ds].groupby('year').count())
                    ds_mmean = isdata[varstr]
                    ds_mval = valdata[varstr]

                    nb_mcs = 2 # number mcs per year minimum to calc year max from, other years are excluded
                    ds_val = np.sum(ds_mval>=nb_mcs)  # number of years where enough MCSs exist

                    if ds_val < valmin:
                        ds_mean = np.nan
                        ds_std = np.nan
                    else:
                        #ipdb.set_trace()
                        ds_mean = np.nanmean(ds_mmean[ds_mval>=nb_mcs])
                        ds_std = np.nanstd(ds_mmean[ds_mval >= nb_mcs])

                except IndexError:
                    ds_mean = np.nan
                    ds_val = np.nan
                    ds_std = np.nan


                (outdic[varstrr])[isq, issh] = ds_mean
                (outdic[varstrr+'_val'])[isq, issh] = ds_val
                (outdic[varstrr + '_std'])[isq, issh] = ds_std

    outdic['xbins'] = xbins
    outdic['ybins'] = ybins

    return outdic



def basic_1d_binning(xvar, xbins):
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
    outdic['nb'] = []


    for issh, shl in enumerate(xbins[0:-1]):

        poss_ds = np.sum((xvar > shl) & (xvar <= xbins[issh + 1]))

        outdic['nb'].append(poss_ds)



    outdic['xbins'] = (np.round(xbins[0:-1]+((xbins[1::]-xbins[0:-1])/2),2))


    return outdic


def basic_1d_binning_mean(xvar, xbins):
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
    outdic['nb'] = []


    for issh, shl in enumerate(xbins[0:-1]):

        mask = (xvar > shl) & (xvar <= xbins[issh + 1])

        outdic['nb'].append(np.nanmean(xvar[mask]))



    outdic['xbins'] = (np.round(xbins[0:-1]+((xbins[1::]-xbins[0:-1])/2),2))


    return outdic


def var2_binning_mean(xvar, yvar, xbins):
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
    outdic['y'] = []
    outdic['ycount'] = []
    outdic['ystd'] = []


    for issh, shl in enumerate(xbins[0:-1]):

        mask = (xvar > shl) & (xvar <= xbins[issh + 1])

        #outdic['y'].append(np.nanmean(yvar[mask]))
        #print(shl)

        try:
            #outdic['y'].append(np.percentile(yvar[mask],95))
            #outdic['y'].append(np.median(yvar[mask]))
            #ipdb.set_trace()
            outdic['y'].append(np.mean(yvar[mask]))
            outdic['ystd'].append(np.std(yvar[mask]))
            outdic['ycount'].append(np.sum(np.isfinite(yvar[mask])))
            #ipdb.set_trace()

        except:
            # outdic['y'].append(np.nan)
            # outdic['ystd'].append(np.nan)
            # outdic['ycount'].append(np.nan)

            outdic['y'].append(0)
            outdic['ystd'].append(0)
            outdic['ycount'].append(0)

        #ipdb.set_trace()

    outdic['xbins'] = (np.round(xbins[0:-1]+((xbins[1::]-xbins[0:-1])/2),2))


    return outdic


def var2_binning_percentile(xvar, yvar, xbins, p):
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
    outdic['y'] = []


    for issh, shl in enumerate(xbins[0:-1]):

        mask = (xvar > shl) & (xvar <= xbins[issh + 1]) & np.isfinite(yvar)

        outdic['y'].append(np.percentile(yvar[mask], p))

    outdic['xbins'] = (np.round(xbins[0:-1]+((xbins[1::]-xbins[0:-1])/2),2))


    return outdic




def perc_1d_binning(data, xvar, xbins, perc):
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
    outdic['data'] = []


    for issh, shl in enumerate(xbins[0:-1]):

        poss_ds = (xvar > shl) & (xvar <= xbins[issh + 1])
        valdat = (data[poss_ds])[np.isfinite(data[poss_ds])]
        dat = np.percentile(valdat, perc)

        outdic['data'].append(dat)

        #ipdb.set_trace()

    outdic['xbins'] = (np.round(xbins[0:-1]+((xbins[1::]-xbins[0:-1])/2),2))


    return outdic


def var2_binning_threshold(xvar, yvar, xbins, gt=None, lt=None):
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
    outdic['y'] = []


    for issh, shl in enumerate(xbins[0:-1]):

        mask = (xvar > shl) & (xvar <= xbins[issh + 1]) & np.isfinite(yvar)

        if (lt!=None) & (gt!=None):
            outdic['y'].append(np.sum((yvar[mask] < lt) | (yvar[mask] > gt)) / float(np.sum(np.isfinite(yvar[mask]))))
        elif lt != None:
            outdic['y'].append(np.sum(yvar[mask]<lt)/float(np.sum(np.isfinite(yvar[mask]))))
        else:
            outdic['y'].append(np.sum(yvar[mask]>gt)/float(np.sum(np.isfinite(yvar[mask]))))

    outdic['xbins'] = (np.round(xbins[0:-1]+((xbins[1::]-xbins[0:-1])/2),2))


    return outdic


def var2_binning_threshold_CDF(xvar, yvar, xbins, gt=None):
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
    outdic['y'] = []


    for issh, shl in enumerate(xbins):

        mask = (xvar <= shl)
        outdic['y'].append(np.sum(yvar[mask]>=gt)/float(np.sum(yvar>=gt)))


    outdic['xbins'] = xbins


    return outdic