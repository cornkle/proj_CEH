import numpy as np
import xarray as xr
from utils import u_mann_kendall as mk
from scipy import stats
import matplotlib.pyplot as plt
import bottleneck
import scipy
import pdb


def shift_lons(ds, lon_dim='lon'):
    """ Shift longitudes from [0, 360] to [-180, 180] """
    lons = ds[lon_dim].values
    new_lons = np.empty_like(lons)
    mask = lons > 180
    new_lons[mask] = -(360. - lons[mask])
    new_lons[~mask] = lons[~mask]
    ds[lon_dim].values = new_lons
    return ds

def linear_trend_mk(x, eps=0.001, alpha=0.01, nb_missing=None):

    #pf = np.polyfit(np.arange(len(x)), x, 1)
    pf, slope, int, p, ind = mk.test(np.arange(len(x)),x.squeeze().values, eps=eps, alpha=alpha, Ha='upordown')

    # we need to return a dataarray or else xarray's groupby won't be happy


    if nb_missing is not None:
        if np.nansum(x.values==0)>=10:
            p = np.nan
            slope = np.nan

    ds = xr.Dataset()
    ds['slope'] = xr.DataArray(slope,)
    ds['pval'] = xr.DataArray(p, )
    ds['ind'] = xr.DataArray(ind)

    return ds

def linear_trend_lingress(x, nb_missing=None):

    slope, intercept, r, p, std_err = stats.linregress(np.arange(len(x)), x)

    # we need to return a dataarray or else xarray's groupby won't be happy

    if nb_missing is not None:
        if np.nansum(x.values==0)>=10:
            p = np.nan
            slope = np.nan

    ds = xr.Dataset()
    ds['slope'] = xr.DataArray(slope,)
    ds['pval'] = xr.DataArray(p, )
    ds['r'] = xr.DataArray(r,)

    return ds


def flip_lat(ds):
    """
    Flips latitude of dataset. Only works with correct latitude name...
    :param ds:
    :return:
    """
    ds = ds.sel(latitude=slice(None, None, -1))
    return ds



################### correlation computation with quick parallel dask usage
def covariance_gufunc(x, y):
    return ((x - x.mean(axis=-1, keepdims=True))
            * (y - y.mean(axis=-1, keepdims=True))).mean(axis=-1)

def pearson_correlation_gufunc(x, y):
    return covariance_gufunc(x, y) / (x.std(axis=-1) * y.std(axis=-1))

def spearman_correlation_gufunc(x, y):
    x_ranks = bottleneck.rankdata(x, axis=-1)
    y_ranks = bottleneck.rankdata(y, axis=-1)
    return pearson_correlation_gufunc(x_ranks, y_ranks)

def spearman_correlation(x, y, dim):
    return xr.apply_ufunc(
        spearman_correlation_gufunc, x, y,
        input_core_dims=[[dim], [dim]],
        dask='parallelized',
        output_dtypes=[float])

def pearson_correlation(x, y, dim):
    return xr.apply_ufunc(
        pearson_correlation_gufunc, x, y,
        input_core_dims=[[dim], [dim]],
        dask='parallelized',
        output_dtypes=[float])


