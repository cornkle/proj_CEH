import numpy as np
import xarray as xr
from utils import u_mann_kendall as mk, u_arrays as ua
from scipy import stats
import matplotlib.pyplot as plt
import bottleneck
import scipy
import ipdb


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
        if np.nansum(x.values==0)>=nb_missing:
            p = np.nan
            slope = np.nan
            ind = 0

    ds = xr.Dataset()
    ds['slope'] = xr.DataArray(slope,)
    ds['pval'] = xr.DataArray(p, )
    ds['ind'] = xr.DataArray(ind)

    return ds

def linear_trend_lingress(x, nb_missing=None):

    if np.isnan(x).all():
        ds = xr.Dataset()
        ds['slope'] = xr.DataArray(np.nan, )
        ds['pval'] = xr.DataArray(np.nan, )
        ds['r'] = xr.DataArray(np.nan, )
        return ds

    slope, intercept, r, p, std_err = stats.linregress(np.arange(len(x)), x)

    # we need to return a dataarray or else xarray's groupby won't be happy

    if nb_missing is not None:
        if np.nansum(x.values==0)>=nb_missing:
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
    try:
        ds = ds.sel(latitude=slice(None, None, -1))
    except ValueError:
        try:
            ds = ds.sel(lat=slice(None, None, -1))
        except ValueError:
            ds = ds.sel(Latitude=slice(None, None, -1))
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


def cut_box(arr, xpos=None, ypos=None, dist=None):
    """

    :param xpos: x coordinate in domain for kernel centre point
    :param ypos: y coordinate in domain for kernel centre point
    :param arr: numpy array (2d)
    :param dist: distance from kernel centre point to kernel edge (total width = 2*dist+1)
    :return: the kernel of dimensions (2*dist+1, 2*dist+1)
    """

    if dist == None:
        'Distance missing. Please provide distance from kernel centre to edge (number of pixels).'
        return

    if arr.ndim == 2:
        kernel = ua.cut_kernel(arr.values,xpos, ypos,dist)
        if kernel.shape != (dist * 2 + 1, dist * 2 + 1):
            print("Please check kernel dimensions, there is something wrong")
            ipdb.set_trace()
    elif arr.ndim == 3:
        kernel = ua.cut_kernel_3d(arr.values,xpos, ypos,dist)

        if kernel.shape != (arr.shape[0], dist * 2 + 1, dist * 2 + 1):
            print("Please check kernel dimensions, there is something wrong")
            ipdb.set_trace()
    else:
        print('Dimension of array not supported, please check')
        ipdb.set_trace()

    if arr.ndim == 3:
        try:
            levels = arr.level.values
        except AttributeError:
            levels = arr.pressure.values

        return xr.DataArray(kernel, dims=['level','y','x'],
                            coords={'level' : arr.level.values})
    else:
        return xr.DataArray(kernel, dims=['y','x'])
