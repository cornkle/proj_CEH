import pdb
import scipy.spatial.qhull as qhull
import numpy as np
from scipy.interpolate import griddata
import pyproj


def _interp_weights(xyz, uvw, d=None):

    """
    :param xyz: flattened coords of current grid
    :param uvw: flattened coords of target grid
    :param d: number of dimensions of new grid
    :return: triangulisation lookup table, point weights
    """
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def _interpolate(values, vtx, wts, fill_value=np.nan):

    """
    :param values: flattened data values
    :param vtx: lookup table
    :param wts: point weights
    :param fill_value: fill value for extrapolated regions
    :return: interpolated data
    """
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = fill_value
    return ret

def interpolation_weights(x, y, new_x, new_y, irregular_1d=False):

    """
    :param x: current x variables (1 or 2d)
    :param y: current y variables (1 or 2d)
    :param new_x: target x vars
    :param new_y: target y vars
    :keyword irregular_1d = False , set to True if input is non-ordered points in 1d arrays
    :return:  triangulisation lookup table, point weights, 2d shape - inputs for interpolation func
    """

    if (x.ndim == 1) & ~irregular_1d:
        grid_xs, grid_ys = np.meshgrid(x, y)
    else:
        grid_xs = x
        grid_ys = y

    if new_x.ndim == 1:
        new_xs, new_ys = np.meshgrid(new_x, new_y)
    else:
        new_xs = new_x
        new_ys = new_y

    points = np.array((grid_xs.flatten(), grid_ys.flatten())).T
    inter = np.array((np.ravel(new_xs), np.ravel(new_ys))).T

    inds, weights = _interp_weights(points, inter, d=2)

    return inds, weights, new_xs.shape


def interpolate_data(data, inds, weights, shape):

    """
    This routine interpolates only over the 2d plane i.e. spatial interpolation
    :param data: original data, 2d, 3d or 4d (e.g. incl. time steps and pressure levels).

    :param inds: lookup table from weights func
    :param weights: index weights from weights func
    :param shape: 2d shape of plane
    :return: interpolated data, same number of dimensions as input data
    """

    if (data.ndim < 2) | (data.ndim > 4):
        print('Error. Only data with 2 - 4 dimensions allowed.')
        return
    # interpolate 2d arrays
    coll = []
    if data.ndim > 2:
        for d in data:
            if d.ndim == 2:

                d2d = _interpolate(d.flatten(), inds, weights)
                d2d = d2d.reshape(shape)
                coll.append(d2d[None, ...])

            if d.ndim == 3:
                plevs = []

                for pl in d:
                    pdb.set_trace()
                    pl2d = _interpolate(pl.flatten(), inds, weights)
                    pl2d = pl2d.reshape(shape)
                    plevs.append(pl2d[None, ...])
                if len(plevs) > 1:
                    plevs = np.concatenate(plevs, axis=0)
                coll.append(plevs[None, ...])
        if len(coll) > 1:
            coll = np.concatenate(coll, axis=0)
    else:
        d2d = _interpolate(data.flatten(), inds, weights)
        d2d = d2d.reshape(shape)
        coll = d2d


    return coll


def regrid_irregular_quick(x, y, new_x, new_y, data):

    """
    Combines all steps of data interpolation, does not provide weights etc
    Useful for quick interpolation of single 3d to 4d arrays
    :param x: array, current x coordinates
    :param y: array, current y coordinates
    :param new_x: array, new x coordinates
    :param new_y: array, new y coordinates
    :param data: array, input data
    :return:
    """
    inds, weights, shape = interpolation_weights(x, y, new_x, new_y)
    # interpolate 2d arrays
    coll = interpolate_data(data, inds, weights, shape)

    return coll


