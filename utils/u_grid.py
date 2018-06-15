import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
import pdb
import math
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import itertools

proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')

"""Define a grid with salem
lon: list or numpy array of longitudes. Has to contain at least min max longitudes of the grid
lat: list or numpy array of latitudes. Has to contain at least min max latitudes of the grid
proj: projection from pyproj to create the grid. If not set default is mercator:
      proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')
dx: number of meters per pixel, must be integer, e.g. 5000 for 5km pixels
"""

def _interp_weights(xyz, uvw, d=None):
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def _interpolate(values, vtx, wts, fill_value=np.nan):
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = fill_value
    return ret


def make(lon, lat, dx, proj=proj):
    if lon.ndim == 1:
        grid_lons, grid_lats = np.meshgrid(lon, lat)
    else:
        grid_lons = lon
        grid_lats = lat
    # Transform lon, lats to the mercator projection
    x, y = pyproj.transform(salem.wgs84, proj, grid_lons, grid_lats)
    # take the min and max
    xmax, xmin = np.max(x), np.min(x)
    ymax, ymin = np.max(y), np.min(y)
    # Count the number of pixels
    nx, r = divmod(xmax - xmin, dx)
    ny, r = divmod(ymax - ymin, dx)
    # Here one could add + 1 to be sure that the last pixel is always included
    grid = salem.Grid(nxny=(nx, ny), dxdy=(dx, dx), ll_corner=(xmin, ymin), proj=proj)
    return grid


"""Transforms latitudes and longitudes to salem grid coordinates and creates griddata input
Returns:
points: An (y,x) tuple of lons, lats transformed to floating coordinates in grid coordinate system
inter: An (y,x) tuple of the grid integer coordinates

Can be used as input for the scipy.interpolate.griddata function with
inter= points at which to interpolate data (xi)
points= data point coordinates

"""


def griddata_input(lon, lat, grid):
    xi, yi = grid.ij_coordinates

    # Transform lons, lats to grid
    x, y = grid.transform(lon.flatten(), lat.flatten(), crs=salem.wgs84)
    # Convert for griddata input
    points = np.array((y, x)).T
    inter = np.array((np.ravel(yi), np.ravel(xi))).T

    return inter, points



"""
Quick linear regrid function for irregular data that takes lon, lat, data of the same dimensions and interpolates
the data to the salem grid

lon: numpy array of longitudes
lat: numpy array of latitudes
data: numpy array of data values
grid: salem grid

returns: the regridded data, linearly interpolated.
"""


def quick_regrid(lon, lat, data, grid):
    if lon.ndim == 1:
        grid_lons, grid_lats = np.meshgrid(lon, lat)
    else:
        grid_lons = lon
        grid_lats = lat

    inter, points = griddata_input(grid_lons, grid_lats, grid)

    # Interpolate using delaunay triangularization
    data = griddata(points, data.flatten(), inter, method='linear')
    data = data.reshape((grid.ny, grid.nx))

    return data

def regrid_irregular(x, y, new_x, new_y,  data):
    if x.ndim == 1:
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

    # Interpolate using delaunay triangularization
    coll = []
    for d in data:

        d2d = griddata(points, d.flatten(), inter, method='linear')

        d2d = d2d.reshape((new_xs.shape[0], new_xs.shape[1]))
        coll.append(d2d[None,...])
    if len(coll)>1:
        coll = np.concatenate(coll, axis=0)

    return coll

def regrid_irregular_quick(x, y, new_x, new_y,  data):

    if data.ndim < 2:
        print('Error. Data has less than 2 dimensions.')
        return

    if x.ndim == 1:
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
    # delaunay triangularization
    vtx, wts = _interp_weights(points, inter, d=2)
    # interpolate 2d arrays
    coll = []
    for d in data:
        if d.ndim == 2:
            d2d = _interpolate(d.flatten(), vtx, wts)
            d2d = d2d.reshape((new_xs.shape))
            coll.append(d2d[None,...])
        if d.ndim == 3:
            plevs = []
            for pl in d:
                pl2d = _interpolate(pl.flatten(), vtx, wts)
                pl2d = pl2d.reshape((new_xs.shape))
                plevs.append(pl2d[None, ...])
            if len(plevs) > 1:
                plevs = np.concatenate(plevs, axis=0)
            coll.append(plevs[None, ...])

    if len(coll)>1:
        coll = np.concatenate(coll, axis=0)

    return coll

def interpolation_weights(x, y, new_x, new_y):

    if x.ndim == 1:
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

    return inds, weights


def interpolate_data(data, inds, weights):

    if data.ndim < 2 | data.ndim > 4:
        print('Error. Only data with 2 - 4 dimensions allowed.')
        return

    # interpolate 2d arrays
    coll = []
    for d in data:
        if d.ndim == 2:
            d2d = _interpolate(d.flatten(), inds, weights)
            d2d = d2d.reshape((data.shape[-2::]))
            coll.append(d2d[None, ...])
        if d.ndim == 3:
            plevs = []
            for pl in d:
                pl2d = _interpolate(pl.flatten(), inds, weights)
                pl2d = pl2d.reshape((data.shape[-2::]))
                plevs.append(pl2d[None, ...])
            if len(plevs) > 1:
                plevs = np.concatenate(plevs, axis=0)
            coll.append(plevs[None, ...])

    if len(coll) > 1:
        coll = np.concatenate(coll, axis=0)

    return coll




