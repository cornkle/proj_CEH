import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata
import ipdb
import scipy.spatial.qhull as qhull
import xarray as xr



"""Define a grid with salem
lon: list or numpy array of longitudes. Has to contain at least min max longitudes of the grid
lat: list or numpy array of latitudes. Has to contain at least min max latitudes of the grid
proj: projection from pyproj to create the grid. If not set default is mercator:
      proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')
dx: number of meters per pixel, must be integer, e.g. 5000 for 5km pixels
"""


def make(lon, lat, dx, keep_ll=False):
    if lon.ndim == 1:
        grid_lons, grid_lats = np.meshgrid(lon, lat)
    else:
        grid_lons = lon
        grid_lats = lat
    if not keep_ll:
        # Transform lon, lats to the mercator projection
        proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')
        x, y = pyproj.transform(salem.wgs84, proj, grid_lons, grid_lats)
    else:
        proj = salem.wgs84
        x, y = grid_lons, grid_lats

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


def grid_to_latlon(grid):

    coords = grid.ll_coordinates
    return ((coords[0])[0,:], (coords[1])[:,0])


def refactor_da(da, factor):
    grid = da.salem.grid.regrid(factor=factor)
    data = grid.lookup_transform(da)

    times = ['year', 'month', 'time', 'date']
    for vn in times:
        if vn in da.coords:
            time = da[vn]

    lon, lat = grid_to_latlon(grid)


    return xr.DataArray(data, coords=[time, lat, lon],
                          dims=['time', 'latitude', 'longitude'])





