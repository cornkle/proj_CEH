import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata

"""Define a grid with salem
lon: list or numpy array of longitudes. Has to contain at least min max longitudes of the grid
lat: list or numpy array of latitudes. Has to contain at least min max latitudes of the grid
proj: projection from pyproj to create the grid for example:
      proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')
dx: number of meters per pixel, must be integer, e.g. 5000 for 5km pixels
"""


def make(lon, lat, proj, dx):
    # Transform lon, lats to the mercator projection
    x, y = pyproj.transform(salem.wgs84, proj, lon, lat)
    # take the min and max
    xmax, xmin = np.max(x), np.min(x)
    ymax, ymin = np.max(y), np.min(y)
    # Count the number of pixels
    nx, r = divmod(xmax - xmin, dx)
    ny, r = divmod(ymax - ymin, dx)
    # Here one could add + 1 to be sure that the last pixel is always included
    grid = salem.Grid(nxny=(nx, ny), dxdy=(dx, dx), ll_corner=(xmin, ymin), proj=proj)
    return grid


def griddata_input(lon, lat, grid):
    xi, yi = grid.ij_coordinates

    # Transform lons, lats to grid
    x, y = grid.transform(lon.flatten(), lat.flatten(), crs=salem.wgs84)
    # Convert for griddata input
    points = np.array((y, x)).T
    inter = np.array((np.ravel(yi), np.ravel(xi))).T

    return inter, points
