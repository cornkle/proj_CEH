import salem
import pyproj
import numpy as np
from scipy.interpolate import griddata

"""Define a grid with salem
lon: list or numpy array of longitudes. Has to contain at least min max longitudes of the grid
lat: list or numpy array of latitudes. Has to contain at least min max latitudes of the grid
proj: projection from pyproj to create the grid for example:
      proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')
dx: number of kilometer per pixel, must be integer
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



