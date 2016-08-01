# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 14:39:16 2016

@author: cornkle
"""
import salem
import pyproj
import numpy as np
import math
from scipy.interpolate import griddata

proj = pyproj.Proj('+proj=merc +lat_0=0. +lon_0=0.')

def creategrid(min_lon, max_lon, min_lat, max_lat, cell_size_deg, mesh=False):
    #’’’Output grid within geobounds and specifice cell size
    #cell_size_deg should be in decimal degrees’’’

    min_lon = math.floor(min_lon)
    max_lon = math.ceil(max_lon)
    min_lat = math.floor(min_lat)
    max_lat = math.ceil(max_lat)
    
    lon_num = (max_lon - min_lon)/cell_size_deg
    lat_num = (max_lat - min_lat)/cell_size_deg
    
    grid_lons = np.zeros(lon_num) # fill with lon_min
    grid_lats = np.zeros(lat_num) # fill with lon_max
    grid_lons = grid_lons + (np.asarray(range(lon_num))*cell_size_deg)
    grid_lats = grid_lats + (np.asarray(range(lat_num))*cell_size_deg)
    
    grid_lons, grid_lats = np.meshgrid(grid_lons, grid_lats)
    grid_lons = np.ravel(grid_lons)
    grid_lats = np.ravel(grid_lats)
    #if mesh = True:
    # grid_lons = grid_lons
    # grid_lats = grid_lats
    
    return grid_lons, grid_lats
    
    
    
def sGrid(lon,lat,res):
    
    
     # Transform lon, lats to the mercator projection
    x, y = pyproj.transform(salem.wgs84, proj,lon, lat)
    # take the min and max
    xmax, xmin = np.max(x), np.min(x)
    ymax, ymin = np.max(y), np.min(y)
    # Count the number of pixels
    dx = res
    nx, r = divmod(xmax - xmin, dx)
    ny, r = divmod(ymax - ymin, dx)
    # Here one could add + 1 to be sure that the last pixel is always included
    grid = salem.Grid(nxny=(nx, ny), dxdy=(dx, dx), ll_corner=(xmin, ymin), proj=proj) 

    return grid
    
    
def regridData(lon,lat,data, grid):
            
    xi, yi = grid.ij_coordinates

    # Transform lons, lats to grid
    xm, ym = grid.transform( lon.flatten(), lat.flatten(), crs=salem.wgs84)
        
    # Convert for griddata input 
    mpoints = np.array((ym, xm)).T   
    inter = np.array((np.ravel(yi), np.ravel(xi))).T
        
    # Interpolate using delaunay triangularization 
    data = griddata(mpoints, data.flatten(), inter, method='linear')
    data = data.reshape((grid.ny, grid.nx))

    return data    
    
    
def bilinear_interpolation(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)