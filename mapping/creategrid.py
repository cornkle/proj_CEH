# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 14:39:16 2016

@author: cornkle
"""

import numpy as np
import math
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