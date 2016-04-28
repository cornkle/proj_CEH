# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 14:14:32 2016

@author: cornkle
"""

import mpl_toolkits.basemap.pyproj as pyproj 

# manuallyincluded projection, proj4 notation
isn2004=pyproj.Proj("+proj=lcc +lat_1=64.25 +lat_2=65.75 +lat_0=65 +lon_0=-19 +x_0=1700000 +y_0=300000 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1")

wgs84=pyproj.Proj("+init=EPSG:4326") # latlon with wgs84 datum, is a geographic coordinate system, lat lon! 

mercator=pyproj.Proj("+init=EPSG:3857") # spatial reference system mercator - is a projected coodinate sysem, x y!09