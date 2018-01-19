import salem
from salem.utils import get_demo_file
import xarray as xr
import numpy as np
from salem import get_demo_file, open_xr_dataset, GeoTiff, wgs84
from utils import u_grid as ug
import pdb
import math

def read_tiff(file):

    g = GeoTiff(file)
    ls = g.get_vardata()

    return ls

def grid_tiff(file, grid):

    g = GeoTiff(file)
    # Spare memory
    ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
    g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
                 crs=wgs84, margin=10)
    ls = g.get_vardata()
    ls_on_grid = grid.lookup_transform(ls, grid=grid)

    return ls_on_grid


def parallax_corr_msg(slon, slat, plon, plat, height):

    er = 6378.137 # earth radius km
    mh = 35786 #msg satellite height
    tot = er + mh
    lat_diff = np.abs(slat-plat)
    lon_diff = np.abs(slon-plon)
    lax_lat = (height * tot * math.sin(math.radians(lat_diff))) / (er*(tot*math.cos(math.radians(lat_diff))-(er+height)))
    lax_lon = (height * tot * math.sin(math.radians(lon_diff))) / (er * (tot * math.cos(math.radians(lon_diff)) - (er + height)))

    lax_y = er * lax_lat
    lax_x = er * lax_lon

    #print('Degree_lat_lon' , lax_lat, lax_lon)
    #print('Km_y_x', lax_y, lax_x)

    if plon < slon:
        lax_x = lax_x * -1

    if plat < slat:
        lax_y = lax_y * -1

    return (lax_x, lax_y)






