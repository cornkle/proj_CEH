import salem
from salem.utils import get_demo_file
import xarray as xr
import numpy as np
from salem import get_demo_file, open_xr_dataset, GeoTiff, wgs84
from utils import u_grid as ug

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