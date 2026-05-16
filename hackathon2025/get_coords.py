"""Get coordinates"""

import numpy as np
import xarray as xr
import healpy

"""K-scale config file"""

REGIONS = {'africa': [(-15., 20.),(5.,25)], 'samerica': [(-86., -26.),(-30., 14.)], 'sea': [(90., 154.),(-18., 25.)]}  #'africa': [(-20., 55.),(-35., 27.)


def get_nn_lon_lat_index_pt(nside, lons, lats):
    """for subsetting HEALPix ICON at grid point"""
    lons2, lats2 = np.meshgrid(lons, lats)
    out = xr.DataArray(
        healpy.ang2pix(nside, lons2, lats2, nest=True, lonlat=True),
        coords=[("lat", lats), ("lon", lons)],
    )

    return out


def get_nn_lon_lat_index(zoom=8, lat_range=(-30., 30.), lon_range=(-180., 180.)):
    """for subsetting HEALPix ICON out onto regular lat/lon grid"""
    lat_min = lat_range[0]
    lat_max = lat_range[1]
    lon_min = lon_range[0]
    lon_max = lon_range[1]
    lat_diff = lat_max - lat_min
    lon_diff = lon_max - lon_min
    if zoom == 6:  # ~100km -> 1 deg. = 1 pixel
        n = 1
    elif zoom == 7:  # ~30km -> 1 deg. = 3 pixel
        n = 3
    elif zoom == 8:  # ~20km -> 1 deg. = 5 pixel
        n = 5
    elif zoom == 9:  # ~10km -> 1 deg. = 10 pixel
        n = 10
    elif zoom == 10:  # ~5km -> 1 deg. = 20 pixel
        n = 20
    else:
        print("Don't know what to do with %s"%zoom)

    lats = np.linspace(lat_min, lat_max, int(lat_diff*n))
    lons = np.linspace(lon_min, lon_max, int(lon_diff*n))

    lons2, lats2 = np.meshgrid(lons, lats)
    nside = 2**zoom

    out = xr.DataArray(healpy.ang2pix(nside, lons2, lats2, nest=True, lonlat=True),
        coords=[("lat", lats), ("lon", lons)],
    )

    return out