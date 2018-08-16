import numpy as np


def shift_lons(ds, lon_dim='lon'):
    """ Shift longitudes from [0, 360] to [-180, 180] """
    lons = ds[lon_dim].values
    new_lons = np.empty_like(lons)
    mask = lons > 180
    new_lons[mask] = -(360. - lons[mask])
    new_lons[~mask] = lons[~mask]
    ds[lon_dim].values = new_lons
    return ds
