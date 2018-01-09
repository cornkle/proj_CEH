import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import pandas as pd
from utils import u_plot as up
import salem


def era_Tlapse(month, temp, lon, lat):
    file = '/users/global/cornkle/data/ERA-I monthly/t_3d_2004_monthly.nc'

    da = xr.open_dataset(file)
    da = da.sel(longitude=lon, latitude=lat, method='nearest')
    da = da.isel(time=month-1)

    g = 9.80665
    t = da['t']-273.15
    z = da['z']
    zm = z / g  ## geopotential / gravity constant

    X = np.abs(t-temp)
    idx = np.argmin(X.values)
    height = zm.values[idx]

    return height
