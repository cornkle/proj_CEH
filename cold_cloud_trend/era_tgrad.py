import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from utils import u_plot as up



file = '/users/global/cornkle/data/ERA-I monthly/ERA-WA-Monthly-2mTemp.nc'


da = xr.open_dataarray(file)

# define a function to compute a linear trend of a timeseries
def linear_trend(x):

    pf = np.polyfit(np.arange(len(x)), x, 1)
    # we need to return a dataarray or else xarray's groupby won't be happy
    return xr.DataArray(pf[0])

# stack lat and lon into a single dimension called allpoints
stacked = da.stack(allpoints=['latitude','longitude'])
# apply the function over allpoints to calculate the trend at each point
trend = stacked.groupby('allpoints').apply(linear_trend)
# unstack back to lat lon coordinates
trend_unstacked = trend.unstack('allpoints')

up.quick_map(trend_unstacked[0,:,:])