import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
from scipy import stats


def era_Tlapse(month, temp, lon, lat):
    #file = '/users/global/cornkle/data/ERA-I monthly/t_3d_2004_monthly.nc'
    file ='/localscratch/wllf030/cornkle/obs_data/ERA-I/ERA-Int-MonthlyAvg-4D-TUVWZ_Africa.nc'
    da = xr.open_dataset(file)
    da = da.sel(longitude=lon, latitude=lat, method='nearest')
    da = da.isel(month=month-1)
    print('Read ERA data')

    g = 9.80665
    t = da['t']-273.15
    z = da['z']
    zm = z / g  ## geopotential / gravity constant
    ismin = np.argmin(t.values)

    gradient, intercept, r_value, p_value, std_err = stats.linregress(zm[ismin:ismin+2], t[ismin:ismin+2])
    t[0:ismin + 1] = gradient * zm[0:ismin + 1] + intercept ## linear tropopause correction (after tmin, temp rises!)

    X = np.abs(t-temp)
    idx = np.argmin(X.values)
    height = zm.values[idx]

    #plt.plot(zm, t)
    # plt.plot(zm[0:ismin+1], gradient*zm[0:ismin+1]+intercept)
    return height
