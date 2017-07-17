import pandas as pd
import numpy as np
import xarray as xr
import pdb
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from utils import u_plot as up
import salem


def slp():
    file = '/localscratch/wllf030/cornkle/obs_data/ERA-I/ERA-Int-Monthly-UVSLP.nc'


    da = xr.open_dataset(file)

    u = da['u10']
    v = da['v10']

    # ds.isel(time=ds['time.month']==8)
    u = u.isel(time=((u['time.month']>=6) & (u['time.month']<=9)))
    v = v.isel(time=((v['time.month']>=6) & (v['time.month']<=9)))

    # u = u.isel(month=(u['month']==8))
    # v = v.isel(month=(u['month']==8))


    u = u.mean(dim='time')
    v = v.mean(dim='time')

    u= u.sel(latitude=slice(20,0), longitude=slice(-20,20))
    v= v.sel(latitude=slice(20,0), longitude=slice(-20,20))

    ws  = np.sqrt(u**2+v**2)
    # uu = u[4::7, 4::7]
    # vv = v[4::7, 4::7]

    map = u.salem.get_map()

    f=plt.figure()
    ax = f.add_subplot(111)
    map.set_data(ws)

    map.visualize(ax=ax)
    qu = ax.quiver(u['longitude'].values, u['latitude'].values, u.values, v.values,
                   transform=map.transform(ax=ax), zorder=99)
    #ma# p.append_colorbar(ax=ax)


def hpa850():
    file = '/localscratch/wllf030/cornkle/obs_data/ERA-I/ERA-Int-MonthlyAvg-4D-TUVWZ.nc'


    da = xr.open_dataset(file)

    u = da['u']
    v = da['v']

    # ds.isel(time=ds['time.month']==8)
    u = u.isel(month=((u['month']>=6) & (u['month']<=9)))
    v = v.isel(month=((v['month']>=6) & (v['month']<=9)))


    u = u.mean(dim='month')
    v = v.mean(dim='month')
    print(u['level'])
    u= u.sel(level=900, latitude=slice(20,0), longitude=slice(-20,20))
    v= v.sel(level=900, latitude=slice(20,0), longitude=slice(-20,20))

    ws  = np.sqrt(u**2+v**2)
    # uu = u[4::7, 4::7]
    # vv = v[4::7, 4::7]

    map = u.salem.get_map()

    f=plt.figure()
    ax = f.add_subplot(111)
    map.set_data(ws)

    map.visualize(ax=ax)
    qu = ax.quiver(u['longitude'].values, u['latitude'].values, u.values, v.values,
                   transform=map.transform(ax=ax), zorder=99)