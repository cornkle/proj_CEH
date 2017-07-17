import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import numpy as np
import pdb
from utils import u_arrays as ua



fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
files = ua.locate(".nc", '/users/global/cornkle/MCSfiles/WA15_big_-40_15W-20E_size_zR/')
top = xr.open_dataarray(fpath)
ds = xr.open_dataset(files[130])



f = plt.figure(figsize=(15, 7), dpi=400)
ax = plt.axes(projection=ccrs.PlateCarree())
ds['tc_lag0'].plot.contourf('lon', 'lat', projection=ccrs.PlateCarree(), vmax=-40)
top.plot.contour('lon', 'lat', projection=ccrs.PlateCarree(), vmax=1000, levels=[50,100,250,500,750,1000], cmap='terrain')
ds['p'].plot.contour('lon', 'lat', projection=ccrs.PlateCarree(), vmin=30, cmap='Reds')
ax.coastlines()
# Gridlines
xl = ax.gridlines(draw_labels=True);
xl.xlabels_top = False
xl.ylabels_right = False
# Countries
ax.add_feature(cartopy.feature.BORDERS, linestyle='--');

plt.show()

f = plt.figure(figsize=(15, 7), dpi=400)
ax = plt.axes(projection=ccrs.PlateCarree())
#ds['tc_lag0'].plot.contourf('lon', 'lat', projection=ccrs.PlateCarree(), vmax=50)
ds['p'].plot.contour('lon', 'lat', projection=ccrs.PlateCarree(), vmin=30, cmap='Reds')
top.plot.contourf('lon', 'lat', projection=ccrs.PlateCarree(), vmax=1000, levels=[0,30,50,100,250,500,750,1000], cmap='terrain')
ds['p'].plot.contour('lon', 'lat', projection=ccrs.PlateCarree(), vmin=30, cmap='Reds')
ax.coastlines()
# Gridlines
xl = ax.gridlines(draw_labels=True);
xl.xlabels_top = False
xl.ylabels_right = False
# Countries
ax.add_feature(cartopy.feature.BORDERS, linestyle='--');

plt.show()


