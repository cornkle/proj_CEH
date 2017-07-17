import salem
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import pdb

path = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/'
file2 = path+'blob_map_35km_-75_sum_16-19UTC.nc'
file=path+'blob_map_35km_-75_sum_0-3UTC.nc'
fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
spath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'

ds = xr.open_dataarray(file)
top = xr.open_dataarray(fpath)
ds2 = xr.open_dataarray(file2)

ds.name = '100k'
ds2.name = '30k'

ds = ds.sel(lon=slice(-17.5,20), lat=slice(4.5,20))  # lake chad lon=slice(10,20), lat=slice(10,15)
ds2 = ds2.sel(lon=slice(-17.5,20), lat=slice(4.5,20))   # volta lon=slice(-10,8), lat=slice(4,10)
top = top.sel(lon=slice(-17.5,20), lat=slice(4.5,20))


# ds[ds == 0]=np.nan
# ds2[ds2 == 0] =np.nan

map = ds.salem.get_map(cmap='Greys')
map.set_shapefile(rivers=True)
# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

lakes = salem.read_shapefile(salem.get_demo_file('ne_50m_lakes.shp'), cached=True)
map.set_shapefile(lakes, edgecolor='k', facecolor='none', linewidth=2,)

grid = ds.salem.grid


#f,((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4,2,figsize = (18,15))

f,((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize = (11,6))



map.set_plot_params(levels=[-32,-16, -8,-4,4,8,16,32], cmap='RdBu')  #levels=[-0.6,-0.3, -0.15, 0.15, 0.3,0.6] levels=[-24,-16, -8,-4,4,8,16,20]
#map.set_plot_params(levels=[-100,-75,-50,-40,-30,30,40,50,75,100], cmap='RdBu')
dat = (ds2.values - ds.values)


map.set_data(dat)
map.visualize(ax=ax3, addcbar=True, title='Afternoon - Nighttime') #, title='>90km 1800UTC - 0300UTC'

map.set_plot_params(vmin=0, vmax=90, nlevels=10, cmap='viridis')
dat = (ds2.values*2)
map.set_data(dat)
map.visualize(ax=ax1, addcbar=True, title='1600-1900 UTC, -70C frequency')

map.set_plot_params(vmin=0, vmax=70, nlevels=8, cmap='viridis')
dat = (ds.values*2)
map.set_data(dat)
map.visualize(ax=ax2, addcbar=True, title='0000-0300 UTC, -70C frequency')

zuse = map.set_topography(top, relief_factor=1.4)
map.set_plot_params(vmax=1000, cmap='topo')
map.set_data(zuse)
map.visualize(ax=ax4, addcbar=True, title='Topography') #, title='Topography'


plt.tight_layout()


plt.savefig(spath+'/overview.png', dpi=300)


