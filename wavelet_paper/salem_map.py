import salem
from salem.utils import get_demo_file
import xarray as xr
import matplotlib.pyplot as plt
import pdb
import numpy as np


file = '/users/global/cornkle/MCSfiles/blob_map_100abs_sum.nc'
file2 = '/users/global/cornkle/MCSfiles/blob_map_30no.nc'
file3 = '/users/global/cornkle/MCSfiles/blob_map_noscthresh_sum.nc'
fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'

ds = xr.open_dataarray(file)
top = xr.open_dataarray(fpath)
ds2 = xr.open_dataarray(file2)
ds3 = xr.open_dataarray(file3)

ds = ds.sel(lon=slice(-10,8), lat=slice(4,10))  # lake dings lon=slice(10,20), lat=slice(10,15)
ds2 = ds2.sel(lon=slice(-10,8), lat=slice(4,10))
ds3 = ds3.sel(lon=slice(-10,8), lat=slice(4,10))
top = top.sel(lon=slice(-10,8), lat=slice(4,10))

ds.name = '100k'
ds2.name = '30k'

map = ds.salem.get_map(cmap='viridis')
#map.set_shapefile(oceans=True)
map.set_shapefile(rivers=True)
# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

lakes = salem.read_shapefile(salem.get_demo_file('ne_50m_lakes.shp'), cached=True)
map.set_shapefile(lakes, edgecolor='k', facecolor='none', linewidth=2,)

map.set_points(-1.32, 12.22)
map.set_points(2.1, 13.61)

grid = ds.salem.grid
grid50 = grid.regrid(factor=0.1)

grid50 = grid50.to_dataset()

#ds = ds.sum(dim='time')
#ds3 = ds3.sum(dim='time')

# ds_min, lut = grid50.salem.lookup_transform(ds, method=np.min, return_lut=True)
# ds_max = grid50.salem.lookup_transform(ds, method=np.max, lut=lut)
# ds2_min = grid50.salem.lookup_transform(ds2, method=np.min, lut=lut)
# ds2_max = grid50.salem.lookup_transform(ds2, method=np.max, lut=lut)
#
# pdb.set_trace()
#
# ds_min = ds.salem.transform(ds_min)
# ds_max = ds.salem.transform(ds_max)
# ds2_min = ds.salem.transform(ds2_min)
# ds2_max = ds.salem.transform(ds2_max)
#
# ds = (ds-ds_min.min(dim=['lon'])) / (ds_max.max(dim=['lon']) - ds_min.min(dim=['lon']))
#
# ds2 = (ds2-ds2_min.min(dim=['lon'])) / (ds2_max.max(dim=['lon']) - ds2_min.min(dim=['lon']))


ds = (ds-ds.min()) / (ds.max()- ds.min())  # dim=['lon']
ds2 = (ds2-ds2.min()) / (ds2.max()- ds2.min())
ds3 = (ds3-ds3.min()) / (ds3.max()- ds3.min())
#
# ds = ds - ds.mean(dim='lon')
# ds2 = ds2 - ds2.mean(dim='lon')
# ds3 = ds3 - ds3.mean(dim='lon')



f,((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3,2,figsize = (12,8))


map.set_data(ds)
map.set_plot_params(vmin=0., vmax=0.8, cmap='viridis')

map.visualize(ax=ax1, title='>100km')

map.set_data(ds2)
map.visualize(ax=ax2, title='>30km')

map.set_plot_params(vmin=-0.2, vmax=0.2, cmap='RdBu')
map.set_data(ds- ds.mean(dim='lon'))
map.visualize(ax=ax3, title='>100km-lonmean')

map.set_data(ds2- ds2.mean(dim='lon'))
map.visualize(ax=ax4, title='<30km-lonmean')

map.set_data(ds-ds3)
map.visualize(ax=ax5, title='<30km - >100km')

map.set_data((ds- ds.mean(dim='lon'))-(ds3- ds3.mean(dim='lon')))
map.visualize(ax=ax6, title='<30km-mean - >100km-mean')

plt.tight_layout()
plt.show()


f= plt.figure()
z = map.set_topography(top, relief_factor=1.4)
map.set_plot_params(vmax=1000, cmap='topo')
map.set_data(z)
map.set_contour(z)
map.visualize()
