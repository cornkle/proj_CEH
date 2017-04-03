import salem
from salem.utils import get_demo_file
import xarray as xr
import matplotlib.pyplot as plt
import pdb
import numpy as np


file = '/users/global/cornkle/MCSfiles/blob_map_June.nc'
fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'

ds = xr.open_dataarray(file)
top = xr.open_dataarray(fpath)
map = ds.salem.get_map(cmap='viridis')

grid = ds.salem.grid

print(ds.time.values[0], ds.time.values[-1])
ds = ds.sum(dim='time')

# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

lakes = salem.read_shapefile(salem.get_demo_file('ne_50m_lakes.shp'),
                              cached=True)

mask_default = grid.region_of_interest(shape=oceans)
mask_default_lakes = grid.region_of_interest(shape=lakes)


map.set_data(ds)
#map.set_contour(z)
map.set_shapefile(oceans=True)
map.set_shapefile(rivers=True)

map.set_shapefile(lakes, edgecolor='k', facecolor='none', linewidth=2,)
map.set_plot_params(vmax=30)
map.set_points(-1.32, 12.22)
map.set_points(2.1, 13.61)
f= plt.figure()
map.visualize()


z = map.set_topography(top, relief_factor=1.4)
f= plt.figure()
map.set_plot_params(vmax=1000, cmap='topo')
map.set_data(z)
map.set_contour(z)
map.visualize()