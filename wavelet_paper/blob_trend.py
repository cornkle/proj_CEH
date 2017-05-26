import salem
from salem.utils import get_demo_file
import xarray as xr
import matplotlib.pyplot as plt
import pdb
import numpy as np
from functools import partial
from salem import get_demo_file, open_xr_dataset, GeoTiff, wgs84

file = '/users/global/cornkle/MCSfiles/blob_maps_-75/blob_map_90km_18UTC.nc'
file2 = '/users/global/cornkle/MCSfiles/blob_maps_-75/blob_map_30km_18UTC.nc'

fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
lst = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/MOD11C1_M_LSTDA_2015-06-01_rgb_3600x1800.tif'


ds = xr.open_dataarray(file)
top = xr.open_dataarray(fpath)
ds2 = xr.open_dataarray(file2)

years = np.unique(ds['time.year'])

ds = ds.sel(lon=slice(-12,10), lat=slice(11,16))  # lake chad lon=slice(10,20), lat=slice(10,15)
ds2 = ds2.sel(lon=slice(-12,10), lat=slice(11,16))   # volta lon=slice(-10,8), lat=slice(4,10)
#

print(np.unique(ds['time.year']))
dsum = ds.groupby('time.year').sum(['time', 'lon', 'lat'])
dsum2 = ds2.groupby('time.year').sum(['time', 'lon', 'lat'])

f = plt.figure()

#plt.plot(years,dsum, label='90k')
plt.plot(years,dsum2/np.sum(dsum2)-dsum/np.sum(dsum), label='30k')

plt.legend()

# top = top.sel(lon=slice(-10,-2), lat=slice(4.5,7.4))

# lon=slice(-3,3), lat=slice(6,9.5)  volta
# lon=slice(-10,-2), lat=slice(4.5,7.4) tai park
# lon=slice(-2,2), lat=slice(11,13) ouaga
#niamey = lon=slice(-3,6), lat=slice(11,18)

ds.name = '100k'
ds2.name = '30k'

# ds[ds == 0]=np.nan
# ds2[ds2 == 0] =np.nan
# #
# # perc = ds.quantile(0.95)
# # perc2 =ds2.quantile(0.95)
#
# map = ds.salem.get_map(cmap='viridis')
# #map.set_shapefile(oceans=True)
# map.set_shapefile(rivers=True)
# # read the ocean shapefile (data from http://www.naturalearthdata.com)
# oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
#                               cached=True)
#
# lakes = salem.read_shapefile(salem.get_demo_file('ne_50m_rivers_lake_centerlines.shp'), cached=True)
# map.set_shapefile(lakes, edgecolor='k', facecolor='none', linewidth=2,)
#
# grid = ds.salem.grid
#
# mask_lakes = grid.region_of_interest(shape=lakes)
#
# ds = ds - ds.mean(dim='lon')
# ds2 = ds2 - ds2.mean(dim='lon')



