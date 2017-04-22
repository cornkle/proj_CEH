import salem
from salem.utils import get_demo_file
import xarray as xr
import matplotlib.pyplot as plt
import pdb
import numpy as np
from functools import partial
from salem import get_demo_file, open_xr_dataset, GeoTiff, wgs84

file = '/users/global/cornkle/MCSfiles/blob_map_90km_sum_18UTC.nc'
file2 = '/users/global/cornkle/MCSfiles/blob_map_30km_sum_18UTC.nc'

fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
lst = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/MOD11C1_M_LSTDA_2015-06-01_rgb_3600x1800.tif'


ds = xr.open_dataarray(file)
top = xr.open_dataarray(fpath)
ds2 = xr.open_dataarray(file2)



# ds = ds.sel(lon=slice(-10,-2), lat=slice(4.5,7.4))  # lake chad lon=slice(10,20), lat=slice(10,15)
# ds2 = ds2.sel(lon=slice(-10,-2), lat=slice(4.5,7.4))   # volta lon=slice(-10,8), lat=slice(4,10)
#
# top = top.sel(lon=slice(-10,-2), lat=slice(4.5,7.4))

# lon=slice(-3,3), lat=slice(6,9.5)  volta
# lon=slice(-10,-2), lat=slice(4.5,7.4) tai park
# lon=slice(-2,2), lat=slice(11,13) ouaga
#niamey = lon=slice(-3,6), lat=slice(11,18)

ds.name = '100k'
ds2.name = '30k'

ds[ds == 0]=np.nan
ds2[ds2 == 0] =np.nan
#
# perc = ds.quantile(0.95)
# perc2 =ds2.quantile(0.95)

map = ds.salem.get_map(cmap='viridis')
#map.set_shapefile(oceans=True)
map.set_shapefile(rivers=True)
# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

lakes = salem.read_shapefile(salem.get_demo_file('ne_50m_rivers_lake_centerlines.shp'), cached=True)
map.set_shapefile(lakes, edgecolor='k', facecolor='none', linewidth=2,)

# srtm = open_xr_dataset(get_demo_file('hef_srtm.tif'))
# srtm_on_ds = ds.salem.lookup_transform(srtm)

grid = ds.salem.grid
# grid50 = grid.regrid(factor=0.1)
#
# srtm_on_ds = ds.salem.lookup_transform(top)
#
# g = GeoTiff(lst)
# # Spare memory
# ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
# g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
#              crs=wgs84, margin=10)
# ls = g.get_vardata()
# ls = np.array(ls, dtype=float)
# ls[ls >200] = np.nan
#
# ls = grid.map_gridded_data(ls, g.grid)
#
# lst_on_ds = ds.salem.lookup_transform(ls, grid=grid)
#
#
mask_lakes = grid.region_of_interest(shape=lakes)
#
# larr = ds.copy()
# larr.values = mask_lakes


# ds = (ds-1) / (perc- 1)  # dim=['lon']
# ds2 = (ds2-1) / (perc2- 1)


# ds = (ds-ds.min()) / (ds.max()- ds.min())  # dim=['lon']
# ds2 = (ds2-ds2.min()) / (ds2.max()- ds2.min())

#
ds = ds - ds.mean(dim='lon')
ds2 = ds2 - ds2.mean(dim='lon')

# ds = ds.where(ds<=1)
# ds2 = ds2.where(ds2<=1)


f,((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3,2,figsize = (12,8))


map.set_data(ds)
map.set_plot_params(vmin=0., vmax=20, cmap='viridis')

map.visualize(ax=ax1, title='>100km 2100UTC')

map.set_data(ds2)
map.visualize(ax=ax2, title='<30km 1800UTC')

map.set_plot_params(vmin=-0.4, vmax=0.4, cmap='RdBu')
map.set_data(ds- ds.mean(dim='lon'))
map.visualize(ax=ax3, title='>100km-lonmean')

map.set_data(ds2- ds2.mean(dim='lon'))
map.visualize(ax=ax4, title='<30km-lonmean')

z = map.set_topography(top, relief_factor=1.4)
map.set_plot_params(vmax=1000, cmap='topo')
map.set_data(z)
map.visualize(ax=ax5, title='Topography')

zuse = map.set_topography(top, relief_factor=1.4)
map.set_plot_params(vmax=1000, cmap='topo')
map.set_data(zuse)
# map.set_points(-7,6.3, color='black')
# map.set_points(-7.45,6.1, color='black')
# map.set_points(-7.3,5.18, color='black')
# map.set_points(-5.8,7.17, color='black')
map.set_points(-1.53,12.26, color='black')
map.set_points(-1.53,12.5, color='black')
map.set_points(-1.11,12.06, color='black')
map.set_points(-1.5,11.69, color='black')
map.set_points(2.1, 13.61, color='black')
map.visualize(ax=ax6, title='Topography')
#
# map.set_data((ds- ds.mean(dim='lon'))-(ds3- ds3.mean(dim='lon')))
# map.visualize(ax=ax6, title='<100 - mean sys frequency')

plt.tight_layout()
plt.show()

#
# f= plt.figure()
# dummy= map.set_topography()
# map.set_plot_params(vmin=-0.5, vmax=0.5, cmap='RdBu')
# map.set_data((ds2- ds2.mean(dim='lon'))-(ds- ds.mean(dim='lon')))
# map.visualize()
#
# lats = slice(4.5,6.5)
# topo = srtm_on_ds.sel(lat=lats).mean(dim='lat')
# temp = lst_on_ds.sel(lat=lats).mean(dim='lat')
#
# f=plt.figure()
# dsp = ds.sel(lat=lats).mean(dim='lat')
# dsp2 = ds2.sel(lat=lats).mean(dim='lat')
# plt.plot(ds.lon, (dsp-dsp.min())/(dsp.max()-dsp.min()), color='blue', label='>100km')
# plt.plot(ds.lon, (dsp2-dsp2.min())/(dsp2.max()-dsp2.min()), color='orange', label='<30km')
# plt.plot(ds.lon, larr.sel(lat=lats).mean(dim='lat'))
# plt.plot(ds.lon, (topo-topo.min())/(topo.max()-topo.min()), color='green', label='topo')
# plt.plot(ds.lon, (temp-temp.min())/(temp.max()-temp.min()), color='red', label='temp')

plt.legend()
