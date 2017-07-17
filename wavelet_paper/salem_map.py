import salem
from salem.utils import get_demo_file
import xarray as xr
import matplotlib.pyplot as plt
import pdb
import numpy as np
from functools import partial
from salem import get_demo_file, open_xr_dataset, GeoTiff, wgs84
from scipy.stats.stats import pearsonr
from utils import u_grid as ug

path = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/'

# file2 = path+'blob_map_30km_sum_18UTC.nc'
# file = path+'blob_map_30km_sum_3UTC.nc'
file2 = path+'blob_map_35km_-75_sum_16-17UTC.nc'

file=path+'blob_map_35km_-75_sum_0-3UTC.nc' #blob_map_35km_-75_sum_0-3UTC.nc'

fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
#lst = '/users/global/cornkle/data/LandCover/evergreen_trees.tif'
lst = '/users/global/cornkle/data/LandCover/evergreen_trees.tif'
vegfra = '/users/global/cornkle/data/MODIS/vegfra/Average.tif'

tfile = '/users/global/cornkle/data/MODIS/LST/acp_0130MYD_h17v08_regridded.nc'


ds = xr.open_dataarray(file)
top = xr.open_dataarray(fpath)
ds = xr.open_dataarray(file)
ds2 = xr.open_dataarray(file2)
t = xr.open_dataset(tfile)


#t = ds.salem.lookup_transform(tg, grid=bgrid)

# tai park
#coord = [-8.55,-6,5,8,5.9,6.3, -7.6, -7.4]
coord = [-8.55,-6.5,5,6.7,5.8,6.25, -7.6, -7.4]
name='Tai park'
#lake volta
# coord=[-3,3,6,9.5,6.4,8,7.9,8]
# name='volta'
#kumasi
# coord = [-3.0,-1.1,5.7,7.2,6.6,6.8,-1.8,-1.5]
# name='kumasi'
#Ouaga
# coord=[-2,-0.8,12,12.8,12.3,12.4,-1.6,-1.4]
# name = 'Ouaga'
#Mali Niger wetland, topo trigger, wetland decay
# coord=[-6.8,-5.3,13.8,14.9, 14.2,14.6, -4.9, -4.3]
# name='Mopti region'
#Tamale
# coord = [-1.7,-0.1,8.5, 10.5,9.3, 9.6, -1.1,-0.78]
# name='tamale'
#Burkina national park
# name='burkina_nat_park'
# coord = [0.5, 3, 11, 12.5, 11.4,11.7, 2, 3]
#nazingA
# coord = [-2.5,-0.7,10.8,11.5,11,11.5,-2.2,-1.1]
# name='nazinga'
# coord=[2,4,12.45,15,14,14.5,2.3,4] #13.6,15.2, 15.4,12.5  # weird anticorrelation
# name = 'Bonkoukou'
# coord=[-14,-13,12.5,16,14.8,15.5,-15,-14] #13.6,15.2, 15.4,12.5  # weird anticorrelation
# name = 'east senegal'



ds = ds.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
ds2 = ds2.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
top = top.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
t = t.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))

pdb.set_trace()

ds.name = '0-3UTC'
ds2.name = '16-18UTC'

map = ds.salem.get_map(cmap='viridis')
#map.set_shapefile(oceans=True)
map.set_shapefile(rivers=True)
# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

river = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/rivers/ne_10m_rivers_lake_centerlines.shp', cached=True)
lakes = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/lakes/ne_10m_lakes.shp', cached=True)
map.set_shapefile(lakes, edgecolor='k', facecolor='none', linewidth=2,)

srtm = open_xr_dataset(get_demo_file('hef_srtm.tif'))
srtm_on_ds = ds.salem.lookup_transform(srtm)

grid = ds.salem.grid
srtm_on_ds = ds.salem.lookup_transform(top)
#t_on_ds = ds.salem.lookup_transform(t, grid=grid)
#
g = GeoTiff(lst)
# Spare memory
ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
             crs=wgs84, margin=10)
ls = g.get_vardata()
ls = np.array(ls, dtype=float)
ls[ls >200] = np.nan
ls = grid.map_gridded_data(ls, g.grid)
lst_on_ds = ds.salem.lookup_transform(ls, grid=grid)

g = GeoTiff(vegfra)
# Spare memory
ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
             crs=wgs84, margin=10)
ls = g.get_vardata()
ls = np.array(ls, dtype=float)
ls[ls >100] = np.nan
ls = grid.map_gridded_data(ls, g.grid)
vegfra_on_ds = ds.salem.lookup_transform(ls, grid=grid)


mask_river = grid.region_of_interest(shape=river)
mask_lakes = grid.region_of_interest(shape=lakes, roi=mask_river, all_touched=True)
mask_river[mask_lakes]= 1
larr = ds.copy()
larr.values = mask_river
ds2f = ds2.values.flatten()
dsf = ds.values.flatten()
print(pearsonr(ds2.values.flatten(), lst_on_ds.values.flatten()))

print(pearsonr((ds2f-ds2f.min())/(ds2f.max()-ds2f.min()), (dsf-dsf.min())/(dsf.max()-dsf.min())))

f,((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize = (12,8))
#(ax5, ax6))

map.set_plot_params(vmin=0., vmax=30, nlevels=9, cmap='viridis')

map.set_data(ds.values)

map.visualize(ax=ax2, title='0-3UTC')
map.set_plot_params(vmin=0., vmax=30, nlevels=9, cmap='viridis')
map.set_data(ds2.values)
map.visualize(ax=ax1, title='16-17UTC')


map.set_plot_params( vmax=100, vmin=0, cmap='viridis', nlevels=11)#( vmax=22, vmin=19, cmap='viridis')#( vmax=100, vmin=50, cmap='viridis')
#map.set_data(t_on_ds-273.15)
map.set_data(lst_on_ds)
map.visualize(ax=ax3, title='Vegetation fraction')

z = map.set_topography(top, relief_factor=1.4)
map.set_plot_params(vmax=700, vmin=0, cmap='topo')
map.set_data(z)
map.visualize(ax=ax4, title='Topography')

#ax4.set_axis_off()

# zuse = map.set_topography(top, relief_factor=1.4)
# map.set_plot_params(vmax=700, cmap='topo')
# map.set_data(zuse)
# # map.set_points(-7,6.3, color='black')
# # map.set_points(-7.45,6.1, color='black')
# # map.set_points(-7.3,5.18, color='black')
# # map.set_points(-5.8,7.17, color='black')
# map.set_points(-1.53,12.26, color='black')
# map.set_points(-1.53,12.5, color='black')
# map.set_points(-1.11,12.06, color='black')
# map.set_points(-1.5,11.69, color='black')
# map.set_points(2.1, 13.61, color='black')
# map.visualize(ax=ax6, title='Topography')


plt.tight_layout()
plt.savefig('/users/global/cornkle/VERA/plots/leeds_june_2017/map_'+name+'.png', dpi=300)

#
lats = slice(coord[4], coord[5])
topo = srtm_on_ds.sel(lat=lats).mean(dim='lat')
temp = lst_on_ds.sel(lat=lats).mean(dim='lat')
veg = vegfra_on_ds.sel(lat=lats).mean(dim='lat')
dsp = ds.sel(lat=lats).mean(dim='lat')
dsp2 = ds2.sel(lat=lats).mean(dim='lat')
tt = t.sel(lat=lats).mean(dim='lat')

f=plt.figure(figsize=(11,4))
ax = f.add_subplot(111)

plt.plot(ds.lon, (dsp-dsp.min())/(dsp.max()-dsp.min()), color='darkblue', label='0-3UTC', marker='o')
ax.plot(ds.lon, (dsp2-dsp2.min())/(dsp2.max()-dsp2.min()), color='orangered', label='16-17UTC', marker='o')
ax.plot(ds.lon, (tt-tt.min())/(tt.max()-tt.min()), color='orange', label='LST')
ax.plot(ds.lon, larr.sel(lat=lats).mean(dim='lat'), label='River', linestyle='dotted')
ax1 = ax.twinx()
ax1.plot(ds.lon, topo, color='brown', label='Topo', linestyle='dotted')

#ax.axvline(0,ymin=0, ymax=1)
#ax.plot(ds.lon,(veg-veg.min())/(veg.max()-veg.min()), color='seagreen', label='Vegetation fraction')
ax.plot(ds.lon,  (temp-temp.min())/(temp.max()-temp.min()), color='seagreen', label='Evergreen trees')
ax.set_xlabel('Longitude')
ax.set_ylabel('Normalised value range (-)')
ax.legend()
ax1.legend()
ax1.set_title('Sub-cloud features <-75C, <35km | MCSs > 15000km2')

print(pearsonr(dsp.values.flatten(), dsp2.values.flatten()))

lons = slice(coord[6], coord[7])#(-7.6, -7.4)
topo2 = srtm_on_ds.sel(lon=lons).mean(dim='lon')
temp2 = lst_on_ds.sel(lon=lons).mean(dim='lon')
veg2 = vegfra_on_ds.sel(lon=lons).mean(dim='lon')
dsp2 = ds.sel(lon=lons).mean(dim='lon')
dsp22 = ds2.sel(lon=lons).mean(dim='lon')
#tt2 = t.sel(lat=lats).mean(dim='lon')
#
# f=plt.figure()
# ax = f.add_subplot(111)
#
# plt.plot(ds.lat, (dsp2-dsp2.min())/(dsp2.max()-dsp2.min()), color='red', label='0-3UTC', marker='o')
# ax.plot(ds.lat, (dsp22-dsp22.min())/(dsp22.max()-dsp22.min()), color='blue', label='16-17UTC', marker='o')
# #ax.plot(ds.lat, (tt2-tt2.min())/(tt2.max()-tt2.min()), color='orange', label='LST')
# ax.plot(ds.lat, larr.sel(lon=lons).mean(dim='lon'), label='River', linestyle='dotted')
# ax1 = ax.twinx()
# ax1.plot(ds.lat, topo2, color='brown', label='Topo', linestyle='dotted')
#
# #ax.axvline(0,ymin=0, ymax=1)
# #ax.plot(ds.lon,(veg-veg.min())/(veg.max()-veg.min()), color='green', label='Vegetation fraction')
# ax.plot(ds.lat,  (temp2-temp2.min())/(temp2.max()-temp2.min()), color='green', label='Evergreen trees')
# ax.set_xlabel('Latitude')
# ax.set_ylabel('Normalised value range (-)')
# ax.legend()
# ax1.legend()



plt.savefig('/users/global/cornkle/VERA/plots/leeds_june_2017/cross_'+name+'.png', dpi=300)
plt.show()
