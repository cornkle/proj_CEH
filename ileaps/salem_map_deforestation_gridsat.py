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
import os
import shapely.geometry as shpg

path = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/circles/'
hours = '18'
mon = 'MAM'
file18 = '/users/global/cornkle/data/OBS/gridsat/gridsat_netcdf/circles/single/'+mon+'-69_18UTC/*.nc'

fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
vegfra = '/users/global/cornkle/data/LandCover/evergreen_trees.tif'
lst = '/users/global/cornkle/data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_loss_10N_010W.tif'


if not os.path.isfile(path+'ds_-69_sum_'+mon+'_'+hours+'UTC_2006.nc'):

    #ds3 = xr.open_mfdataset(file3)
    ds18 = xr.open_mfdataset(file18)
    print('Starting to write nc files')

    ds18_hist = ds18['power'][(ds18['time.hour']== 18) & (ds18['time.year']>=2006) & (ds18['time.year']<=2009)].sum(dim='time')
    ds18_hist.to_netcdf(path + 'ds_-69_sum_'+hours+'UTC_2006.nc')

    ds18_present = ds18['power'][(ds18['time.hour']== 18) & (ds18['time.year']>=2012)& (ds18['time.year']<=2015)].sum(dim='time')
    ds18_present.to_netcdf(path + 'ds_-69_sum_'+hours+'UTC_2012.nc')
else:
    ds18_hist = xr.open_dataarray(path+'ds_-69_sum_'+hours+'UTC_2006.nc')
    ds18_present = xr.open_dataarray(path + 'ds_-69_sum_'+hours+'UTC_2012.nc')


top = xr.open_dataarray(fpath)

tdummy = xr.open_mfdataset('/users/global/cornkle/data/MODIS/LST_MOD11C3/clim/2006-2009_*_day.nc')
tdummy2 = xr.open_mfdataset('/users/global/cornkle/data/MODIS/LST_MOD11C3/clim/2012-2015_*_day.nc')

# coord = [-8.55,-5,5,7.7,5.8,6.25, -7.6, -7.4]
# name='Tai park'

coord = [-8,-5.5,6,8.5,7,7.3, -7.6, -7.4]
name='Tai park'

ds18_hist = ds18_hist.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
ds18_present = ds18_present.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
top = top.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
tdummy = tdummy.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]), month=5)
tdummy2 = tdummy2.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]), month=5)
t2 = tdummy2['lst']
t = tdummy['lst']

ds18_hist.name = '2004-2008'
ds18_present.name = '2011-2015'

map = ds18_present.salem.get_map(cmap='viridis')
#map.set_shapefile(oceans=True)
map.set_shapefile(rivers=True)
# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

river = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/rivers/ne_10m_rivers_lake_centerlines.shp', cached=True)
lakes = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/lakes/ne_10m_lakes.shp', cached=True)
map.set_shapefile(lakes, edgecolor='k', facecolor='none', linewidth=2,)


srtm_on_ds = ds18_present.salem.lookup_transform(top)
t_on_ds = ds18_present.salem.transform(t)
t2_on_ds = ds18_present.salem.transform(t2)

grid = ds18_present.salem.grid
#deforestation
g = GeoTiff(lst)
ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
             crs=wgs84, margin=10)
ls = g.get_vardata()
ls = np.array(ls, dtype=float)
lst_on_ds = ds18_present.salem.lookup_transform(ls, grid=g.grid)

#evergreen forest
g = GeoTiff(vegfra)

ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
             crs=wgs84, margin=10)
ls = g.get_vardata()
ls = np.array(ls, dtype=float)
ls[ls >100] = np.nan
ls = grid.map_gridded_data(ls, g.grid)
vegfra_on_ds = ds18_present.salem.lookup_transform(ls, grid=grid)

mask_river = grid.region_of_interest(shape=river)
mask_lakes = grid.region_of_interest(shape=lakes, roi=mask_river, all_touched=True)
mask_river[mask_lakes]= 1
larr = ds18_present.copy()
larr.values = mask_river

ds2f = ds18_present.values.flatten()
dsf = ds18_hist.values.flatten()
print(pearsonr(ds18_present.values.flatten(), lst_on_ds.values.flatten()))
print(pearsonr((ds2f-ds2f.min())/(ds2f.max()-ds2f.min()), (dsf-dsf.min())/(dsf.max()-dsf.min())))


f,((ax1, ax2, ax3) , (ax4, ax5, ax6)) = plt.subplots(2,3,figsize = (14,8))

map.set_plot_params(vmin=1., vmax=10, nlevels=10, cmap='viridis')
map.set_data(ds18_hist.values)
geom = shpg.LineString(((coord[0], coord[4]), (coord[1], coord[4])))
map.set_geometry(geom, zorder=99, color='r')
geom = shpg.LineString(((coord[0], coord[5]), (coord[1], coord[5])))
map.set_geometry(geom, zorder=99, color='r')
map.visualize(ax=ax2, title=ds18_hist.name)

map.set_plot_params(vmin=1., vmax=10, nlevels=10, cmap='viridis')
map.set_data(ds18_present.values)
map.visualize(ax=ax1, title=ds18_present.name)

map.set_plot_params( vmax=1, vmin=0, cmap='viridis', nlevels=11)
map.set_data(lst_on_ds)
map.visualize(ax=ax3, title='Deforestation fraction')

map.set_plot_params( vmax=100, vmin=0, cmap='viridis', nlevels=11)
map.set_data(vegfra_on_ds)
map.visualize(ax=ax4, title='Evergreen forest 2004')

z = map.set_topography(top, relief_factor=1.4)

map.set_plot_params(vmax=28, vmin=23, cmap='jet')
map.set_data(t_on_ds.values-273.15)
map.visualize(ax=ax5, title='Temperature')

map.set_plot_params(vmax=700, vmin=0, cmap='topo')
map.set_data(z)
map.visualize(ax=ax6, title='Topography')
plt.tight_layout()

# plt.savefig('/users/global/cornkle/VERA/plots/map_'+name+'.png', dpi=300)

lats = slice(coord[4], coord[5])
topo = srtm_on_ds.sel(lat=lats).mean(dim='lat')
temp = lst_on_ds.sel(lat=lats).mean(dim='lat')
veg = vegfra_on_ds.sel(lat=lats).mean(dim='lat')
dsp = ds18_hist.sel(lat=lats).mean(dim='lat')
dsp2 = ds18_present.sel(lat=lats).mean(dim='lat')
tt = t_on_ds.sel(lat=lats).mean(dim='lat')
tt2 = t2_on_ds.sel(lat=lats).mean(dim='lat')

f=plt.figure(figsize=(11,4))
ax = f.add_subplot(111)

# plt.plot(ds18_present.lon, (dsp-dsp.min())/(dsp.max()-dsp.min()), color='darkblue', label=ds18_hist.name, marker='o')
# ax.plot(ds18_present.lon, (dsp2-dsp2.min())/(dsp2.max()-dsp2.min()), color='orangered', label=ds18_present.name, marker='o')
ax.plot(ds18_present.lon, dsp, color='darkblue', label=ds18_hist.name, marker='o')
ax.plot(ds18_present.lon, dsp2, color='orangered', label=ds18_present.name, marker='o')
# ax.plot(ds18_present.lon, (tt-tt.min())/(tt.max()-tt.min()), color='red', label='Temp2007')
# ax.plot(ds18_present.lon, (tt2-tt2.min())/(tt2.max()-tt2.min()), color='orange', label='Temp2011')
geom = shpg.LineString(((coord[0], coord[4]), (coord[1], coord[4])))
map.set_geometry(geom, zorder=99, color='r')

geom = shpg.LineString(((coord[0], coord[5]), (coord[1], coord[5])))
map.set_geometry(geom, zorder=99, color='r')
#ax.plot(ds.lon, (tt-tt.min())/(tt.max()-tt.min()), color='orange', label='LST')
#ax.plot(ds18_present.lon, larr.sel(lat=lats).mean(dim='lat'), label='River', linestyle='dotted')
ax1 = ax.twinx()
ax1.plot(ds18_present.lon, tt-273.15, color='red', label='Temp2007')
ax1.plot(ds18_present.lon, tt2-273.15, color='orange', label='Temp2011')
#ax1.plot(ds18_present.lon, topo, color='brown', label='Topo', linestyle='dotted')

#ax.axvline(0,ymin=0, ymax=1)
#ax.plot(ds18_present.lon,(veg-veg.min())/(veg.max()-veg.min()), color='seagreen', label='Evergreen trees')
ax.plot(ds18_present.lon,  (temp-temp.min())/(temp.max()-temp.min()), color='grey', label='Deforestation')
ax.set_xlabel('Longitude')
ax.set_ylabel('Normalised value range (-)')
ax.legend()
ax1.legend()
ax1.set_title('Sub-cloud features <-75C, <35km | MCSs > 15000km2')

print(pearsonr(dsp.values.flatten(), dsp2.values.flatten()))