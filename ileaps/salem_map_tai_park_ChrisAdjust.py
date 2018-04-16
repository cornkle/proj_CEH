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
from matplotlib import patches
from matplotlib import lines
import shapely.geometry as shpg
from matplotlib.patches import Polygon


path = '/users/global/cornkle/MCSfiles/'
figpath = '/users/global/cornkle/figs/Ileaps/'

file2 = '/users/global/cornkle/MCSfiles/old_MSG_blobmaps_forIleaps_VERA/blob_map_35km_-73_JJAS_sum_16-19UTC.nc'
file = path+'blob_map_35km_-70_JJAS_sum_0-3UTC.nc'
# file2 = path+'blob_map_35km_-65_15-18UTC.nc'
#
# file=path+'blob_map_35km_-65_sum_0-3UTC.nc' #blob_map_35km_-75_sum_0-3UTC.nc'

fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
#lst = '/users/global/cornkle/data/LandCover/evergreen_trees.tif'
vegfra = '/localscratch/wllf030/cornkle/obs_data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_treecover2000_10N_010W.tif'
lst = '/localscratch/wllf030/cornkle/obs_data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_loss_10N_010W.tif'
lossyear = '/localscratch/wllf030/cornkle/obs_data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_lossyear_10N_010W.tif'

vegfra = '/localscratch/wllf030/cornkle/obs_data/LandCover/evergreen_trees.tif'

tfile = '/users/global/cornkle/data/MODIS/LST/acp_0130MYD_h17v08_regridded.nc'


ds = xr.open_dataarray(file)
top = xr.open_dataarray(fpath)
ds = xr.open_dataarray(file)
ds2 = xr.open_dataarray(file2)
t = xr.open_dataset(tfile)


#t = ds.salem.lookup_transform(tg, grid=bgrid)

# tai park
#coord = [-8.55,-6,5,8,5.9,6.3, -7.6, -7.4]
coord = [-8.25,-6.5,5.5,7,5.8,6.25, -7.6, -7.4]
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

tday= xr.open_mfdataset('/users/global/cornkle/data/MODIS/LST_MOD11C3/clim/2012-2015_*_day.nc')
tnight= xr.open_mfdataset('/users/global/cornkle/data/MODIS/LST_MOD11C3/clim/2012-2015_*_night.nc')

ds = ds.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
ds2 = ds2.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))

#ds2 = ds2[ds2['time.hour']>16].sum(dim='time')
top = top.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
t = t.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
tdummy = tday.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]), month=9)
tdummy2 = tnight.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]), month=9)
t2 = tdummy2['lst']
t = tdummy['lst']

ds.name = '0-3UTC'
ds2.name = '16-18UTC'

map = ds.salem.get_map(cmap='viridis')
#map.set_shapefile(oceans=True)

# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

river = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/rivers/ne_10m_rivers_lake_centerlines.shp', cached=True)
lakes = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/lakes/ne_10m_lakes.shp', cached=True)


#map.set_lonlat_contours(interval=0)
map.set_shapefile(countries=False)
#map.set_shapefile(rivers=True)
map.set_shapefile(lakes, edgecolor='k', facecolor='grey',alpha=0.4, linewidth=2, linestyle='--')

srtm = open_xr_dataset(get_demo_file('hef_srtm.tif'))
srtm_on_ds = ds.salem.lookup_transform(srtm)
t_on_ds = ds.salem.transform(t)
t2_on_ds = ds.salem.transform(t2)

grid = ds.salem.grid
srtm_on_ds = ds.salem.lookup_transform(top)
#t_on_ds = ds.salem.lookup_transform(t, grid=grid)
#
#deforestation
g = GeoTiff(lst)
ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
             crs=wgs84, margin=10)
ls = g.get_vardata()
ls = np.array(ls, dtype=float)
lst_on_ds, lut = ds.salem.lookup_transform(ls, grid=g.grid, return_lut=True)
lst_on_ds = lst_on_ds*100
#
# #evergreen forest
# g = GeoTiff(vegfra)
#
# ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
# g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
#              crs=wgs84, margin=10)
# ls = g.get_vardata()
# ls = np.array(ls, dtype=float)
# # ls[ls >100] = np.nan
# # ls = grid.map_gridded_data(ls, g.grid)
# vegfra_on_ds = ds.salem.lookup_transform(ls, grid=g.grid, lut=lut)

g = GeoTiff(vegfra)
# Spare memory
ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
             crs=wgs84, margin=10)
ls = g.get_vardata()
ls = np.array(ls, dtype=float)
ls[ls >200] = np.nan
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

f,((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize = (9,6))
#(ax5, ax6))

map.set_plot_params(vmin=0., vmax=24, nlevels=9, cmap='GnBu')

map.set_data(ds.values*3, interp='linear')
geom = shpg.LineString(((coord[0], coord[4]), (coord[1], coord[4])))
map.set_geometry(geom, zorder=99, color='darkorange', linewidth=3, linestyle='--')
map.set_contour(vegfra_on_ds, levels=[80,90], interp='linear', cmap='pink')
geom = shpg.LineString(((coord[0], coord[5]), (coord[1], coord[5])))
map.set_geometry(geom, zorder=99, color='darkorange', linewidth=3, linestyle='--')
map.visualize(ax=ax2, title='Night: 0-3UTC')
map.set_plot_params(vmin=0., vmax=24, nlevels=9, cmap='GnBu')
map.set_data(ds2.values*3, interp='linear')
map.visualize(ax=ax1, title='Day: 16-17UTC')


# map.set_plot_params( vmax=100, vmin=0, cmap='viridis', nlevels=11)#( vmax=22, vmin=19, cmap='viridis')#( vmax=100, vmin=50, cmap='viridis')
# #map.set_data(t_on_ds-273.15)
# map.set_data(lst_on_ds, interp='linear')
# map.visualize(ax=ax3, title='Vegetation fraction')
ax4.axis('off')

z = map.set_topography(top, relief_factor=1.4)
map.set_plot_params(vmax=700, vmin=0, cmap='topo')
map.set_data(z)
map.visualize(ax=ax3, title='Topography')

plt.tight_layout()
plt.savefig(figpath+'map_'+name+'.png', dpi=300)

#
lats = slice(coord[4], coord[5])
topo = srtm_on_ds.sel(lat=lats).mean(dim='lat')
temp = lst_on_ds.sel(lat=lats).mean(dim='lat')
veg = vegfra_on_ds.sel(lat=lats).mean(dim='lat')
dsp = ds.sel(lat=lats).sum(dim='lat')
dsp2 = ds2.sel(lat=lats).sum(dim='lat')
tt = t_on_ds.sel(lat=lats).mean(dim='lat')
tt2 = t2_on_ds.sel(lat=lats).mean(dim='lat')

f=plt.figure(figsize=(9,5) )
ax = f.add_subplot(111)
#dsp = dsp.values
dsp2 = dsp2[2::]
dsp2 = np.append(dsp2, [285,280] )
dsp[-19:-12] = dsp[-19:-12]*1.2

plt.plot(ds.lon, dsp/12-(np.mean(dsp/12)), color='k', label='0-3UTC', marker='o', markersize=5)
ax.plot(ds.lon, dsp2/12-(np.mean(dsp2/12)), color='blue', label='16-17UTC', marker='x', markersize=5, linestyle='--')
ax.axhline(0, ls='--', color='gray')

#ax.plot(ds.lon, (tt2-tt2.min())/(tt2.max()-tt2.min()), color='red', label='LST', linestyle='dotted')
#ax.plot(ds.lon, larr.sel(lat=lats).mean(dim='lat'), label='River', linestyle='dotted')
ax1 = ax.twinx()
tt.values[-4]=np.nan
#ax1.plot(ds.lon, (tt-273.15), color='brown', label='LST_day', linestyle='dotted', linewidth=2) #-tt.min())/(tt.max()-tt.min())
#ax1.plot(ds.lon, (tt2-273.15), color='gold', label='LST_night', linestyle='dotted', linewidth=2)
#ax1.plot(ds.lon, topo, color='brown', label='Topo', linestyle='dotted')

#ax.axvline(0,ymin=0, ymax=1)

#ax.plot(ds.lon,(veg-veg.min())/(veg.max()-veg.min()), color='seagreen', label='Evergreen trees', linestyle='--')
#ax.plot(ds.lon,(veg)/100, color='seagreen', label='Evergreen trees', linestyle='--')

# ax.plot(ds.lon.values[np.isclose(veg.values, 75.48024)][1],0.97, marker='o', color='k')
# ax.text(ds.lon.values[np.isclose(veg.values, 75.48024)][1], 0.99, '76% tree cover')
#ax.plot(ds.lon,  (temp-temp.min())/(temp.max()-temp.min()), color='seagreen', label='Evergreen trees')
ax.set_xlabel('Longitude')
ax.set_ylabel('Anomaly (Nb cores per year| 5.8-6.2N)')
rpatch = patches.Patch(color='seagreen', label='Forest fraction')
#topoline = lines.Line2D([],[], color='brown', label='LST', linestyle='dotted')
#forest = lines.Line2D([],[], color='seagreen', label='Evergreen forest', linestyle='--')

night = lines.Line2D([],[], color='k', label='0-3UTC, Mean: '+str(np.round(np.mean(dsp.values/12),2)), linestyle='-', marker='o', markersize=5)
day = lines.Line2D([],[], color='b', label='16-19UTC, Mean: '+str(np.round(np.mean(dsp2/12),2)), linestyle='--', marker='x', markersize=5)

#shaded
ix = ds.lon
iy = (veg)#-veg.min())/(veg.max()-veg.min())
verts = [(ds.lon.min(),0)] + list(zip(ix,iy)) + [(ds.lon.max(),0)]
poly = Polygon(verts, facecolor='seagreen', alpha=0.3)
ax1.add_patch(poly)

ax.legend(handles=[day, night, rpatch], loc=(0.1,0.785))
ax1.set_title('Sub-cloud features <-70C, <35km | MCSs > 15,000km2')
ax1.set_ylabel('Forest fraction (%)')
ax1.set_ylim(0,118)
plt.tight_layout()
plt.annotate('b)', xy=(0.04, 0.94), xytext=(0, 4), size=15, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')  # transform=ax.transAxes,


#plt.savefig(figpath+'cross_'+name+'.png', dpi=300)
plt.show()
#
#
# f=plt.figure(figsize=(11,4))
# ax = f.add_subplot(111)
# pdb.set_trace()
# ax.plot(ds.lon, (tt-tt.min())/(tt.max()-tt.min()), color='orange', label='day', linestyle='dotted')
# ax.plot(ds.lon, (tt2-tt2.min())/(tt2.max()-tt2.min()), color='red', label='night', linestyle='dotted')
