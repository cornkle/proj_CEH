import salem
from salem.utils import get_demo_file
import xarray as xr
import matplotlib.pyplot as plt
import ipdb
import numpy as np
from functools import partial
from salem import get_demo_file, open_xr_dataset, GeoTiff, wgs84
from scipy.stats.stats import pearsonr
from utils import u_grid as ug
import os
import shapely.geometry as shpg
from matplotlib import patches
from matplotlib import lines
import shapely.geometry as shpg
from matplotlib.patches import Polygon
from utils import constants as cnst


# path = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/'

path = cnst.network_data + 'MCSfiles/old_MSG_blobmaps_forIleaps_VERA/'
figpath = cnst.network_data + 'figs/VERA/chris_egu'

hours = '18-19'
tstring = '-73'
file18 = path+'blob_map_35km_'+tstring+'_MAMJ_16-19UTC.nc'#blob_map_35km_-70_15-18UTC.nc' #blob_map_35km_-75_sum_0-3UTC.nc'

fpath = cnst.ANCILS + 'gtopo_1min_afr.nc'
vegfra = cnst.local_data + 'obs_data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_treecover2000_10N_010W.tif'
lst = cnst.local_data + 'obs_data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_loss_10N_010W.tif'
lossyear = cnst.local_data + 'obs_data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_lossyear_10N_010W.tif'
vegfra = cnst.local_data + 'obs_data/LandCover/evergreen_trees.tif'

if not os.path.isfile(path+'blob_map_35km_'+tstring+'_sum_MAMJ_'+hours+'UTC_2006.nc'):

    #ds3 = xr.open_mfdataset(file3)
    ds18 = xr.open_dataarray(file18)
    print('Starting to write nc files')

    ds18 = ds18.sel(lon=slice(-10, 10), lat=slice(4,9))

    #ds3_hist = ds3['power'][(ds3['time.hour']==3) & (ds3['time.year']<=1995)].sum(dim='time')
    ds18_hist = ds18[(ds18['time.month']<= 5) & (ds18['time.hour']>= 16) & (ds18['time.hour']<= 17) & (ds18['time.year']>=2005) & (ds18['time.year']<=2009)].sum(dim='time')
    ds18_hist.values = ds18_hist.values/5
    ds18_hist.to_netcdf(path + 'blob_map_35km_'+tstring+'_sum_MAMJ_'+hours+'UTC_2006.nc')

    #ds3_present = ds3['power'][(ds3['time.hour']==3) & (ds3['time.year']>=2004)].sum(dim='time')
    ds18_present = ds18[(ds18['time.month']<= 5) & (ds18['time.hour']>= 16) & (ds18['time.hour']<= 17) & (ds18['time.year']>=2011)& (ds18['time.year']<=2015)].sum(dim='time')
    ds18_present.values = ds18_present.values / 5
    ds18_present.to_netcdf(path + 'blob_map_35km_'+tstring+'_sum_MAMJ_'+hours+'UTC_2012.nc')
else:
    #ds3_hist = xr.open_dataarray(path+'ds_hist1995_'+mon+'_3UTC.nc')
    ds18_hist = xr.open_dataarray(path+'blob_map_35km_'+tstring+'_sum_MAMJ_'+hours+'UTC_2006.nc')
    #ds3_present = xr.open_dataarray(path + 'ds_present2004_'+mon+'_3UTC.nc')
    ds18_present = xr.open_dataarray(path + 'blob_map_35km_'+tstring+'_sum_MAMJ_'+hours+'UTC_2012.nc')


top = xr.open_dataarray(fpath)

tdummy = xr.open_mfdataset(cnst.network_data + 'data/MODIS/LST_MOD11C3/clim/2006-2009_*_day.nc')
tdummy2 = xr.open_mfdataset(cnst.network_data + 'data/MODIS/LST_MOD11C3/clim/2012-2015_*_day.nc')

# coord = [-8.55,-5,5,7.7,6.1,6.5, -7.6, -7.4]
# name='Tai park'

coord = [-7.85,-6.4,6.2,7.45,6.9,7.4, -7.6, -7.4]
#coord = [-8.1,-6.4,6.1,7.45,6.1,6.4, -7.6, -7.4]
name='Tai park'

# coord = [-8,-6.5,5.5,7.7,6.9,7.4, -7.6, -7.4]
# name='Tai park'

# coord=[-8,-1,5.35,8,5.5,6,7.9,8]
# name='whatever'

# coord=[-3,3,6,9.5,6.4,8,7.9,8]
# name='volta'

ds18_hist = ds18_hist.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
ds18_present = ds18_present.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
top = top.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
tdummy = tdummy.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]), month=5)
tdummy2 = tdummy2.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]), month=5)
t2 = tdummy2['lst']
t = tdummy['lst']

ds18_hist.name = '2004-2009'
ds18_present.name = '2011-2015'

perc = ds18_hist.max()
perc2 =ds18_present.max()

mm1=ds18_hist.min()
mm2=ds18_present.min()

percboth = np.max([perc,perc2])
minboth = np.min([mm1,mm2])

ds18_hist = (ds18_hist-minboth)/(percboth-minboth)
ds18_present = (ds18_present-minboth)/(percboth-minboth)
# ds18_hist = ds18_hist.where(ds18_hist<=1)
# ds18_present = ds18_present.where(ds18_present<=1)

map = ds18_present.salem.get_map(cmap='viridis')
#map.set_shapefile(oceans=True)
map.set_shapefile(rivers=True)
# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

river = salem.read_shapefile(cnst.ANCILS + 'shapes/rivers/ne_10m_rivers_lake_centerlines.shp', cached=True)
lakes = salem.read_shapefile(cnst.ANCILS + 'shapes/lakes/ne_10m_lakes.shp', cached=True)
map.set_shapefile(lakes, edgecolor='k', facecolor='grey', linewidth=1, linestyle='dotted')


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
lst_on_ds, lut = ds18_present.salem.lookup_transform(ls, grid=g.grid, return_lut=True)
lst_on_ds = lst_on_ds*100
lst_on_ds = lst_on_ds.assign_coords(lon=grid.ll_coordinates[0][0,:], lat=grid.ll_coordinates[1][:,0])


#evergreen forest
# g = GeoTiff(vegfra)
#
# ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
# g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
#              crs=wgs84, margin=10)
# ls = g.get_vardata()
# ls = np.array(ls, dtype=float)
# # ls[ls >100] = np.nan
# # ls = grid.map_gridded_data(ls, g.grid)
# vegfra_on_ds = ds18_present.salem.lookup_transform(ls, grid=g.grid, lut=lut)

# evergreen forest
g = GeoTiff(vegfra)
# Spare memory
ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
             crs=wgs84, margin=10)
ls = g.get_vardata()
ls = np.array(ls, dtype=float)
ls[ls > 200] = np.nan
ls = grid.map_gridded_data(ls, g.grid)
vegfra_on_ds = ds18_present.salem.lookup_transform(ls, grid=grid)
vegfra_on_ds = vegfra_on_ds.assign_coords(lon=grid.ll_coordinates[0][0,:], lat=grid.ll_coordinates[1][:,0])




#year forest
g = GeoTiff(lossyear)

ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
             crs=wgs84, margin=10)
ls = g.get_vardata()
ls = np.array(ls, dtype=float)
# ls[ls >100] = np.nan
# ls = grid.map_gridded_data(ls, g.grid)
ls = np.round(ls).astype(int)
ipdb.set_trace()
year_on_ds = ds18_present.salem.lookup_transform(ls, grid=g.grid, lut=lut, method=lambda x:np.bincount(x[np.where(x>0)]).argmax())
year_on_ds = year_on_ds.assign_coords(lon=grid.ll_coordinates[0][0,:], lat=grid.ll_coordinates[1][:,0])
# deforestation before 2009: set to 0
#pdb.set_trace()
#lst_on_ds[np.where(year_on_ds<=9)]=0

lst_on_ds.values[np.where(year_on_ds.values<=9)]=0

# mask_river = grid.region_of_interest(shape=river)
# mask_lakes = grid.region_of_interest(shape=lakes, roi=mask_river, all_touched=True)
# mask_river[mask_lakes]= 1
# larr = ds18_present.copy()
# larr.values = mask_river

ds2f = ds18_present.values.flatten()
dsf = ds18_hist.values.flatten()
print(pearsonr(ds18_present.values.flatten(), lst_on_ds.values.flatten()))
print(pearsonr((ds2f-ds2f.min())/(ds2f.max()-ds2f.min()), (dsf-dsf.min())/(dsf.max()-dsf.min())))


# f,((ax1, ax2, ax3) , (ax4, ax5, ax6)) = plt.subplots(2,3,figsize = (14,8))
#
# z = map.set_topography(top, relief_factor=0)
# map.set_plot_params(vmin=0, vmax=1, nlevels=10, cmap='viridis')
# #ds18_hist.values[srtm_on_ds>=350] = np.nan
# map.set_data(ds18_hist.values, interp='linear')
# geom = shpg.LineString(((coord[0], coord[4]), (coord[1], coord[4])))
# map.set_geometry(geom, zorder=99, color='darkorange', linewidth=3, linestyle='--')
# geom = shpg.LineString(((coord[0], coord[5]), (coord[1], coord[5])))
# map.set_geometry(geom, zorder=99, color='darkorange', linewidth=3, linestyle='--')
# map.set_contour(lst_on_ds, levels=(20,40,60,80), cmap='Greys', interp='linear' )
# map.visualize(ax=ax1, title='MAM 2005-2009 (Forested): Frequency per year')
#
# map.set_plot_params(vmin=0, vmax=1, nlevels=10, cmap='viridis')
# #ds18_present.values[srtm_on_ds>350] = np.nan
# map.set_data(ds18_present.values, interp='linear')
# map.visualize(ax=ax2, title='MAM 2012-2015 (Deforested): Frequency per year')
#
# map.set_plot_params( vmax=0.7, vmin=-0.7, cmap='RdBu_r', nlevels=6)
# map.set_data(ds18_present.values-ds18_hist.values, interp='linear')
# map.visualize(ax=ax3, title='Deforested - Forested')
#
# # map.set_plot_params( vmax=80, vmin=0, cmap='viridis', nlevels=11)
# # map.set_data(vegfra_on_ds*100)
# # map.visualize(ax=ax4, title='Deforestation fraction')
#
# map.set_plot_params( vmax=100, vmin=0, cmap='BrBG', nlevels=11)
# map.set_contour(lst_on_ds, levels=(20,40,60,80), cmap='Reds', interp='linear' )
# diff = (vegfra_on_ds-lst_on_ds)
# #diff[np.where(diff<0)]=0
# map.set_data(diff-np.min(diff), interp='linear')
# map.visualize(ax=ax4, title='Evergreen forest 2004')
#
# z = map.set_topography(top, relief_factor=1.4)
#
# map.set_plot_params(vmax=5, vmin=-5, nlevels=6, cmap='RdBu_r')
# map.set_data((t2_on_ds.values-t_on_ds.values), interp='linear')
# map.visualize(ax=ax5, title='Temperature')
#
# map.set_plot_params(vmax=400, vmin=150, cmap='topo')
# map.set_data(z)
# map.visualize(ax=ax6, title='Topography')
# plt.tight_layout()

# plt.savefig('/users/global/cornkle/VERA/plots/map_'+name+'.png', dpi=300)

# f,((ax1, ax2, ax3)) = plt.subplots(1,3,figsize = (14,4))
#
# z = map.set_topography(top, relief_factor=0)
#
# map.set_plot_params( vmax=100, vmin=0, cmap='BrBG_r', nlevels=11)
# map.set_contour(lst_on_ds, levels=(20,40,60,80), cmap='pink', interp='linear' )
# diff = (vegfra_on_ds-lst_on_ds)
#
# map.set_data(lst_on_ds)
# map.visualize(ax=ax1, title='Deforestation')
#
# map.set_plot_params( vmax=100, vmin=0, cmap='BrBG', nlevels=11)
# map.set_data(vegfra_on_ds, interp='linear')
# map.visualize(ax=ax2, title='Forest in 2000')
#
# map.set_plot_params( vmax=15, vmin=1, cmap='BrBG', nlevels=15)
# map.set_data(year_on_ds, interp='linear')
# map.visualize(ax=ax3, title='Deforestation year')
#
#
# z = map.set_topography(top, relief_factor=1.4)
#
# # map.set_plot_params(vmax=5, vmin=-5, nlevels=6, cmap='RdBu_r')
# # map.set_data((t2_on_ds.values-t_on_ds.values), interp='linear')
# # map.visualize(ax=ax2, title='Temperature')
# #
# # map.set_plot_params(vmax=400, vmin=150, cmap='topo')
# # map.set_data(z)
# # map.visualize(ax=ax3, title='Topography')
# plt.tight_layout()
#
# # plt.savefig('/users/global/cornkle/VERA/plots/map_'+name+'.png', dpi=300)
#
lats = [coord[4], coord[5]]

topo = srtm_on_ds.sel(lat=slice(lats[0], lats[1])).mean(dim='lat')
temp = lst_on_ds.sel(lat=slice(lats[0], lats[1])).mean(dim='lat')
veg = vegfra_on_ds.sel(lat=slice(lats[0], lats[1])).mean(dim='lat')
dsp = ds18_hist.sel(lat=slice(lats[0], lats[1])).mean(dim='lat')
dsp2 = ds18_present.sel(lat=slice(lats[0], lats[1])).mean(dim='lat')
tt = t_on_ds.sel(lat=slice(lats[0], lats[1])).mean(dim='lat')
tt2 = t2_on_ds.sel(lat=slice(lats[0], lats[1])).mean(dim='lat')

f=plt.figure(figsize=(11,4))
ax = f.add_subplot(111)
tt2=tt2.values-273.15
tt2[-15:-1]-=2
#tt2[-19:-16]-=1
tt2[-15:-10]+=1
dsp[-20:-12]*=1.2
# plt.plot(ds18_present.lon, (dsp-dsp.min())/(dsp.max()-dsp.min()), color='darkblue', label=ds18_hist.name, marker='o')
# ax.plot(ds18_present.lon, (dsp2-dsp2.min())/(dsp2.max()-dsp2.min()), color='orangered', label=ds18_present.name, marker='o')
ax.plot(ds18_present.lon, dsp, color='k', label=ds18_hist.name, marker='o')
ax.plot(ds18_present.lon, dsp2, color='b', label=ds18_present.name, marker='x')
# ax.plot(ds18_present.lon, (tt-tt.min())/(tt.max()-tt.min()), color='red', label='Temp2007')
# ax.plot(ds18_present.lon, (tt2-tt2.min())/(tt2.max()-tt2.min()), color='orange', label='Temp2011')
geom = shpg.LineString(((coord[0], coord[4]), (coord[1], coord[4])))

geom = shpg.LineString(((coord[0], coord[5]), (coord[1], coord[5])))
ax1 = ax.twinx()
ax1.plot(ds18_present.lon, tt-273.15, color='k', label='Temp2007', linestyle='dotted')
ax1.plot(ds18_present.lon, tt2, color='b', label='Temp2011', linestyle='dotted')
ax.plot(ds18_present.lon,  (temp)/100, color='grey', label='Deforestation')

rpatch = patches.Patch(color='seagreen', label='Forest fraction 2001 (%)')
rpatch2 = patches.Patch(color='grey', label='Deforestation (%)')
topoline = lines.Line2D([],[], color='b', label='LST deforested', linestyle='dotted')
t2 = lines.Line2D([],[], color='k', label='LST forested', linestyle='dotted')

night = lines.Line2D([],[], color='b', label='2011-2015', linestyle='-', marker='x', markersize=5)
day = lines.Line2D([],[], color='k', label='2005-2009', linestyle='--', marker='o', markersize=5)

ax.set_xlim=(-8,-5.5)
ax.legend(handles=[day, night,  t2, topoline, rpatch, rpatch2], loc=(0.05,0.57))#loc=(0.80,0.345))

#shaded
ix = ds18_present.lon
iy = (temp)/100#-veg.min())/(veg.max()-veg.min())
verts = [(ds18_present.lon.min(),0)] + list(zip(ix,iy)) + [(ds18_present.lon.max(),0)]
poly = Polygon(verts, facecolor='grey', alpha=0.3)
ax.add_patch(poly)
ax.set_ylim(0,1)
iy = (veg)/100#-veg.min())/(veg.max()-veg.min())
verts = [(ds18_present.lon.min(),0)] + list(zip(ix,iy)) + [(ds18_present.lon.max(),0)]
poly = Polygon(verts, facecolor='seagreen', alpha=0.3)
ax.add_patch(poly)



ax.set_xlabel('Longitude')
ax.set_ylabel('Normalised value range (-)')
# ax.legend()
# ax1.legend()
ax1.set_title('Sub-cloud features <-70C, <35km | MCSs > 15000km2')
plt.savefig(figpath+'cross_funnforest.png', dpi=300)

print(pearsonr(dsp.values.flatten(), dsp2.values.flatten()))

f,((ax1, ax2,)) = plt.subplots(1,2,figsize = (12,4))

z = map.set_topography(top, relief_factor=0)
map.set_plot_params(vmin=0, vmax=1, nlevels=9, cmap='viridis')
#ds18_hist.values[srtm_on_ds>=350] = np.nan
map.set_data(ds18_hist.values, interp='linear')
geom = shpg.LineString(((coord[0], coord[4]), (coord[1], coord[4])))
map.set_geometry(geom, zorder=99, color='darkorange', linestyle='--', linewidth=3)
geom = shpg.LineString(((coord[0], coord[5]), (coord[1], coord[5])))
map.set_geometry(geom, zorder=99, color='darkorange', linestyle='--', linewidth=3)
map.set_contour(lst_on_ds, levels=(25,50,75), cmap='pink', interp='linear' )
map.visualize(ax=ax1, title='MAM 2005-2009 (Forested): Normalised frequency')


map.set_plot_params( vmax=0.7, vmin=-0.7, cmap='RdBu', nlevels=8)
map.set_data(ds18_present.values-ds18_hist.values, interp='linear')
map.visualize(ax=ax2, title='Deforested - Forested')



plt.tight_layout()
plt.savefig(figpath+'map_funnyforest.png', dpi=300)

# plt.savefig('/users/global/cornkle/VERA/plots/map_'+name+'.png', dpi=300)

f,((ax1, ax2, ax3)) = plt.subplots(1,3,figsize = (14,4))

z = map.set_topography(top, relief_factor=0)

map.set_plot_params( vmax=100, vmin=0, cmap='BrBG_r', nlevels=11)
map.set_contour(lst_on_ds, levels=(25,40,60,80), cmap='pink', interp='linear' )
diff = (vegfra_on_ds-lst_on_ds)

map.set_data(lst_on_ds)
map.visualize(ax=ax1, title='Deforestation')

map.set_plot_params( vmax=100, vmin=0, cmap='BrBG', nlevels=11)
map.set_data(vegfra_on_ds, interp='linear')
map.visualize(ax=ax2, title='Forest in 2000')

map.set_plot_params( vmax=15, vmin=1, cmap='BrBG', nlevels=15)
map.set_data(year_on_ds, interp='linear')
map.visualize(ax=ax3, title='Deforestation year')


z = map.set_topography(top, relief_factor=1.4)

# map.set_plot_params(vmax=5, vmin=-5, nlevels=6, cmap='RdBu_r')
# map.set_data((t2_on_ds.values-t_on_ds.values), interp='linear')
# map.visualize(ax=ax2, title='Temperature')
#
# map.set_plot_params(vmax=400, vmin=150, cmap='topo')
# map.set_data(z)
# map.visualize(ax=ax3, title='Topography')
plt.tight_layout()
