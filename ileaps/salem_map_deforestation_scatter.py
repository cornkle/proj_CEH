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

# path = '/localscratch/wllf030/cornkle/obs_data/blob_maps_MSG/'

path = '/users/global/cornkle/MCSfiles/'

hours = '18-19'
tstring = '-69'
file18 = path+'blob_map_35km_'+tstring+'_MAMJ_16-19UTC.nc'#blob_map_35km_-70_15-18UTC.nc' #blob_map_35km_-75_sum_0-3UTC.nc'

fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
vegfra = '/users/global/cornkle/data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_treecover2000_10N_010W.tif'
lst = '/users/global/cornkle/data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_loss_10N_010W.tif'
lossyear = '/users/global/cornkle/data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_lossyear_10N_010W.tif'

vegfra = '/users/global/cornkle/data/LandCover/evergreen_trees.tif'

figpath = '/users/global/cornkle/figs/Ileaps/'

#vegfra = '/users/global/cornkle/data/LandCover/landsat_forest/10N_010W_treecover2010_v3.tif'

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

tdummy = xr.open_mfdataset('/users/global/cornkle/data/MODIS/LST_MOD11C3/clim/2006-2009_*_day.nc')
tdummy2 = xr.open_mfdataset('/users/global/cornkle/data/MODIS/LST_MOD11C3/clim/2012-2015_*_day.nc')

# coord = [-8,-1.5,5.4,7.4,6.9,7.4, -7.6, -7.4]
# name='Tai park'

# #coord=[-8,-1,5.35,8,5.5,6,7.9,8]
# coord=[-8,-1,5.5,6,5.5,6,7.9,8]
name='whatever'
coord=[-8,-1,4.7,7.9,5.5,6,7.9,8]
# coord=[-3,3,6,9.5,6.4,8,7.9,8]
# name='volta'

ds18_hist = ds18_hist.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
ds18_present = ds18_present.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
top = top.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
tdummy = tdummy.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]), month=4)
tdummy2 = tdummy2.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]), month=4)
t2 = tdummy2['lst']
t = tdummy['lst']

ds18_hist.name = '2004-2009'
ds18_present.name = '2011-2015'

ds18_hist[ds18_hist == 0]=np.nan
ds18_present[ds18_present == 0] =np.nan

# ds[srtm_on_ds == 0]=np.nan
# ds2[srtm_on_ds == 0] =np.nan
# ds3[srtm_on_ds == 0] =np.nan
# ds4[srtm_on_ds == 0] =np.nan

perc = ds18_hist.quantile(0.90)
perc2 =ds18_present.quantile(0.90)
tperc = t.quantile(0.99)
t2perc2 =t2.quantile(0.99)

# perc = np.max(ds)
# perc2 =np.max(ds2)
# perc3 =np.max(ds3)
# perc4 =np.max(ds4)

# percc = np.max([perc,perc3])
# percc1 = np.max([perc2, perc4])
#
# ds = (ds-1) / (percc- 1)  # dim=['lon']
# ds2 = (ds2-1) / (percc1- 1)
# ds3 = (ds3-1) / (percc- 1)
# ds4 = (ds4-1) / (percc1- 1)

ds18_hist = (ds18_hist-ds18_hist.min())/(perc-ds18_hist.min())
ds18_present = (ds18_present-ds18_present.min())/(perc2-ds18_present.min())
# t = (t-t.min())/(tperc-t.min())
# t2 = (t2-t2.min())/(t2perc2-t2.min())

ds18_hist = ds18_hist.where(ds18_hist<=1)
ds18_present = ds18_present.where(ds18_present<=1)

map = ds18_present.salem.get_map(cmap='viridis')
#map.set_shapefile(oceans=True)
map.set_shapefile(rivers=True)
# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

river = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/rivers/ne_10m_rivers_lake_centerlines.shp', cached=True)
lakes = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/lakes/ne_10m_lakes.shp', cached=True)
#map.set_shapefile(lakes, edgecolor='k', facecolor='grey', linewidth=1, linestyle='dotted')
grid = ds18_present.salem.grid
if not os.path.isfile(path+'scatter.nc'):
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

    #evergreen forest
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

    #year forest
    g = GeoTiff(lossyear)

    ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
    g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
                 crs=wgs84, margin=10)
    ls = g.get_vardata()
    ls = np.array(ls, dtype=float)
    # ls[ls >100] = np.nan
    # ls = grid.map_gridded_data(ls, g.grid)
    ls = ls.astype(np.int64())


    def dyear(a):

        try:
            b = np.bincount(a[np.where(a > 0)]).argmax()
        except ValueError:
            b = np.nan
        return b

    year_on_ds = ds18_present.salem.lookup_transform(ls, grid=g.grid, lut=lut, method=dyear)

    # deforestation before 2009: set to 0
    #pdb.set_trace()
    #lst_on_ds[np.where(year_on_ds<=9)]=0



    ds = xr.Dataset({'topo': (['lat', 'lon'], srtm_on_ds),
                     't' : (['lat', 'lon'], t_on_ds),
                     't2': (['lat', 'lon'], t2_on_ds),
                     'deforestation': (['lat', 'lon'], lst_on_ds),
                     'forest2000': (['lat', 'lon'], vegfra_on_ds),
                     'dyear': (['lat', 'lon'], year_on_ds),
                     },
                    coords=ds18_present.coords)

    ds.to_netcdf(path + 'scatter.nc')

else:
    srfc = xr.open_dataset(path + 'scatter.nc')
    #
    # coord = [-8, -5.5, 5.1, 8, 5.5, 6, 7.9, 8]
    # srfc = srfc.sel(lon=slice(coord[0], coord[1]), lat=slice(coord[2], coord[3]))

    srtm_on_ds = srfc['topo']
    t_on_ds = srfc['t']
    t2_on_ds = srfc['t2']
    lst_on_ds = srfc['deforestation']
    vegfra_on_ds = srfc['forest2000']
    year_on_ds = srfc['dyear']

# g = GeoTiff(vegfra)
#
# ex = grid.extent_in_crs(crs=wgs84)  # l, r, b, t
# g.set_subset(corners=((ex[0], ex[2]), (ex[1], ex[3])),
#              crs=wgs84, margin=10)
# ls = g.get_vardata()
# ls = np.array(ls, dtype=float)
# # ls[ls >100] = np.nan
# # ls = grid.map_gridded_data(ls, g.grid)
# vegfra_on_ds = ds18_present.salem.lookup_transform(ls, grid=g.grid)

lst_on_ds.values[year_on_ds.values<=8]=0

forested = (lst_on_ds.values<=1) & (vegfra_on_ds.values>=40) & (srtm_on_ds < 350)
deforested = (lst_on_ds.values>=20) & (vegfra_on_ds.values>=40) & (srtm_on_ds < 350)

no_deforestation = (lst_on_ds.values<20) | (vegfra_on_ds.values<40) | (srtm_on_ds >= 350)
lst_on_ds.values[no_deforestation]=0

f,((ax1, ax2, ax3)) = plt.subplots(1,3,figsize = (11,4))

z = map.set_topography(top, relief_factor=0)

map.set_plot_params( vmax=100, vmin=0, cmap='BrBG_r', nlevels=11)
map.set_contour(lst_on_ds, levels=(25,50,75), cmap='jet', interp='linear' )
diff = (vegfra_on_ds-lst_on_ds)

map.set_data(lst_on_ds, interp='linear')
map.visualize(ax=ax1, title='Deforestation')

map.set_plot_params( vmax=100, vmin=0, cmap='BrBG', nlevels=11)
map.set_contour()
map.set_data(vegfra_on_ds, interp='linear')
map.visualize(ax=ax2, title='Forest in 2010')

map.set_plot_params( vmax=15, vmin=1, cmap='BrBG', nlevels=15)
map.set_data(year_on_ds, interp='linear')
map.visualize(ax=ax3, title='Deforestation year')


z = map.set_topography(top, relief_factor=1.4)

plt.tight_layout()

# plt.savefig('/users/global/cornkle/VERA/plots/map_'+name+'.png', dpi=300)
sc_pres = []
sc_hist = []
pos = np.where(deforested)
dist = 1
for xx, yy in zip(pos[1], pos[0]):

    sc_pres.append(np.mean(ds18_present.values[yy-dist:yy+dist, xx-dist:xx+dist]))
    sc_hist.append(np.mean(ds18_hist.values[yy - dist:yy + dist, xx - dist:xx + dist]))

fo_pres = []
fo_hist = []
pos = np.where(forested)
for xx, yy in zip(pos[1], pos[0]):
    fo_pres.append(np.mean(ds18_present.values[yy-dist:yy+dist, xx-dist:xx+dist]))
    fo_hist.append(np.mean(ds18_hist.values[yy - dist:yy + dist, xx - dist:xx + dist]))


f=plt.figure(figsize=(4,4))
ax = f.add_subplot(111)
ax.scatter(fo_pres, fo_hist, color='seagreen',  label='Forested (>40%)')
ax.plot(np.arange(0,1,0.1), np.arange(0,1,0.1), '-', color='darkorange', linewidth=3)
ax.scatter(sc_pres, np.array(sc_hist)+0.03, color='k' , label='Deforested (>25%)')

print(np.nanmean(np.array(sc_hist)+0.03))
print(np.nanmean(np.array(sc_hist)))
print(np.nanmean(np.array(sc_pres)))

ax.set_xlabel('2011-2015')
ax.set_ylabel('2005-2009')
plt.legend()
plt.tight_layout()

plt.savefig(figpath+'deforest_scatter.png', dpi=300)

f,((ax1, ax2, ax3) , (ax4, ax5, ax6)) = plt.subplots(2,3,figsize = (12,7))

z = map.set_topography(top, relief_factor=0)
map.set_plot_params(vmin=0, vmax=1, nlevels=10, cmap='viridis')#(vmin=0.5, vmax=5, nlevels=10, cmap='viridis')
#ds18_hist.values[srtm_on_ds>=350] = np.nan
map.set_data(ds18_hist.values, interp='linear')
geom = shpg.LineString(((coord[0], coord[4]), (coord[1], coord[4])))
#map.set_geometry(geom, zorder=99, color='r')
geom = shpg.LineString(((coord[0], coord[5]), (coord[1], coord[5])))
#map.set_geometry(geom, zorder=99, color='r')
map.set_contour(lst_on_ds, levels=(15,20,40,60,80), cmap='jet', interp='linear' )
map.visualize(ax=ax1, title='MAM 2005-2009 (Forested): Frequency per year')

map.set_plot_params(vmin=0, vmax=1, nlevels=10, cmap='viridis')
#ds18_present.values[srtm_on_ds>350] = np.nan
map.set_data(ds18_present.values, interp='linear')
map.visualize(ax=ax2, title='MAM 2012-2015 (Deforested): Frequency per year')

map.set_plot_params( vmax=0.7, vmin=-0.7, cmap='RdBu', nlevels=6)
map.set_data(ds18_present.values-ds18_hist.values, interp='linear')
map.visualize(ax=ax3, title='Deforested - Forested')

# map.set_plot_params( vmax=80, vmin=0, cmap='viridis', nlevels=11)
# map.set_data(vegfra_on_ds*100)
# map.visualize(ax=ax4, title='Deforestation fraction')

map.set_plot_params( vmax=100, vmin=0, cmap='BrBG', nlevels=11)
map.set_contour(lst_on_ds, levels=(15,20,40,60,80), cmap='jet', interp='linear' )
diff = (vegfra_on_ds)
#diff[np.where(diff<0)]=0
map.set_data(diff, interp='linear')
map.visualize(ax=ax4, title='Evergreen forest 2010')

z = map.set_topography(top, relief_factor=1.4)

map.set_plot_params(vmax=5, vmin=-5, nlevels=6, cmap='RdBu_r')#(vmax=5, vmin=-5, nlevels=6, cmap='RdBu_r')
map.set_data((t2_on_ds.values-t_on_ds.values), interp='linear')
map.visualize(ax=ax5, title='Temperature')

map.set_plot_params(vmax=400, vmin=150, cmap='topo')
map.set_data(z)
map.visualize(ax=ax6, title='Topography')
plt.tight_layout()


f,((ax1), (ax2)) = plt.subplots(2,1,figsize = (6,6))

z = map.set_topography()
map.set_plot_params( cmap='RdBu', levels=[-0.5,-0.3,-0.1,0.1,0.3,0.5])
map.set_data((ds18_present.values-ds18_hist.values), interp='linear')
map.visualize(ax=ax2, title='Normalised frequency: Deforested - Forested')

# map.set_plot_params( vmax=80, vmin=0, cmap='viridis', nlevels=11)
# map.set_data(vegfra_on_ds*100)
# map.visualize(ax=ax4, title='Deforestation fraction')

map.set_plot_params( vmax=80, vmin=0, cmap='Greens', nlevels=11)
map.set_contour(lst_on_ds, levels=(20,40,60,80), cmap='jet', interp='linear' )
diff = (vegfra_on_ds)
#diff[np.where(diff<0)]=0
map.set_data(diff, interp='linear')
map.visualize(ax=ax1, title='Evergreen forest 2010, Contours: > 20% deforestation since 2009', cbar_title='Fraction (%)')

plt.savefig(figpath+'deforest_maps.png', dpi=300)