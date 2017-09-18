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


#coord = [-8.5,-1.5,4.8,7.8,6.9,7.4, -7.6, -7.4]
name='Tai park'

coord=[-8,-1,4.8,8,5.5,6,7.9,8]
# coord=[-8,-1,5.5,6,5.5,6,7.9,8]
# name='whatever'

# coord=[-3,3,6,9.5,6.4,8,7.9,8]
# name='volta'

ds18_hist = ds18_hist.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
ds18_present = ds18_present.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))
top = top.sel(lon=slice(coord[0],coord[1]), lat=slice(coord[2],coord[3]))

ds18_hist.name = '2004-2009'
ds18_present.name = '2011-2015'



map = ds18_present.salem.get_map(cmap='viridis')
#map.set_shapefile(oceans=True)
map.set_shapefile(rivers=True)
# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

river = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/rivers/ne_10m_rivers_lake_centerlines.shp', cached=True)
lakes = salem.read_shapefile('/users/global/cornkle/data/pythonWorkspace/proj_CEH/shapes/lakes/ne_10m_lakes.shp', cached=True)
map.set_shapefile(lakes, edgecolor='k', facecolor='lightblue', linewidth=1, linestyle='dotted')

if not os.path.isfile(path+'forest.nc'):

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
            b = np.bincount(a[np.where(a>0)]).argmax()
        except ValueError:
            b = np.nan
        return b

    year_on_ds = ds18_present.salem.lookup_transform(ls, grid=g.grid, lut=lut, method=dyear)
    print('Passed year_on_ds')
    # deforestation before 2009: set to 0
    #pdb.set_trace()
    #lst_on_ds[np.where(year_on_ds<=9)]=0



    ds = xr.Dataset({
                     'deforestation': (['lat', 'lon'], lst_on_ds),
                     'forest2000': (['lat', 'lon'], vegfra_on_ds),
                     'dyear': (['lat', 'lon'], year_on_ds),
                     },
                    coords=ds18_present.coords)

    ds.to_netcdf(path + 'forest.nc')

else:
    srfc = xr.open_dataset(path + 'forest.nc')
    #
    # coord = [-8, -5.5, 5.1, 8, 5.5, 6, 7.9, 8]
    # srfc = srfc.sel(lon=slice(coord[0], coord[1]), lat=slice(coord[2], coord[3]))

    lst_on_ds = srfc['deforestation']
    vegfra_on_ds = srfc['forest2000']
    year_on_ds = srfc['dyear']

lst_on_ds.values[year_on_ds.values<=9]=0

f,((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize = (8,4))

z = map.set_topography(top, relief_factor=0)

map.set_plot_params( vmax=100, vmin=0, cmap='BrBG_r', nlevels=11)
map.set_contour(lst_on_ds, levels=(15,20,40,60,80), cmap='jet', interp='linear' )
diff = (vegfra_on_ds-lst_on_ds)

map.set_data(lst_on_ds, interp='linear')
map.visualize(ax=ax1, title='Deforestation')

map.set_plot_params( vmax=100, vmin=0, cmap='BrBG', nlevels=11)
map.set_data(vegfra_on_ds)
map.visualize(ax=ax2, title='Forest in 2000')

map.set_plot_params( vmax=15, vmin=1, cmap='BrBG', nlevels=15)
map.set_data(year_on_ds, interp='linear')
map.visualize(ax=ax3, title='Deforestation year')


ax4.axis('off')