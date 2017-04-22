import salem
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np


file = '/users/global/cornkle/MCSfiles/blob_map_90km_sum_18UTC.nc'
file2 = '/users/global/cornkle/MCSfiles/blob_map_30km_sum_18UTC.nc'
file3 = '/users/global/cornkle/MCSfiles/blob_map_90km_sum_3UTC.nc'
file4 = '/users/global/cornkle/MCSfiles/blob_map_30km_sum_3UTC.nc'
fpath = '/users/global/cornkle/data/pythonWorkspace/proj_CEH/topo/gtopo_1min_afr.nc'
spath = '/users/global/cornkle/C_paper/wavelet/figs/paper/'

ds = xr.open_dataarray(file)
top = xr.open_dataarray(fpath)
ds2 = xr.open_dataarray(file2)
ds3 = xr.open_dataarray(file3)
ds4 = xr.open_dataarray(file4)



ds.name = '100k'
ds2.name = '30k'

ds = ds.sel(lon=slice(-17,30), lat=slice(4.5,20))  # lake chad lon=slice(10,20), lat=slice(10,15)
ds2 = ds2.sel(lon=slice(-17,30), lat=slice(4.5,20))   # volta lon=slice(-10,8), lat=slice(4,10)
ds3 = ds3.sel(lon=slice(-17,30), lat=slice(4.5,20))
ds4 = ds4.sel(lon=slice(-17,30), lat=slice(4.5,20))
top = top.sel(lon=slice(-17,30), lat=slice(4.5,20))

ds[ds == 0]=np.nan
ds2[ds2 == 0] =np.nan
ds3[ds3 == 0] =np.nan
ds4[ds4 == 0] =np.nan

perc = ds.quantile(0.95)
perc2 =ds2.quantile(0.95)
perc3 =ds3.quantile(0.95)
perc4 =ds4.quantile(0.95)

ds = (ds-1) / (perc- 1)  # dim=['lon']
ds2 = (ds2-1) / (perc2- 1)
ds3 = (ds3-1) / (perc3- 1)
ds4 = (ds4-1) / (perc4- 1)

ds = ds.where(ds<=1)
ds2 = ds2.where(ds2<=1)
ds3 = ds3.where(ds2<=1)
ds4 = ds4.where(ds2<=1)

# ds[ds.where(ds<=1)]=-999
# ds2[ds2 > 1]=-999
# ds3[ds3 > 1]=-999
# ds4[ds4 > 1]=-999

map = ds.salem.get_map(cmap='Greys')
map.set_shapefile(rivers=True)
# read the ocean shapefile (data from http://www.naturalearthdata.com)
oceans = salem.read_shapefile(salem.get_demo_file('ne_50m_ocean.shp'),
                              cached=True)

lakes = salem.read_shapefile(salem.get_demo_file('ne_50m_lakes.shp'), cached=True)
map.set_shapefile(lakes, edgecolor='k', facecolor='none', linewidth=2,)

grid = ds.salem.grid


f,((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4,2,figsize = (18,15))


map.set_data(ds)
map.set_plot_params(vmin=0., vmax=0.7, cmap='Greys')

map.visualize(ax=ax1, title='>90km 1800UTC', addcbar=False)

map.set_data(ds2)
map.visualize(ax=ax2, title='<30km 1800UTC', addcbar=False)

kw1 = map.get_colorbarbase_kwargs()



map.set_data(ds3)
map.visualize(ax=ax3, title='>90km 0300UTC')

map.set_data(ds4)
map.visualize(ax=ax4, title='<30km 0300UTC')




map.set_plot_params(levels=[-0.6,-0.4,-0.2,0.2,0.4,0.6], cmap='RdBu')
map.set_data(ds - ds3)
map.visualize(ax=ax5, title='>90km 1800UTC - 0300UTC', addcbar=False)

map.set_data(ds2-ds4)
map.visualize(ax=ax6, title='<30km 1800UTC - 0300UTC', addcbar=False)


zuse = map.set_topography(top, relief_factor=1.4)
map.set_plot_params(vmax=1500, cmap='topo')
map.set_data(zuse)
map.visualize(ax=ax7, title='<30km 0300UTC')

ax8.axis('off')


fsiz = 18
x = 0.02
plt.annotate('a)', xy=(x, 0.98), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
plt.annotate('b)', xy=(x, 0.73), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
plt.annotate('c)', xy=(x, 0.49), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')
plt.annotate('d)', xy=(x, 0.25), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
                 textcoords='offset points')


plt.tight_layout()

f.subplots_adjust(right=0.94)
cax = f.add_axes([0.96,0.5,0.025,0.4])
mpl.colorbar.ColorbarBase(cax, **kw1)

plt.savefig(spath+'/scales_map.png', dpi=400)

plt.close('all')
