import salem
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import pdb


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

ds = ds.sel(lon=slice(-17,20), lat=slice(4.5,20))  # lake chad lon=slice(10,20), lat=slice(10,15)
ds2 = ds2.sel(lon=slice(-17,20), lat=slice(4.5,20))   # volta lon=slice(-10,8), lat=slice(4,10)
ds3 = ds3.sel(lon=slice(-17,20), lat=slice(4.5,20))
ds4 = ds4.sel(lon=slice(-17,20), lat=slice(4.5,20))
top = top.sel(lon=slice(-17,20), lat=slice(4.5,20))

ds[ds == 0]=np.nan
ds2[ds2 == 0] =np.nan
ds3[ds3 == 0] =np.nan
ds4[ds4 == 0] =np.nan

perc = ds.quantile(0.99)
perc2 =ds2.quantile(0.99)
perc3 =ds3.quantile(0.99)
perc4 =ds4.quantile(0.99)

# perc = np.max(ds)
# perc2 =np.max(ds2)
# perc3 =np.max(ds3)
# perc4 =np.max(ds4)

percc = np.max([perc,perc3])
percc1 = np.max([perc2, perc4])
#
# ds = (ds-1) / (percc- 1)  # dim=['lon']
# ds2 = (ds2-1) / (percc1- 1)
# ds3 = (ds3-1) / (percc- 1)
# ds4 = (ds4-1) / (percc1- 1)

ds = (ds-1) / (perc- 1)  # dim=['lon']
ds2 = (ds2-1) / (perc2- 1)
ds3 = (ds3-1) / (perc3- 1)
ds4 = (ds4-1) / (perc4- 1)

# ds = ds.where(ds<=1)
# ds2 = ds2.where(ds2<=1)
# ds3 = ds3.where(ds3<=1)
# ds4 = ds4.where(ds4<=1)

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


#f,((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4,2,figsize = (18,15))

f,((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize = (11,6))


# map.set_data(ds)
# map.set_plot_params(vmin=0., vmax=0.7, cmap='Blues')
#
# map.visualize(ax=ax1, title='>90km 1800UTC', addcbar=False)
#
# map.set_data(ds2)
# map.visualize(ax=ax2, title='<30km 1800UTC', addcbar=False)
#
# kw1 = map.get_colorbarbase_kwargs()
#
#
#
# map.set_data(ds3)
# map.visualize(ax=ax3, title='>90km 0300UTC', addcbar=False)
#
# map.set_data(ds4)
# map.visualize(ax=ax4, title='<30km 0300UTC', addcbar=False)




map.set_plot_params(levels=[-0.6,-0.4,-0.2,0.2,0.4,0.6], cmap='RdBu')
dat = ds.values - ds3.values
print(np.std(dat[np.isfinite(dat)]))
dat[dat<0] = dat[dat<0]-0.03

map.set_data(dat)
map.visualize(ax=ax1, title='>90km 1800UTC - 0300UTC', addcbar=False)

dat = ds2.values - ds4.values

dat[dat>0] = dat[dat>0]+0.03
dat[dat<0] = dat[dat<0]-0.02
print(np.std(dat[np.isfinite(dat)]))
map.set_data(dat)
map.set_lonlat_contours(add_ytick_labels=False)
map.visualize(ax=ax2, title='<35km 1800UTC - 0300UTC', addcbar=False)
kw2 = map.get_colorbarbase_kwargs()

for tick in ax1.xaxis.get_major_ticks():
                tick.label.set_fontsize(9)
for tick in ax2.xaxis.get_major_ticks():
                tick.label.set_fontsize(9)
for tick in ax1.yaxis.get_major_ticks():
                tick.label.set_fontsize(9)
for tick in ax2.yaxis.get_major_ticks():
                tick.label.set_fontsize(9)


zuse = map.set_topography(top, relief_factor=1.4)
map.set_lonlat_contours(add_ytick_labels=True)
map.set_plot_params(vmax=2000, cmap='topo')
map.set_data(zuse)
map.visualize(ax=ax3, title='Topography', addcbar=False)
for tick in ax3.yaxis.get_major_ticks():
                tick.label.set_fontsize(9)
for tick in ax3.xaxis.get_major_ticks():
                tick.label.set_fontsize(9)
kw = map.get_colorbarbase_kwargs()

ax4.axis('off')


# fsiz = 18
# x = 0.02
# plt.annotate('a)', xy=(x, 0.98), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
#                  textcoords='offset points')
# plt.annotate('b)', xy=(x, 0.73), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
#                  textcoords='offset points')
# plt.annotate('c)', xy=(x, 0.49), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
#                  textcoords='offset points')
# plt.annotate('d)', xy=(x, 0.25), xytext=(0, 4), size=fsiz, xycoords=('figure fraction', 'figure fraction'),
#                  textcoords='offset points')


plt.tight_layout()

f.subplots_adjust(right=0.89)
# cax = f.add_axes([0.95,0.55,0.025,0.4])
# salem.DataLevels(ds, nlevels=8)
# #kw1=salem.DataLevels.colorbarbase(cax, **kw1)
# mpl.colorbar.ColorbarBase(cax, **kw1)

cax = f.add_axes([0.90,0.57,0.017,0.35])
dl=salem.DataLevels(ds, nlevels=8)
#kw1=salem.DataLevels.colorbarbase(cax, **kw1)
cb1 = mpl.colorbar.ColorbarBase(cax, **kw2, label = 'Nocturnal   |   Afternoon')
#dl.set_extend(extend='both')
#cb1.set_ticklabels(['']*6)
#f.colorbar(cax).set_yticklabels(['','','','','',''])


cax = f.add_axes([0.46,0.14,0.017,0.30])
mpl.colorbar.ColorbarBase(cax, **kw, label='m')

plt.savefig(spath+'/scales_map.png', dpi=300)

plt.close('all')
