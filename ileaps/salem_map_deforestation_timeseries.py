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
path = '/prj/vera/cores/'
figpath = cnst.network_data + 'figs/VERA/chris_egu/'

hours = '18-19'
tstring = '-73'
#file18 = path+'blob_map_35km_'+tstring+'_MAMJ_16-19UTC.nc'#blob_map_35km_-70_15-18UTC.nc' #blob_map_35km_-75_sum_0-3UTC.nc'

fpath = cnst.ANCILS + 'gtopo_1min_afr.nc'
vegfra = cnst.local_data + 'obs_data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_treecover2000_10N_010W.tif'
lst = cnst.local_data + 'obs_data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_loss_10N_010W.tif'
lossyear = cnst.local_data + 'obs_data/LandCover/landsat_forest/Hansen_GFC-2015-v1.3_lossyear_10N_010W.tif'
vegfra = cnst.local_data + 'obs_data/LandCover/evergreen_trees.tif'

file_out = figpath+'blob_map_35km_hist_MAM_16-17UTC_2006.nc'
file_in = path+'coresPower_*'
ds18 = xr.open_mfdataset(file_in)

ds18 = ds18['blobs'][(ds18['time.month']<= 5) & (ds18['time.hour']>= 16) & (ds18['time.hour']<= 18)]
ds18.values = (ds18.values > 0).astype(int)
ts = ds18.sel(lon=slice(-7, -6.9), lat=slice(6.9,7.2)).groupby('time.year').mean()
plt.figure()
plt.plot(np.unique(ds18['time.year'].values), ts)
