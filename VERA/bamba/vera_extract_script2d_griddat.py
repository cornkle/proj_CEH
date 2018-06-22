import xarray as xr
import pdb
import numpy as np
import glob
import os
from utils import u_interpolate as uint


#veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'

### 2d vars , xmh*.pc*.nc files
folder = '/scratch/ssf/'
out = '/users/global/cornkle/w2018_bamba/clean_files_raw/'

box = [130, 400, 5, 100 ]  # x1, x2, y1, y2 Tai park box

flist = glob.glob(folder+'current.pc*.nc')

dummy = xr.open_dataset(flist[0])
#dummy = dummy.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))

lats = dummy.latitude_t[:, 0].values
lons = dummy.longitude_t[0, :].values
latmin = dummy.latitude_t.min().values
latmax = dummy.latitude_t.max().values
lonmin = dummy.longitude_t.min().values
lonmax = dummy.longitude_t.max().values
dist = np.round(np.float(np.mean(lats[1::] - lats[0:-1])), decimals=4)

lat_regular = np.arange(latmin, latmax, dist)
lon_regular = np.arange(lonmin, lonmax, dist)

# # interpolation weights for new grid
# inds, weights, shape = uint.interpolation_weights(lons, lats, lon_regular, lat_regular)

varsdat = dummy
varsdat = varsdat.isel(TH1_MN=5, TH1_MN_rad=5)

vkeys = varsdat.keys()
# for key in vkeys:
key = 'SH'
data = varsdat[key].values

regridded = uint.griddata_comparison(data,dummy.longitude_t.values, dummy.latitude_t.values, lon_regular, lat_regular, isll=True)
ds = xr.Dataset(coords={'lat': lat_regular, 'lon': lon_regular})
ds['SH']= (('lat', 'lon'), regridded)
ds.to_netcdf(out+'test.nc')




