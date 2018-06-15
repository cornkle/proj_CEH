import xarray as xr
import pdb
import matplotlib.pyplot as plt
import numpy as np
from utils import constants as cnst, u_arrays, u_grid
import salem
import glob



veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'
topo = cnst.WA_TOPO_1MIN
### 2d vars , xmh*.pc*.nc files
folder = '/users/global/cornkle/w2018_bamba/linked/'

flist = glob.glob(folder+'*.nc')

name_dict ={
    'STASH_m01s01i201_2' : 'sw_in_net',
   'STASH_m01s01i235_2' : 'sw_in',
   'STASH_m01s02i201_2' : 'lw_in_net',
   'STASH_m01s02i207_2' : 'lw_in',
   'STASH_m01s03i217_2' : 'SH',
   'STASH_m01s03i234_2' : 'LH',
   'STASH_m01s05i205_2' : 'ConvRain',
   'STASH_m01s05i216_2' : 'TotalRain',
   'STASH_m01s04i203_2' : 'LargeScaleRain',
    }


box = [325, 475, 93, 235 ]  # x1, x2, y1, y2
for f in flist:
    vegdat = xr.open_dataarray(veg)
    varsdat = xr.open_dataset(f)

    vegdat = vegdat.isel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))
    varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))

    varsdat = varsdat[list(name_dict.keys())]

    varsdat.rename(name_dict, inplace=True)

    varsdat.coords['grid_latitude_t'] = varsdat.latitude_t.values.mean(axis=1)
    varsdat.coords['grid_longitude_t'] = varsdat.longitude_t.values.mean(axis=0)

    varsdat['forest_frac'] = (('pseudo_dim', 'grid_latitude_t', 'grid_longitude_t'), vegdat.values[0,:,:][None,...])


    newname = f.replace('pc', 'small')
    varsdat.to_netcdf(newname)

