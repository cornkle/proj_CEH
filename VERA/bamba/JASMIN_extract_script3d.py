import xarray as xr
import pdb
import numpy as np
import glob
import os
from utils import u_interpolate as uint


#veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'

### 2d vars , xmh*.pc*.nc files
folder = '/scratch/ssf/'
out = '/users/global/cornkle/w2018_bamba/linked/small/'

atmo_bstream = '.pb'
current = 'xmhkga'
past = 'xmhkha'
timex = [current, past]

dict_cstream ={
    'STASH_m01s30i203' : 'w_pl',
    'STASH_m01s30i201' : 'u_pl',
    'STASH_m01s30i202' : 'v_pl',
    'STASH_m01s30i204' : 'T_pl',
    'STASH_m01s30i205' : 'q_pl',
    'STASH_m01s30i206' : 'rh_pl',
    'STASH_m01s30i208' : 'omega_pl'
    }

units = {
         'w_pl' : 'm s-1',
         'u_pl' : 'm s-1',
         'v_pl' : 'm s-1',
         'T_pl' : 'K',
         'q_pl' : 'kg kg-1',
         'rh_pl' : '%',
         'omega_pl': 'Pa s-1'
         }

box = [320, 500, 93, 250 ]  # x1, x2, y1, y2 Tai park box

flist = glob.glob(folder+'xmhk*.nc')

dummy = xr.open_dataset(flist[0])
dummy = dummy.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))

lats = dummy.latitude_t[:, 0]
lons = dummy.longitude_t[0, :]
latmin = dummy.latitude_t.min()
latmax = dummy.latitude_t.max()
lonmin = dummy.longitude_t.min()
lonmax = dummy.longitude_t.max()
dist = np.round(np.float(np.mean(lats[1::].values - lats[0:-1].values)), decimals=4)

lat_regular = np.arange(latmin+5*dist, latmax - 5*dist, dist)
lon_regular = np.arange(lonmin+5*dist, lonmax - 5*dist, dist)

# interpolation weights for new grid
inds, weights, shape = uint.interpolation_weights(lons, lats, lon_regular, lat_regular)

for tt in timex:
    flist = glob.glob(folder + tt + atmo_bstream +'*.nc')

    for f in flist[0:2]:

        fname = f.split(os.sep)[-1]
        varsdat = xr.open_dataset(f)
        varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]),
                               P_ECMWF=slice(950, 550), grid_longitude_uv=slice(box[0], box[1]),
                               grid_latitude_uv=slice(box[2], box[3]))

        time = varsdat.TH1_MN.values
        plevels = varsdat.P_ECMWF.values
        ds = xr.Dataset(coords={'time': time, 'lat': lat_regular, 'lon': lon_regular, 'plevels' : plevels} )

        varsdat = varsdat[list(dict_cstream.keys())]
        varsdat.rename(dict_cstream, inplace=True)

        vkeys = varsdat.keys()
        for key in vkeys:

            if ('longitude' in key) | ('latitude' in key) | ('TH1' in key) :
                continue
            print('Doing variable', key)
            data = varsdat[key].values
            regridded = uint.interpolate_data(data, inds, weights, shape)
            unit = units[key]
            ds[key] = (('time', 'plevels', 'lat', 'lon'), regridded)
            ds[key].attrs['unit'] = unit

        # varsdat['forest_frac'] = (('pseudo_dim', 'grid_latitude_t', 'grid_longitude_t'), vegdat.values[0,:,:][None,...])

        if current in fname:
            fname = fname.replace(current, 'current')
        if past in fname:
            fname = fname.replace(past, 'past')

        ds.to_netcdf(out+fname)
        ds.close()



