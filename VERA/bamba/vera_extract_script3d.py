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
flux_cstream = '.pc'
current = 'xmhkga'
past = 'xmhkha'

keep_vars2d = ['STASH_m01s00i024', 'STASH_m01s03i237', 'STASH_m01s03i236', 'STASH_m01s16i222', 'STASH_m01s03i226' , 'STASH_m01s03i225'] # LST, spec hum 1.5m, T 1.5m, p at sl,10m v wind, 10m u wind
keep_vars3d = ['STASH_m01s30i203', 'STASH_m01s30i201', 'STASH_m01s30i202', 'STASH_m01s30i204', 'STASH_m01s30i205', 'STASH_m01s30i206', 'STASH_m01s30i208' ] # u on P, v on P, T on P, spec H on P, rel H on P, omega on P

vars2d_names = ['lst', 'q2', 'T2', 'slp', 'v10', 'u10']
vars3d_names = ['w_pl', 'u_pl', 'v_pl', 'T_pl', 'q_pl', 'rh_pl', 'omega_pl']
joined_names = vars2d_names+vars3d_names
dict_bstream = {}
for id, n in enumerate(keep_vars2d+keep_vars3d):
    dict_bstream[n] = joined_names[id]

dict_cstream ={
   'STASH_m01s01i201_2' : 'sw_net',
   'STASH_m01s01i235_2' : 'sw_in',
   'STASH_m01s02i201_2' : 'lw_net',
   'STASH_m01s02i207_2' : 'lw_in',
   'STASH_m01s03i217_2' : 'SH',
   'STASH_m01s03i234_2' : 'LH',
   'STASH_m01s05i205_2' : 'ConvRain',
   'STASH_m01s05i216_2' : 'TotalRain',
   'STASH_m01s04i203_2' : 'LargeScaleRain',
    }

units = {'sw_net': 'W m-2',
         'sw_in': 'W m-2',
         'lw_net': 'W m-2',
         'lw_in': 'W m-2',
         'SH': 'W m-2',
         'LH': 'W m-2',
         'ConvRain': 'kg m-2s-1',
         'TotalRain': 'kg m-2s-1',
         'LargeScaleRain': 'kg m-2s-1',
         'lst' : 'K' ,
         'q2' : 'kg kg-1',
         'T2' : 'K',
         'slp' : 'Pa',
         'v10' : 'm s-1',
         'u10' : 'm s-1',
         'w_pl' : 'm s-1',
         'u_pl' : 'm s-1',
         'v_pl' : 'm s-1',
         'T_pl' : 'K',
         'q_pl' : 'kg kg-1',
         'rh_pl' : '%',
         'omega_pl': 'Pa s-1'
         }

box = [325, 475, 93, 235 ]  # x1, x2, y1, y2 Tai park box

flist = glob.glob(folder+'xmhk*.nc')

# vegdat = xr.open_dataarray(veg)
# vegdat = vegdat.isel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))

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

for f in flist[0:2]:

    fname = f.split(os.sep)[-1]
    varsdat = xr.open_dataset(f)
    varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))
    time = varsdat.TH1_MN.values

    ds = xr.Dataset(coords={'time': time, 'lat': lat_regular, 'lon': lon_regular} )
    if atmo_bstream in fname:
        varsdat = varsdat[list(dict_bstream.keys())]
        varsdat.rename(dict_bstream, inplace=True)
    if flux_cstream in fname:
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
        ds[key] = (('time', 'lat', 'lon'), regridded)
        ds[key].attrs['unit'] = unit

    # varsdat['forest_frac'] = (('pseudo_dim', 'grid_latitude_t', 'grid_longitude_t'), vegdat.values[0,:,:][None,...])

    if current in fname:
        fname.replace(current, 'current')
    if past in fname:
        fname.replace(past, 'past')

    ds.to_netcdf(out+fname)
    pdb.set_trace()


