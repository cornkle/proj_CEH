import xarray as xr
import pdb
import numpy as np
from scipy.interpolate import griddata
import glob
import os


def regrid(lon, lat, new_lon, new_lat, data):

    if new_lon.ndim == 1:
        grid_lons, grid_lats = np.meshgrid(new_lon, new_lat)
    else:
        grid_lons = new_lon
        grid_lats = new_lat

    points = np.array((grid_lons.flatten(), grid_lats.flatten())).T
    inter = np.array((np.ravel(lon), np.ravel(lat))).T

    # Interpolate using delaunay triangularization
    data = griddata(points, data.flatten(), inter, method='linear')
    data = data.reshape((grid_lons.shape[0], grid_lons.shape[1]))

    return data


veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'

### 2d vars , xmh*.pc*.nc files
folder = '/users/global/cornkle/w2018_bamba/linked/'

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

cleanup = {



}
box = [325, 475, 93, 235 ]  # x1, x2, y1, y2 Tai park box

flist = glob.glob(folder+'*.nc')

vegdat = xr.open_dataarray(veg)
vegdat = vegdat.isel(longitude=slice(box[0], box[1]), latitude=slice(box[2], box[3]))

for f in flist:

    fname = f.split(os.sep)[-1]
    varsdat = xr.open_dataset(f)
    varsdat = varsdat.isel(grid_longitude_t=slice(box[0], box[1]), grid_latitude_t=slice(box[2], box[3]))

    if atmo_bstream in fname:
        varsdat = varsdat[list(dict_bstream.keys())]
        varsdat.rename(dict_bstream, inplace=True)
    if flux_cstream in fname:
        varsdat = varsdat[list(dict_cstream.keys())]
        varsdat.rename(dict_cstream, inplace=True)


    #only defendable for small box
    varsdat.coords['grid_latitude_t'] = varsdat.latitude_t.values.mean(axis=1)
    varsdat.coords['grid_longitude_t'] = varsdat.longitude_t.values.mean(axis=0)

    varsdat['forest_frac'] = (('pseudo_dim', 'grid_latitude_t', 'grid_longitude_t'), vegdat.values[0,:,:][None,...])


    newname = f.replace('pc', 'small')
    varsdat.to_netcdf(newname)

