import xarray as xr
import pdb
import numpy as np
import glob
import os
from utils import u_interpolate as uint


#veg = '/users/global/cornkle/w2018_bamba/qrparm.cci.4km.nc'

### 2d vars , xmh*.pc*.nc files
folder = '/scratch/ssf/'
out = '/users/global/cornkle/w2018_bamba/regrid_raw/'

atmo_bstream = '.pb'
current = 'current'
past = 'past'
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

flist = glob.glob(folder+'current*.pb*.nc')

dummy = xr.open_dataset(flist[0])
box = [0, 690, 0, 300 ]

dummy = dummy.isel(grid_latitude_t=slice(box[2],box[3]), grid_longitude_t=slice(box[0],box[1]),
                   grid_latitude_uv=slice(box[2],box[3]), grid_longitude_uv=slice(box[0],box[1]))

lats = dummy.latitude_t[:, 0].values
lons = dummy.longitude_t[0, :].values
latmin = dummy.latitude_t.min().values
latmax = dummy.latitude_t.max().values
lonmin = dummy.longitude_t.min().values
lonmax = dummy.longitude_t.max().values
dist = np.round(np.float(np.mean(lats[1::] - lats[0:-1])), decimals=4)

lat_regular = np.arange(latmin+10*dist, latmax - 10*dist, dist)
lon_regular = np.arange(lonmin+10*dist, lonmax - 10*dist, dist)

inds_uv, weights_uv, shr = uint.interpolation_weights(dummy.longitude_uv.values, dummy.latitude_uv.values, lon_regular, lat_regular)

for tt in timex:
    flist = glob.glob(folder + tt+'*.pb*.nc')
    dates = []
    for f in flist:
        strs = f.split('.')
        dates.append(strs[1][2:10])
    dates = np.unique(dates)

    for d in dates:

        ds = xr.Dataset(coords={'lat': lat_regular, 'lon': lon_regular})

        for f in flist:

            fname = f.split(os.sep)[-1]
            varsdat = xr.open_dataset(f)


            if 'time' not in ds.keys():
                time = varsdat.TH1.values
                ds['time'] = time


            varsdat = varsdat[list(units.keys())]

            varsdat = varsdat.isel(grid_longitude_uv=slice(box[0], box[1]), grid_latitude_uv=slice(box[2], box[3]))
            varsdat = varsdat.sel(P_ECMWF=slice(975,200))

            vkeys = varsdat.keys()
            for key in vkeys:

                if ('longitude' in key) | ('latitude' in key) | ('TH1' in key) | ('height' in key) | ('P_ECMWF' in key) :
                    continue
                print('Doing variable', key)

                data = varsdat[key].squeeze().values

                times = []
                for dat in data:
                    regridded = uint.interpolate_data(dat, inds_uv, weights_uv, shr)

                    times.append(regridded[None, ...])
                regridded = np.concatenate(times, axis=0)

                unit = units[key]
                ds[key] = (('time', 'plevel', 'lat', 'lon'), regridded)
                ds[key].attrs['unit'] = unit

            if current in fname:
                fname = fname.replace(current, 'current')
            if past in fname:
                fname = fname.replace(past, 'past')


        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in ds.data_vars}
        ds.to_netcdf(out+tt+'_3d_'+d+'.nc', format='NETCDF4', encoding=encoding)
        ds.close()


